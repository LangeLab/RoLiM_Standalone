import os
import re

import pymysql.cursors
import pandas as pd
import numpy as np

from src import plots, utils, extraction, sequence

from src.globals import (
    COMPOUND_RESIDUES,
    SITES, MEROPS_POSITIONS, 
    THREE_LETTER_CODES, SINGLE_LETTER_CODES
)

def connect_to_merops_database():
    '''
    Connect to the MEROPS database.
    Parameters:
        None

    Returns:
        merops_connection -- 
    '''

    # Connect to MEROPS database.
    merops_connection = pymysql.connect(
        host='localhost',
        user='zigzag',
        password='MyStrongP@ss123', # Use your own password.
        db='merops',
        cursorclass=pymysql.cursors.DictCursor
    )

    return merops_connection

def partition_proteases(rows,fields,columns):
    '''

    Parameters:
        rows -- 
        fields --
        columns -- 
    Returns:
        partitioned_results -- 
    '''

    # Partition result set by protease.
    proteases = {}
    for row in rows:
        df_row = []
        for field in fields:
            df_row.append(row[field])
        protease = row['Protease']
        if protease in proteases.keys():
            proteases[protease].append(df_row)
        else:
            proteases[protease] = [df_row]

    # Generate Pandas DataFrame for each protease and append
    # to list of protease substrate DataFrames.
    partitioned_results = {}    
    for protease,df_rows in proteases.items():
        partitioned_results[protease] = pd.DataFrame(df_rows,columns=columns)

    return partitioned_results


def retrieve_substrates(
        names=[],
        ids=[],
        organisms=['Homo sapiens']
    ):
    '''
    Connects to MEROPS, retrieves protease substrate sequences,
        and returns a Pandas DataFrame containing the sequences
        for each protease in the Substrate_search table.
    
    Parameters:
        names -- List. Default empty. Optional list of protease
                        names used to filter results.
        ids -- List. Default empty. Optional list of substrate IDs
                        used to filter results.
        organisms -- List. Default ['Homo sapiens']. Optional list
                            of organisms used to filter results.
                            Homo sapiens must be overridden with empty
                            list or alternative organisms.

    Returns:
        protease_substrates -- Dict. Maps protease name to Pandas 
                                    DataFrame for each protease with
                                    associated substrutes in MEROPS.
    '''

    # Connect to MEROPS database.
    connection = connect_to_merops_database()

    # Get all substrates from MEROPS Substrate_search table.
    try:
        with connection.cursor() as cursor:
            conditional_string = ""
            if names != []:
                conditional_string += (
                    " AND Substrate_search.Protease IN ('"
                    + "','".join(names)
                    + "')"
                )
            if ids != []:
                conditional_string += (
                    " AND Substrate_search.cleavageID IN ('"
                    + "','".join(ids)
                    + "')"
                )
            if organisms != []:
                conditional_string += (
                    " AND Substrate_search.organism IN ('"
                    + "','".join(organisms)
                    + "')"
                ) 
            
            for site in SITES:
                conditional_string += (
                    " AND %s" % site
                    + " IN ('"
                    + "' ,'".join(THREE_LETTER_CODES)
                    + "')"
                )
            
            query = (
                "SELECT Substrate_search.Protease, {0} FROM Substrate_search"
                + " INNER JOIN (SELECT code FROM activity_status WHERE organism"
                + " = 'human' AND status = 'active') as human_active"
                + " ON Substrate_search.code = human_active.code WHERE cleavage_type"
                + " NOT IN ('synthetic','theoretical'){1}"
            ).format(", Substrate_search.".join(SITES), conditional_string)
            cursor.execute(query)
            all_substrates = cursor.fetchall()
    
    finally:
        connection.close()
    
    protease_substrates = partition_proteases(
        all_substrates, SITES, MEROPS_POSITIONS
    )

    return protease_substrates  

def retrieve_vectorized_substrates(
        names=[]
    ):
    '''
    Retrieves vectorized substrates from vectorized_substrates
        table in merops MySQL database.

    Parameters:
        names -- List. Default empty. Optional list of protease
                        names used to filter results.

    Returns:
        vectorized_substrates -- Dict.
    '''

    conditional_string = ""
    if names != []:
        conditional_string += " WHERE Protease in ('%s')" % ("', '".join(names))

    connection = connect_to_merops_database()
    try:
        with connection.cursor() as cursor:
            query = "SELECT * FROM vectorized_substrates%s" % (conditional_string)
            cursor.execute(query)
            all_substrates = cursor.fetchall()
    finally:
        connection.close()

    fields = []
    for site in SITES:
        for residue in SINGLE_LETTER_CODES:
            fields.append(site+"_"+residue)

    positions = []
    for position in MEROPS_POSITIONS:
        for residue in SINGLE_LETTER_CODES:
            positions.append(position+"_"+residue)

    vectorized_substrates = partition_proteases(
        all_substrates, fields, positions
    )

    return vectorized_substrates

def retrieve_protease_patterns(
        names=[], 
        enable_compound_residues=True
    ):
    '''
    Retrieve protease substrate patterns from
        protease_patterns table in merops database.
        Query optionally filtered by protease name.

    Parameters:
        names -- List. Default empty.
    Returns:
        protease_patterns -- Dict,
    '''

    connection = connect_to_merops_database()

    if enable_compound_residues:
        protease_pattern_table = 'protease_patterns'
    else:
        protease_pattern_table = 'protease_patterns_no_cr'

    if names != []:
        conditional_string = " WHERE Protease in ('%s')" % ("', '".join(names))
    else:
        conditional_string = ""

    protease_patterns = {}

    try:
        with connection.cursor() as cursor:
            query = "SELECT Protease, pattern FROM {0}{1}".format(
                protease_pattern_table,
                conditional_string
            )
            cursor.execute(query)
            while True:
                row = cursor.fetchone()
                if row:
                    if row['Protease'] in protease_patterns:
                        protease_patterns[row['Protease']].append(row['pattern'])
                    else:
                        protease_patterns[row['Protease']] = [row['pattern']]
                else:
                    break
    finally:
        connection.close()

    return protease_patterns

def insert_vectorized_substrate(substrate):
    '''

    Parameters:
        substrate --
    
    Returns:
        None
    '''

    columns = ['Protease']
    for site in SITES:
        for residue in SINGLE_LETTER_CODES:
            columns.append(site+"_"+residue)

    connection = connect_to_merops_database()

    try:
        with connection.cursor() as cursor:
            query = (
                "INSERT INTO vectorized_substrates ({0}) VALUES ('{1}}',{2})"
            ).format(", ".join(columns), substrate[0], ", ".join(substrate[1:]))
            cursor.execute(query)
        connection.commit()
    finally:
        connection.close()

def vectorize_substrates(substrates):
    '''
    Pulls substrates from merops database Substrate_search
        table, converts to single-letter amino acid encoding,
        inserts substrates into vectorized_substrates table.
    
    Parameters:
        substrates --

    Returns:
        None
    '''

    vector_template = pd.Series(
        np.zeros(len(THREE_LETTER_CODES), dtype=np.int8),
        index=THREE_LETTER_CODES
    )

    for protease,sequences in substrates.items():
        for i,sequence in sequences.iterrows():
            skip_sequence = False
            row = [protease.strip("'\"")]
            for residue in sequence:
                positional_vector = vector_template.copy()
                if residue in positional_vector.index:
                    positional_vector[residue] = 1
                else:
                    skip_sequence = True
                    break
                row += positional_vector.tolist()
            if skip_sequence:
                print('Skipping sequence')
                continue
            row = [str(i) for i in row]
            print(row)
            insert_vectorized_substrate(row)

def create_vectorized_substrates_table():
    '''
    Creates table for vectorized substrates containing
        a positional column for each residue in backgroun
        data.

    Parameters:
        None

    Returns:
        None
    '''

    # Generate column name and data type string.
    columns = [
        "id int NOT NULL AUTO_INCREMENT",
        "Protease varchar(100) NOT NULL",
    ]
    for site in SITES:
        for residue in SINGLE_LETTER_CODES:
            columns.append(site + "_" + residue + " int NOT NULL")
    columns.append("PRIMARY KEY (id)")

    # Open merops database connection and create vectorized_substrates table.
    connection = connect_to_merops_database()
    try:
        with connection.cursor() as cursor:
            query = (
                "CREATE TABLE vectorized_substrates ({}) "
            ).format(", ".join(columns))
            cursor.execute(query)
        connection.commit()
    finally:
        connection.close() 
    print(query)

def retrieve_protease_family_code(protease):
    '''

    Parameters:
        protease -- String.

    Returns:
        protease_family_code -- String.
    '''

    connection = connect_to_merops_database()

    try:
        with connection.cursor() as cursor:
            query = (
                "SELECT DISTINCT code FROM Substrate_search"
                + " WHERE Protease = '{}'"
            ).format(protease)
            cursor.execute(query)
            protease_code = cursor.fetchone()["code"]
    finally:
        connection.close()

    return protease_code[:protease_code.find('.')]

def create_protease_patterns_table(enable_compound_residues=True):
    '''
    Creates protease_patterns table in merops database.

    Parameters:
        None
    
    Returns:
        None
    '''

    connection = connect_to_merops_database()

    try:
        with connection.cursor() as cursor:
            if enable_compound_residues:
                protease_pattern_table = 'protease_patterns'
            else:
                protease_pattern_table = 'protease_patterns_no_cr'
            query = (
                "CREATE TABLE {} (id int NOT NULL AUTO_INCREMENT,"
                + " Protease varchar(100) NOT NULL, pattern varchar(64) NOT NULL,"
                + " PRIMARY KEY (id))"
            ).format(protease_pattern_table)
            cursor.execute(query)
        connection.commit()
    finally:
        connection.close()

def insert_protease_patterns(protease_patterns, enable_compound_residues=True):
    '''

    Parameters:
        protease_patterns -- Dict.

    Returns:
        None
    '''

    for protease,patterns in protease_patterns.items():
        for pattern in patterns:
            connection = connect_to_merops_database()
            try:
                with connection.cursor() as cursor:
                    if enable_compound_residues:
                        protease_pattern_table = 'protease_patterns'
                    else:
                        protease_pattern_table = 'protease_patterns_no_cr'
                    query = (
                        "INSERT INTO {0} (Protease, pattern)"
                        + " VALUES ('{1}', '{2}')"
                    ).format(protease_pattern_table, protease, pattern)
                    cursor.execute(query)
                connection.commit()
            finally:
                connection.close()

def swap_protease_patterns(target):
    '''

    Parameters:

    Returns:

    '''

    connection = connect_to_merops_database()

    try:
        with connection.cursor() as cursor:
            cursor.execute("DROP TABLE protease_patterns")
            cursor.execute("CREATE TABLE protease_patterns LIKE protease_patterns_06")
            cursor.execute("INSERT INTO protease_patterns SELECT * FROM %s" % target)
        connection.commit()
    finally:
        connection.close()

def generate_merops_heatmap(
        pattern_container, 
        position_labels
    ):
    """
    Generate heatmap for each pattern in MEROPS database with detected
        patterns.

    Parameters:
        pattern_container -- PatternContainer instance.
        position_labels -- String.

    Returns:
        None
    """
    protease_pattern_heatmap_title = (
        pattern_container.title + ' - Protease Pattern Matches (Percent Positions Matched)'
    )
    output_prefix = re.sub(r'\W+', ' ', pattern_container.title).strip().replace(" ", "_")
    # Set output path for absolute frequency protease pattern heat map.
    protease_pattern_heatmap_output_path = (
        pattern_container.output_directory
        + '/figures/'
        + output_prefix
        + '_protease_pattern_heatmap.svg'
    )
    if pattern_container.background.compound_residues is not None:
        enable_compound_residues = True
    else:
        enable_compound_residues = False
    protease_patterns = retrieve_protease_patterns(
        enable_compound_residues=enable_compound_residues
    )
    protease_labels = plots.generate_protease_labels(protease_patterns)
    non_exact_scoring_matrix = plots.generate_non_exact_protease_pattern_matrix(
        pattern_container,
        protease_patterns,
        protease_labels
    )
    # Generate absolute frequency protease pattern heatmap.
    plots.generate_protease_pattern_heatmap(
        protease_pattern_heatmap_title,
        pattern_container,
        non_exact_scoring_matrix,
        protease_labels,
        position_labels,
        protease_pattern_heatmap_output_path
    )

    # Export heatmap data as a tab-separated table.
    pattern_labels = plots.generate_pattern_labels(position_labels, pattern_container)
    protease_pattern_table_output_path = (
        pattern_container.output_directory
        + '/figures/'
        + output_prefix
        + '_protease_pattern_heatmap.txt'
    )
    plots.export_protease_heatmap_table(
        non_exact_scoring_matrix,
        pattern_labels,
        protease_labels,
        protease_pattern_table_output_path
    )

def extract_protease_substrate_patterns(
        background=None,
        percentage_frequency_cutoff=0.05,
        max_depth=None,
        p_value_cutoff=0.001,
        minimum_occurrences=2,
        fold_change_cutoff=1,
        multiple_testing_correction=True,
        positional_weighting=True,
        allow_compound_residue_decomposition=True,
        enable_compound_residues=True,
        position_specific=True,
        width=8,
        center=True
    ):
    '''
    Run pattern extraction on MEROPS protease substrate data sets.
    
    Parameters:
        

    Returns:
        None
    '''

    # Replace existing protease_patterns table with empty table.
    connection = connect_to_merops_database()

    if enable_compound_residues:
        protease_pattern_table = 'protease_patterns'
    else:
        protease_pattern_table = 'protease_patterns_no_cr'

    with connection.cursor() as cursor:
        cursor.execute("SHOW TABLES LIKE '{}'".format(protease_pattern_table))
        if cursor.fetchone():
            cursor.execute("DROP TABLE {}".format(protease_pattern_table))
            connection.commit()
    create_protease_patterns_table(
        enable_compound_residues=enable_compound_residues
    )

    if enable_compound_residues:
        compound_residues = COMPOUND_RESIDUES
    else:
        compound_residues = None

    # If no background is specified, use default SwissProt background.
    if background == None:
        # Load context sequences from fasta.
        context = utils.parse_fasta(
            # os.path.join('media', 'defaults', 'uniprot.fasta')
            # TODO: Set a default fasta and use path to it.
        )
        # Generate new Background instance.
        background = sequence.Background(
            context['sequence'].tolist(),
            compound_residues=compound_residues,
            position_specific=position_specific,
            width=width,
            center=center
        )

    # Get substrates for all MEROPS proteases and loop through the sets.
    pattern_containers = []
    protease_substrates = retrieve_substrates()
    for protease, substrates in protease_substrates.items():
        print('{}\n'.format(protease))
        # Initialize output directory name.
        if enable_compound_residues:
            protease_output_directory = os.path.join(
                'media', 'merops', 'proteases', protease.replace(" ", "_")
            )
        else:
            protease_output_directory = os.path.join(
                'media', 'merops', 'proteases_no_cr', protease.replace(" ", "_")
            )
        # Prepare substrate sets for pattern extraction.
        single_letter_substrates = utils.convert_encoding(
            substrates,
            1
        )
        single_letter_substrates.drop_duplicates(inplace=True)
        substrate_tensor = sequence.vectorize_sequences(
            single_letter_substrates,
            background
        )
        protease_substrate_sample = sequence.Sample(
            sequence_df=single_letter_substrates,
            sequence_tensor=substrate_tensor
        )
        
        min_occurrences = (
            percentage_frequency_cutoff
            * protease_substrate_sample.sequence_tensor.shape[0]
        )
        if min_occurrences < 2:
            min_occurrences = 2
        
        # Run pattern extraction on substrate set.
        patterns = extraction.PatternContainer(
            protease_substrate_sample,
            background,
            protease,
            protease_output_directory,
            minimum_occurrences=min_occurrences
        )

        # Post-process extracted patterns and save pattern outputs.
        patterns.prune_patterns()

        if len(patterns.pattern_list) > 0:
            # Prepare output directories.
            try:
                os.makedirs(protease_output_directory)
            except FileExistsError:
                pass
            patterns.generate_pattern_outputs()

            # Generate clustermap for protease substrate set. 
            clustermap_title = (
                protease
                + ' Substrates -'
                + ' Mean Sequence-Pattern Positional Substitution Probability'
            )
            clustermap_output_directory = protease_output_directory + '/figures'
            try:
                os.makedirs(clustermap_output_directory)
            except:
                pass
            clustermap_output_path = clustermap_output_directory + '/clustermap.svg'

            position_labels = plots.generate_position_labels(
                protease_substrate_sample.sequence_df
            )
            pattern_labels = [
                label[:(label.find('{') - 4)] for label in plots.generate_pattern_labels(
                    position_labels,
                    patterns
                )
            ]
            pattern_similarity_matrix = plots.calculate_pattern_similarity_matrix(
                protease_substrate_sample.sequence_df,
                patterns.pattern_list,
                pattern_labels,
                plots.SUBSTITUTION_MATRIX
            )

            if np.any(pattern_similarity_matrix.to_numpy().astype(np.bool)):
                try:
                    plots.generate_sequence_clustermap(
                        clustermap_title,
                        pattern_similarity_matrix,
                        clustermap_output_path
                    )
                except ValueError:
                    print('clustermap for {} failed.'.format(protease))

            clusters = [
                ','.join(pattern.character_pattern().tolist())
                for pattern in patterns.pattern_list
            ]

            substrate_patterns = {protease: clusters}
            insert_protease_patterns(
                substrate_patterns,
                enable_compound_residues=enable_compound_residues
            )
            pattern_containers.append(patterns)

    for pattern_container in pattern_containers:
        generate_merops_heatmap(pattern_container, position_labels)