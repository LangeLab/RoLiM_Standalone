import os
import re
import csv
import json
import numpy as np
import pandas as pd
from collections import Counter
from typing import List, Dict, Optional, Union

from src import utils
from src.globals import CompoundResidue, Sample, ExtendedSequences
from src.globals import (
    COMPOUND_RESIDUES, SWISSPROT_ACCESSION_PATTERN, PRECOMPUTED_FILES
)

def detect_delimiter(file_path):
    # Detect delimiter from file extension (supports comma or tab).
    file_extension = file_path[-file_path[::-1].find('.'):].lower()
    if file_extension == 'csv':
        delimiter = ','
    else:
        delimiter = '\t'

    return delimiter

def generate_positions(center, num_positions):
    """
    Generate position labels for sequence dataframe columns or series
        indexes.

    Parameters:
        center -- Boolean. If true, column labels will count away from
                        center (leftward direction = non-prime,
                        rightward direction = prime). If false, column
                        labels will count from left to right beginning
                        with 1.
        num_positions -- Int. Total number of position labels.

    Returns:
        positions -- List. Contains position labels.
    """
    if center == True:
        positions = []

        p1_pos = int(num_positions / 2)

        for i in range(num_positions):
            if (i<p1_pos):
                position_label = 'p' + str(p1_pos-i)
                positions.append(position_label)
            else:
                position_label = 'p' + str(abs(p1_pos-(i+1))) + "'"
                positions.append(position_label)
    else:
        positions = ['p{}'.format(i) for i in range(1, num_positions + 1)]

    return positions

def sequences_to_df(
        sequences,
        center=True,
        redundancy_level='sequence'
    ):
    """
    Split sequence strings to Pandas DataFrame with one column per
        position.

    Parameters:
        sequences -- List-like object. Contains one sequence per row.
        center -- Boolean. If true, column labels will count away from
                        center (leftward direction = non-prime,
                        rightward direction = prime). If false, column
                        labels will count from left to right beginning
                        with 1.
    Returns:
        sequence_df -- Pandas DataFrame. Contains input data set split
                        to columns by position and accompanying fields
                        as specified above.
    """
    split_sequences = [list(sequence.rstrip().upper()) for sequence in sequences]
    cols = generate_positions(center, len(split_sequences[0]))
    sequence_df = pd.DataFrame(split_sequences, columns=cols)
    if redundancy_level == 'sequence':
        sequence_df.drop_duplicates(inplace=True)

    return sequence_df

def get_residue_index(residue, residue_dict):
    """

    Parameters:
        residue -- String.
        residue_dict -- Dict.

    Returns:
        residue_index -- Int.
    """

    try:
        residue_index = residue_dict[residue]
    except KeyError:
        residue_index = -1

    return residue_index

def vectorize_sequences(
        data,
        background,
        empty_position_value=0
    ):
    """

    Parameters:
        data -- Pandas DataFrame.
        background -- Background instance.
        empty_position_value -- Int.

    Returns:
        sequence_tensor -- Numpy array. Dtype=np.int8
    """

    residue_dict = {residue: i for i, residue in enumerate(background.ordered_residues)}
    # Construct tensor from residue indices.
    index_array = data.map(
        lambda x: get_residue_index(x, residue_dict)
    ).to_numpy(dtype=np.int8)
    sequence_tensor = np.eye(len(background.ordered_residues), dtype=np.int8)[index_array]
    # Remove unknown residues from tensor.
    sequence_tensor[index_array == -1] = 0
    # Add compound residues to tensor.
    if background.compound_residues is not None:
        sequence_tensor[:, :, len(background.alphabet):] = np.inner(
            sequence_tensor[:, :, :len(background.alphabet)],
            background.compound_residue_matrix
        )

    return sequence_tensor

def empty_position_vector(length, empty_position_value=0):
    '''

    Parameters:
        length -- Int. Length of alphabet.
        empty_position_value -- Int. Expects value of 0 or 1. If 0,
                                matching amino acids coded as 1. If 0,
                                matching amino acids coded as 1.

    Returns:
        empty_position -- 1D Numpy Array.
    '''

    if empty_position_value == 1:
        empty_position = np.ones(length, dtype='int8')
    elif empty_position_value == 0:
        empty_position = np.zeros(length, dtype='int8')

    return empty_position

def get_pattern_constituents(pattern_matrix, background):
    """
    Return pattern featuring all original residues, including compound
        residues, plus compound residue constituents.

    Parameters:
        pattern_matrix -- 2D Numpy Array.
        background -- Background instance.

    Returns:
        constituent_pattern -- 2D Numpy Array.
    """

    try:
        # Extend compound residues with blank.
        constituents = np.concatenate(
            (
                background.compound_residue_matrix,
                np.zeros(
                    (1, background.compound_residue_matrix.shape[1]), dtype=np.int8
                )
            ),
            axis=0
        )
    except AttributeError as e:
        raise e
    else:
        # Generate index array from pattern matrix.
        compound_residue_index = np.argmax(pattern_matrix, axis=1) - len(background.alphabet)

        # Map simple positions to blank in compound residue matrix.
        compound_residue_index[
            compound_residue_index < 0] = len(background.compound_residue_matrix)
        
        # Get compound residue constituents using index array.
        positional_residue_constituents = constituents[compound_residue_index]

        # Combine constituents with simple fixed residues.
        fixed_position_residues = np.logical_or(
            pattern_matrix[:, :len(background.alphabet)],
            positional_residue_constituents)

        constituent_pattern = np.concatenate(
            (
                fixed_position_residues,
                pattern_matrix[:, len(background.alphabet):]
            ),
            axis=1
        )

        return constituent_pattern

def vectorize_pattern(
        pattern,
        background,
        empty_position_value=0,
        add_constituents=True,
        add_compound_residues=True
    ):
    '''
    Takes pattern of Pandas Series type, converts single
        letter encoding to binary encoding. Returns vectorized
        pattern as Numpy array.

    Parameters:
        pattern -- Pandas Series.
        background -- Instance of Background class.
        empty_position_value -- Int. Expects value of 0 or 1. If 0,
                                matching amino acids coded as 1. If 0,
                                matching amino acids coded as 1.
    
    Returns:
        pattern_matrix -- Numpy Array.
    '''

    pattern_positions = []
    for position in pattern:
        residue = position.strip('[]').upper()
        
        position_vector = pd.Series(
            empty_position_vector(
                len(background.ordered_residues),
                empty_position_value
            ),
            index=background.ordered_residues
        )

        # Encode residue with non-null value.
        if residue in background.ordered_residues:
            position_vector[residue] = int(not empty_position_value)

        # Add positional residues to pattern.
        positional_residues = position_vector.tolist()
        pattern_positions += positional_residues

    # Convert pattern to numpy array.
    pattern_matrix = np.array(pattern_positions, dtype='int8').reshape(
        len(pattern),
        len(background.ordered_residues)
    )

    try:
        if add_constituents:
            # Add constituents to positions containing compound residue
            constituent_pattern_matrix = get_pattern_constituents(pattern_matrix, background)
        else:
            background.compound_residue_matrix
    except AttributeError:
        pass
    else:
        if add_compound_residues:
            # Add compound residues matching constituents.
            compound_residues = np.inner(
                pattern_matrix[:, :len(background.alphabet)],
                background.compound_residue_matrix
            )
            compound_residue_pattern_matrix = pattern_matrix.copy()
            compound_residue_pattern_matrix[:, len(background.alphabet):] = compound_residues

            # Combine compound residue and constituent pattern
            pattern_matrix = np.logical_or(
                constituent_pattern_matrix,
                compound_residue_pattern_matrix
            ).astype(np.int8)
        elif add_constituents:
            pattern_matrix = constituent_pattern_matrix

    return pattern_matrix

def load_prealigned_file(
        prealigned_file_path,
        background,
        center=True,
        redundancy_level='none',
        title=''
    ):
    """
    Top-level helper function to load pre-aligned sequences from text
        file.
    
    Parameters:
        prealigned_file_path -- String. Path to prealigned file.
        background -- Background instance.
        center -- Boolean. Count sequence positions away from center.

    Returns:
        sample -- Sample objects as values.
    """
    
    delimiter = detect_delimiter(prealigned_file_path)
    prealigned_sequences = pd.read_csv(
        prealigned_file_path,
        delimiter=delimiter,
        header=0
    )

    if len(prealigned_sequences.columns) == 1:
        prealigned_sequences.columns = ['aligned_sequence']
        sequence_df = sequences_to_df(
            prealigned_sequences['aligned_sequence'].tolist(),
            center=center,
            redundancy_level=redundancy_level
        )
        sequence_tensor = vectorize_sequences(sequence_df, background)
        samples = {
            str(title): Sample(
                sequence_df=sequence_df,
                sequence_tensor=sequence_tensor,
                original_sequences=prealigned_sequences
            )
        }
    else:
        prealigned_sequences.rename(
            columns={
                prealigned_sequences.columns[0]: 'input_sample_name',
                prealigned_sequences.columns[1]: 'aligned_sequence'
            },
            inplace=True
        )
        prealigned_samples = dict(tuple(prealigned_sequences.groupby('input_sample_name')))
        samples = {}
        for sample_name, prealigned_sequences in prealigned_samples.items():
            sequence_df = sequences_to_df(
                prealigned_sequences['aligned_sequence'].tolist(),
                center=center, redundancy_level=redundancy_level
            )
            sequence_tensor = vectorize_sequences(sequence_df, background)
            samples[str(sample_name)] = Sample(
                sequence_df=sequence_df,
                sequence_tensor=sequence_tensor,
                original_sequences=prealigned_sequences
            )
    
    return samples

def merge_sequences(matches):
    """
    If multiple possible sequence extensions are found in proteome,
        combine matches to single sequence with ambiguous positions
        coded as 'X'.

    Parameters:
        matches -- Set.

    Returns:
        merged_sequence -- String.
    """

    # Convert matches to data frame for disambiguation.
    matches = pd.DataFrame([list(match) for match in sorted(list(matches))])
    # Merge unique matches with 'X' in ambiguous positions.
    merged_sequence = ''.join((matches.loc[0, i] if unambiguous else 'X'
                    for i, unambiguous
                    in enumerate(matches.nunique(axis=0) == 1)))

    return merged_sequence

def extend_sequence(
        match, 
        context_element, 
        width, 
        non_prime_width, 
        terminal
    ):
    """
    Extend sequence in X-terminal direction.

    Parameters:
        match -- Match object.
        width -- Int.
        non_prime_width -- Int.
        terminal -- String.

    Returns:
        extended_sequence -- String.
    """
    
    if terminal == 'n':
        match_index = match.start(0)
    elif terminal == 'c':
        match_index = match.end(0)

    start = (
        match_index - non_prime_width if match_index >= non_prime_width else 0
    )
    end = match_index + width - non_prime_width
    non_prime = context_element[start:match_index].rjust(non_prime_width, '-')
    prime = context_element[match_index:end].ljust(width - non_prime_width, '-')
    extended_sequence = non_prime + prime

    return extended_sequence

def match_segment_to_context(
        segment,
        context_element,
        width,
        non_prime_width,
        terminal
    ):
    """
    Regular expression match of segment string to context element
        string. Return matches.

    Parameters:
        segment -- String.
        context_element -- String.
        width -- Int.
        non_prime_width -- Int.
        terminal -- String. 'n', 'c', or both.

    Returns:
        matches -- Set.
    """

    matches = ExtendedSequences(n_term=set(), c_term=set())

    # Set of unique non-prime matches in context sequence.
    for match in re.finditer(segment, context_element):
        if terminal in ['n', 'N']:
            matches.n_term.add(
                extend_sequence(match, context_element, width, non_prime_width, 'n')
            )
        elif terminal in ['c', 'C']:
            matches.c_term.add(
                extend_sequence(match, context_element, width, non_prime_width, 'c')
            )
        elif terminal == 'both':
            matches.n_term.add(
                extend_sequence(match, context_element, width, non_prime_width, 'n')
            )
            matches.c_term.add(
                extend_sequence(match, context_element, width, non_prime_width, 'c')
            )

    return matches

def get_all_ids_from_context(context, precomputed):
    if precomputed is not None:
        if os.path.basename(precomputed) in PRECOMPUTED_FILES['swissprot_human']:
            context_ids = context['swissprot_id'].tolist()
        else:
            context_ids = context['id'].tolist()
    else:
        context_ids = context['id'].tolist()

    return context_ids

def check_context_id_types(context_id, context_id_types, context):
    for context_id_type in context_id_types:
        context_sequences = context[
            context[context_id_type] == context_id
        ]['sequence'].tolist()
        num_sequences = len(context_sequences)
        if num_sequences > 0:
            break
            
    return context_sequences, num_sequences

def get_context_sequences(context_id, context):
    parsed_id = utils.extract_accid(context_id)
    if parsed_id:
        context_id_types = ['swissprot_id', 'id']
        context_sequences, num_sequences = check_context_id_types(
            parsed_id,
            context_id_types,
            context
        )    
    else:
        context_id_types = ['id']
        context_sequences, num_sequences = check_context_id_types(
            context_id,
            context_id_types,
            context
        )

    return context_sequences, num_sequences


def align_sequences(
        context,
        sequences,
        width=15,
        terminal='N',
        redundancy_level='none',
        first_protein_only=True,
        original_row_merge='all',
        original_sequences=None,
        # TODO: Check this now False
        require_context_id=False,
        precomputed=None
    ):
    """
    Take unaligned sequences and context. Map unaligned sequences to
        context. Truncated prime segment to half width. Extend sequence
        by half width in non-prime direction based on context. Replace
        ambiguous residues in non-prime segment with "X". Fill missing
        positions with "-".

    Parameters:
        context -- Pandas DataFrame. Expected columns: [
                        'id' (ID from supplied FASTA file),
                        'sequence' (sequence strings),
                    ]
        sequences -- Pandas DataFrame. Expected columns: [
                        'context_id' (IDs matching context
                                        elements, sep=';'),
                        'sequence' (sequence strings),
                    ]
        width -- Int. Number of positions in aligned sequence output.
        terminal -- String. 'n' or 'N' for extension of N-term sequences,
                        'c' or 'C' for extension of C-term sequences,
                        'both' for extension of intact peptides.

    Returns:
        aligned_sequences -- List. Aligned sequences of specified width
                                containing prime segment from sequences
                                and non-prime segment from matching
                                context.
    """

    non_prime_width = width // 2

    # Remove exact peptides and protein ID duplicates unless disabled.
    if redundancy_level != 'none':
        sequences.drop_duplicates(inplace=True) 
    
    # print('Aligning {} sequences to {} context elements'.format(len(sequences), len(context)))

    # Generate aligned sequence by mapping sequence segments to context.
    aligned_sequences = []
    for i, sequence in sequences.iterrows():
        # Select pre-mapped context elements if available.
        try:
            context_ids = sequence['context_id'].replace(' ', '').split(';')
        except:
            # Exclude sequence if no IDs are required but not provided.
            if require_context_id:
                continue
            # Otherwise, try matching against all context sequences.
            else:
                context_ids = get_all_ids_from_context(context, precomputed)
                context_elements = context['sequence'].tolist()
                first_protein_only = False
        else:            
            # Match against all context sequences if ID required and
            # IDs column is present but blank.
            if (len(max(context_ids, key=len)) == 0) and not require_context_id:
                context_ids = get_all_ids_from_context(context, precomputed)
                context_elements = context['sequence'].tolist()
                first_protein_only = False
            # Otherwise exclude sequence.
            elif len(max(context_ids, key=len)) == 0:
                continue
            # Match against provided IDs.
            else:
                swissprot_ids = []
                context_elements = []
                # Get context sequences using provided ID.
                for context_id in context_ids:
                    # Special case for precoumpted human Swiss-Prot.
                    if precomputed is not None:
                        if os.path.basename(precomputed) in PRECOMPUTED_FILES['swissprot_human']:
                            parsed_id = utils.extract_accid(context_id)
                            if SWISSPROT_ACCESSION_PATTERN.match(parsed_id):
                                swissprot_sequences = context[
                                    context['swissprot_id'] == parsed_id
                                ]['sequence'].tolist()
                                num_sequences = len(swissprot_sequences)
                                if  num_sequences == 1:
                                    swissprot_ids.append(parsed_id)
                                    context_elements += swissprot_sequences
                                    if first_protein_only:
                                        break
                                elif num_sequences > 1:
                                    raise AssertionError(
                                        '{} was found more than one time in the proteome.'.format(
                                            parsed_id
                                        )
                                    )
                    # If not precomputed, use contexts
                    else:
                        # Get context sequences using provided ID.
                        context_sequences, num_sequences = get_context_sequences(context_id, context)
                        if num_sequences == 1:
                            context_elements += context_sequences
                            if first_protein_only:
                                break
                        elif num_sequences > 1:
                            raise AssertionError(
                                '{} was found more than one time in the proteome.'.format(
                                    context_id
                                )
                            )
                if len(swissprot_ids) > 0:
                    context_ids = swissprot_ids
        
        # print('Sequence: ', sequence['sequence'])
        # print('Aligning sequence {} to {} context elements'.format(i, len(context_elements)))
        # print("Context IDs: ", context_ids)
        # print("Context Elements: ", context_elements) # context_elements are protein sequences
        # print("Width: ", width, "Non-prime width: ", non_prime_width, "Terminal: ", terminal)

        # Regular expression matching of sequence to context.
        extended_sequences = {}
        for j, context_element in enumerate(context_elements):
            matches = match_segment_to_context(
                segment = sequence['sequence'].rstrip().upper(),
                context_element = context_element,
                width = width,
                non_prime_width = non_prime_width,
                terminal = terminal
            )
            if matches.n_term or matches.c_term:
                context_id = context_ids[j]
                extended_sequences[context_id] = ExtendedSequences(
                    n_term=set(),
                    c_term=set()
                )
                for n_term_match in sorted(list(matches.n_term)):
                    extended_sequences[context_id].n_term.add(n_term_match)
                for c_term_match in sorted(list(matches.c_term)):
                    extended_sequences[context_id].c_term.add(c_term_match)
        
        # Add all aligned instances of peptide if redundancy_level is 'none'.
        if original_row_merge == 'none':
            for context_id, context_element_matches in extended_sequences.items():
                for extended_sequence in sorted(list(context_element_matches.n_term)):
                    aligned_sequences.append([i, extended_sequence, context_id])
                for extended_sequence in sorted(list(context_element_matches.c_term)):
                    aligned_sequences.append([i, extended_sequence, context_id])
        # Merge aligned instances of peptide at protein level otherwise.
        elif original_row_merge == 'protein':
            for context_id, context_element_matches in extended_sequences.items():
                if context_element_matches.n_term:
                    aligned_sequences.append(
                        [i, merge_sequences(context_element_matches.n_term), context_id]
                    )
                if context_element_matches.c_term:
                    aligned_sequences.append(
                        [i, merge_sequences(context_element_matches.c_term), context_id]
                    )
        elif original_row_merge == 'all':
            context_ids = []
            all_extended_sequences = ExtendedSequences(
                n_term=set(),
                c_term=set()
            )
            for context_id, context_element_matches in extended_sequences.items():
                for extended_sequence in context_element_matches.n_term:
                    all_extended_sequences.n_term.add(extended_sequence)
                for extended_sequence in context_element_matches.c_term:
                    all_extended_sequences.c_term.add(extended_sequence)
                context_ids.append(context_id)
            context_id = '; '.join(context_ids)
            if all_extended_sequences.n_term:
                aligned_sequences.append(
                    [i, merge_sequences(all_extended_sequences.n_term), context_id]
                )
            if all_extended_sequences.c_term:
                aligned_sequences.append(
                    [i, merge_sequences(all_extended_sequences.c_term), context_id]
                )

    # Generate aligned sequence data frame for redundancy processing.
    aligned_sequences_df = pd.DataFrame(
        aligned_sequences,
        columns=['original_row', 'extended_sequence', 'context_id']
    )
    # Redundancy handling.
    if redundancy_level == 'protein':
        aligned_sequences_df.drop_duplicates(
            subset=['extended_sequence', 'context_id'],
            inplace=True
        )
    elif redundancy_level == 'sequence':
        aligned_sequences_df.drop_duplicates(
            subset=['extended_sequence'],
            inplace=True
        )

    # Update tabular output DataFrame.
    if original_sequences is not None:
        aligned_sequence_list = []
        matched_context_id_list = []
        for i, original_row in original_sequences.iterrows():
            aligned_row_results = aligned_sequences_df[
                aligned_sequences_df['original_row'] == i
            ]
            aligned_sequence = '; '.join(aligned_row_results['extended_sequence'].tolist())
            aligned_sequence_list.append(
                aligned_sequence if aligned_sequence != '' else np.nan
            )
            matched_context_id = '; '.join(aligned_row_results['context_id'].tolist())
            matched_context_id_list.append(
                matched_context_id if matched_context_id != '' else np.nan

            )
        original_sequences['aligned_sequence'] = aligned_sequence_list
        original_sequences['matched_context_id'] = matched_context_id_list

    aligned_sequences = aligned_sequences_df['extended_sequence'].tolist()

    return aligned_sequences

def peptides_to_sample(
        peptides,
        context,
        background,
        center=True,
        width=15,
        terminal='N',
        redundancy_level='none',
        first_protein_only=True,
        original_row_merge='all',
        require_context_id=True
    ):
    """
    Align peptides and return data frame of positional residues.

    Parameters:
        peptides --
        context --
        width --
        terminal --

    Returns:
        sample -- Sample instance.
    """

    original_sequences = peptides.copy()
    if len(original_sequences.columns) == 1:
        original_sequences.columns = [
            'input_sequence',
        ]
    elif len(original_sequences.columns) == 2:
        original_sequences.columns = [
            'input_sequence',
            'input_context_id',
        ]
    else:
        original_sequences.rename(
            columns={
                original_sequences.columns[0]: 'input_sample_name',
                original_sequences.columns[1]: 'input_sequence',
                original_sequences.columns[2]: 'input_context_id',
            },
            inplace=True
        )

    aligned_sequences = align_sequences(
        context,
        peptides,
        width=width,
        terminal=terminal,
        redundancy_level=redundancy_level,
        first_protein_only=first_protein_only,
        original_row_merge=original_row_merge,
        original_sequences=original_sequences,
        require_context_id=require_context_id,
        precomputed=background.precomputed
    )
    
    sequence_df = sequences_to_df(
        aligned_sequences,
        center=center,
        redundancy_level=redundancy_level
    )
    sequence_tensor = vectorize_sequences(sequence_df, background)
    sample = Sample(
        sequence_df=sequence_df,
        sequence_tensor=sequence_tensor,
        original_sequences=original_sequences
    )

    return sample

def import_peptide_list(
        peptide_list_file,
        delimiter='\t',
    ):
    """
    Import peptide list and row-matched protein IDs from text file.

    Parameters:
        peptide_list_file -- File-like object. Path to peptide list.
        delimiter -- String. Peptide list file delimiter.
        require_context_id -- Boolean. Discard rows without context ID.

    Returns:
        peptide_list -- Pandas DataFrame.
    """

    # Import peptide list from file.
    peptide_list = pd.read_csv(peptide_list_file, delimiter=delimiter, header=0)

    # Assign proper headers based on assumed data format.
    if len(peptide_list.columns) == 1:
        peptide_list.columns = ['sequence']
    elif len(peptide_list.columns) == 2:
        peptide_list.columns = ['sequence', 'context_id']
    else:
        peptide_list.rename(
            columns={
                peptide_list.columns[0]: 'sample_name',
                peptide_list.columns[1]: 'sequence',
                peptide_list.columns[2]: 'context_id'
            },
            inplace=True
        )

    return peptide_list

def load_peptide_list_file(
        peptide_list_path,
        context,
        background,
        center=True,
        width=15,
        terminal='N',
        require_context_id=False,
        redundancy_level='none',
        first_protein_only=True,
        original_row_merge='all',
        title='',
        verbose=False
    ):
    """
    Top-level helper function to load and extend peptides from text
        file.

    Parameters:
        peptide_list_path -- String. Path to peptide list.
        context -- Pandas DataFrame. Context proteome data frame.
        terminal -- String. Defines extension direction.
        require_context_id -- Boolean.
        redundancy_level -- String. Defines level of repeated sequence
                                redundancy elimination.
        first_protein_only -- Boolean. Only use the first context ID
                                for each row in the foreground data set
                                during alignment and extension.
        original_row_merge -- String. Defines behavior for merging
                                multiple proteome matches from each
                                row in the foreground data set during 
                                alignment and extension.
        title -- String. Analysis title.

    Returns:
        sample -- Sample instance.
    """

    # Detect delimiter from file extension (supports comma or tab).
    delimiter = detect_delimiter(peptide_list_path)
    peptide_list = import_peptide_list( peptide_list_path, delimiter=delimiter )
    
    try:
        sample_peptides = dict(tuple(peptide_list.groupby('sample_name')))
    except:
        sample_peptides = {title: peptide_list}

    samples = {}
    for sample_name, peptides in sample_peptides.items():
        if verbose:
            print('\t - Processing sample {} with {} Peptides'.format(sample_name, len(peptides)))
        samples[str(sample_name)] = peptides_to_sample(
            peptides,
            context,
            background,
            center=center,
            width=width,
            terminal=terminal,
            redundancy_level=redundancy_level,
            first_protein_only=first_protein_only,
            original_row_merge=original_row_merge,
            require_context_id=require_context_id
        )

    return samples

class Background:
    """
    Methods and data structures for background frequency data.
    
    Attributes:
        sequences (List[str]): List of background sequences.
        background_dict (Dict[str, float]): Percentage frequencies mapped to single letter codes.
        background_vector (np.ndarray): 1D Numpy array of percentage frequencies.
    """

    def __init__(
            self,
            sequences: Optional[List[str]] = None,
            background_dict: Optional[Dict[str, float]] = None,
            compound_residues: Optional[Dict[str, CompoundResidue]] = COMPOUND_RESIDUES,
            position_specific: bool = False,
            width: int = 15,
            regular_expression: Optional[str] = None,
            initial_background_size_limit: Optional[int] = None,
            fast: bool = False,
            center: bool = False,
            precomputed: str = None
        ):
        """
        Default background mode. Initialize from list of sequences.

        Args:
            sequences (List[str]): Background sequence strings.
            background_dict (Dict[str, float]): Background frequency dictionary.
            compound_residues (Dict[str, CompoundResidue]): Compound residues.
            position_specific (bool): Whether to use position-specific background.
            width (int): Width of the sequences.
            regular_expression (str): Regular expression for sequence matching.
            initial_background_size_limit (int): Initial background size limit.
            fast (bool): Whether to use fast mode.
            center (bool): Whether to center the sequences.
            precomputed (str): Path to precomputed background data.

        Returns:
            None
        """
        self.precomputed = precomputed
        self._position_specific = position_specific
        self._fast = fast

        self._background_sequences = sequences
        if sequences:
            self._background_dict = self._calculate_percentage_frequencies()
        elif background_dict:
            self._background_dict = background_dict
        else:
            raise ValueError('No background provided.')

        self._alphabet = sorted(self._background_dict.keys())
        self._compound_residues = compound_residues
        if compound_residues:
            self._background_dict.update(self._generate_compound_residues(compound_residues))

        self._background_vector = self._vectorize_background()

        if precomputed:
            self._background_df = pd.read_csv(precomputed, sep=',', header=0)
        elif position_specific:
            if fast:
                self._generate_background_tensor_efficient(sequences, width, regular_expression, initial_background_size_limit)
            else:
                self._generate_background_df(sequences, width, regular_expression, initial_background_size_limit, center)

    def _generate_background_df(
            self,
            sequences: List[str],
            width: int,
            regular_expression: Optional[str],
            initial_background_size_limit: Optional[int],
            center: bool
        ) -> None:
        """
        Generates Pandas DataFrame containing aligned background sequences from reference data set.

        Args:
            sequences (List[str]): List of sequences.
            width (int): Width of the sequences.
            regular_expression (str): Regular expression for sequence matching.
            initial_background_size_limit (int): Initial background size limit.
            center (bool): Whether to center the sequences.

        Returns:
            None
        """
        self._background_df = sequences_to_df(
            self.get_all_centered(sequences, width, regular_expression),
            center=center
        )

    def _generate_background_tensor(
            self,
            sequences: List[str],
            width: int,
            regular_expression: Optional[str],
            initial_background_size_limit: Optional[int]
        ) -> None:
        """
        Generates background tensor as 3D NumPy array.

        Args:
            sequences (List[str]): List of sequences.
            width (int): Width of the sequences.
            regular_expression (str): Regular expression for sequence matching.
            initial_background_size_limit (int): Initial background size limit.

        Returns:
            None
        """
        if not initial_background_size_limit:
            initial_background_size_limit = len(sequences)

        self._background_tensor = vectorize_sequences(
            sequences_to_df(
                self.get_all_centered(sequences, width, regular_expression)
            ),
            self
        )

    def _generate_background_tensor_efficient(
            self,
            sequences: List[str],
            width: int,
            regular_expression: Optional[str],
            initial_background_size_limit: Optional[int]
        ) -> None:
        """
        NumPy-based background tensor generation solution.

        Args:
            sequences (List[str]): List of sequences.
            width (int): Width of the sequences.
            regular_expression (str): Regular expression for sequence matching.
            initial_background_size_limit (int): Initial background size limit.

        Returns:
            None
        """
        if not initial_background_size_limit:
            initial_background_size_limit = len(sequences)

        residue_dict = {residue: i for i, residue in enumerate(self._ordered_residues)}
        self._background_tensor = sequences_to_df(
            self.get_all_centered(sequences, width, regular_expression)
        ).applymap(lambda x: residue_dict[x]).to_numpy(dtype=np.int8)

        if self._compound_residues:
            self._background_tensor[:, :, len(self._alphabet):] = np.inner(
                self._background_tensor[:, :, :len(self._alphabet)],
                self.compound_residue_matrix
            )

    @classmethod
    def from_csv(cls, background_path: str) -> 'Background':
        """
        Alternative background mode. Initialize from list of background frequencies.

        Args:
            background_path (str): Path to background frequency CSV file.

        Returns:
            Background: Background instance.
        """
        background_dict = {}
        with open(background_path) as background_file:
            reader = csv.reader(background_file, delimiter=',')
            for row in reader:
                background_dict[row[0]] = float(row[1])
        
        return cls(None, background_dict=background_dict)

    @classmethod
    def from_json(cls, background_path: str) -> 'Background':
        """
        Generate background frequencies dictionary from JSON file containing percentage frequencies mapped to single letter codes.

        Args:
            background_path (str): Path to background frequency JSON file.

        Returns:
            Background: Background instance.
        """
        with open(background_path) as background_file:
            background_dict = json.load(background_file)

        return cls(None, background_dict=background_dict)

    @property
    def fast(self) -> bool:
        return self._fast
    
    @property
    def background_sequences(self) -> Optional[List[str]]:
        return self._background_sequences

    @property
    def compound_residues(self) -> Optional[Dict[str, CompoundResidue]]:
        return self._compound_residues

    @property
    def compound_residue_codes(self) -> List[str]:
        return self._compound_residue_codes

    @property
    def compound_residue_frequencies(self) -> Dict[str, float]:
        return self._compound_residue_frequencies

    @property
    def compound_residue_matrix(self) -> pd.DataFrame:
        return self._compound_residue_matrix

    @property
    def ordered_residues(self) -> List[str]:
        return self._ordered_residues
    
    @background_sequences.setter
    def background_sequences(self, sequences: List[str]) -> None:
        self._background_sequences = sequences
        self._background_dict = self._calculate_percentage_frequencies()
        self._background_vector = self._vectorize_background()

    @property
    def background_dict(self) -> Dict[str, float]:
        return self._background_dict

    @property
    def background_vector(self) -> np.ndarray:
        return self._background_vector

    @property
    def alphabet(self) -> List[str]:
        return self._alphabet

    @property
    def background_tensor(self) -> np.ndarray:
        return self._background_tensor

    @property
    def background_df(self) -> pd.DataFrame:
        return self._background_df

    @property
    def position_specific(self) -> bool:
        return self._position_specific
    
    def _calculate_percentage_frequencies(self, remove_ambiguous_residues: bool = True) -> Dict[str, float]:
        """
        Calculate background frequencies from list of sequences.

        Args:
            remove_ambiguous_residues (bool): Whether to remove ambiguous residues.

        Returns:
            Dict[str, float]: Letter code frequencies.
        """
        absolute_frequencies = Counter()
        for sequence in self.background_sequences:
            absolute_frequencies += Counter(sequence)

        if remove_ambiguous_residues:
            absolute_frequencies.pop('X', None)

        total = sum(absolute_frequencies.values())
        percentage_frequencies = {residue: (absolute_frequency / total) * 100
                                  for residue, absolute_frequency in absolute_frequencies.items()}

        return percentage_frequencies

    def _generate_compound_residues(self, compound_residues: Dict[str, CompoundResidue]) -> Dict[str, float]:
        """
        Load compound residues from dictionary of CompoundResidue objects. Calculate expected compound residue frequencies
        as sum of expected constituent single residue frequencies. Return dictionary of expected frequencies, which is used
        to update background frequency dictionary.

        Args:
            compound_residues (Dict[str, CompoundResidue]): Dictionary of CompoundResidue objects.

        Returns:
            Dict[str, float]: Expected compound residue frequencies.
        """
        self._compound_residue_codes = sorted(compound_residues.keys())
        self._compound_residue_frequencies = {}
        compound_residue_vectors = []

        for compound_residue in self._compound_residue_codes:
            compound_residue_vector = pd.Series(np.zeros(len(self.alphabet), dtype='int8'), index=self.alphabet)
            compound_residue_frequency = 0

            for residue in compound_residues[compound_residue].residues:
                compound_residue_vector[residue] = 1
                compound_residue_frequency += self._background_dict[residue]

            self._compound_residue_frequencies[compound_residue] = compound_residue_frequency
            compound_residue_vectors.append(compound_residue_vector)

        self._compound_residue_matrix = pd.DataFrame(compound_residue_vectors, index=self._compound_residue_codes)

        return self._compound_residue_frequencies

    def _vectorize_background(self) -> np.ndarray:
        """
        Convert background frequencies dictionary to Numpy array.

        Returns:
            np.ndarray: 1D Numpy array of percentage frequencies.
        """
        try:
            self._ordered_residues = self.alphabet + self._compound_residue_codes
        except AttributeError:
            self._ordered_residues = self.alphabet

        background_frequency_vector = np.array(
            [self._background_dict[residue] for residue in self._ordered_residues],
            dtype=np.float64
        )

        return background_frequency_vector

    def get_all_centered(
            self,
            sequences: List[str],
            width: int,
            regular_expression: Optional[str] = None
        ) -> List[str]:
        """
        Get all residue-centered sequences of specified width from context data set.

        Args:
            sequences (List[str]): List of sequences.
            width (int): Width of the sequences.
            regular_expression (str): Regular expression for sequence matching.

        Returns:
            List[str]: List of centered sequences.
        """
        if not regular_expression:
            regular_expression = r'(?=(.{{{0}}}))'.format(width)
        elif len(regular_expression) == 1 and regular_expression.isalpha():
            flanking_width = width // 2
            residue_cases = regular_expression.lower() + regular_expression.upper()
            regular_expression = r'(?=(.{{{0}}}[{1}].{{{0}}}))'.format(flanking_width, residue_cases)

        centered_sequences = []
        for sequence in sequences:
            centered_sequences += re.findall(regular_expression, sequence)

        return centered_sequences

    def to_csv(self, output_path: str) -> None:
        """Save background_df to CSV."""
        self._background_df.to_csv(output_path, header=True, index=False)