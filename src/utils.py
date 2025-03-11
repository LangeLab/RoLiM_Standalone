import re
import time
from tqdm import tqdm

import numpy as np
import pandas as pd

from .globals import (
    SINGLE_LETTER_CODES, SWISSPROT_ACCESSION_PATTERN,
    COMPOUND_RESIDUES, CompoundResidue,
    encode_3to1, encode_1to3
)

def getTime() -> float:
    """
    Get the current time for timer

    Returns:
        float: The current time in seconds.
    """
    return time.time()

def prettyTimer(
        seconds: float
    ) -> str:
    """
    Better way to show elapsed time

    Args:
        seconds (float): The number of seconds to convert to a pretty format.

    Returns:
        str: The elapsed time in a pretty format.
    
    Examples:
        >>> prettyTimer(100)
        '00h:01m:40s'

        >>> prettyTimer(1000)
        '00h:16m:40s'

        >>> prettyTimer(10000)
        '02h:46m:40s'
    """
    m, s = divmod(seconds, 60)
    h, m = divmod(m, 60)
    return "%02dh:%02dm:%02ds" % (h, m, s)

def extract_accid(
        identifier: str
    ) -> str:
    """
    Extracts the accession number from a given identifier string 
    using a global regex pattern for Swiss-Prot (UniProt) accession numbers.

    Args:
        identifier (str): The identifier string to parse.

    Returns:
        str: The parsed Swiss-Prot accession number, or None if not found.

    Examples:
        >>> extract_accid("P12345")
        'P12345'

        >>> extract_accid(">Some text P12345-1 more text")
        'P12345-1'

        >>> extract_accid("No accession number here")
        None
    """

    try:
        swissprot_accession_number = [
            identifier[match.start():match.end()]
            for match in re.finditer( SWISSPROT_ACCESSION_PATTERN, identifier )
        ][0]
    except IndexError:
        swissprot_accession_number = None

    return swissprot_accession_number

def parse_fasta(
        fasta_path: str,
        swissprot: bool = True
    ) -> pd.DataFrame:
    """
    Parses a fasta file and returns a pandas DataFrame containing the sequences
    and their IDs. If swissprot is set to True, the Swiss-Prot accession number
    is also extracted and added to the DataFrame.
    
    Args:
        fasta_path (str): Path to fasta file.
        swissprot (bool): If True, parse Swiss-Prot accession numbers. Default is True.

    Returns:
        pd.DataFrame: DataFrame containing the fasta sequences with columns 
        'id', 'sequence', and 'swissprot_id'.
    """

    with open(fasta_path) as fasta:
        sequences = []
        for line in fasta:
            # Skip leading comments (not common any more).
            if line[0] == ';':
                continue
            elif line[0] == '>':
                # Add sequence and ID to sequences if not first ID.
                try:
                    sequences.append([sequence_id, sequence, swissprot_id])
                except NameError:
                    pass
                finally:
                    sequence_id = line.rstrip()[1:]
                    sequence = ''
                    swissprot_id = extract_accid(sequence_id) if swissprot else None
            else:
                # Add characters to sequence if sequence has ID.
                try:
                    sequence += line.rstrip().upper()
                except NameError:
                    continue
        sequences.append([sequence_id, sequence, swissprot_id])

    cols = ['id', 'sequence', 'swissprot_id']
    fasta_df = pd.DataFrame(sequences, columns=cols)

    return fasta_df

def load_compound_residues(
        compound_residues_path: str,
    ) -> dict:
    """    
    Parses the compound residues file and returns a dictionary containing the
    compound residue specification to be used for grouping amino acids. If no
    path is provided, the default compound residues variable are used.

    Args:
        compound_residues_path (str): Path to the compound residues file.

    Returns:
        dict: Dictionary containing the compound residue specification.
    """
    
    if compound_residues_path is None:
        return COMPOUND_RESIDUES
    
    compound_residues = {}
    with open(compound_residues_path, 'r') as f:
        for i, line in enumerate(f):
            description, residues = line.strip().split('\t')
            compound_residues[str(i + 1)] = CompoundResidue(
                description=description,
                residues=[residue.strip() for residue in residues.split(';')]
            )

    return compound_residues

def record_compound_residues_to_file(
        compound_residues: dict,
        current_path: str
    ) -> None:
    """
    Writes the compound residues dictionary to a file for record keeping.
    """
    # Define the filepath
    file_path = current_path + 'compound_residues_definition.txt'
    # Open the file and write the data
    with open(file_path, 'w') as f:
        # Write header
        f.write('Compound Residue Code\tDescription\tConstituents\n')
        # Write data rows
        for code, specifications in compound_residues.items():
            f.write(
                f"{code}\t{specifications.description}\t{';'.join(specifications.residues)}\n"
            )

def convert_encoding(
        sequences: pd.DataFrame,
        encode: int
    ) -> pd.DataFrame:
    """
    Converts amino acid from single-letter to three-letter encoding and vice versa.
    Takes and returns a Pandas DataFrame.

    Args:
        sequences (pd.DataFrame): Contains sequences in three-letter code.
        encode (int): Either 1 or 3 based on desired output encoding.

    Returns:
        pd.DataFrame: Contains sequences from input data converted from three-letter
                    to single-letter encoding or vice versa.

    Examples:
        >>> df = pd.DataFrame([['Ala', 'Arg'], ['Asn', 'Asp']])
        >>> convert_encoding(df, 1)
        0  1
        0  A  R
        1  N  D

        >>> df = pd.DataFrame([['A', 'R'], ['N', 'D']])
        >>> convert_encoding(df, 3)
            0    1
        0  Ala  Arg
        1  Asn  Asp
    """

    # TODO: Consider future optimization by moving out of pandas applymap
    # TODO: Checks if the input for 3 or 1 is valid
    # Uses pre-defined dictionaries from globals.py
    if encode == 3:
        converted = sequences.applymap(lambda x: encode_1to3[x])
    elif encode == 1:
        converted = sequences.applymap(lambda x: encode_3to1[x.capitalize()])
    else:
        print('Please enter a valid output encoding format')

    return converted

def generate_positions(
        center: bool,
        num_positions: int
    ) -> list:
    """
    Generates position labels for sequences.

    Args:
        center (bool): If true, column labels will count away from the center 
                       (leftward direction = non-prime, rightward direction = prime).
                       If false, column labels will count from left to right beginning with 1.
        num_positions (int): Total number of position labels to generate.

    Returns:
        list: Contains position labels.

    Examples:
        >>> generate_positions(True, 8)
        ['p4', 'p3', 'p2', 'p1', "p0'", "p1'", "p2'", "p3'"]

        >>> generate_positions(False, 8)
        ['p1', 'p2', 'p3', 'p4', 'p5', 'p6', 'p7', 'p8']
    """

    if center:
        positions = []
        p1_pos = int(num_positions / 2)

        for i in range(num_positions):
            if i < p1_pos:
                position_label = 'p' + str(p1_pos - i)
                positions.append(position_label)
            else:
                position_label = 'p' + str(abs(p1_pos - (i + 1))) + "'"
                positions.append(position_label)
    else:
        positions = ['p{}'.format(i) for i in range(1, num_positions + 1)]

    return positions

def load_sequences(
        data: str,
        center: bool
    ) -> pd.DataFrame:
    """
    Generate and return a pandas dataframe containing input sequences.
    Additional columns include sequence ID, cluster ID, alignmentOffset 
    (+ shifts sequence left, - shifts sequence right), and peptide score.

    Args:
        data (str): Path to input data set in text file form.
        center (bool): If true, column labels will count away from the center 
                       (leftward direction = non-prime, rightward direction = prime).
                       If false, column labels will count from left to right beginning with 1.

    Returns:
        pd.DataFrame: Contains input data set split to columns by position and 
                      accompanying fields as specified above.

    Examples:
        >>> load_sequences('path/to/data.txt', True)
        # Returns a DataFrame with sequences split into columns by position
    """

    split_sequences = []

    with open(data) as sequences:
        for sequence in sequences:
            split_sequences.append(list(sequence.rstrip().upper()))

    cols = generate_positions(center, len(split_sequences[0]))

    sequence_df = pd.DataFrame(split_sequences, columns=cols)

    return sequence_df

def filter_missing_residues(
        data: np.ndarray,
        condition: str
    ) -> np.ndarray:
    """
    Filters vectorized sequences based on missing residues.

    Args:
        data (np.ndarray): 3D Numpy Array.
        condition (str): If 'any', remove sequences containing any missing residues.
                         If 'all', remove sequences missing all residues.

    Returns:
        np.ndarray: Contains vectorized sequences filtered on missing residues and 
                    user-specified filtering condition.

    Examples:
        >>> data = np.array([[[1, 2], [0, 0]], [[3, 4], [5, 6]]])
        >>> filter_missing_residues(data, 'any')
        array([[[3, 4],
                [5, 6]]])

        >>> filter_missing_residues(data, 'all')
        array([[[1, 2],
                [0, 0]],
               [[3, 4],
                [5, 6]]])
    """

    if condition == 'any':
        filtered_data = data[tuple(np.where(np.any(np.all(
            data == 0, axis=2), axis=1) == False))]
    elif condition == 'all':
        filtered_data = data[tuple(np.where(np.all(np.all(
            data == 0, axis=2), axis=1) == False))]

    return filtered_data

def filter_unknown_residues(
        data: pd.DataFrame,
        condition: str
    ) -> pd.DataFrame:
    """
    Filters character sequences based on unknown residue codes using
    a list of expected residue codes from the background data set.

    Args:
        data (pd.DataFrame): Contains input data set split to columns by position 
                             and accompanying fields as specified above.
        condition (str): If 'any', remove sequences containing any missing residues.
                         If 'all', remove sequences missing all residues.

    Returns:
        pd.DataFrame: Contains input data filtered based on expected residue codes 
                      and user-specified filtering condition.

    Examples:
        >>> df = pd.DataFrame([['A', 'R'], ['N', 'X'], ['D', 'E']])
        >>> filter_unknown_residues(df, 'any')
           0  1
        0  A  R
        2  D  E

        >>> filter_unknown_residues(df, 'all')
           0  1
        0  A  R
        1  N  X
        2  D  E
    """

    if condition == 'any':
        filtered_data = data.iloc[
            data[data.isin(SINGLE_LETTER_CODES)].dropna(how='any').index
        ]
    elif condition == 'all':
        filtered_data = data.iloc[
            data[data.isin(SINGLE_LETTER_CODES)].dropna(how='all').index
        ]

    return filtered_data

def import_substitution_matrix(
        substitution_matrix_file: str
    ) -> pd.DataFrame:
    """
    Import substitution matrix to pandas dataframe. Takes a text file as
    an argument and returns a pandas dataframe.

    Args:
        substitution_matrix_file (str): Path to text file containing amino acid 
                                        substitution probabilities.

    Returns:
        pd.DataFrame: Symmetrical matrix containing pairwise amino acid 
                      substitution probabilities.

    Examples:
        >>> df = import_substitution_matrix('path/to/substitution_matrix.txt')
        >>> df.head()
             A    R    N    D    C    Q    E    G    H    I  ...    P    S    T    W    Y    V    B    Z    X    *
        A  4.0 -1.0 -2.0 -2.0  0.0 -1.0 -1.0  0.0 -2.0 -1.0  ...  -1.0  1.0  0.0 -3.0 -2.0  0.0 -2.0 -1.0  0.0 -4.0
        R -1.0  5.0  0.0 -2.0 -3.0  1.0  0.0 -2.0  0.0 -3.0  ...  -2.0 -1.0 -1.0 -3.0 -2.0 -3.0 -1.0  0.0 -1.0 -4.0
        N -2.0  0.0  6.0  1.0 -3.0  0.0  0.0  0.0  1.0 -3.0  ...  -2.0  1.0  0.0 -4.0 -2.0 -3.0  3.0  0.0 -1.0 -4.0
        D -2.0 -2.0  1.0  6.0 -3.0  0.0  2.0 -1.0 -1.0 -3.0  ...  -1.0  0.0 -1.0 -4.0 -3.0 -3.0  4.0  1.0 -1.0 -4.0
        C  0.0 -3.0 -3.0 -3.0  9.0 -3.0 -4.0 -3.0 -3.0 -1.0  ...  -3.0 -1.0 -1.0 -2.0 -2.0 -1.0 -3.0 -3.0 -2.0 -4.0
    """

    with open(substitution_matrix_file) as submat:
        matrix = []
        for i, line in enumerate(submat):
            if i == 0:
                headers = line.split()
                continue

            line_values = line.split()
            line_values.pop(0)

            for j, value in enumerate(line_values):
                line_values[j] = float(value)

            matrix.append(line_values)
        
    substitution_matrix_df = pd.DataFrame(matrix, index=headers, columns=headers)

    return substitution_matrix_df
