import re
import csv
import json
import numpy as np
import pandas as pd
from collections import Counter
from typing import List, Dict, Optional, Union

from .globals import COMPOUND_RESIDUES, CompoundResidue, Sample

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

def sequences_to_df(sequences,
                    center=True,
                    redundancy_level='sequence'):
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
    index_array = data.applymap(
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

def load_prealigned_file(prealigned_file_path,
                            background,
                            center=True,
                            redundancy_level='none',
                            title=''):
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
            precomputed: str = ''
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