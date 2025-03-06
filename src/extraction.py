import os
import re

import numpy as np
import pandas as pd
import scipy.stats as st
from scipy import special

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

def pattern_sequence_match(pattern_matrix, background, sequence_tensor):
    """
    Matches pattern to all sequences in a sample. Returns a list of the same
        length as the sample, containing an index-matched 0 for non-matching
        sequences or 1 for matching sequences.

    Parameters:
        pattern -- Pattern instance.

    Returns:
        pattern_sequence_matches -- List.
    """
        
    # Get constituent pattern.
    try:
        constituent_pattern = get_pattern_constituents(
            pattern_matrix,
            background
        )
    except AttributeError:
        constituent_pattern = pattern_matrix

    # Intersect pattern with sample sequence tensor.
    positional_matches = np.sum(
        np.any(
            np.logical_and(
                sequence_tensor,
                constituent_pattern
            ),
            axis=2
        ).astype(np.int8),
        axis=1
    )
    pattern_sequence_matches  = np.zeros_like(positional_matches)
    pattern_sequence_matches[
        positional_matches == np.sum(np.any(constituent_pattern, axis=1).astype(np.int8))
    ] = 1

    return pattern_sequence_matches

def pattern_to_regular_expression(character_pattern):
    '''
    Converts series-style patterns to regular expressions.

    Parameters:
        patterns -- List.
    
    Returns:
        regular_expressions -- List.
    '''

    regular_expression = r"(?=("
    separation = 0
    for position in character_pattern.tolist():
        stripped_position = position.replace('[', '').replace(']', '')
        if stripped_position != '.':
            if separation != 0:
                regular_expression += (
                    r"(.){"
                    + re.escape(str(separation))+ r"}" 
                    + position
                )
                separation = 0
            else:
                regular_expression += position
        else:
            separation += 1
    regular_expression += r"))"
    
    return regular_expression

