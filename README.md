# RoLiM_Standalone

RoLiM (Robust Likelihood Method) is a standalone tool for the analysis of residue enrichment in protein sequences. The tool is designed to identify positions and residues that are significantly enriched in a foreground dataset compared to a background dataset. RoLiM uses a robust likelihood method to calculate the enrichment scores and p-values for each position and residue in the sequences. The analysis results include the enriched positions and residues, their statistical significance, and other relevant information.

## Installation

To install the RoLiM standalone tool, follow these steps:

| WIP: Installation instructions will be provided here. |

## Example Usage

The following example contains all the required and optional parameters for the main function of the RoLiM standalone tool. The example demonstrates how to set up an analysis with the specified parameters.

```bash
python -m src.run_rolim \
    --analysis_name "prealigned_test" \
    --foreground_format "prealigned" \
    --foreground_filename "data/prealigned_text_file.txt" \
    --context_format "fasta" \
    --context_filename "data/uniprot.fasta" \
    --width 15 \
    --output_path "data" \
    --precomputed "" \
    --p_value_cutoff 0.001 \
    --fold_change_cutoff 1.0 \
    --correction_method "fdr_bh" \
    --min_occurrence 20 \
    --extension_direction "N" \
    --center_sequences \
    --positional_weighting \
    --enable_compound_grouping \
    --position_specific \
    --redundancy_level "none" \
    --first_protein_only \
    --original_row_merge "all" \
    --generate_logo_maps \
    --verbose
```

This can be run from the main directory and don't need to be in the src directory. The `--verbose` flag is optional and can be used to display additional information during the analysis. The example demonstrates how to set up an analysis with prealigned sequences as the foreground input, fasta sequences as the context input, and various optional parameters to customize the analysis.

The following is a non-prealigned example:

```bash
python -m src.run_rolim \
    --analysis_name "peptideList_test" \
    --foreground_format "peptide_list" \
    --foreground_filename "data/text_file_peptide_list.txt" \
    --context_format "fasta" \
    --context_filename "data/uniprot.fasta" \
    --width 8 \
    --output_path "data" \
    --precomputed "" \
    --p_value_cutoff 0.001 \
    --fold_change_cutoff 1.0 \
    --correction_method "fdr_bh" \
    --min_occurrence 20 \
    --extension_direction "N" \
    --center_sequences \
    --positional_weighting \
    --enable_compound_grouping \
    --position_specific \
    --redundancy_level "none" \
    --first_protein_only \
    --original_row_merge "all" \
    --generate_logo_maps \
    --verbose
```

## Detailed Description of Parameters

This section describes the required and optional parameters for the main function, including their types, requirements, and validations.

### I. Required Parameters

* **A. `analysis_name`**
    * **Description:** The name of the analysis. Used for creating the output directory at the specified output path.
    * **Type:** String
    * **Requirements:** Required
    * **Validations:** Must be a string.
    * **Documentation:** The analysis name is used to create a directory in the output path where the results of the analysis will be stored. This directory will contain the results of the enrichment analysis, including the enriched positions and residues, as well as other relevant information saved.

* **B. `foreground_format`**
    * **Description:** The format of the foreground input.
    * **Type:** String
    * **Requirements:** Required
    * **Validations:** Must be one of: `prealigned`, `fasta`, `peptide_list`.
    * **Documentation:** The foreground format specifies the type of input data provided for the analysis. The supported formats include prealigned sequences, fasta files, and peptide lists. (The examples are provided in the GitHub readme). The format determines how the input data is processed and analyzed in the enrichment analysis.

* **C. `foreground_filename`**
    * **Description:** The foreground input as a file path or string representing a path.
    * **Type:** String
    * **Requirements:** Required
    * **Validations:** Must be a valid file path or a string that can be resolved to a valid file path.
    * **Documentation:** The foreground filename specifies the path to the input data file or the string representing the path to the input data. The input data contains the sequences or peptides that will be analyzed for enrichment. The format of the input data is determined by the `foreground_format` parameter.

* **D. `context_format`**
    * **Description:** The format of the context input (background or reference sequence).
    * **Type:** String
    * **Requirements:** Required
    * **Validations:** Must be one of: `sequences`, `fasta`.
    * **Documentation:** The context format specifies the type of input data provided for the background or reference sequence. The supported formats include a list of sequences or a fasta file containing the reference sequences. The context data is used to calculate the background frequencies for the enrichment analysis.

* **E. `context_filename`**
    * **Description:** The context input as a file path or string representing a path.
    * **Type:** String
    * **Requirements:** Required
    * **Validations:** Must be a valid file path or a string that can be resolved to a valid file path.
    * **Documentation:** The context filename specifies the path to the input data file or the string representing the path to the input data. The input data contains the reference sequences used to calculate the background frequencies for the enrichment analysis. The format of the input data is determined by the `context_format` parameter.

### II. Optional Parameters

* **A. `output_path`**
    * **Description:** The path to the output directory.
    * **Type:** String
    * **Default:** Current Working Directory
    * **Requirements:** Optional
    * **Validations:** Must be a string representing a valid directory path. If not provided, the current working directory is used.
    * **Documentation:** The output path specifies the directory where the result folder of the enrichment analysis will be saved. If the output path is not provided, the results will be saved in the current working directory.

* **B. `p_value_cutoff`**
    * **Description:** The p-value cutoff for the enrichment analysis.
    * **Type:** Float
    * **Default:** `0.001`
    * **Requirements:** Optional
    * **Validations:** Must be a float between 0 and 1 (inclusive).
    * **Documentation:** The p-value cutoff is used to filter the enriched positions and residues based on their statistical significance. Positions and residues with p-values below the cutoff are considered statistically significant and are included in the results.

* **C. `fold_change_cutoff`**
    * **Description:** The fold change cutoff for the enrichment analysis.
    * **Type:** Float
    * **Default:** `1.0`
    * **Requirements:** Optional
    * **Validations:** Must be a float greater than 0.
    * **Documentation:** The fold change cutoff is used to filter the enriched positions and residues based on their fold change values. Positions and residues with fold changes above the cutoff are considered enriched and are included in the results.

* **D. `correction_method`**
    * **Description:** The method used to correct for multiple testing.
    * **Type:** String
    * **Default:** `fdr_bh`
    * **Requirements:** Optional
    * **Validations:** Must be one of: `bonferroni`, `sidak`, `holm`, `fdr_bh`, `qvalue`.
    * **Documentation:** The correction method is used to adjust the p-values for multiple testing in the enrichment analysis. The available methods include Bonferroni, Sidak, Holm, Benjamini-Hochberg (fdr_bh), and Q-value. The correction method helps control the false discovery rate in the analysis.

* **E. `min_occurrence`**
    * **Description:** The minimum frequency of a position/residue pair in the foreground data set required for enrichment analysis.
    * **Type:** Integer
    * **Default:** `20`
    * **Requirements:** Optional
    * **Validations:** Must be an integer greater than 1.
    * **Documentation:** The minimum occurrence parameter specifies the minimum frequency of a position/residue pair in the foreground data set required for the pair to be considered enriched. Positions and residues with frequencies below this threshold are not included in the enrichment results.

* **F. `max_depth`**
    * **Description:**  (Seems not used, not sure what it is)
    * **Type:**  (Unknown - likely can be blank/None)
    * **Default:** (Likely blank/None)
    * **Requirements:** Optional
    * **Validations:** (Likely no validation if unused)
    * **Documentation:** (Unknown option)

* **G. `extend_sequences`**
    * **Description:** (Unknown option)
    * **Type:** Boolean
    * **Default:** `False`
    * **Requirements:** Optional
    * **Validations:** Must be a boolean (True/False or 1/0).
    * **Documentation:** (Unknown option)

* **H. `extension_direction`**
    * **Description:** Specifies the direction in which the sequences are extended.
    * **Type:** String
    * **Default:** `N`
    * **Requirements:** Optional
    * **Validations:** Must be one of: `N`, `C`, `both`.
    * **Documentation:** The direction in which the sequences are extended. By default, Rolim supports the alignment and extension of unaligned foreground datasets using the selected context dataset. This parameter defines the direction in which the sequences are extended. 'N' extends the sequences to the N-terminus, 'C' extends the sequences to the C-terminus, and 'both' extends the sequences in both directions (produces a separate sequence for each extension direction).

* **I. `width`**
    * **Description:** The width of the sequence window used for the enrichment analysis.
    * **Type:** Integer
    * **Default:** `8`
    * **Requirements:** Optional
    * **Validations:** Must be an integer greater than 2. Must be 8 for the enhanced analysis with Merops.
    * **Documentation:** The width of the sequence window used for the enrichment analysis. If peptides are supplied and extension/alignment is enabled, this is the final length of each extended and aligned sequence. If pre-aligned sequences are supplied, each supplied sequence must be of this length.

* **J. `center_sequences`**
    * **Description:** Enables the centering of the sequence position number.
    * **Type:** Boolean
    * **Default:** `True`
    * **Requirements:** Optional
    * **Validations:** Must be a boolean (True/False or 1/0).
    * **Documentation:** Centers the sequence position number. If True, the sequence position number is centered in the sequence window (p4, p3, p2, p1, p1', p2', p3', p4'). If False, the sequence position number is at the N-terminus of the sequence window (p1, p2, p3, p4, p5, p6, p7, p8). The advanced analysis with Merops requires the sequence position number to be centered.

* **K. `positional_weighting`**
    * **Description:** Enables positional weighting calculation
    * **Type:** Boolean
    * **Default:** `True`
    * **Requirements:** Optional
    * **Validations:** Must be a boolean (True/False or 1/0).
    * **Documentation:** Enables the optional positional weighting term in positional residue enrichment calculation. Positional weight is calculated as (1/#distinct residues in position).

* **L. `enable_compound_grouping`**
    * **Description:** Enables compound residue grouping.
    * **Type:** Boolean
    * **Default:** `False`
    * **Requirements:** Optional
    * **Validations:** Must be a boolean (True/False or 1/0).
    * **Documentation:** Enables aggregation of single amino acids into groups of biochemically and/or structurally similar amino acids, which may be cumulatively enriched in a position.

* **M. `compound_residues`**
    * **Description:** Path to the compound residues file.
    * **Type:** String
    * **Default:** `None`
    * **Requirements:** Optional
    * **Validations:** Must be a string representing a valid file path. The file needs to be structured specifically to be parsed (check for style in GitHub docs).
    * **Documentation:** Path to a file containing compound residues. If `enable_compound_grouping` is enabled and this is None, the default compound residues (as defined in `globals.py`) are used. The file should be structured in a specific format to be parsed correctly (details in the GitHub documentation).

* **N. `position_specific`**
    * **Description:** Enables position-specific analysis.
    * **Type:** Boolean
    * **Default:** `True`
    * **Requirements:** Optional
    * **Validations:** Must be a boolean (True/False or 1/0).
    * **Documentation:** Enables position-specific analysis. Enables background frequency calculation from a complete, position-specific background derived from the context dataset used for an analysis. When disabled, the background frequencies are averaged across all positions of the context dataset and dynamically updated when position/residue pairs are eliminated from the foreground dataset.

* **O. `redundancy_level`**
    * **Description:** Select sequence redundancy elimination level.
    * **Type:** String
    * **Default:** `none`
    * **Requirements:** Optional
    * **Validations:** Must be one of: `none`, `sequence`, `protein`.
    * **Documentation:** Selects the level of redundancy elimination in the foreground dataset. 'none' disables redundancy elimination, 'sequence' eliminates redundant sequences from the foreground dataset, and 'protein' eliminates redundant proteins (with the same id) from the foreground dataset.


