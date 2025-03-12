# Convenience script to run the ROLIM analysis with the peptide list test data
# Allows better modifiability of the parameters, detailed in the README.md
# To run the script, execute the following command in the terminal:
# chmod +x run.sh && ./run.sh

python -m src.run_rolim \
    --analysis_name "250304_host_EVOSEP30" \
    --foreground_format "peptide_list" \
    --foreground_filename "data/2025-03-04_host_EVOSEP30k_high_noncanon_cleaned.txt" \
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
    --redundancy_level "protein" \
    --first_protein_only \
    --original_row_merge "all" \
    --generate_logo_maps \
    --verbose