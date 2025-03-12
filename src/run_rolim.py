import os
import re
import argparse
import pandas as pd
from src import utils, sequence, extraction # Custom modules

def parse_arguments() -> argparse.Namespace:
    """
    Parses command-line arguments and performs validations.

    Returns:
        argparse.Namespace: The parsed arguments.

    Raises:
        ValueError: If any of the argument validations fail.
    """

    parser = argparse.ArgumentParser(
        description="Positional Residue Enrichment Analysis"
    )

    # Required arguments
    required = parser.add_argument_group('Required Arguments')
    required.add_argument(
        "-a", "--analysis_name", 
        required=True, 
        type=str, 
        help="The name of the analysis."
    )
    required.add_argument(
        "-ff", "--foreground_format", 
        required=True, 
        type=str, 
        choices=['prealigned', 'fasta', 'peptide_list'], 
        help="The format of the foreground input."
    ) # TODO: Need to add more options as shown in original
    required.add_argument(
        "-fn", "--foreground_filename", 
        required=True, 
        type=str, 
        help="The foreground input filename."
    )
    required.add_argument(
        "-cf", "--context_format", 
        required=True, 
        type=str, 
        choices=['sequences', 'fasta'], 
        help="The format of the context input."
    )
    required.add_argument(
        "-cn", "--context_filename", 
        required=True, 
        type=str, 
        help="The context input filename."
    )

    # Optional arguments
    optional = parser.add_argument_group('Optional Arguments')
    optional.add_argument(
        "-o", "--output_path", 
        type=str, 
        default=os.getcwd(), 
        help="The path to the output directory (default: current directory)."
    )
    optional.add_argument(
        "-pc", "--precomputed",
        type=str,
        default=None,
        help="The path to the precomputed file for various uses. (default: None)"
    )
    optional.add_argument(
        "-p", "--p_value_cutoff", 
        type=float, 
        default=0.001, 
        help="The p-value cutoff for the enrichment analysis (default: 0.001)."
    )
    optional.add_argument(
        "-fc", "--fold_change_cutoff", 
        type=float, 
        default=1.0, 
        help="The fold change cutoff for the enrichment analysis (default: 1.0)."
    )
    optional.add_argument(
        "-cm", "--correction_method", 
        type=str, 
        default='fdr_bh', 
        choices=['bonferroni', 'sidak', 'holm', 'fdr_bh', 'qvalue'], 
        help="The method used to correct for multiple testing (default: fdr_bh)."
    )
    optional.add_argument(
        "-mo", "--min_occurrence", 
        type=int, 
        default=20, 
        help="The minimum frequency of a position/residue pair (default: 20)."
    )
    optional.add_argument(
        "--max_depth", 
        help=argparse.SUPPRESS
    )  # Hidden option
    optional.add_argument(
        "--extend_sequences", 
        action='store_true', 
        default=False, 
        help=argparse.SUPPRESS
    )  # Hidden option
    optional.add_argument(
        "-ed", "--extension_direction", 
        type=str, 
        default='N', 
        choices=['N', 'C', 'both'], 
        help="The direction in which the sequences are extended (default: N)."
    )
    optional.add_argument(
        "-w", "--width", 
        type=int, 
        default=15, 
        help="The width of the sequence window (default: 15)."
    )  # Default changed to 15
    optional.add_argument(
        "--center_sequences", 
        action='store_true', 
        default=True, 
        help="Centers the sequence position number (default: True)."
    )
    optional.add_argument(
        "--positional_weighting", 
        action='store_true', 
        default=True, 
        help="Enables positional weighting calculation (default: True)."
    )
    optional.add_argument(
        "--enable_compound_grouping", 
        action='store_true', 
        default=False, 
        help="Enables compound grouping (default: False)."
    )
    optional.add_argument(
        "--compound_residues", 
        type=str, 
        default=None, 
        help="Path to the compound residues file."
    )
    optional.add_argument(
        "--position_specific", 
        action='store_true', 
        default=True, 
        help="Enables position specific analysis (default: True)."
    )
    optional.add_argument(
        "-rl", "--redundancy_level", 
        type=str, 
        default='none', 
        choices=['none', 'sequence', 'protein'], 
        help="Select sequence redundancy elimination level (default: none)."
    )
    optional.add_argument(
        "--first_protein_only", 
        action='store_true', 
        default=True, 
        help="Select first protein only (default: True)."
    )
    optional.add_argument(
        "--original_row_merge", 
        type=str, 
        default='all', 
        choices=['all', 'none', 'protein'], 
        help="Select original row merge (default: all)."
    )
    optional.add_argument(
        "--generate_logo_maps",
        action='store_true',
        default=False,
        help="Generate logo maps (default: False)."
    )
    optional.add_argument(
        "-v", "--verbose", 
        action='store_true', 
        default=False, 
        help="Prints detailed information during the analysis to the command line."
    )

    args = parser.parse_args()

    # Validations
    if not os.path.exists(args.output_path):
        raise ValueError("Output path does not exist.")

    if not 0 <= args.p_value_cutoff <= 1:
        raise ValueError("P-value cutoff must be between 0 and 1.")

    if args.fold_change_cutoff <= 0:
        raise ValueError("Fold change cutoff must be greater than 0.")

    if args.min_occurrence <= 1:
        raise ValueError("Minimum occurrence must be greater than 1.")

    if args.width <= 2:
        raise ValueError("Width must be greater than 2.")
    
    if args.width == 8 and not args.center_sequences:
        raise ValueError("For width 8, center_sequences must be True.")

    if args.compound_residues and not os.path.exists(args.compound_residues):
        raise ValueError("Compound residues file does not exist.")

    return args


# Run the script
if __name__ == "__main__":
    st = utils.getTime()
    # Parse command-line arguments
    try:
        args = parse_arguments()
        verbose = args.verbose
        if verbose:
            print("Arguments parsed successfully")
            print("Starting the analysis...\n")
            # print(args)
    except ValueError as e:
        print(f"Error: {e}")
    except FileNotFoundError as e:
        print(f"Error: {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

    # Define paths

    ## Main output path
    output_path = os.path.join(args.output_path, args.analysis_name)
    ## Subdirectories
    ### Records (keeping record of used data, analysis parameters, etc.)
    record_path = os.path.join(output_path, 'records')
    if not os.path.exists(record_path):
        os.makedirs(record_path)
    ### Results (results of the analysis)
    results_path = os.path.join(output_path, 'results')
    if not os.path.exists(results_path):
        os.makedirs(results_path)
    
    # Run the complete analysis   

    ## Step 1: Load the Context Data
    if verbose:
        print("1 - Loading Context Data")
    ### Load the reference data from sequences
    if args.context_format == 'fasta':
        context = utils.parse_fasta(args.context_filename, swissprot=True)
    elif args.context_format == 'sequences':
        raise NotImplementedError("Sequence format not implemented yet.")
    
    ## Step 2: Load compound residues
    if verbose:
        print("2 - Loading Compound Residues")
    ### Load the compound residues either from the default or user-defined file
    compound_residues = utils.load_compound_residues(args.compound_residues)
    ### Save the compound residues definition as a file
    utils.record_compound_residues_to_file(compound_residues, record_path)
    
    ## Step 3: Establish background instance
    if verbose:
        print("3 - Establishing Background Instance")
    # Using the sequence.Background class to establish the 
    #   background instance for context(reference/background/fasta)
    background = sequence.Background(
        context['sequence'].tolist(),
        background_dict=None, # TODO: Not implemented yet
        compound_residues=compound_residues,
        position_specific=args.position_specific,
        width=args.width,
        center=args.center_sequences,
        initial_background_size_limit=None, # TODO: Not implemented yet
        # TODO: The fast=T is actually slower, need to investigate
        fast=False,
        # TODO: Adding the usage of pre-computed background
        precomputed=args.precomputed,
    )

    ## Step 4: Establish the foreground instance
    if verbose:
        print("4 - Establishing Foreground Instance")
    ### Check the format of the foreground data
    if args.foreground_format == 'prealigned':
        samples = sequence.load_prealigned_file(
            prealigned_file_path=args.foreground_filename,
            background=background,
            center=args.center_sequences,
            redundancy_level=args.redundancy_level,
            title=args.analysis_name,
        )
    elif args.foreground_format == 'fasta':
        raise NotImplementedError("Fasta format not implemented yet.")
    elif args.foreground_format == 'peptide_list':
        samples = sequence.load_peptide_list_file(
            args.foreground_filename,
            context=context,
            background=background,
            center=args.center_sequences,
            width=args.width,
            terminal=args.extension_direction,
            require_context_id=True,
            redundancy_level=args.redundancy_level,
            first_protein_only=args.first_protein_only,
            original_row_merge=args.original_row_merge,
            title=args.analysis_name,
            verbose=verbose,
        )
    else:
        raise ValueError("Invalid foreground format.")
    
    ## Step 5: Run pattern extraction on the foreground
    if verbose:
        print("5 - Running Pattern Extraction on Foreground")
    
    ### Split by Unique Sample conditions
    sample_output_paths = []
    for sample_name in samples.keys():
        # TODO: Following is likely cause issues...
        sample_directory = re.sub(r'\W+', ' ', sample_name).strip().replace(" ", "_")
        sample_output_path = os.path.join(results_path, sample_directory)
        sample_output_paths.append((sample_name, sample_output_path))
    
    ### Running the analysis for each sample
    all_pattern_containers = []
    summary_tables = []
    for sample_name, sample_output_path in sample_output_paths:
        if verbose: 
            print(f"\t - Extracting patterns for group {sample_name}...")
        patterns = extraction.PatternContainer(
            samples[sample_name],
            background,
            str(sample_name),
            sample_output_path,
            initial_pattern=None,
            initial_removed_positional_residues = None,
            max_depth=args.max_depth,
            p_value_cutoff=args.p_value_cutoff,
            minimum_occurrences=args.min_occurrence,
            fold_change_cutoff=args.fold_change_cutoff,
            multiple_testing_correction=args.correction_method,
            positional_weighting=args.positional_weighting,
            # Att: Likely Not used...
            allow_compound_residue_decomposition=args.enable_compound_grouping,
            set_reduction = True,
            verbose=verbose,
        )
        if (args.width == 8) and args.center_sequences:
            summary_tables.append(
                patterns.post_processing(
                    proteolysis_data=True,
                    cluster_sequences=True,
                    logo_maps=args.generate_logo_maps, 
                )
            )
        else:
            summary_tables.append(
                patterns.post_processing(
                    proteolysis_data=False,
                    cluster_sequences=False,
                    logo_maps=args.generate_logo_maps,
                )
            )

        all_pattern_containers.append(patterns)

    ## Step 6: Generate the summary table
    # TODO: Extremely slow, need to investigate...
    if verbose:
        print("6 - Generating Summary Table")
    # Add sequence summary table to summary output directory.
    if len(summary_tables) > 1:
        unique_pattern_strings = []
        unique_pattern_list = []
        for pattern_container in all_pattern_containers:
            for pattern in pattern_container.pattern_list:
                pattern_string = ''.join(pattern.character_pattern())
                if pattern_string not in unique_pattern_strings:
                    unique_pattern_strings.append(pattern_string)
                    unique_pattern_list.append(pattern)
        all_original_sequences = pd.concat(
            [sample.original_sequences for sample in samples.values()]
        )
        summary_table = extraction.PatternContainer.generate_summary_table(
            all_original_sequences,
            unique_pattern_list
        )
        summary_table.to_csv(
            os.path.join(results_path, f'{args.analysis_name}_summary_table.txt'),
            index=False, header=True, sep='\t'
        )
    
    ## Step 7: Record the analysis parameters
    if verbose:
        print("7 - Recording Analysis Parameters")
    elapsed = utils.getTime() - st
    if verbose:
        print(f"Elapsed time: {utils.prettyTimer(elapsed)}")
    # Call the log function
    utils.log_analysis_parameters(
        args, os.path.join(record_path, "analysis_parameters.log"), elapsed
    )

