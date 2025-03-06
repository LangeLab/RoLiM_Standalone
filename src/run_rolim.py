import os
import argparse
import utils # Custom module

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
    )
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
    
    # Parse command-line arguments
    try:
        args = parse_arguments()
        verbose = args.verbose
        if verbose:
            print("Arguments parsed successfully:")
            print(args)
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
    utils.record_compound_residues_to_file(compound_residues, args.output_path)
    ## Step 3: Establish background instance