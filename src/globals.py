import re
from collections import namedtuple

## Define the global variables

# TODO: This is a placeholder would need to make use of precomputed files
PRECOMPUTED_FILES = {
    'MEROPS': 'data/merops_precomputed.csv',
    'PHOSPHO': 'data/phospho_precomputed.csv',
    'SWISSPROT': 'data/swissprot_precomputed.csv',
    'UNIPROT': 'data/uniprot_precomputed.csv',
    'swissprot_human': 'data/swissprot_human_precomputed.csv',
}

# regex pattern for swissprot (uniprot) accession
SWISSPROT_ACCESSION_PATTERN = re.compile(
    r"[OPQ][0-9][A-Z0-9]{3}[0-9](-\d+)?|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}(-\d+)?"
) # Added a way to handle isoform addons

# Define the amino acids single and three letter codes
SINGLE_LETTER_CODES = [
    'A','R','N','D','B','C','E','Q','Z','G','H','I','L','K',
    'M','F','P','S','T','W','Y','V',
]
THREE_LETTER_CODES = [
    'Ala','Arg','Asn','Asp','Asx','Cys','Glu','Gln','Glx',
    'Gly','His','Ile','Leu','Lys','Met','Phe','Pro','Ser',
    'Thr','Trp','Tyr','Val',
]
# Dictionary for single:three
encode_1to3 = {
    single:THREE_LETTER_CODES[i]
    for i, single in enumerate(SINGLE_LETTER_CODES)
}
# Dictionary for three:single
encode_3to1 = {
    three:single
    for single, three in encode_1to3.items()
}

# Site and position variables
SITES = [
    "Site_P4", "Site_P3", "Site_P2", "Site_P1", "Site_P1prime",
    "Site_P2prime", "Site_P3prime", "Site_P4prime",
]

MEROPS_POSITIONS = [
    "p4","p3","p2","p1","p1'","p2'","p3'","p4'",
]

PHOSPHO_POSITIONS = [
    "p1","p2","p3","p4","p5","p6","p7",
    "p8","p9","p10","p11","p12","p13",
]

POSITIONS = MEROPS_POSITIONS
IGNORE_POSITIONS = ['p7']

# Catalytic types from MEROPS for label in plots
CATALYTIC_TYPES = {
    'A':'Aspartic', 'C':'Cysteine','G':'Glutamic','M':'Metallo',
    'N':'Asparagine','P':'Mixed','S':'Serine','T':'Threonine',
    'U':'Unknown',
}
# Pyro-glu definitions
PYRO_GLU = [
    'Glu->pyro-Glu',
    'Gln->pyro-Glu',
    '2 Oxidation (M),Gln->pyro-Glu',
    'Oxidation (M),Gln->pyro-Glu',
    'Oxidation (M),Glu->pyro-Glu',
]

# define namedTuples for the data
CompoundResidue = namedtuple(
    'CompoundResidue',
    ['description', 'residues']
)

ExtendedSequences = namedtuple(
    'ExtendedSequences',
    ['n_term', 'c_term']
)

Sample = namedtuple(
    'Sample',
    ['sequence_df', 'sequence_tensor', 'original_sequences'],
    defaults=[None]
)

ParentStats = namedtuple(
    'ParentStats',
    ['pattern_id', 'size', 'bonferroni_m', 'expected_frequency']
)

PatternFrequency = namedtuple(
    'PatternFrequency',
    ['absolute', 'percentage']
)

PatternSimilarityMatrix = namedtuple(
    'PatternSimilarityMatrix',
    ['matrix', 'num_clustered', 'num_unclustered']
)

## Default compound residues
COMPOUND_RESIDUES = {
    '1': CompoundResidue(
        description='Nonpolar, aliphatic R groups',
        residues=[ 'G', 'A', 'V', 'L', 'M', 'I', ]
    ),
    '2': CompoundResidue(
        description='Aromatic R groups',
        residues=[ 'F', 'Y', 'W', ]
    ),
    '3': CompoundResidue(
        description='Polar, uncharged R groups',
        residues=[ 'S', 'T', 'C', 'P', 'N', 'Q', ]
    ),
    '4': CompoundResidue(
        description='Positively charged R groups',
        residues=[ 'K', 'R', 'H', ]
    ),
    '5': CompoundResidue(
        description='Negatively charged R groups',
        residues=[ 'D', 'E', ]
    )
}
