ATTRIBUTES = ['in_central_chain', 'KR_ep_target', 'KR_red_target',
              'latest_elongation_o', 'latest_elongation_methyl', 'DH_target',
              'ER_target', 'domain_type', 'unknown_index', 'chiral_c_ep',
              'n_atom_nmeth', 'leaving_oh_o', 'leaving_oh_h', 'c2_acid',
              'terminal_c', 'terminal_o']

ALL_PKS_ELONGATION_UNITS = ['malonylcoa', 'methylmalonylcoa',
                            'methoxymalonylacp', 'ethylmalonylcoa', 'pk']

AMINOACID_ONE_LETTER_TO_NAME = {"A": "alanine",
                                "R": "arginine",
                                "N": "asparagine",
                                "D": "aspartic acid",
                                "C": "cysteine",
                                "E": "glutamic acid",
                                "Q": "glutamine",
                                "G": "glycine",
                                "H": "histidine",
                                "I": "isoleucine",
                                "L": "leucine",
                                "K": "lysine",
                                "M": "methionine",
                                "F": "phenylalanine",
                                "P": "proline",
                                "S": "serine",
                                "T": "threonine",
                                "W": "tryptophan",
                                "Y": "tyrosine",
                                "V": "valine"}
PRENYL_TRANSFERASE_SUBSTRATES_TO_SMILES = {
    "DIMETHYLALLYL": r"CC=C(C)C",
    "GERANYL": r"C\C=C(/C)CCC=C(C)C",
    "FARNESYL": r"C\C=C(/C)CC\C=C(/C)CCC=C(C)C",
    "GERANYLGERANYL": r"CC(C)=CCC\C(C)=C\CC\C=C(/C)CCC=C(C)C",
    "SQUALENE": r"CC(=CCC/C(=C/CC/C(=C/CC/C=C(/CC/C=C(/CCC=C(C)C)\C)\C)/C)/C)C",
    "PHYTOENE": r"CC(=CCC/C(=C/CC/C(=C/CC/C(=C/C=C\C=C(/C)\CC/C=C(\C)/CC/C=C(\C)/CCC=C(C)C)/C)/C)/C)C"

}