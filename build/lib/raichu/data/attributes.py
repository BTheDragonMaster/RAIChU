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

AMINOACID_ONE_LETTER_TO_SMILES = {"A": "N[C@@]([H])(C)C(=O)", "C": "N[C@@]([H])(CS)C(=O)",
                                  "D": "N[C@@]([H])(CC(=O)O)C(=O)", "E": "N[C@@]([H])(CCC(=O)O)C(=O)",
                                  "F": "N[C@@]([H])(Cc1ccccc1)C(=O)", "G": "NCC(=O)", "H": "N[C@@]([H])(CC1=CN=C-N1)C(=O)",
                                  "I": "N[C@@]([H])([C@]([H])(CC)C)C(=O)", "K": "N[C@@]([H])(CCCCN)C(=O)",
                                  "L": "N[C@@]([H])(CC(C)C)C(=O)", "M": "N[C@@]([H])(CCSC)C(=O)",
                                  "N": "N[C@@]([H])(CC(=O)N)C(=O)", "P": "N1[C@@]([H])(CCC1)C(=O)",
                                  "Q": "N[C@@]([H])(CCC(=O)N)C(=O)", "R": "N[C@@]([H])(CCCNC(=N)N)C(=O)",
                                  "S": "N[C@@]([H])(CO)C(=O)", "T": "N[C@@]([H])([C@]([H])(O)C)C(=O)",
                                  "V": "N[C@@]([H])(C(C)C)C(=O)", "W": "N[C@@]([H])(CC(=CN2)C1=C2C=CC=C1)C(=O)",
                                  "Y": "N[C@@]([H])(Cc1ccc(O)cc1)C(=O)"}

PRENYL_TRANSFERASE_SUBSTRATES_TO_SMILES = {
    "DIMETHYLALLYL": r"CC=C(C)C",
    "GERANYL": r"C\C=C(/C)CCC=C(C)C",
    "FARNESYL": r"C\C=C(/C)CC\C=C(/C)CCC=C(C)C",
    "GERANYLGERANYL": r"CC(C)=CCC\C(C)=C\CC\C=C(/C)CCC=C(C)C",
    "SQUALENE": r"CC(=CCC/C(=C/CC/C(=C/CC/C=C(/CC/C=C(/CCC=C(C)C)\C)\C)/C)/C)C",
    "PHYTOENE": r"CC(=CCC/C(=C/CC/C(=C/CC/C(=C/C=C\C=C(/C)\CC/C=C(\C)/CC/C=C(\C)/CCC=C(C)C)/C)/C)/C)C"
}

SUGARS_TO_SMILES = {

}