from pikachu.reactions.basic_reactions import condensation
from pikachu.reactions.functional_groups import find_bonds, BondDefiner, GroupDefiner, find_atoms
from pikachu.general import draw_structure
from pikachu.smiles.smiles import Smiles
from pk_attach_to_domain import attach_to_domain_nrp
from raichu_drawer import Drawer
from find_central_peptide_chain import find_central_chain_nrp

LEAVING_OH_BOND = BondDefiner('Leaving -OH group bond', 'C(=O)(O)CN', 0, 2)
N_AMINO_ACID = GroupDefiner('Nitrogen atom amino acid', 'NCC(=O)O', 0)

def condensation_nrps(amino_acid, nrp_intermediate):
    """
    Returns the NRPS condensation product as a PIKAChU Structure object

    nrp_intermediate: PIKAChU Structure object of the NRP intermediate
    amino_acid: PIKAChU Structure object of the amino acid
    """
    # Define the bond attached to the -OH leaving group
    # nrp_intermediate.set_connectivities()
    # nrp_intermediate.find_cycles()
    # nrp_intermediate.refresh_structure()
    # nrp_intermediate.set_atom_neighbours()
    found_bonds = find_bonds(LEAVING_OH_BOND, nrp_intermediate)
    print(found_bonds, 'found bonds')
    assert len(found_bonds) == 1
    oh_bond = found_bonds[0]



    # Define the bond attached to the -H leaving group
    n_atoms_aa = find_atoms(N_AMINO_ACID, amino_acid)
    assert len(n_atoms_aa) == 1
    n_atom = n_atoms_aa[0]
    for bond in n_atom.bonds:
        for neighbour in bond.neighbours:
            if neighbour.type == 'H':
                h_bond = bond
                break

    # Carry out condensation reaction using build-in PIKAChU function
    condensation_product = condensation(nrp_intermediate, amino_acid, oh_bond, h_bond)[0]
    condensation_product.refresh_structure()
    condensation_product.set_connectivities()
    condensation_product.set_atom_neighbours()
    condensation_product.find_cycles()
    print('ja')
    return condensation_product

def make_nrp(list_amino_acids):
    """
    Returns the NRP product as Structure object

    list_amino_acids: [strings] of amino acid names
    """
    # Parse list SMILES amino acids attached to PCP to dict name -> Structure
    lines_aa = open('PARAS_smiles.txt', 'r', encoding='utf8').readlines()
    dict_aa_structure = {}

    for line in lines_aa:
        line = line.strip()
        name, smiles = line.split()
        name = name.upper()
        # structure_aa = Smiles(smiles).smiles_to_structure()
        # structure_aa.find_cycles()
        dict_aa_structure[name] = smiles

    # Make every name in amino acid list all caps
    for i in range(len(list_amino_acids)):
        list_amino_acids[i] = list_amino_acids[i].upper()

    # Take amino acid Structure object from dict and add to growing NRP chain
    nrp_chain_intermediate = Smiles(dict_aa_structure[list_amino_acids[0]]).smiles_to_structure()
    list_amino_acids = list_amino_acids[1:]
    for amino_acid_name in list_amino_acids:
        amino_acid_struct = Smiles(dict_aa_structure[amino_acid_name]).smiles_to_structure()
        print(amino_acid_name)
        nrp_chain_intermediate = condensation_nrps(amino_acid_struct, nrp_chain_intermediate)

    # Refresh chain intermediate
    nrp_chain_intermediate.refresh_structure()
    nrp_chain_intermediate.set_connectivities()
    nrp_chain_intermediate.set_atom_neighbours()
    nrp_chain_intermediate.find_cycles()

    return nrp_chain_intermediate


if __name__ == "__main__":
    test_peptide = make_nrp(['alanine', 'valine', 'tyrosine', 'citrulline',  'threonine', 'cysteine', 'norcoronamicacid', '(2S,3R)-2-amino-3-hydroxy-4-(4-nitrophenyl)butanoate'])
    attached_test_peptide = attach_to_domain_nrp(test_peptide, 'PCP')
    Drawer(attached_test_peptide)
    test_peptide2 = make_nrp(['d-threonine', 'valine', 'cysteine'])
    attached_test_peptide2 = attach_to_domain_nrp(test_peptide2, 'PCP')
    Drawer(attached_test_peptide2)
