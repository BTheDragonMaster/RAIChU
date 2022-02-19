from pikachu.reactions.basic_reactions import condensation
from pikachu.reactions.functional_groups import find_bonds, BondDefiner, GroupDefiner, find_atoms
from pikachu.general import draw_structure
from pikachu.smiles.smiles import Smiles
from pk_attach_to_domain import attach_to_domain_nrp
from pikachu.reactions.basic_reactions import combine_structures
from raichu_drawer import Drawer
from random import randint
from pikachu.chem.structure import Structure

THIOESTERBOND = BondDefiner('thioester_bond', 'SC(C)=O', 0, 1)
THIOESTER_CARBON = GroupDefiner('thioester carbon', 'SC(C)=O', 1)
LEAVING_OH_BOND = BondDefiner('Leaving -OH group bond', 'C(=O)(O)CN', 0, 2)
N_AMINO_ACID = GroupDefiner('Nitrogen atom amino acid', 'NCC(=O)O', 0)
C1_AMINO_ACID = GroupDefiner('C1 atom amino acid', 'NCC(=O)O', 1)
C2_AMINO_ACID = GroupDefiner('C2 atom amino acid', 'NCC(=O)O', 2)

def condensation_nrps(amino_acid, nrp_intermediate):
    """
    Returns the NRPS condensation product as a PIKAChU Structure object

    nrp_intermediate: PIKAChU Structure object of the NRP intermediate
    amino_acid: PIKAChU Structure object of the amino acid
    """
    # If this is the first elongation reaction, determine central peptide
    if not any(hasattr(atom, 'in_central_chain') for atom in nrp_intermediate.graph):
        n_atoms_aa = find_atoms(N_AMINO_ACID, nrp_intermediate)
        c1_atoms_aa = find_atoms(C1_AMINO_ACID, nrp_intermediate)
        c2_atoms_aa = find_atoms(C2_AMINO_ACID, nrp_intermediate)
        assert len(n_atoms_aa) == 1
        assert len(c1_atoms_aa) == 1
        assert len(c2_atoms_aa) == 1
        for atom in nrp_intermediate.graph:
            if atom == n_atoms_aa[0]:
                atom.in_central_chain = True
            elif atom == c1_atoms_aa[0]:
                atom.in_central_chain = True
            elif atom == c2_atoms_aa[0]:
                atom.in_central_chain = True
            else:
                atom.in_central_chain = False

    # Check if the intermediate is a thioester intermediate, if so: convert
    is_thioester = False
    found_bonds_thioester = find_bonds(THIOESTERBOND, nrp_intermediate)
    if len(found_bonds_thioester) > 0:
        is_thioester = True
        nrp_intermediate, oh_bond = sulphur_to_hydroxyl(nrp_intermediate)


    # Define reaction targets
    if not is_thioester:
        found_bonds = find_bonds(LEAVING_OH_BOND, nrp_intermediate)
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

    #Determine atoms in amino acid that end up in central peptide chain
    n_atoms_aa = find_atoms(N_AMINO_ACID, amino_acid)
    c1_atoms_aa = find_atoms(C1_AMINO_ACID, amino_acid)
    c2_atoms_aa = find_atoms(C2_AMINO_ACID, amino_acid)
    assert len(n_atoms_aa) == 1
    assert len(c1_atoms_aa) == 1
    assert len(c2_atoms_aa) == 1
    for atom in amino_acid.graph:
        if atom == n_atoms_aa[0]:
            atom.in_central_chain = True
        elif atom == c1_atoms_aa[0]:
            atom.in_central_chain = True
        elif atom == c2_atoms_aa[0]:
            atom.in_central_chain = True
        else:
            atom.in_central_chain = False

    # Carry out condensation reaction using build-in PIKAChU function
    condensation_product = condensation(nrp_intermediate, amino_acid, oh_bond, h_bond)[0]
    condensation_product.refresh_structure()
    condensation_product.set_connectivities()
    condensation_product.set_atom_neighbours()
    condensation_product.find_cycles()

    return condensation_product


def sulphur_to_hydroxyl(thioester_structure):
    """
    Identifies and removes the sulphur (and attached atoms/domains) in the
    thioester structure and replaces it with a hydroxyl group, creating a
    carboxylic acid group. Returns the product of this reaction, also as a
    PIKAChU Structure object.

    thioester_structure: PIKAChU Structure object of a thioester (R-C(S)=O)
    """
    found_bonds_thioester = find_bonds(THIOESTERBOND, thioester_structure)
    found_carbon_thioester = find_atoms(THIOESTER_CARBON, thioester_structure)
    assert len(found_carbon_thioester)
    assert len(found_bonds_thioester) == 1
    carbon_thioester = found_carbon_thioester[0]
    sh_bond = found_bonds_thioester[0]
    thioester_structure.break_bond(sh_bond)
    one, two = thioester_structure.split_disconnected_structures()
    if len(one.graph) < 3:
        residual = one
        thioester_structure = two
    else:
        residual = two
        thioester_structure = one




    methanol = Smiles('CO').smiles_to_structure()
    for atom in methanol.graph:
        if atom.type == 'C':
            carbon_methanol = atom
            for neighbour in carbon_methanol.neighbours:
                if neighbour.type == 'O':
                    oxygen_hydroxyl = neighbour
    methanol.break_bond(methanol.bond_lookup[carbon_methanol][oxygen_hydroxyl])
    one, two = methanol.split_disconnected_structures()
    if len(one.graph) == 4:
        residual = one
        hydroxyl = two
    else:
        residual = two
        hydroxyl = one


    # Add hydroxyl group to the carbon atom of the intermediate to create a
    # carboxylic acid group
    combined = combine_structures((hydroxyl, thioester_structure))
    combined.get_connectivities()
    combined.set_connectivities()
    combined.set_atom_neighbours()
    combined.make_bond_lookup()
    for atom in combined.graph:
        atom_neighbours = []
        atom_neighbour_types = []
        for neighbour in atom.neighbours:
            atom_neighbour_types.append(neighbour.type)
            atom_neighbours.append(neighbour)
        if atom.type == 'C' and len(atom_neighbours) == 2 and atom_neighbour_types.count('O') == 1 and atom_neighbour_types.count('C') == 1:
            carbon_thioester = atom
        elif atom.type == 'O' and len(atom_neighbours) == 1 and atom_neighbour_types.count('H') == 1:
            oxygen_hydroxyl = atom
    next_bond_nr = combined.find_next_atom_nr()
    combined.make_bond(carbon_thioester, oxygen_hydroxyl, next_bond_nr)

    combined.set_atoms()
    combined.make_bond_lookup()
    combined.set_connectivities()
    combined.set_atom_neighbours()
    combined.find_cycles()
    for bond_nr, bond in combined.bonds.items():
        bond.set_bond_summary()
    oh_bond = combined.bond_lookup[oxygen_hydroxyl][carbon_thioester]



    return combined, oh_bond

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

    # Make amino acid Structure object from dict and add to growing NRP chain
    nrp_chain_intermediate = Smiles(dict_aa_structure[list_amino_acids[0]]).smiles_to_structure()

    # Check if it is a amino acid
    if not (len(nrp_chain_intermediate.find_substructures(Smiles('CN').smiles_to_structure())) > 0\
            and len(nrp_chain_intermediate.find_substructures(Smiles('C(O)=O').smiles_to_structure())) > 0):
        raise ValueError(f'The starter structure: {list_amino_acids[0]}, is not an amino acid')

    # Determine cental peptide chain atoms in first amino acid
    n_atoms_aa = find_atoms(N_AMINO_ACID, nrp_chain_intermediate)
    c1_atoms_aa = find_atoms(C1_AMINO_ACID, nrp_chain_intermediate)
    c2_atoms_aa = find_atoms(C2_AMINO_ACID, nrp_chain_intermediate)
    assert len(n_atoms_aa) == 1
    assert len(c1_atoms_aa) == 1
    assert len(c2_atoms_aa) == 1
    for atom in nrp_chain_intermediate.graph:
        if atom == n_atoms_aa[0]:
            atom.in_central_chain = True
        elif atom == c1_atoms_aa[0]:
            atom.in_central_chain = True
        elif atom == c2_atoms_aa[0]:
            atom.in_central_chain = True
        else:
            atom.in_central_chain = False

    list_amino_acids = list_amino_acids[1:]
    for amino_acid_name in list_amino_acids:
        amino_acid_struct = Smiles(dict_aa_structure[amino_acid_name]).smiles_to_structure()
        if not (len(amino_acid_struct.find_substructures(
                Smiles('CN').smiles_to_structure())) > 0
                and len(amino_acid_struct.find_substructures(
                    Smiles('C(O)=O').smiles_to_structure())) > 0):
            raise ValueError(
                f'The structure: {list_amino_acids[0]}, is not an amino acid')
        print(amino_acid_name)
        n_atoms_aa = find_atoms(N_AMINO_ACID, amino_acid_struct)
        c1_atoms_aa = find_atoms(C1_AMINO_ACID, amino_acid_struct)
        c2_atoms_aa = find_atoms(C2_AMINO_ACID, amino_acid_struct)
        assert len(n_atoms_aa) == 1
        assert len(c1_atoms_aa) == 1
        assert len(c2_atoms_aa) == 1
        for atom in amino_acid_struct.graph:
            if atom == n_atoms_aa[0]:
                atom.in_central_chain = True
            elif atom == c1_atoms_aa[0]:
                atom.in_central_chain = True
            elif atom == c2_atoms_aa[0]:
                atom.in_central_chain = True
            else:
                atom.in_central_chain = False
        nrp_chain_intermediate = condensation_nrps(amino_acid_struct, nrp_chain_intermediate)

    # Refresh chain intermediate
    nrp_chain_intermediate.refresh_structure()
    nrp_chain_intermediate.set_connectivities()
    nrp_chain_intermediate.set_atom_neighbours()
    nrp_chain_intermediate.find_cycles()

    return nrp_chain_intermediate


if __name__ == "__main__":
    # nrp = make_nrp(['valine', 'glycine', 'cysteine'])
    # att = attach_to_domain_nrp(nrp, 'PCP')
    # Drawer(att)
    # test_peptide = make_nrp(['alanine', 'valine', 'tyrosine', 'citrulline',  'threonine', 'cysteine', 'norcoronamicacid', '(2S,3R)-2-amino-3-hydroxy-4-(4-nitrophenyl)butanoate'])
    # attached_test_peptide = attach_to_domain_nrp(test_peptide, 'PCP')
    # Drawer(attached_test_peptide)
    # test_peptide2 = make_nrp(['d-threonine', 'valine', 'cysteine'])
    # attached_test_peptide2 = attach_to_domain_nrp(test_peptide2, 'PCP')
    # Drawer(attached_test_peptide2)
    peptide = make_nrp(['alanine', 'glycolicacid', 'proline'])
    # Drawer(peptide)
    attached = attach_to_domain_nrp(peptide, 'PCP')
    print(attached.graph)
    Drawer(attached)



