from pikachu.reactions.basic_reactions import condensation
from pikachu.reactions.functional_groups import find_bonds, BondDefiner, GroupDefiner, find_atoms
from pikachu.smiles.smiles import Smiles
from pikachu.reactions.basic_reactions import combine_structures
from class_domain import ATTRIBUTES


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
    # Initialize annotations in all atoms in amino acid
    for atom in amino_acid.graph:
        if not hasattr(atom.annotations, 'in_central_chain'):
            for attribute in ATTRIBUTES:
                atom.annotations.add_annotation(attribute, False)

    # If this is the first elongation reaction, determine central peptide
    if not any(atom.annotations.in_central_chain for atom in nrp_intermediate.graph):
        n_atoms_aa = find_atoms(N_AMINO_ACID, nrp_intermediate)
        c1_atoms_aa = find_atoms(C1_AMINO_ACID, nrp_intermediate)
        c2_atoms_aa = find_atoms(C2_AMINO_ACID, nrp_intermediate)
        assert len(n_atoms_aa) == 1
        assert len(c1_atoms_aa) == 1
        assert len(c2_atoms_aa) == 1
        for atom in nrp_intermediate.graph:
            if atom == n_atoms_aa[0]:
                atom.annotations.in_central_chain = True
            elif atom == c1_atoms_aa[0]:
                atom.annotations.in_central_chain = True
            elif atom == c2_atoms_aa[0]:
                atom.annotations.in_central_chain = True
            else:
                atom.annotations.in_central_chain = False

    # Check if the intermediate is a thioester intermediate, if so: convert
    is_thioester = False
    found_bonds_thioester = find_bonds(THIOESTERBOND, nrp_intermediate)
    if len(found_bonds_thioester) > 0:
        is_thioester = True
        nrp_intermediate, oh_bond = sulphur_to_hydroxyl(nrp_intermediate)
        for atom in nrp_intermediate.graph:
            if not hasattr(atom.annotations, 'in_central_chain'):
                for attribute in ATTRIBUTES:
                    atom.annotations.add_annotation(attribute, False)


    # Define reaction targets
    if not is_thioester:
        found_bonds = find_bonds(LEAVING_OH_BOND, nrp_intermediate)
        assert len(found_bonds) == 1
        oh_bond = found_bonds[0]

    # If the amino acid unit contains an unknown moiety, add a number to the
    # '*' atom to display in the structure drawing
    pk_chain_atom_types = []
    for atom in nrp_intermediate.graph:
        pk_chain_atom_types.append(atom.type)
    nr_unknown_atoms = pk_chain_atom_types.count('*')
    for atom in amino_acid.graph:
        if atom.type == '*':
            atom.unknown_index = nr_unknown_atoms + 1

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
            atom.annotations.in_central_chain = True
        elif atom == c1_atoms_aa[0]:
            atom.annotations.in_central_chain = True
        elif atom == c2_atoms_aa[0]:
            atom.annotations.in_central_chain = True
        else:
            atom.annotations.in_central_chain = False

    # Carry out condensation reaction using build-in PIKAChU function
    condensation_product = condensation(nrp_intermediate, amino_acid, oh_bond, h_bond)[0]

    # Refresh condensation product
    condensation_product.refresh_structure()
    condensation_product.set_connectivities()
    condensation_product.set_atom_neighbours()
    condensation_product.find_cycles()

    # Initialize annotations for atoms that don't have any yet
    for atom in condensation_product.graph:
        if not hasattr(atom.annotations, 'in_central_chain'):
            for attribute in ATTRIBUTES:
                atom.annotations.add_annotation(attribute, False)

    return condensation_product


def sulphur_to_hydroxyl(thioester_structure):
    """
    Identifies and removes the sulphur (and attached atoms/domains) in the
    thioester structure and replaces it with a hydroxyl group, creating a
    carboxylic acid group. Returns the product of this reaction, also as a
    PIKAChU Structure object.

    thioester_structure: PIKAChU Structure object of a thioester (R-C(S)=O)
    """
    # Find thioester bond in the input structure
    found_bonds_thioester = find_bonds(THIOESTERBOND, thioester_structure)
    found_carbon_thioester = find_atoms(THIOESTER_CARBON, thioester_structure)
    assert len(found_carbon_thioester)
    assert len(found_bonds_thioester) == 1

    carbon_thioester = found_carbon_thioester[0]
    sh_bond = found_bonds_thioester[0]
    thioester_structure.break_bond(sh_bond)

    # Break S-C bond in thioester
    one, two = thioester_structure.split_disconnected_structures()
    if len(one.graph) < 3:
        residual = one
        thioester_structure = two
    else:
        residual = two
        thioester_structure = one

    # Create Structure object hydroxyl group
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

    # Refresh combined Structure object (no new bond yet)
    combined.get_connectivities()
    combined.set_connectivities()
    combined.set_atom_neighbours()
    combined.make_bond_lookup()
    next_bond_nr = combined.find_next_bond_nr()

    # Define carbon and oxygen atom of the new hydroxyl group, make bond
    # between these atoms
    for atom in combined.graph:
        atom_neighbours = []
        atom_neighbour_types = []
        for neighbour in atom.neighbours:
            atom_neighbour_types.append(neighbour.type)
            atom_neighbours.append(neighbour)
        if atom.type == 'C' and len(atom_neighbours) == 2 and\
                atom_neighbour_types.count('O') == 1 and\
                atom_neighbour_types.count('C') == 1:
            carbon_thioester = atom
        elif atom.type == 'O' and len(atom_neighbours) == 1 and\
                atom_neighbour_types.count('H') == 1:
            oxygen_hydroxyl = atom
    combined.make_bond(carbon_thioester, oxygen_hydroxyl, next_bond_nr)


    combined.set_atoms()
    combined.make_bond_lookup()
    combined.set_connectivities()
    combined.set_atom_neighbours()
    combined.find_cycles()
    for bond_nr, bond in combined.bonds.items():
        bond.set_bond_summary()
    oh_bond = combined.bond_lookup[oxygen_hydroxyl][carbon_thioester]

    for atom in combined.graph:
        if not hasattr(atom.annotations, 'in_central_chain'):
            for attribute in ATTRIBUTES:
                atom.annotations.add_annotation(attribute, False)


    return combined, oh_bond







