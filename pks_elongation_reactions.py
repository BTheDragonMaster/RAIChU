from pikachu.general import read_smiles
from pikachu.reactions.functional_groups import BondDefiner, combine_structures
from central_atoms_pk_starter import find_central_atoms_pk_starter
from attributes import ATTRIBUTES

COABOND = BondDefiner('CoA_bond', 'CC(NCCC(NCCSC)=O)=O', 8, 9)
THIOESTERBOND = BondDefiner('thioester_bond', 'SC(C)=O', 0, 1)


PKS_SUBUNIT_TO_MONOMER = {'malonylcoa': ['CC=O', 0, 1],
                          'methylmalonylcoa': ['O=CCC', 2, 1],
                          'methoxymalonylacp': ['O=CCOC', 2, 1],
                          'ethylmalonylcoa': ['O=CCCC', 2, 1],
                          'pk': ['O=CC*', 2, 1]}


def pks_elongation(pk_chain, elongation_monomer):
    """
    Returns a Structure object of the PK chain after a single elongation step
    elongation step.

    structure: PIKAChU Structure object of the RC(=O)S PK intermediate before
    the elongation step
    elongation_monomer: ['SMILES_elongation_unit', index_c_to_c, index_c_to_s']
    """

    # If this is the first elongation reaction on the starter unit, define
    # central chain atoms in the starter unit
    if not any(atom.annotations.in_central_chain for atom in pk_chain.graph):
        pk_chain = find_central_atoms_pk_starter(pk_chain)

    # Reset atom colours to black
    for atom in pk_chain.graph:
        atom.draw.colour = 'black'

    # Defining the structure of the elongation units
    elongation_monomer_struct = read_smiles(elongation_monomer[0])

    # Add annotation attributes to elongation monomer
    elongation_monomer_struct.add_attributes(ATTRIBUTES, boolean=True)

    for atom in elongation_monomer_struct.graph:
        if atom.nr == elongation_monomer[1]:
            # C0 needs to be attached to the C atom of the C-S bond in the PK chain
            c_to_pkchain = atom

            c_to_pkchain.annotations.in_central_chain = True
            if elongation_monomer[0] == 'O=CCC':
                c_to_pkchain.chiral = 'clockwise'
            for atom in c_to_pkchain.neighbours:
                if atom.type == 'H':
                    h_to_remove = atom
                    continue

        elif atom.nr == elongation_monomer[2]:
            # C1 needs to be attached to the S atom of the C-S bond in the PK chain
            c_to_s = atom
            c_to_s.annotations.in_central_chain = True
            for atom in c_to_s.neighbours:
                if atom.type == 'H':
                    h_to_remove2 = atom
                    continue

    # If the elongation unit contains an unknown moiety, add a number to the
    # '*' atom to display in the structure drawing
    pk_chain_atom_types = []
    for atom in pk_chain.graph:
        pk_chain_atom_types.append(atom.type)
    nr_unknown_atoms = pk_chain_atom_types.count('*')
    if elongation_monomer[0] == 'O=CC*':
        for atom in elongation_monomer_struct.graph:
            if atom.type == '*':
                atom.annotations.unknown_index = nr_unknown_atoms + 1


    # Remove the Hs in the malonylunit in order to add it to the pk chain
    elongation_monomer_struct.remove_atom(h_to_remove)
    elongation_monomer_struct.set_atom_neighbours()
    elongation_monomer_struct.remove_atom(h_to_remove2)
    for atom in elongation_monomer_struct.graph:
        atom.set_connectivity()
    elongation_monomer_struct.set_atom_neighbours()

    # Find thioester bond in growing chain, define Cs to attach new unit
    pk_chain.set_connectivities()
    pk_chain.refresh_structure()
    thioester_bonds = find_bonds(pk_chain, THIOESTERBOND)
    assert len(thioester_bonds) == 1
    for bond in thioester_bonds:
        if bond.atom_1.type == 'S':
            s_pkchain = bond.atom_1
            c_pkchain = bond.atom_2
        elif bond.atom_2.type == 'S':
            s_pkchain = bond.atom_2
            c_pkchain = bond.atom_1

    # Breaking the thioester bond in the PK chain
    for atom in pk_chain.graph:
        atom.hybridise()
    for bond in thioester_bonds[:]:
        pk_chain.break_bond(bond)
    pk_chain.get_connectivities()
    pk_chain.set_connectivities()
    pk_chain.set_atom_neighbours()
    pk_chain.make_bond_lookup()

    # Refresh malonylunit

    elongation_monomer_struct.refresh_structure()

    # Combining structures PK chain and elongation unit into one Structure object
    pk_chain_and_malonyl = (elongation_monomer_struct, pk_chain)

    combined = combine_structures(pk_chain_and_malonyl)

    combined.get_connectivities()
    combined.set_connectivities()
    combined.set_atom_neighbours()
    combined.make_bond_lookup()

    # Adding the bonds to form the new molecule after the single elongation step
    new_bond_nr = combined.find_next_bond_nr()
    combined.make_bond(c_to_pkchain, c_pkchain, new_bond_nr)
    new_bond_nr = combined.find_next_bond_nr()
    combined.make_bond(c_to_s, s_pkchain, new_bond_nr)
    combined.get_connectivities()
    combined.set_connectivities()
    combined.set_atom_neighbours()
    combined.make_bond_lookup()
    for bond_nr, bond in combined.bonds.items():
        bond.set_bond_summary()

    return combined


def find_bonds(structure, bond_type):
    """
    Returns a list of Pikachu.bond objects of the bonds of the indicated type
    in the structure of interest

    structure: Pikachu.structure object of the molecule of interest
    bond_type: BondDefiner object indicating the bond type that is searched
    in the indicated structure
    """
    for atom in structure.graph:
        atom.get_connectivity()
    locations = structure.find_substructures(bond_type.structure)
    bonds = []
    for match in locations:
        atom_1 = match.atoms[bond_type.atom_1]
        atom_2 = match.atoms[bond_type.atom_2]
        bond = structure.bond_lookup[atom_1][atom_2]
        bonds.append(bond)

    return bonds


# def add_malonylunit(pk_chain):
#     """
#     Returns a Structure object of the PK chain after a single malonyl-CoA
#     elongation step
#
#     pk_chain: Structure object of the RC(=O)S PK intermediate before the
#     elongation step
#     """
#     elongation_product = pks_elongation(pk_chain, ['CC=O', 0, 1])
#
#     return elongation_product
#
#
# def add_methylmalonylunit(pk_chain):
#     """
#     Returns a Structure object of the PK chain after a single methylmalonyl-CoA
#     elongation step
#
#     pk_chain: Structure object of the RC(=O)S PK intermediate before the
#     elongation step
#     """
#     elongation_product = pks_elongation(pk_chain, ['O=CCC', 2, 1])
#
#     return elongation_product
#
#
# def add_methoxymalonylunit(pk_chain):
#     """
#     Returns a Structure object of the PK chain after a single
#     methoxymalonyl-ACP elongation step
#
#     pk_chain: Structure object of the RC(=O)S PK intermediate before the
#     elongation step
#     """
#     elongation_product = pks_elongation(pk_chain, ['O=CCOC', 2, 1])
#
#     return elongation_product
#
#
# def add_ethylmalonylunit(pk_chain):
#     """
#     Returns a Structure object of the PK chain after a single
#     ethylmalonyl-CoÄ elongation step
#
#     pk_chain: Structure object of the RC(=O)S PK intermediate before the
#     elongation step
#     """
#     elongation_product = pks_elongation(pk_chain, ['O=CCCC', 2, 1])
#
#     return elongation_product
#
#
# def add_pkunit(pk_chain):
#     """
#     Returns a Structure object of the PK chain after a single elongation step
#     using the 'unknown wildcard polyketide elongation unit' pk
#
#     pk_chain: Structure object of the RC(=O)S PK intermediate before the
#     elongation step
#     """
#     elongation_product = pks_elongation(pk_chain, ['O=CC*', 2, 1])
#
#     return elongation_product



