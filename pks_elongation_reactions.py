from raichu_drawer import *
from pikachu.smiles.smiles import *
from pikachu.chem.structure import *
from pikachu.reactions.functional_groups import BondDefiner
from central_atoms_pk_starter import find_central_atoms_pk_starter
#from pikachu.drawing.drawing import *

COABOND = BondDefiner('CoA_bond', 'CC(NCCC(NCCSC)=O)=O', 8, 9)
THIOESTERBOND = BondDefiner('thioester_bond', 'SC(C)=O', 0, 1)


def combine_structures(structures):
    """
    Returns a single combined Structure object from multiple input Structures

    Structures: tuple of to be combined Structure objects
    """
    new_graph = {}
    new_bonds = {}
    struct1, struct2 = structures
    bonds1 = struct1.bonds
    bonds2 = struct2.bonds
    new_bonds_1 = {}

    #make sure you dont add atoms with the same atom.nr
    atom_nrs = []
    for atom in struct2.graph:
        atom_nrs.append(atom.nr)
    index = 0
    for atom in struct1.graph:
        while index in atom_nrs:
            if index not in atom_nrs:
                atom.nr = index
                break
            index += 1
        else:
            atom.nr = index
            for other_atom in struct1.graph:
                if atom in other_atom.neighbours:
                    atom.nr = index
            atom_nrs.append(index)
            index += 1

    #refresh bond_lookup and structure
    struct1.make_bond_lookup()
    struct1.refresh_structure()

    # make sure you add no double bond nrs
    start = 0
    bonds = []
    bond_nrs = []
    for atom in struct2.graph:
        for bond in atom.bonds:
            if bond.nr not in bond_nrs:
                bond_nrs.append(bond.nr)
    new_nrs = []

    #Bond objects are identified only by Bond.nr, reset each Bond.nr first to
    #nr that is not used in both to-be combined structures, to make sure
    #no bonds are skipped in the loop
    for atom in struct1.graph:
        for bond in atom.bonds[:]:
            bond.nr = 9999999
    for atom in struct1.graph:
        for bond in atom.bonds[:]:
            if bond not in bonds:
                while start in bond_nrs:
                    if start not in bond_nrs and start not in new_nrs:
                        bond.nr = start
                        bonds.append(bond)
                        bond_nrs.append(start)
                        new_nrs.append(start)
                        break
                    start += 1
                else:
                    if start not in new_nrs:
                        bond.nr = start
                        bonds.append(bond)
                        bond_nrs.append(start)
                        new_nrs.append(start)
                        start += 1

    # Refresh second structure
    struct2.get_connectivities()
    struct2.set_connectivities()
    struct2.set_atom_neighbours()
    struct2.make_bond_lookup()

    #added 29-11-21: bond nrs are not refreshed as keys in structure.bonds dict
    struct1.bonds = {}
    for bond_nr, bond in bonds1.items():
        updated_bond_nr = bond.nr
        struct1.bonds[updated_bond_nr] = bond
    bonds1 = struct1.bonds

    #In order to prevent aliasing issues, make all bonds new Bond instances for
    #the bonds in structure 1:
    for nr, bond in bonds1.items():
        atom1 = bond.atom_1
        atom2 = bond.atom_2
        bond_type = bond.type
        bond_summary = bond.bond_summary
        bond_electrons = bond.electrons
        bond_nr = bond.nr
        new_bond = Bond(atom1, atom2, bond_type, bond_nr)
        new_bond.bond_summary = bond_summary
        new_bonds_1[bond_nr] = new_bond
        new_bond.electrons = bond_electrons

    #Make new structure object of the combined structure
    bonds2.update(new_bonds_1)
    for structure in structures:
        new_graph.update(structure.graph)
        new_bonds.update(structure.bonds)
    new_structure = Structure(new_graph, bonds2)
    new_structure.make_bond_lookup()
    new_structure.set_connectivities()
    new_structure.refresh_structure()
    new_structure.find_cycles()
    for bond_nr, bond in new_structure.bonds.items():
        bond.set_bond_summary()
    return new_structure

def pks_elongation(pk_chain, elongation_monomer):
    """
    Returns a Structure object of the PK chain after a single elongation step
    elongation step.

    structure: PIKAChU Structure object of the RC(=O)S PK intermediate before
    the elongation step
    elongation_monomer: ['SMILES_elongation_unit', index_c_to_c, index_c_to_s']
    """
    # If this is is the first elongation reaction on the starter unit, define
    # central chain atoms in the starter unit
    pk_chain = find_central_atoms_pk_starter(pk_chain)

    # Reset atom colours to black
    for atom in pk_chain.graph:
        atom.draw.colour = 'black'

    # Defining the structure of the elongation units
    elongation_monomer_struct = Smiles(elongation_monomer[0]).smiles_to_structure()
    for atom in elongation_monomer_struct.graph:
        if atom.nr == elongation_monomer[1]:
            # C0 needs to be attached to the C atom of the C-S bond in the PK chain
            c_to_pkchain = atom
            c_to_pkchain.in_central_chain = True
            for atom in c_to_pkchain.neighbours:
                if atom.type == 'H':
                    h_to_remove = atom
                    continue

        elif atom.nr == elongation_monomer[2]:
            # C1 needs to be attached to the S atom of the C-S bond in the PK chain
            c_to_s = atom
            c_to_s.in_central_chain = True
            for atom in c_to_s.neighbours:
                if atom.type == 'H':
                    h_to_remove2 = atom
                    continue

    #If the elongation unit contains an unknown moiety, add a number to the
    #'*' atom to display in the structure drawing
    pk_chain_atom_types = []
    for atom in pk_chain.graph:
        pk_chain_atom_types.append(atom.type)
    nr_unknown_atoms = pk_chain_atom_types.count('*')
    if elongation_monomer[0] == 'O=CC*':
        for atom in elongation_monomer_struct.graph:
            if atom.type == '*':
                atom.unknown_index = nr_unknown_atoms + 1


    #Remove the Hs in the malonylunit in order to add it to the pk chain
    elongation_monomer_struct.remove_atom(h_to_remove)
    elongation_monomer_struct.set_atom_neighbours()
    elongation_monomer_struct.remove_atom(h_to_remove2)
    for atom in elongation_monomer_struct.graph:
        atom.set_connectivity()
    elongation_monomer_struct.set_atom_neighbours()

    #find thioester bond in growing chain, define Cs to attach new unit
    thioester_bonds = find_bonds(pk_chain, THIOESTERBOND)
    for bond in thioester_bonds:
        if bond.atom_1.type == 'S':
            s_pkchain = bond.atom_1
            c_pkchain = bond.atom_2
        elif bond.atom_2.type == 'S':
            s_pkchain = bond.atom_2
            c_pkchain = bond.atom_1

    #Breaking the thioester bond in the PK chain
    for bond in thioester_bonds[:]:
        pk_chain.break_bond(bond)
    pk_chain.get_connectivities()
    pk_chain.set_connectivities()
    pk_chain.set_atom_neighbours()
    pk_chain.make_bond_lookup()

    #refresh malonylunit
    elongation_monomer_struct.refresh_structure()

    # Combining structures PK chain and elongation unit into one Structure object
    pk_chain_and_malonyl = (elongation_monomer_struct, pk_chain)
    combined = combine_structures(pk_chain_and_malonyl)
    combined.get_connectivities()
    combined.set_connectivities()
    combined.set_atom_neighbours()
    combined.make_bond_lookup()

    #Adding the bonds to form the new molecule after the single elongation step
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


def remove_coa(structure_with_coa):
    """
    Returns a Pikachu.structure object of the molecule where the CoA has been
    removed

    structure_with_coa: Pikachu.structure object of the molecule with CoA
    still attached
    """
    coa_bonds = find_bonds(structure_with_coa, COABOND)
    for bond in coa_bonds[:]:
        structure_with_coa.break_bond(bond)
    split = structure_with_coa.split_disconnected_structures()
    coa, remainder = split
    return remainder


def add_malonylunit(pk_chain):
    """
    Returns a Structure object of the PK chain after a single malonyl-CoA
    elongation step

    pk_chain: Structure object of the RC(=O)S PK intermediate before the
    elongation step
    """
    elongation_product = pks_elongation(pk_chain, ['CC=O', 0, 1])

    return elongation_product

def add_methylmalonylunit(pk_chain):
    """
    Returns a Structure object of the PK chain after a single methylmalonyl-CoA
    elongation step

    pk_chain: Structure object of the RC(=O)S PK intermediate before the
    elongation step
    """
    elongation_product = pks_elongation(pk_chain, ['O=CCC', 2, 1])

    return elongation_product


def add_methoxymalonylunit(pk_chain):
    """
    Returns a Structure object of the PK chain after a single
    methoxymalonyl-ACP elongation step

    pk_chain: Structure object of the RC(=O)S PK intermediate before the
    elongation step
    """
    elongation_product = pks_elongation(pk_chain, ['O=CCOC', 2, 1])

    return elongation_product

def add_ethylmalonylunit(pk_chain):
    """
    Returns a Structure object of the PK chain after a single
    ethylmalonyl-CoÄ elongation step

    pk_chain: Structure object of the RC(=O)S PK intermediate before the
    elongation step
    """
    elongation_product = pks_elongation(pk_chain, ['O=CCCC', 2, 1])

    return elongation_product

def add_pkunit(pk_chain):
    """
    Returns a Structure object of the PK chain after a single elongation step
    using the 'unknown wildcard polyketide elongation unit' pk

    pk_chain: Structure object of the RC(=O)S PK intermediate before the
    elongation step
    """
    elongation_product = pks_elongation(pk_chain, ['O=CC*', 2, 1])

    return elongation_product


def repeat_elongations(pk_chain, times):
    for i in range(times-1):
        pk_chain = add_methylmalonylunit(pk_chain)
    return pk_chain

if __name__ == "__main__":
    #Load structure starter unit malonyl-CoA
    #malonylcoa = Smiles('CC(C)(COP(=O)(O)OP(=O)(O)OCC1C(C(C(O1)N2C=NC3=C(N=CN=C32)N)O)OP(=O)(O)O)C(C(=O)NCCC(=O)NCCSC(=O)CC(=O)O)O').smiles_to_structure()
    #Drawer(malonylcoa)
    #Remove CoA-group from malonyl-CoA
    #malonyl = remove_coa(malonylcoa)


    malonyl = Smiles('SC(CC)=O').smiles_to_structure()
    #Drawer(malonyl)

    #Performing a single elongation step on the malonyl using a malonyl-CoA
    new_chain = add_malonylunit(malonyl)
    Drawer(new_chain)

