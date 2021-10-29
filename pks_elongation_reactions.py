from pikachu.smiles.smiles import *
from pikachu.chem.structure import *
#from pikachu.drawing.drawing import *
from the_drawer import *



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


class BondDefiner():
    """
    Class to store an object that defines a bond type

    Attributes
    ----------
    name: str, name of the bond type
    smiles: Smiles object indicating the structure of the bond
    atom_1: int, index of atom 1 in the bond in the smiles string
    atom_2: int, index of atom 2 in the bond in the smiles string
    """
    def __init__(self, name, smiles, atom_nr_1, atom_nr_2):
        self.name = name
        self.smiles = Smiles(smiles)
        self.structure = self.smiles.smiles_to_structure()
        self.find_atoms(atom_nr_1, atom_nr_2)
        self.bond = self.structure.bond_lookup[self.atom_1][self.atom_2]

    def __repr__(self):
        return self.name

    def find_atoms(self, atom_nr_1, atom_nr_2):
        self.atom_1 = None
        self.atom_2 = None
        for atom in self.structure.graph:
            if atom.nr == atom_nr_1:
                self.atom_1 = atom
            elif atom.nr == atom_nr_2:
                self.atom_2 = atom
            if self.atom_1 and self.atom_2:
                break
        if not self.atom_1 or not self.atom_2:
            raise Exception("Can't find atoms adjacent to bond.")

COABOND = BondDefiner('CoA_bond', 'CC(NCCC(NCCSC)=O)=O', 8, 9)
THIOESTERBOND = BondDefiner('thioester_bond', 'SC(C)=O', 0, 1)


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
    # Reset atom colours to black
    for atom in pk_chain.graph:
        atom.draw.colour = 'black'
    # Defining the structure of the elongation units when malonyl-CoA is
    # used as the elongation unit during PK synthesis
    malonylunit = Smiles('CC=O').smiles_to_structure()
    for atom in malonylunit.graph:
        if atom.nr == 0:
            # C0 needs to be attached to the C atom of the C-S bond in the PK chain
            c_to_pkchain = atom
            c_to_pkchain.central_chain = True
            for atom in c_to_pkchain.neighbours:
                if atom.type == 'H':
                    h_to_remove = atom
                    continue
        elif atom.nr == 1:
            # C1 needs to be attached to the S atom of the C-S bond in the PK chain
            c_to_s = atom
            c_to_s.central_chain = True
            for atom in c_to_s.neighbours:
                if atom.type == 'H':
                    h_to_remove2 = atom
                    continue


    #Remove the Hs in the malonylunit in order to add it to the pk chain
    malonylunit.remove_atom(h_to_remove)
    malonylunit.set_atom_neighbours()
    malonylunit.remove_atom(h_to_remove2)
    for atom in malonylunit.graph:
        atom.set_connectivity()
    malonylunit.set_atom_neighbours()

    #make sure you dont add atoms with the same atom.nr
    atom_nrs = []
    for atom in pk_chain.graph:
        atom_nrs.append(atom.nr)
    start = 0
    for atom in malonylunit.graph:
        while start in atom_nrs:
            if start not in atom_nrs:
                atom.nr = start
                break
            start += 1
        else:
            atom.nr = start
            for other_atom in malonylunit.graph:
                if atom in other_atom.neighbours:
                    atom.nr = start
            atom_nrs.append(start)
            start += 1

    #refresh bond_lookup
    malonylunit.make_bond_lookup()

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
    malonylunit.refresh_structure()

    # make sure you add no doule bond nrs
    start = 0
    bonds = []
    bond_nrs = []
    for atom in pk_chain.graph:
        for bond in atom.bonds:
            if bond.nr not in bond_nrs:
                bond_nrs.append(bond.nr)
    for atom in malonylunit.graph:
        for bond in atom.bonds:
            if bond not in bonds:
                while start in bond_nrs:
                    if start not in bond_nrs:
                        bond.nr = start
                        bonds.append(bond)
                        break
                    start += 1
                else:
                    bond.nr = start
                    bonds.append(bond)
                    bond_nrs.append(start)
                    start += 1
    pk_chain.get_connectivities()
    pk_chain.set_connectivities()
    pk_chain.set_atom_neighbours()
    pk_chain.make_bond_lookup()

    # Combining structures PK chain and elongation unit into one Structure object
    pk_chain_and_malonyl = (malonylunit, pk_chain)
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

def add_methylmalonylunit(pk_chain):
    """
    Returns a Structure object of the PK chain after a single methylmalonyl-CoA
    elongation step

    pk_chain: Structure object of the RC(=O)S PK intermediate before the
    elongation step
    """
    # Reset atom colours to black
    for atom in pk_chain.graph:
        atom.draw.colour = 'black'
    # Defining the structure of the elongation units when methylmalonyl-CoA is
    # used as the elongation unit during PK synthesis
    methylmalonylunit = Smiles('O=CCC').smiles_to_structure()
    for atom in methylmalonylunit.graph:
        if atom.nr == 2:
            # C0 needs to be attached to the C atom of the C-S bond in the PK chain
            c_to_pkchain = atom
            c_to_pkchain.central_chain = True
            c_to_pkchain.chiral = 'counterclockwise'
            for atom in c_to_pkchain.neighbours:
                if atom.type == 'H':
                    h_to_remove = atom
                    continue
        elif atom.nr == 1:
            # C1 needs to be attached to the S atom of the C-S bond in the PK chain
            c_to_s = atom
            c_to_s.central_chain = True
            for atom in c_to_s.neighbours:
                if atom.type == 'H':
                    h_to_remove2 = atom
                    continue


    #Remove the Hs in the malonylunit in order to add it to the pk chain
    methylmalonylunit.remove_atom(h_to_remove)
    methylmalonylunit.set_atom_neighbours()
    methylmalonylunit.remove_atom(h_to_remove2)
    for atom in methylmalonylunit.graph:
        atom.set_connectivity()
    methylmalonylunit.set_atom_neighbours()

    #make sure you dont add atoms with the same atom.nr
    atom_nrs = []
    for atom in pk_chain.graph:
        atom_nrs.append(atom.nr)
    start = 0
    for atom in methylmalonylunit.graph:
        while start in atom_nrs:
            if start not in atom_nrs:
                atom.nr = start
                break
            start += 1
        else:
            atom.nr = start
            for other_atom in methylmalonylunit.graph:
                if atom in other_atom.neighbours:
                    atom.nr = start
            atom_nrs.append(start)
            start += 1

    #refresh bond_lookup
    methylmalonylunit.make_bond_lookup()

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

    #refresh methylmalonylunit
    methylmalonylunit.refresh_structure()

    # make sure you add no doule bond nrs
    start = 0
    bonds = []
    bond_nrs = []
    for atom in pk_chain.graph:
        for bond in atom.bonds:
            if bond.nr not in bond_nrs:
                bond_nrs.append(bond.nr)
    for atom in methylmalonylunit.graph:
        for bond in atom.bonds:
            if bond not in bonds:
                while start in bond_nrs:
                    if start not in bond_nrs:
                        bond.nr = start
                        bonds.append(bond)
                        break
                    start += 1
                else:
                    bond.nr = start
                    bonds.append(bond)
                    bond_nrs.append(start)
                    start += 1
    pk_chain.get_connectivities()
    pk_chain.set_connectivities()
    pk_chain.set_atom_neighbours()
    pk_chain.make_bond_lookup()

    # Combining structures PK chain and elongation unit into one Structure object
    pk_chain_and_malonyl = (methylmalonylunit, pk_chain)
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
    combined.find_cycles()
    for bond_nr, bond in combined.bonds.items():
        bond.set_bond_summary()

    return combined


def add_methoxymalonylunit(pk_chain):
    """
    Returns a Structure object of the PK chain after a single
    methoxymalonyl-ACP elongation step

    pk_chain: Structure object of the RC(=O)S PK intermediate before the
    elongation step
    """
    # Reset atom colours to black
    for atom in pk_chain.graph:
        atom.draw.colour = 'black'
    # Defining the structure of the elongation units when methylmalonyl-CoA is
    # used as the elongation unit during PK synthesis
    methoxymalonylunit = Smiles('O=CCOC').smiles_to_structure()
    for atom in methoxymalonylunit.graph:
        if atom.nr == 2:
            # C0 needs to be attached to the C atom of the C-S bond in the PK chain
            c_to_pkchain = atom
            c_to_pkchain.central_chain = True
            c_to_pkchain.chiral = 'counterclockwise'
            for atom in c_to_pkchain.neighbours:
                if atom.type == 'H':
                    h_to_remove = atom
                    continue
        elif atom.nr == 1:
            # C1 needs to be attached to the S atom of the C-S bond in the PK chain
            c_to_s = atom
            c_to_s.central_chain = True
            for atom in c_to_s.neighbours:
                if atom.type == 'H':
                    h_to_remove2 = atom
                    continue


    #Remove the Hs in the malonylunit in order to add it to the pk chain
    methoxymalonylunit.remove_atom(h_to_remove)
    methoxymalonylunit.set_atom_neighbours()
    methoxymalonylunit.remove_atom(h_to_remove2)
    for atom in methoxymalonylunit.graph:
        atom.set_connectivity()
    methoxymalonylunit.set_atom_neighbours()

    #make sure you dont add atoms with the same atom.nr
    atom_nrs = []
    for atom in pk_chain.graph:
        atom_nrs.append(atom.nr)
    start = 0
    for atom in methoxymalonylunit.graph:
        while start in atom_nrs:
            if start not in atom_nrs:
                atom.nr = start
                break
            start += 1
        else:
            atom.nr = start
            for other_atom in methoxymalonylunit.graph:
                if atom in other_atom.neighbours:
                    atom.nr = start
            atom_nrs.append(start)
            start += 1

    #refresh bond_lookup
    methoxymalonylunit.make_bond_lookup()

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

    #refresh methylmalonylunit
    methoxymalonylunit.refresh_structure()

    # make sure you add no double bond nrs
    start = 0
    bonds = []
    bond_nrs = []
    for atom in pk_chain.graph:
        for bond in atom.bonds:
            if bond.nr not in bond_nrs:
                bond_nrs.append(bond.nr)
    for atom in methoxymalonylunit.graph:
        for bond in atom.bonds:
            if bond not in bonds:
                while start in bond_nrs:
                    if start not in bond_nrs:
                        bond.nr = start
                        bonds.append(bond)
                        break
                    start += 1
                else:
                    bond.nr = start
                    bonds.append(bond)
                    bond_nrs.append(start)
                    start += 1
    pk_chain.get_connectivities()
    pk_chain.set_connectivities()
    pk_chain.set_atom_neighbours()
    pk_chain.make_bond_lookup()

    # Combining structures PK chain and elongation unit into one Structure object
    pk_chain_and_malonyl = (methoxymalonylunit, pk_chain)
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
    combined.find_cycles()
    for bond_nr, bond in combined.bonds.items():
        bond.set_bond_summary()

    return combined



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


    malonyl = Smiles('SC(CC(=O)O)=O').smiles_to_structure()
    #Drawer(malonyl)

    #Performing a single elongation step on the malonyl using a malonyl-CoA
    new_chain = add_malonylunit(malonyl)
    #Drawer(new_chain)

    #Performing a second elongation step on the malonyl using malonyl-CoA
    new_chain2 = add_malonylunit(new_chain)
    #Drawer(new_chain2)

    #Perform a third elongation step on the chain intermediate using malonyl-CoA
    new_chain3 = add_malonylunit(new_chain2)
    #Drawer(new_chain3)

    #Perform a fourth elongation step on the chain intermediate using methylmalonyl-CoA
    new_chain4 = add_methylmalonylunit(new_chain3)
    #Drawer(new_chain4)

    #Test methylmalonyl-CoA elongation repetition function:
    new_chain5 = repeat_elongations(new_chain4, 10)
    Drawer(new_chain5)



    #Perform methoxymalonyl-ACP elongation reactoin
    new_chain6 = add_methoxymalonylunit(new_chain5)
    Drawer(new_chain6)


