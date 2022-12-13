from typing import Tuple

from pikachu.reactions.functional_groups import find_atoms, find_bonds, combine_structures
from raichu.reactions.general import initialise_atom_attributes
from pikachu.general import read_smiles
from raichu.data.attributes import ATTRIBUTES

def hydroxylation(target_atom, structure):
    """
    Returns the hydroxylated structure thats hydroxylated at the target atom.

    structure: PIKAChU Structure object
    target_atom:  PIKAChU atom object
    """

    hydroxyl_group = read_smiles('o')
    hydroxyl_group.add_attributes(ATTRIBUTES, boolean=True)
    oxygen = hydroxyl_group.atoms[0]
    hydrogen_1 = hydroxyl_group.atoms[1]
    bond_1 = oxygen.get_bond(hydrogen_1)
    hydrogen_2 = target_atom.get_neighbour('H')
    if not hydrogen_2:
        raise Exception("Can't oxidate this atom!")

    bond_2 = target_atom.get_bond(hydrogen_2)

    hydroxylated_structure = combine_structures([structure, hydroxyl_group])

    hydroxylated_structure.break_bond(bond_1)
    hydroxylated_structure.break_bond(bond_2)

    hydroxylated_structure.make_bond(oxygen, target_atom, hydroxylated_structure.find_next_bond_nr())
    hydroxylated_structure.make_bond(hydrogen_1, hydrogen_2,hydroxylated_structure.find_next_bond_nr())

    structures = hydroxylated_structure.split_disconnected_structures()

    for s in structures:
        if oxygen.nr in s.atoms:
            return s


def methylation(target_atom, structure):
    """
    Returns the structure thats methylated at the target atom.

    structure: PIKAChU Structure object
    target_atom:  PIKAChU atom object
    """


    methyl_group = read_smiles('C')
    methyl_group.add_attributes(ATTRIBUTES, boolean=True)
    carbon = methyl_group.atoms[0]
    hydrogen_1 = methyl_group.atoms[1]
    bond_1 = carbon.get_bond(hydrogen_1)

    hydrogen_2 = target_atom.get_neighbour('H')
    if not hydrogen_2:
        raise Exception("Can't methylate this atom!")

    bond_2 = target_atom.get_bond(hydrogen_2)

    methylated_structure = combine_structures([structure, methyl_group])

    methylated_structure.break_bond(bond_1)
    methylated_structure.break_bond(bond_2)

    methylated_structure.make_bond(carbon, target_atom, methylated_structure.find_next_bond_nr())
    methylated_structure.make_bond(hydrogen_1, hydrogen_2, methylated_structure.find_next_bond_nr())

    structures = methylated_structure.split_disconnected_structures()

    for s in structures:
        if carbon.nr in s.atoms:
            return s
