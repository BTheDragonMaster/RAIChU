from typing import Tuple

from pikachu.reactions.functional_groups import find_atoms, find_bonds, combine_structures, BondDefiner, GroupDefiner
from raichu.reactions.general import initialise_atom_attributes
from pikachu.general import read_smiles
from raichu.data.attributes import ATTRIBUTES

SH_BOND = BondDefiner('recent_elongation', 'SC(C)=O', 0, 1)
CO_BOND = BondDefiner('recent_elongation', 'CO', 0, 1)
N_AMINO = GroupDefiner('N_amino', 'CN', 1)
O_OH = GroupDefiner('O_oh', 'CO', 1)
C_CH = GroupDefiner("C_Ch", "CC", 0)
O_BETAPROPRIOLACTONE = GroupDefiner('o_betapropriolactone', 'SC(CCO)=O', 4)


def find_o_betapropriolactone(polyketide):
    """
    Finds and returns the oxygen atom (PIKAChU atom object) in the -OH group
    that shouldn't be used by the thioesterase_circular_product function, as
    this will create a beta-propriolactone compound, which does not occur in
    polyketide synthesis, if present. Otherwise the function returns None

    polyketide: PIKAChU structure object of a polyketide
    """
    o_propriolactone = find_atoms(O_BETAPROPRIOLACTONE, polyketide)
    if len(o_propriolactone) == 0:
        return None
    elif len(o_propriolactone) == 1:
        return o_propriolactone[0]
    else:
        raise ValueError('Error: this molecule is not a polyketide, as the \
        carbon in the beta ketone/hydroxyl group is bound to an additional \
         oxygen atom')


def find_all_o_n_atoms_for_cyclization(chain_intermediate):
    """Performs all thioesterase reactions on the input chain_intermediate
     using all internal amino and -OH groups except for the -OH group that
     leads to the formation of a beta-propriolactone compound, which does not
     occur in polyketide synthesis. Returns a list of PIKAChU Structure objects
     of all possible thioesterase products.

     chain_intermediate: PIKAChU Structure object of a polyketide/NRP
    """
    # Perform first thioesterase reaction, generating linear polyketide/NRP
    chain_intermediate.refresh_structure()
    for atom in chain_intermediate.graph:
        atom.hybridise()
    chain_intermediate_copy = chain_intermediate.deepcopy()
    # Find OH groups in polyketide/NRP, perform cyclization for each -OH group
    chain_intermediate_copy.refresh_structure()
    o_oh_atoms = find_atoms(O_OH, chain_intermediate_copy)
    o_oh_atoms_filtered = []
    for atom in o_oh_atoms:
        if atom not in o_oh_atoms_filtered and any(neighbour.type == 'H' for neighbour in atom.neighbours):
            o_oh_atoms_filtered.append(atom)

    # Find amino gruops in polyketide/NRP, perform cyclization for each group
    chain_intermediate_copy = chain_intermediate.deepcopy()
    chain_intermediate_copy.refresh_structure()
    chain_intermediate.set_connectivities()
    chain_intermediate.set_atom_neighbours()
    amino_n_atoms_filtered = []
    for atom in chain_intermediate_copy.graph:
        if atom.type == 'N':
            n_neighbour_types = []
            for neighbour in atom.neighbours:
                n_neighbour_types.append(neighbour.type)
            if atom not in amino_n_atoms_filtered and n_neighbour_types.count('H') == 2 and n_neighbour_types.count('C') == 1:
                amino_n_atoms_filtered.append(atom)

    # Define -OH group that should not be used to carry out the thioesterase
    # reaction (distance -S and internal -OH group)
    o_not_to_use = find_o_betapropriolactone(chain_intermediate)
    return o_oh_atoms_filtered + amino_n_atoms_filtered


def find_atoms_for_tailoring(chain_intermediate, atom_type):
    """Atoms that can be tailored

     chain_intermediate: PIKAChU Structure object of a polyketide/NRP
     atom_type: type of atom to search for (e.g. C)
    """
    # Perform first thioesterase reaction, generating linear polyketide/NRP
    chain_intermediate.refresh_structure()
    for atom in chain_intermediate.graph:
        atom.hybridise()
    chain_intermediate_copy = chain_intermediate.deepcopy()
    chain_intermediate_copy.refresh_structure()
    atoms_filtered = []
    for atom in chain_intermediate_copy.graph:
        if atom.type == atom_type:
            neighbour_types = []
            for neighbour in atom.neighbours:
                neighbour_types.append(neighbour.type)
            if atom not in atoms_filtered and neighbour_types.count('H') >= 1:
                atoms_filtered.append(atom)

    return atoms_filtered


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

    hydroxylated_structure.make_bond(
        oxygen, target_atom, hydroxylated_structure.find_next_bond_nr())
    hydroxylated_structure.make_bond(
        hydrogen_1, hydrogen_2, hydroxylated_structure.find_next_bond_nr())

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

    methylated_structure.make_bond(
        carbon, target_atom, methylated_structure.find_next_bond_nr())
    methylated_structure.make_bond(
        hydrogen_1, hydrogen_2, methylated_structure.find_next_bond_nr())

    structures = methylated_structure.split_disconnected_structures()

    for s in structures:
        if carbon.nr in s.atoms:
            return s
