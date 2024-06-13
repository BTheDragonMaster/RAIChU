from pikachu.reactions.functional_groups import (
    find_atoms,
    find_bonds,
    combine_structures,
    GroupDefiner,
)
from pikachu.reactions.basic_reactions import internal_condensation
from pikachu.general import read_smiles
from raichu.data.attributes import ATTRIBUTES
from raichu.data.molecular_moieties import PYROPHOSPHATE_BOND
from raichu.reactions.general import initialise_atom_attributes


def reductive_bond_breakage(atom1, atom2, structure):
    """
    Returns the product with a broken bond as a PIKAChU Structure object. 
    If the bond breakage leads to disconnected structures, 
    the structure containing the first atom will be kept.

    structure: PIKAChU Structure object of the intermediate
    atom1: atom on side of the bond to be kept
    atom2. atom on the other side of the bond
    """
    structure = structure.kekulise()
    atom1 = structure.get_atom(atom1)
    atom2 = structure.get_atom(atom2)
    bond = atom1.get_bond(atom2)
    assert bond
    assert bond.type == "single"

    structure.break_bond(bond)
    structure.add_atom("H", [atom1])
    structure.add_atom("H", [atom2])
    initialise_atom_attributes(structure)
    structures = structure.split_disconnected_structures()
    for structure in structures:
        if atom1 in structure.graph:
            structure.refresh_structure(find_cycles=True)
            return structure


def excise_from_structure(atom1, atom2, structure):
    """
    Returns the product without the group between atoms as a PIKAChU Structure object

    structure: PIKAChU Structure object of the intermediate
    atom1: last atom in structure on one side
    atom2. last atom in structure on other side
    """
    # find exact bonds  to the atoms between:

    traveling_atom = atom1
    target = atom2
    visited = []
    dead_ends = []
    atom_before = None
    while traveling_atom != target:
        if all(atom in visited for atom in traveling_atom.neighbours):
            dead_ends.append(traveling_atom)

            traveling_options = [
                atom for atom in traveling_atom.neighbours if atom not in dead_ends
            ]
            if len(traveling_options) > 0:
                traveling_atom = traveling_options[0]
            else:
                raise ValueError(f"Atoms are not connected")
        else:
            neighbour = [
                atom for atom in traveling_atom.neighbours if atom not in visited
            ][0]
            visited.append(neighbour)
            atom_before = traveling_atom
            traveling_atom = neighbour
    atom_before_atom2 = atom_before
    bond_to_break_2 = atom2.get_bond(atom_before_atom2)

    visited = []
    dead_ends = []
    atom_before = None
    traveling_atom = atom2
    target = atom1
    while traveling_atom != target:
        if all(atom in visited for atom in traveling_atom.neighbours):
            dead_ends.append(traveling_atom)

            traveling_options = [
                atom for atom in traveling_atom.neighbours if atom not in dead_ends
            ]
            if len(traveling_options) > 0:
                traveling_atom = traveling_options[0]
            else:
                raise ValueError(f"Atoms are not connected")
        else:
            neighbour = [
                atom for atom in traveling_atom.neighbours if atom not in visited
            ][0]
            visited.append(neighbour)
            atom_before = traveling_atom
            traveling_atom = neighbour
    atom_before_atom1 = atom_before
    bond_to_break_1 = atom1.get_bond(atom_before_atom1)
    structure.break_bond(bond_to_break_1)
    structure.break_bond(bond_to_break_2)
    structure.make_bond(atom1, atom2, structure.find_next_bond_nr())
    structures = structure.split_disconnected_structures()
    for structure in structures:
        if atom1 in structure.graph:
            return structure


def change_chirality(carbon, structure):
    carbon = structure.get_atom(carbon)
    if carbon.chiral == "clockwise":
        carbon.chiral = "counterclockwise"
    elif carbon.chiral == "counterclockwise":
        carbon.chiral = "clockwise"
    initialise_atom_attributes(structure)
    structure.refresh_structure(find_cycles=True)
    return structure


def remove_atom(atom, structure):
    """
    Returns the product without the group attached to the atom as a PIKAChU Structure object

    structure: PIKAChU Structure object of the intermediate
    atom: atom to be removed
    """
    atom = structure.get_atom(atom)
    carbon = atom.get_neighbour("C")
    bond = atom.get_bond(carbon)
    assert bond
    structure.break_bond(bond)
    structure_1, structure_2 = structure.split_disconnected_structures()
    if carbon in structure_1.graph:
        structure = structure_1
    else:
        structure = structure_2
    structure.add_atom("H", [carbon])
    initialise_atom_attributes(structure)
    structure.refresh_structure(find_cycles=True)
    return structure


def remove_group_at_bond(bond, structure):
    """
    Returns the product without the grup attached to the bond as a PIKAChU Structure object

    structure: PIKAChU Structure object of the intermediate
    bond: PIKAChU bond to be broken
    """
    carbon = bond.get_neighbour("C")
    structure.break_bond(bond)
    structure_1, structure_2 = structure.split_disconnected_structures()
    if carbon in structure_1.graph:
        structure = structure_1
    else:
        structure = structure_2
    structure.add_atom("H", [carbon])
    initialise_atom_attributes(structure)
    structure.refresh_structure(find_cycles=True)
    return structure


def dephosphorylation(structure):
    """
    Returns the dephosphorylated product as a PIKAChU Structure object

    structure: PIKAChU Structure object of the intermediate
    """
    pyrophosphate_bonds = find_bonds(PYROPHOSPHATE_BOND, structure)
    if len(pyrophosphate_bonds) > 0:
        pyrophosphate_bond = pyrophosphate_bonds[0]
        carbon = pyrophosphate_bond.get_neighbour("C")
        structure.break_bond(pyrophosphate_bond)
        structure_1, structure_2 = structure.split_disconnected_structures()
        if carbon in structure_1.graph:
            structure = structure_1
        else:
            structure = structure_2
        structure.add_atom("H", [carbon])
        initialise_atom_attributes(structure)
    structure.refresh_structure(find_cycles=True)
    return structure


def double_bond_shift(
    structure,
    old_double_bond_atom1,
    old_double_bond_atom2,
    new_double_bond_atom1,
    new_double_bond_atom2,
):
    """
    Returns the reduced product as a PIKAChU Structure object

    old_double_bond_atom1: C atom1 in old double bond
    old_double_bond_atom2: C atom2 in old double bond
    new_double_bond_atom1: C atom1 in new double bond
    new_double_bond_atom2: C atom2 in new double bond
    structure: PIKAChU Structure object of the intermediate
    """
    new_double_bond = structure.bond_lookup[new_double_bond_atom1][
        new_double_bond_atom2
    ]
    assert new_double_bond
    structure = double_bond_reduction(
        old_double_bond_atom1, old_double_bond_atom2, structure
    )
    structure = single_bond_oxidation(
        new_double_bond_atom1, new_double_bond_atom2, structure
    )

    return structure


def single_bond_oxidation(atom1, atom2, structure):
    """
    Returns the oxidized product as a PIKAChU Structure object

    atom1: C atom1 in single bond
    atom2: C atom2 in single bond
    structure: PIKAChU Structure object of the intermediate
    """

    # remove h atoms
    for selected_atom in [atom1, atom2]:
        bond_to_break = None
        hydrogen = None

        for atom in structure.graph:
            if atom == selected_atom:
                for neighbour in atom.neighbours:
                    if neighbour.type == "H":
                        hydrogen = neighbour
                        for bond in neighbour.bonds:
                            bond_to_break = bond
                        break

        assert bond_to_break and hydrogen
        structure.break_bond(bond_to_break)
    atom1 = structure.get_atom(atom1)
    atom2 = structure.get_atom(atom2)
    single_bond = structure.bond_lookup[atom1][atom2]
    single_bond.make_double()
    structures = structure.split_disconnected_structures()
    for structure in structures:
        if atom1 in structure.graph:
            initialise_atom_attributes(structure)
            structure.refresh_structure(find_cycles=True)
            return structure


def double_bond_reduction(atom1, atom2, structure):
    """
    Returns the reduced product as a PIKAChU Structure object

    atom1: C atom1 in double bond
    atom2: C atom2 in double bond
    structure: PIKAChU Structure object of the intermediate
    """
    double_bond = atom1.get_bond(atom2)
    assert double_bond
    if double_bond.type == "aromatic":
        structure = structure.kekulise()
        atom1 = structure.get_atom(atom1)
        double_bond = atom1.get_bond(atom2)
        double_bond.make_single()
        double_bond.set_bond_summary()
        for neighbour in double_bond.neighbours:
            structure.add_atom("H", [neighbour])
        structure.refresh_structure(find_cycles=True)
        structure.aromatic_cycles = structure.find_aromatic_cycles()
        structure.aromatic_systems = structure.find_aromatic_systems()
    if double_bond.type == "double":
        double_bond.make_single()
        double_bond.set_bond_summary()
        for neighbour in double_bond.neighbours:
            structure.add_atom("H", [neighbour])
    initialise_atom_attributes(structure)
    structure.refresh_structure(find_cycles=True)
    return structure


def epoxidation(atom1, atom2, structure):
    """
    Returns the epoxidated product as a PIKAChU Structure object

    atom1: C atom1 to be epoxidated
    atom2: C atom2 to be epoxidated, needs to be neighbour of atom1
    structure: PIKAChU Structure object of the intermediate
    """
    assert atom1.type == "C"
    assert atom2.type == "C"
    double_bond_reduction(atom1, atom2, structure)
    structure = addition(atom1, "O", structure)
    initialise_atom_attributes(structure)
    structure.refresh_structure()
    # find O atom
    oxygen = atom1.get_neighbour("O")
    structure = oxidative_bond_formation(oxygen, atom2, structure)
    initialise_atom_attributes(structure)
    structure.refresh_structure()
    return structure


def cyclisation(structure, atom1, atom2):
    """Performs cyclisation

    atom1: PIKAChU atom to be used in cyclisation that will remain in final product
    atom2: PIKAChU atom to be used in cyclisation that will be removed from final product -> needs to be oxygen
    structure: PIKAChU Structure object to perform cyclization on
    """
    if atom2.type != "O":
        atom1, atom2 = atom2, atom1
    cyclisation_site_1 = structure.get_atom(atom1)
    h_atom = cyclisation_site_1.get_neighbour("H")
    assert h_atom
    h_bond = cyclisation_site_1.get_bond(h_atom)
    assert h_bond
    cyclisation_site_2 = structure.get_atom(atom2)
    c_atom = cyclisation_site_2.get_neighbour("C")
    assert c_atom
    c_bond = cyclisation_site_2.get_bond(c_atom)
    assert c_bond

    cyclic_product, water = internal_condensation(structure, c_bond, h_bond)

    return cyclic_product


def oxidative_bond_formation(atom1, atom2, structure):
    """Performs cyclisation

    atom1: PIKAChU atom to be used in cyclisation
    atom2: PIKAChU atom to be used in cyclisation
    structure: PIKAChU Structure object to perform cyclization on
    """

    cyclisation_site_1 = structure.get_atom(atom1)
    h_atom = cyclisation_site_1.get_neighbour("H")
    assert h_atom
    h_bond = cyclisation_site_1.get_bond(h_atom)
    assert h_bond
    cyclisation_site_2 = structure.get_atom(atom2)
    h_atom_2 = cyclisation_site_2.get_neighbour("H")
    assert h_atom_2
    h_bond_2 = cyclisation_site_2.get_bond(h_atom_2)
    assert h_bond_2
    structure.break_bond(h_bond)
    structure.break_bond(h_bond_2)

    # Create the bonds
    structure.make_bond(
        cyclisation_site_1, cyclisation_site_2, structure.find_next_bond_nr()
    )
    structure.make_bond(h_atom, h_atom_2, structure.find_next_bond_nr())
    # remove from chirality dict

    cyclisation_site_1 = structure.get_atom(cyclisation_site_1)
    cyclisation_site_2 = structure.get_atom(cyclisation_site_2)

    # update chiral_dict -remove with new pikachu update
    bonds_atom_1 = cyclisation_site_1.bonds[:]
    bonds_atom_2 = cyclisation_site_2.bonds[:]

    for bond in bonds_atom_1:
        if h_atom in bond.chiral_dict:
            bond.chiral_dict[cyclisation_site_2] = bond.chiral_dict[h_atom]
            del bond.chiral_dict[h_atom]
        for atom, atoms_and_chirality in bond.chiral_dict.items():
            if h_atom in atoms_and_chirality:
                atoms_and_chirality[cyclisation_site_2] = atoms_and_chirality[h_atom]
                del atoms_and_chirality[h_atom]

    for bond in bonds_atom_2:
        if h_atom_2 in bond.chiral_dict:
            bond.chiral_dict[cyclisation_site_1] = bond.chiral_dict[h_atom_2]
            del bond.chiral_dict[h_atom_2]
        for atom, atoms_and_chirality in bond.chiral_dict.items():
            if h_atom_2 in atoms_and_chirality:
                atoms_and_chirality[cyclisation_site_1] = atoms_and_chirality[h_atom_2]
                del atoms_and_chirality[h_atom_2]

    # Put the h2 and the product into different Structure instances
    structure.refresh_structure(find_cycles=False)
    structures = structure.split_disconnected_structures()
    product = None
    initialise_atom_attributes(structure)

    # Find out which of the structures is your product and which is your h2

    for structure in structures:
        if cyclisation_site_1 in structure.graph:
            structure.refresh_structure(find_cycles=False)
            initialise_atom_attributes(structure)
            product = structure
    return product


def proteolytic_cleavage(bond, structure, structure_to_keep: str = "end"):
    """Performs proteolytic cleavage

    bond: exact PIKAChU bond object to cleave
    structure: PIKAChU Structure object to perform cleavage on
    structure_to_keep: determines if the leading or the following peptide should be kept ("leader" or "follower")
    """
    carbon = bond.get_neighbour("C")
    assert carbon
    nitrogen = bond.get_neighbour("N")
    assert nitrogen
    structure.break_bond(bond)
    oxygen = structure.add_atom("O", [carbon])
    initialise_atom_attributes(structure)
    structure.add_atom("H", [oxygen])
    structure.add_atom("H", [nitrogen])

    initialise_atom_attributes(structure)

    structures = structure.split_disconnected_structures()
    for structure in structures:
        if structure_to_keep == "start":
            if carbon in structure.graph:
                structure.refresh_structure()
                return structure
        if structure_to_keep == "end":
            if nitrogen in structure.graph:
                structure.refresh_structure()
                return structure


def find_atoms_for_tailoring(chain_intermediate, atom_type):
    """
    Returns atoms that can be tailored

     chain_intermediate: PIKAChU Structure object of a polyketide/NRP
     atom_type: type of atom to search for (e.g. C)
    """

    for atom in chain_intermediate.graph:
        atom.hybridise()
    atoms_filtered = []
    for atom in chain_intermediate.graph:
        if atom.type == atom_type:
            neighbour_types = []
            for neighbour in atom.neighbours:
                neighbour_types.append(neighbour.type)
            if atom not in atoms_filtered and neighbour_types.count("H") >= 1:
                atoms_filtered.append(atom)

    return atoms_filtered


def addition(target_atom, structure_to_add, structure):
    """
    Returns the structure that has a group added at the target atom.

    structure_to_add: SMILES of structure, atom to be added needs to be at first position
    structure: PIKAChU Structure object
    target_atom:  PIKAChU atom object
    """
    FIRST_ATOM = GroupDefiner("first_atom", structure_to_add, 0)
    group_to_add = read_smiles(structure_to_add)
    group_to_add.add_attributes(ATTRIBUTES, boolean=True)
    atom_to_add = find_atoms(FIRST_ATOM, group_to_add)[0]
    target_atom = structure.get_atom(target_atom)
    hydrogen_1 = atom_to_add.get_neighbour("H")
    bond_1 = atom_to_add.get_bond(hydrogen_1)
    hydrogen_2 = target_atom.get_neighbour("H")
    if not hydrogen_2:
        raise Exception("Can't add to this atom!")

    bond_2 = target_atom.get_bond(hydrogen_2)

    combined_structure = combine_structures([group_to_add, structure])
    combined_structure.break_bond(bond_1)
    combined_structure.break_bond(bond_2)

    combined_structure.make_bond(
        atom_to_add, target_atom, combined_structure.find_next_bond_nr()
    )
    combined_structure.make_bond(
        hydrogen_1, hydrogen_2, combined_structure.find_next_bond_nr()
    )
    structures = combined_structure.split_disconnected_structures()
    for structure in structures:
        if atom_to_add.nr in structure.atoms:
            structure.refresh_structure()
            initialise_atom_attributes(structure)
            return structure


def cyclodehydration(structure, atom1, oxygen):
    initialise_atom_attributes(structure)
    carbon = oxygen.get_neighbour("C")
    nitrogen = carbon.get_neighbour("N")
    structure = double_bond_reduction(carbon, oxygen, structure)
    structure = oxidative_bond_formation(carbon, atom1, structure)
    structure = remove_atom(oxygen, structure)
    structure = single_bond_oxidation(carbon, nitrogen, structure)
    structure.refresh_structure()
    initialise_atom_attributes(structure)
    return structure


def macrolactam_formation(structure, o_oh_carboxyl):
    """Performs the macrolactam formation on the input chain_intermediate
    using the -OH group defined by the input O-atom participating in that
    internal -OH group and the N-terminus, returns the circular product as PIKAChU Structure
    object.

    chain_intermediate: PIKAChU Structure object of a polyketide
    o_oh_carboxyl: PIKAChU Atom object of the O-atom in the -OH group that the function should use to perform the
    macrolactame reaction.
    """
    cyclisation_site = structure.get_atom(o_oh_carboxyl)
    c_atom = cyclisation_site.get_neighbour("C")
    assert c_atom
    oh_bond = cyclisation_site.get_bond(c_atom)
    assert oh_bond
    terminal_nitrogen = None
    for atom in structure.graph:
        if (
            #atom.annotations.in_central_chain
            atom.type == "N"
            and [neighbour.type for neighbour in atom.neighbours].count("H") == 2
        ):
            terminal_nitrogen = atom
    assert terminal_nitrogen
    h_atom = terminal_nitrogen.get_neighbour("H")
    h_bond = terminal_nitrogen.get_bond(h_atom)
    cyclic_product, water = internal_condensation(structure, oh_bond, h_bond)
    cyclic_product.refresh_structure()
    initialise_atom_attributes(cyclic_product)
    return cyclic_product
