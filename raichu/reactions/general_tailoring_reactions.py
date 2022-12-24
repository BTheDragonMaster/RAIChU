from typing import Tuple

from pikachu.reactions.functional_groups import find_atoms, find_bonds, combine_structures, BondDefiner, GroupDefiner
from pikachu.reactions.basic_reactions import internal_condensation, condensation
from pikachu.general import read_smiles
from raichu.data.attributes import ATTRIBUTES
from raichu.data.molecular_moieties import CO_BOND, N_AMINO, O_OH, O_BETAPROPRIOLACTONE, SH_BOND, C_CH, C_TERMINAL_OH, N_AMINO_ACID, B_N_AMINO_ACID
from raichu.central_chain_detection.label_central_chain import label_nrp_central_chain
from raichu.reactions.general import label_rest_groups, initialise_atom_attributes, reset_nrp_annotations



def double_bond_reduction(atom1, atom2, structure):
    """
    Returns the reduced product as a PIKAChU Structure object
    
    atom1: C atom1 in double bond
    atom2: C atom2 in double bond
    structure: PIKAChU Structure object of the ripp intermediate
    """
    assert atom1.type == 'C'
    assert atom2.type == 'C'
    cc_bond = atom1.get_bond(atom2)
    assert cc_bond
    if cc_bond.type=="double" :
        cc_bond.make_single()
        for neighbour in cc_bond.neighbours:
            structure.add_atom('H', [neighbour])
    return structure

def epoxidation(atom1, atom2, structure):
    """
    Returns the epoxidated product as a PIKAChU Structure object
    
    atom1: C atom1 to be epoxidated
    atom2: C atom2 to be epoxidated, needs to be neighbour of atom1
    structure: PIKAChU Structure object of the ripp intermediate
    """
    assert atom1.type == 'C'
    assert atom2.type == 'C'
    cc_bond = atom1.get_bond(atom2)
    assert cc_bond
    if cc_bond.type=="double" or cc_bond.type=="aromatic":
        cc_bond.make_single()
        for neighbour in cc_bond.neighbours:
            structure.add_atom('H', [neighbour])
    structure = hydroxylation(atom1, structure)
    initialise_atom_attributes(structure)
    structure.refresh_structure()
    # find O atom
    oxygen = atom1.get_neighbour('O')
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
    if atom2.type != 'O':
        atom1, atom2 = atom2, atom1
    cyclisation_site_1 = structure.get_atom(atom1)
    h_atom = cyclisation_site_1.get_neighbour('H')
    assert h_atom
    h_bond = cyclisation_site_1.get_bond(h_atom)
    assert h_bond
    cyclisation_site_2 = structure.get_atom(atom2)
    c_atom = cyclisation_site_2.get_neighbour('C')
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
    h_atom = cyclisation_site_1.get_neighbour('H')
    assert h_atom
    h_bond = cyclisation_site_1.get_bond(h_atom)
    assert h_bond
    cyclisation_site_2 = structure.get_atom(atom2)
    h_atom_2 = cyclisation_site_2.get_neighbour('H')
    assert h_atom_2
    h_bond_2 = cyclisation_site_2.get_bond(h_atom_2)
    assert h_bond_2

    structure.break_bond(h_bond)
    structure.break_bond(h_bond_2)

    # Create the bonds

    structure.make_bond(cyclisation_site_1, cyclisation_site_2, structure.find_next_bond_nr())
    structure.make_bond(h_atom, h_atom_2, structure.find_next_bond_nr())

    # Put the h2 and the product into different Structure instances

    structures = structure.split_disconnected_structures()

    h2 = None
    product = None

    # Find out which of the structures is your product and which is your h2

    for structure in structures:
        if h_atom in structure.graph:
            h2 = structure
        elif cyclisation_site_1 in structure.graph:
            structure.refresh_structure(find_cycles=True)
            product = structure

    return product


def proteolytic_cleavage(bond, structure, structure_to_keep: str = "follower"):
    """Performs proteolytic cleavage

     bond: exact PIKAChU bond object to cleave
     structure: PIKAChU Structure object to perform cleavage on
     structure_to_keep: determines if the leading or the following peptide should be kept ("leader" or "follower")
    """
    carbon = bond.get_neighbour('C')
    nitrogen = bond.get_neighbour('N')
    structure.break_bond(bond)
    oxygen = structure.add_atom('O', [carbon])
    initialise_atom_attributes(structure)
    structure.add_atom('H', [oxygen])
    structure.add_atom('H', [nitrogen])

    initialise_atom_attributes(structure)

    structures = structure.split_disconnected_structures()
    for structure in structures:
        if structure_to_keep == "leader":
            if carbon in structure.graph:
                structure.refresh_structure()
                return structure
        if structure_to_keep == "follower":
            if nitrogen in structure.graph:
                structure.refresh_structure()
                return structure


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
