from typing import Tuple

from pikachu.reactions.functional_groups import find_atoms, find_bonds, combine_structures, BondDefiner, GroupDefiner
from pikachu.reactions.basic_reactions import internal_condensation, condensation
from pikachu.general import read_smiles
from raichu.data.attributes import ATTRIBUTES
from raichu.data.molecular_moieties import CO_BOND, N_AMINO, O_OH, O_BETAPROPRIOLACTONE, SH_BOND, C_CH, C_TERMINAL_OH, N_AMINO_ACID, B_N_AMINO_ACID
from raichu.central_chain_detection.label_central_chain import label_nrp_central_chain
from raichu.reactions.general import label_rest_groups, initialise_atom_attributes, reset_nrp_annotations

def ribosomal_elongation(amino_acid, chain_intermediate, amino_acid_number=None):
    """
    Returns the ribosomal condensation product as a PIKAChU Structure object
    
    amino_acid_number: Type and number of amino acid in peptide (e.g. A14)
    amino_acid: PIKAChU Structure object of the amino acid
    chain_intermediate: PIKAChU Structure object of the ripp intermediate
    """

    # Initialize annotations in all atoms in amino acid
    initialise_atom_attributes(amino_acid)
    for atom in amino_acid.atoms.values():
        atom.amino_acid = amino_acid_number
    amino_acid.atoms
    oh_bond = find_bonds(C_TERMINAL_OH, chain_intermediate)
    assert len(oh_bond)==1
    oh_bond = oh_bond[0]

    # Label rest groups
    label_rest_groups(amino_acid, chain_intermediate)


    # Find amino group of amino acid or beta-amino acid
    n_atoms_aa = find_atoms(N_AMINO_ACID, amino_acid)
    if not n_atoms_aa:
        n_atoms_aa = find_atoms(B_N_AMINO_ACID, amino_acid)

    assert len(n_atoms_aa) == 1

    n_atom = n_atoms_aa[0]
    h_bond = None

    for bond in n_atom.bonds:
        for neighbour in bond.neighbours:
            if neighbour.type == 'H':
                h_bond = bond
                break

    assert h_bond

    # Label atoms in amino acid that end up in central peptide chain
    label_nrp_central_chain(amino_acid)

    # Reset chiral_c_ep and n_atom_nmeth AtomAnnotations attribute for all atoms in the intermediate;
    # This way, only the newly incorporated amino acid is labelled for tailoring
    reset_nrp_annotations(chain_intermediate)

    # Carry out condensation reaction using build-in PIKAChU function
    condensation_product = condensation(chain_intermediate, amino_acid, oh_bond, h_bond)[0]

    # Refresh condensation product
    condensation_product.refresh_structure(find_cycles=True)

    # Initialize annotations for atoms that don't have any yet
    initialise_atom_attributes(condensation_product)

    return condensation_product


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
    if cyclisation_site_2.type == ('C'):
        h_atom = cyclisation_site_2.get_neighbour('H')
        assert c_atom
        c_bond = cyclisation_site_2.get_bond(h_atom)
        assert c_bond
    else:
        c_atom = cyclisation_site_2.get_neighbour('C')
        assert c_atom
        c_bond = cyclisation_site_2.get_bond(c_atom)
        assert c_bond

    cyclic_product, water = internal_condensation(structure, c_bond, h_bond)

    return cyclic_product


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
