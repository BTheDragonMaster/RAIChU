from modules_to_structure import *
from pikachu.reactions.functional_groups import find_atoms, GroupDefiner
from pikachu.general import draw_structure

SH_BOND = BondDefiner('recent_elongation', 'SC(C)=O', 0, 1)
CO_BOND = BondDefiner('recent_elongation', 'CO', 0, 1)
N_AMINO = GroupDefiner('N_amino', 'CN', 1)
O_OH = GroupDefiner('O_oh', 'CO', 1)
O_BETAPROPRIOLACTONE = GroupDefiner('o_betapropriolactone', 'SC(CCO)=O', 4)



def thioesterase_linear_product(chain_intermediate):
    """Carries out the thioesterase reaction on the polyketide chain
    intermediate, returning the free acid linear polyketide product as a
    PIKAChU Structure object

    chain_intermediate: PIKAChU Structure object of the polyketide chain
    intermediate, either attached to a PKS domain or not
    """
    # Find S-H bond in chain intermediate and break the bond and define atoms
    sh_bonds = find_cs_bond(chain_intermediate)
    for sh_bond in sh_bonds:
        for atom in sh_bond.neighbours:
            if atom.type == 'S':
                sulphur_sh = atom
            elif atom.type == 'C':
                carbon_sh = atom
        chain_intermediate.break_bond(sh_bond)

    # Remove the sulphur and atoms/domains attached from the Structure
    first_part, second_part = chain_intermediate.split_disconnected_structures()
    for atom in chain_intermediate.graph:
        if atom == carbon_sh:
            atom = carbon_sh
    if 2 <= len(first_part.graph) <= 3:
        sh = first_part
        chain_intermediate = second_part
    else:
        sh = second_part
        chain_intermediate = first_part

    # Add a negatively charged oxygen atom to the carbon atom
    chain_intermediate.add_atom('O', [carbon_sh])
    for atom in chain_intermediate.graph[carbon_sh]:
        if atom.type == 'O' and \
        chain_intermediate.bond_lookup[atom][carbon_sh].type == 'single':
            new_oxygen = atom
            chain_intermediate.add_atom('H', [new_oxygen])
    # Make sure this new oxygen atom has all required atom annotations
    for atom in chain_intermediate.graph:
        if not hasattr(atom.annotations, 'in_central_chain'):
            for attribute in ATTRIBUTES:
                atom.annotations.add_annotation(attribute, False)

    # Refresh structure
    chain_intermediate.find_cycles()
    chain_intermediate.refresh_structure()

    return chain_intermediate


def thioesterase_circular_product(chain_intermediate, o_oh_n_amino):
    """Performs the thioesterase reactions on the input chain_intermediate
     using the -OH group defined by the input O-atom participating in that
     internal -OH group, returns the circular product as PIKAChU Structure
     object.

     chain_intermediate: PIKAChU Structure object of a polyketide
     o_oh_n_amino: PIKAChU Atom object of the O-atom in the -OH group or N-atom
     in the amino group that the function should use to perform the
     thioesterase reaction.
    """
    # Find S-H bond in chain intermediate and break the bond
    sh_bonds = find_cs_bond(chain_intermediate)
    for sh_bond in sh_bonds:
        for atom in sh_bond.neighbours:
            if atom.type == 'S':
                sulphur_sh = atom
            elif atom.type == 'C':
                carbon_sh = atom
        chain_intermediate.break_bond(sh_bond)

    # Remove the sulphur and atoms/domains attached from the Structure
    first_part, second_part = chain_intermediate.split_disconnected_structures()
    for atom in chain_intermediate.graph:
        if atom == carbon_sh:
            atom = carbon_sh
    if 2 <= len(first_part.graph) <= 3:
        sh = first_part
        chain_intermediate = second_part
    else:
        sh = second_part
        chain_intermediate = first_part

    # Refresh structure and set orbitals for all electrons in the structure
    # that got lost in the deepcopy
    chain_intermediate.refresh_structure()
    for atom in chain_intermediate.graph:
        for orbital in atom.valence_shell.orbitals:
            for electron in orbital.electrons:
                electron.set_orbital(orbital)
    chain_intermediate.find_cycles()

    # Remove H from internal -OH group and refresh structure
    for atom in chain_intermediate.graph:
        if atom == o_oh_n_amino:
            for neighbour in atom.neighbours:
                if neighbour.type == 'H':
                    h_to_remove = neighbour
    chain_intermediate.remove_atom(h_to_remove)
    chain_intermediate.refresh_structure()
    chain_intermediate.set_connectivities()
    for bond_nr, bond in chain_intermediate.bonds.items():
        bond.set_bond_summary()

    # Form bond to create circular product
    next_bond_nr = chain_intermediate.find_next_bond_nr()
    chain_intermediate.add_bond(carbon_sh, o_oh_n_amino, 'single', next_bond_nr)

    return chain_intermediate


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



def thioesterase_all_products(chain_intermediate):
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
    RaichuDrawer(thioesterase_linear_product(chain_intermediate_copy))

    # Find OH groups in polyketide/NRP, perform cyclization for each -OH group
    chain_intermediate_copy = chain_intermediate.deepcopy()
    chain_intermediate_copy.refresh_structure()
    chain_intermediate.set_connectivities()
    chain_intermediate.set_atom_neighbours()
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
    list_product_drawings = []

    # Perform all possible thioesterase reactions leading to the formation of
    # circular products using the internal amino groups, save Structure objects
    # to list
    for n_amino in amino_n_atoms_filtered:
        chain_intermediate.refresh_structure()
        chain_intermediate_copy = chain_intermediate.deepcopy()
        for atom in chain_intermediate_copy.graph:
            if atom == n_amino:
                atom = n_amino
        product = thioesterase_circular_product(chain_intermediate_copy, n_amino)
        list_product_drawings.append(product)

    # Perform all possible thioesterase reactions leading to the formation of
    # circular products using the internal -OH groups, save Structure objects
    # to list
    for o_oh in o_oh_atoms_filtered:
        chain_intermediate.refresh_structure()
        chain_intermediate_copy = chain_intermediate.deepcopy()
        for atom in chain_intermediate_copy.graph:
            if atom == o_oh:
                atom = o_oh
        if o_oh != o_not_to_use:
            product = thioesterase_circular_product(chain_intermediate_copy, o_oh)
            list_product_drawings.append(product)

    # Draw all products
    for product in list_product_drawings:
        RaichuDrawer(product)

    return list_product_drawings


def find_cs_bond(structure):
    """
    Returns the bond between the sulphur and carbon atom of the polyketide
    chain intermediate

    structure: PIKAChU Structure object
    """
    locations = structure.find_substructures(SH_BOND.structure)
    bonds = []
    for match in locations:
        atom_1 = match.atoms[RECENT_ELONGATION.atom_1]
        atom_2 = match.atoms[RECENT_ELONGATION.atom_2]
        bond = structure.bond_lookup[atom_1][atom_2]
        bonds.append(bond)
    return bonds






