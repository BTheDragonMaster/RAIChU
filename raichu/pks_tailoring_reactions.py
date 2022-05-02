from pikachu.reactions.functional_groups import BondDefiner
from pikachu.chem.chirality import same_chirality

from raichu.attach_to_domain import *

RECENT_ELONGATION = BondDefiner('recent_elongation', 'O=CCC(=O)S', 0, 1)
RECENT_REDUCTION_COH = BondDefiner('recent_reduction_C-OH', 'OCCC(=O)S', 0, 1)
RECENT_REDUCTION_MMAL_CHIRAL_C = GroupDefiner('recent_reduction_mmal_chiral_c', 'CCC(C(=O)S)C', 2)
RECENT_REDUCTION_C = GroupDefiner('recent_reduction_mal', 'OCCC(=O)S', 2)
RECENT_REDUCTION_CC = BondDefiner('recent_reduction_C-C', 'OCCC(=O)S', 1, 2)
RECENT_DEHYDRATION = BondDefiner('recent_dehydration', 'SC(C=CC)=O', 2, 3)
KR_DOMAIN_TYPES = ['A1', 'A2', 'B1', 'B2', 'C1', 'C2']
S_KR = GroupDefiner('C1 atom before KR reaction', 'SC(C)=O', 0)


def carbonyl_to_hydroxyl(double_bond):
    """Alters the double bond in a carbonyl group to a single bond and returns
    this as a PIKAChU Bond object

    double_bond: PIKAChU Bond object, bond between the C and O atom of a
    carbonyl group
    """
    # Assert the correct bond is parsed
    assert double_bond.type == 'double'
    atom_types_in_bond = []
    for atom in double_bond.neighbours:
        atom_types_in_bond.append(atom.type)
    assert 'O' in atom_types_in_bond and 'C' in atom_types_in_bond

    # Define the two electrons in the pi bond (inside the db) that needs
    # to be broken
    electron_1 = None
    electron_2 = None
    electrons_in_db = double_bond.electrons
    for electron in electrons_in_db:
        if electron.atom.type == 'O' and electron.orbital_type == 'p':
            electron_1 = electron
        elif electron.atom.type == 'C' and electron.orbital_type == 'p':
            electron_2 = electron

    assert electron_1 and electron_2

    # Remove the pi electrons from their respective orbital
    orbital_1 = electron_1.orbital
    orbital_2 = electron_2.orbital
    orbital_1.remove_electron(electron_2)
    orbital_2.remove_electron(electron_1)

    # Set bond type to single
    for bond in double_bond.atom_1.bonds:
        if bond == double_bond:
            bond.type = 'single'
    for bond in double_bond.atom_2.bonds:
        if bond == double_bond:
            bond.type = 'single'

    # Remove pi electrons from the Bond between the C and O atom
    for electron in double_bond.electrons[:]:
        if electron == electron_1 and electron.orbital_type == 'p':
            double_bond.electrons.remove(electron)
        elif electron == electron_2 and electron.orbital_type == 'p':
            double_bond.electrons.remove(electron)

    # Change hybridisation of both C and O atoms to sp3
    atom_1, atom_2 = double_bond.neighbours
    atom_1.valence_shell.dehybridise()
    atom_1.valence_shell.hybridise('sp3')
    atom_2.valence_shell.dehybridise()
    atom_2.valence_shell.hybridise('sp3')

    # Change bond_type of Bond instance to single
    double_bond.type = 'single'

    new_single_bond = double_bond
    new_single_bond.set_bond_summary()

    return new_single_bond


def find_betaketon(structure):
    """
    Returns the bond between the carbonyl oxygen and -carbon of the beta ketone
    group in the input structure as a PIKAChU Bond object

    structure: PIKAChU Structure object of the input structure
    """
    locations = structure.find_substructures(RECENT_ELONGATION.structure)
    bonds = []
    for match in locations:
        atom_1 = match.atoms[RECENT_ELONGATION.atom_1]
        atom_2 = match.atoms[RECENT_ELONGATION.atom_2]
        bond = structure.bond_lookup[atom_1][atom_2]
        bonds.append(bond)
    return bonds


def ketoreductase(chain_intermediate, kr_type=None):
    """
    Performs the ketoreductase reaction on the PKS chain intermediate, returns
    the reaction product as a PIKAChU Structure object

    chain_intermediate: PIKAChU Structure object, PKS chain intermediate just
    after a PKS elongation step
    """

    if kr_type:
        assert kr_type in KR_DOMAIN_TYPES

    # Reset all colours to black:
    for atom in chain_intermediate.graph:
        atom.draw.colour = 'black'
    for bond_nr, bond in chain_intermediate.bonds.items():
        bond.set_bond_summary()

    chiral_c = None
    chiral_c_locations = find_atoms(RECENT_REDUCTION_MMAL_CHIRAL_C, chain_intermediate)

    new_single_bond = None

    # Identify beta-ketone bond, identify O- and C-atom participating in bond
    if kr_type in {'A1', 'A2', 'B1', 'B2', 'C2', None}:

        beta_ketone_bonds = find_betaketon(chain_intermediate)
        assert len(beta_ketone_bonds) == 1
        beta_ketone_bond = beta_ketone_bonds[0]

        carbonyl_oxygen = beta_ketone_bond.get_neighbour('O')
        carbonyl_carbon = beta_ketone_bond.get_neighbour('C')
        if not carbonyl_carbon or not carbonyl_oxygen:
            raise Exception('Cannot find atoms in beta ketone bond')

    # Identify beta-ketone bond, identify O- and C-atom participating in bond
        if kr_type in {'A1', 'A2', 'B1', 'B2', None}:

            # Change carbonyl bond to single bond
            new_single_bond = carbonyl_to_hydroxyl(beta_ketone_bond)

            # Add H atom to form hydroxyl group and another H to the C
            chain_intermediate.add_atom('H', [carbonyl_oxygen])
            chain_intermediate.add_atom('H', [carbonyl_carbon])
            for atom in chain_intermediate.graph:
                if not hasattr(atom.annotations, 'in_central_chain'):
                    for attribute in ATTRIBUTES:
                        atom.annotations.add_annotation(attribute, False)

            chain_intermediate.refresh_structure()
            carbonyl_carbon = chain_intermediate.get_atom(carbonyl_carbon)

            if kr_type:
                h_atom = carbonyl_carbon.get_neighbour('H')
                o_atom = carbonyl_carbon.get_neighbour('O')
                c_atoms = carbonyl_carbon.get_neighbours('C')
                bottom_c = None

                top_cs = find_atoms(RECENT_REDUCTION_C, chain_intermediate)
                assert len(top_cs) == 1
                top_c = top_cs[0]
                assert top_c in c_atoms

                for atom in c_atoms:
                    if atom != top_c:
                        bottom_c = atom

                assert bottom_c

                counterclockwise_order = [h_atom, o_atom, top_c, bottom_c]

                if kr_type.startswith('A'):

                    if same_chirality(counterclockwise_order, carbonyl_carbon.neighbours):
                        carbonyl_carbon.chiral = 'clockwise'
                    else:
                        carbonyl_carbon.chiral = 'counterclockwise'

                elif kr_type.startswith('B'):
                    if same_chirality(counterclockwise_order, carbonyl_carbon.neighbours):
                        carbonyl_carbon.chiral = 'counterclockwise'
                    else:
                        carbonyl_carbon.chiral = 'clockwise'
                else:
                    raise ValueError('This type of KR domain is not supported by RAIChU or does not exist')

            # If the type of KR domain is not known, the chirality of the carbonyl carbon is set to None
            else:
                carbonyl_carbon.chiral = None

            # Set bond summary for newly formed bond (cannot do from struct.bonds?)

            for atom in chain_intermediate.graph:
                if atom == carbonyl_carbon:
                    for bond in atom.bonds:
                        for neighbour in bond.neighbours:
                            if neighbour.type == 'O':
                                bond.set_bond_summary()

    # See if the previous elongation step was performed using methylmalonyl-CoA,
    # perform epimerization if required

    if chiral_c_locations:
        chiral_c = chiral_c_locations[0]

    if chiral_c and kr_type:
        if kr_type.endswith('1'):
            chiral_c.chiral = 'clockwise'
        elif kr_type.endswith('2'):
            chiral_c.chiral = 'counterclockwise'
        else:
            raise ValueError('This type of KR domain is not supported\
            by RAIChU or does not exist')

    # Default chirality methyl centre is known
    if chiral_c and not kr_type:
        chiral_c.chiral = 'clockwise'

    # Fix order atom neighbours for chiral methyl centre
    if chiral_c:

        sulphur_locations = find_atoms(S_KR, chain_intermediate)
        assert len(sulphur_locations) == 1
        sulphur = sulphur_locations[0]
        assert sulphur

        c_1 = None
        for neighbour in sulphur.neighbours:
            if neighbour.annotations.in_central_chain and neighbour.type == 'C':
                c_1 = neighbour
        assert c_1

        c_2 = None
        for neighbour in c_1.neighbours:
            if neighbour.annotations.in_central_chain and neighbour != sulphur:
                c_2 = neighbour
        assert c_2
        assert c_2 == chiral_c

        c_3 = None
        for neighbour in c_2.neighbours:
            if neighbour.annotations.in_central_chain and neighbour != c_1:
                c_3 = neighbour
        assert c_3

        hydrogen = None
        first_sidechain_atom = None
        for atom in chiral_c.neighbours:
            if atom != c_1 and atom != c_2 and atom.type != 'H' and \
                    not atom.annotations.in_central_chain:
                first_sidechain_atom = atom
            elif atom.type == 'H':
                hydrogen = atom
        assert hydrogen
        assert first_sidechain_atom

        new_neighbours = [first_sidechain_atom, c_1, c_3, hydrogen]
        chain_intermediate.graph[chiral_c] = new_neighbours

    # Refresh structure :)
    chain_intermediate.set_atom_neighbours()
    chain_intermediate.find_cycles()
    for atom in chain_intermediate.graph:
        atom.set_connectivity()
    for bond_nr, bond in chain_intermediate.bonds.items():
        bond.set_bond_summary()

    # Add colouring to the tailored group
    if not kr_type or kr_type.startswith('A') or kr_type.startswith('B'):
        for atom in new_single_bond.neighbours:
            atom.draw.colour = 'red'

    return chain_intermediate


def dehydratase(chain_intermediate):
    """
    Performs the dehydratase reaction on the PKS chain intermediate, returns
    the reaction product as a PIKAChU Structure object

    chain_intermediate: PIKAChU Structure object, PKS chain intermediate where
    the beta ketone group has been recently reduced by the KR domain
    """
    # Reset atom colours to black
    for atom in chain_intermediate.graph:
        atom.draw.colour = 'black'
    for bond_nr, bond in chain_intermediate.bonds.items():
        bond.set_bond_summary()

    # Find and define the atoms that participate in the bond changes
    co_bond = find_oh_dh(chain_intermediate)

    c1 = None
    c2 = None
    o_oh = None

    for bond in co_bond[:]:
        for neighbour in bond.neighbours:
            if neighbour.type == 'C':
                c2 = neighbour
            elif neighbour.type == 'O':
                o_oh = neighbour
    cc_bond = find_cc_dh(chain_intermediate)

    for bond in cc_bond[:]:
        bond.type = 'double'
        for neighbour in bond.neighbours:
            if neighbour == c2:
                pass
            else:
                c1 = neighbour

    # Remove hydroxyl group from c2
    for bond in co_bond[:]:
        chain_intermediate.break_bond(bond)
    # split = chain_intermediate.split_disconnected_structures()
    # chain_intermediate, oh = split

    # Remove H-atom from c1
    bond_to_break = None
    hydrogen = None
    for atom in chain_intermediate.graph:
        if atom == c1:
            for neighbour in atom.neighbours:
                if neighbour.type == 'H':
                    hydrogen = neighbour
                    for bond in neighbour.bonds:
                        bond_to_break = bond
                    break

    assert bond_to_break and hydrogen

    chain_intermediate.break_bond(bond_to_break)
    # split = chain_intermediate.split_disconnected_structures()
    # for structure in split:
    #     if len(structure.graph) == 1:
    #         h_atom = structure
    #     else:
    #         chain_intermediate = structure

    # Patch. When the bond is broken, the electron is not removed from
    # the orbitals of the C atom.
    for atom in chain_intermediate.graph:
        if atom == c1:
            c1 = atom

    # Make double c1=c2 bond
    form_double_bond(c1, c2, chain_intermediate.bond_lookup[c1][c2])

    # Make bond between removed hydroxyl group and hydrogen atom to form water
    next_bond_nr = chain_intermediate.find_next_bond_nr()
    chain_intermediate.make_bond(hydrogen, o_oh, next_bond_nr)

    # Remove separate water structure from the Structure object
    one, two = chain_intermediate.split_disconnected_structures()
    if len(one.graph) == 3:
        chain_intermediate = two
    else:
        chain_intermediate = one

    # After double bond formation, remove chirality c1 and c2
    for atom in chain_intermediate.graph:
        if atom == c1:
            atom.chiral = None
        elif atom == c2:
            atom.chiral = None

    # Add colouring
    for atom in chain_intermediate.graph:
        if atom == c1:
            for bond in atom.bonds:
                for neighbour in bond.neighbours:
                    if neighbour == c2:
                        for atom1 in bond.neighbours:
                            atom1.draw.colour = 'blue'

    for atom in chain_intermediate.graph:
        if atom == c2:
            for bond in atom.bonds:
                for neighbour in bond.neighbours:
                    if neighbour == c1:
                        for atom1 in bond.neighbours:
                            atom1.draw.colour = 'blue'

    # Refresh structure :)
    for atom in chain_intermediate.graph:
        atom.get_connectivity()
    chain_intermediate.set_atom_neighbours()
    chain_intermediate.set_connectivities()
    chain_intermediate.find_cycles()

    for bond_nr, bond in chain_intermediate.bonds.items():
        bond.set_bond_summary()

    return chain_intermediate


def form_double_bond(atom1, atom2, bond):
    """Forms a double bond between two input carbon atoms which are already
    single bonded

    atom1: Atom object of a carbon atom bonded with atom2 through a single bond
    atom2: Atom object of a carbon atom bonded with atom1 through a single bond
    bond: Bond object of the bond between atom1 and atom2
    """
    # Select the bonding electron in the first atom valence shell
    electron_1 = None
    electron_2 = None
    orbital_1 = None
    orbital_2 = None

    for orbital in atom2.valence_shell.orbitals:
        if len(orbital.electrons) == 1:
            orbital_2 = orbital
            electron_2 = orbital.electrons[0]
    for orbital in atom1.valence_shell.orbitals:
        if len(orbital.electrons) == 1:
            orbital_1 = orbital
            electron_1 = orbital.electrons[0]

    assert electron_1 and electron_2 and orbital_1 and orbital_2

    # Add the lone electron of atom 1 to the orbital of atom 2 and vice versa
    orbital_1.add_electron(electron_2)
    orbital_2.add_electron(electron_1)

    orbital_1.set_bond(bond, 'pi')
    orbital_2.set_bond(bond, 'pi')

    # Change hybridisation
    atom1.valence_shell.dehybridise()
    atom1.valence_shell.hybridise('sp2')
    atom2.valence_shell.dehybridise()
    atom2.valence_shell.hybridise('sp2')

    p_orbital_1 = None
    electrons_in_p = []
    p_orbital_2 = None
    electrons_in_p_2 = []
    electrons_in_sp2 = []
    electrons_in_sp2_2 = []

    # Make sure that the electrons of the new bond are in type p orbitals
    # instead of sp2
    if not orbital_1.orbital_type == 'p':
        electrons_in_sp2 = orbital_1.electrons
        for orbital in atom1.valence_shell.orbitals:
            if orbital.orbital_type == 'p':
                p_orbital_1 = orbital
                electrons_in_p = orbital.electrons
    if not orbital_2.orbital_type == 'p':
        electrons_in_sp2_2 = orbital_1.electrons
        for orbital in atom2.valence_shell.orbitals:
            if orbital.orbital_type == 'p':
                p_orbital_2 = orbital
                electrons_in_p_2 = orbital.electrons

    # Swap electrons in sp and p orbitals if not correct yet
    if not orbital_1.orbital_type == 'p':
        orbital_1.electrons = electrons_in_p
        p_orbital_1.electrons = electrons_in_sp2
    if not orbital_2.orbital_type == 'p':
        orbital_2.electrons = electrons_in_p_2
        p_orbital_2.electrons = electrons_in_sp2_2

    # Set right orbital name for all electrons in all orbitals
    for orbital in atom1.valence_shell.orbitals:
        for electron in orbital.electrons:
            if atom1.valence_shell.atom == electron.atom:
                electron.set_orbital(orbital)
    for orbital in atom2.valence_shell.orbitals:
        for electron in orbital.electrons:
            if atom2.valence_shell.atom == electron.atom:
                electron.set_orbital(orbital)

    newly_double_bond = None

    # Add the new electrons of the pi-bond to the Bond between the two atoms
    for bond in atom1.bonds:
        for neighbour in bond.neighbours:
            if neighbour == atom2:
                newly_double_bond = bond

    assert newly_double_bond

    for orbital in atom1.valence_shell.orbitals:
        if orbital.orbital_type == 'p':
            for electron in orbital.electrons:
                if electron not in newly_double_bond.electrons:
                    newly_double_bond.electrons.append(electron)
    for orbital in atom2.valence_shell.orbitals:
        if orbital.orbital_type == 'p':
            for electron in orbital.electrons:
                if electron not in newly_double_bond.electrons:
                    newly_double_bond.electrons.append(electron)
                    newly_double_bond.set_bond_summary()


def find_oh_dh(structure):
    """
    Returns the Bond between the C and OH group to carry out the DH reaction

    structure: PIKAChU Structure object
    """
    locations = structure.find_substructures(RECENT_REDUCTION_COH.structure)
    structure.refresh_structure()
    structure.get_connectivities()
    structure.set_atom_neighbours()
    bonds = []
    for match in locations:
        atom_1 = match.atoms[RECENT_REDUCTION_COH.atom_1]
        atom_2 = match.atoms[RECENT_REDUCTION_COH.atom_2]
        bond = structure.bond_lookup[atom_1][atom_2]
        bonds.append(bond)

    return bonds


def find_cc_dh(structure):
    """
    Returns the Bond between the two C's that needs to be doubled

    structure: PIKAChU Structure object
    """
    locations = structure.find_substructures(RECENT_REDUCTION_CC.structure)
    bonds = []
    for match in locations:
        atom_1 = match.atoms[RECENT_REDUCTION_CC.atom_1]
        atom_2 = match.atoms[RECENT_REDUCTION_CC.atom_2]
        bond = structure.bond_lookup[atom_1][atom_2]
        bonds.append(bond)

    return bonds


def enoylreductase(chain_intermediate):
    """
    Performs the enoylreductase reaction on the PKS chain intermediate, returns
    the reaction product as a PIKAChU Structure object

    chain_intermediate: PIKAChU Structure object, PKS chain intermediate where
    the beta ketone group has been recently reduced and dehydrated by the KR
    and ER domains, respectively
    """
    # Reset all colours to black:
    for atom in chain_intermediate.graph:
        atom.draw.colour = 'black'
    chain_intermediate.set_atom_neighbours()
    for atom in chain_intermediate.graph:
        atom.set_connectivity()

    # Find double bond, change to single bond
    double_cc_bond = find_double_cc(chain_intermediate)
    atoms_in_double_bond = []
    double_cc_bond[0].make_single()

    atom_1 = None
    atom_2 = None

    for bond in double_cc_bond:
        atom_1 = bond.neighbours[0]
        atom_2 = bond.neighbours[1]
        atoms_in_double_bond.append(atom_1)
        atoms_in_double_bond.append(atom_2)

    assert atom_1 and atom_2

    # Add H-atom to the C atoms participating in the new single bond
    for neighbour in double_cc_bond[0].neighbours:
        chain_intermediate.add_atom('H', [neighbour])

    # Change hybridisation from sp2 to sp3
    atom_1.valence_shell.dehybridise()
    atom_1.valence_shell.hybridise('sp3')
    atom_2.valence_shell.dehybridise()
    atom_2.valence_shell.hybridise('sp3')

    # Give annotation to added H-atom
    for atom in chain_intermediate.graph:
        if not hasattr(atom.annotations, 'in_central_chain'):
            for attribute in ATTRIBUTES:
                atom.annotations.add_annotation(attribute, False)

    # Return chirality C atom in the case that the previous elongation reaction
    # was performed using methylmalonyl-CoA

    for atom in atoms_in_double_bond:
        carbon_neighbours = 0
        for neighbour in atom.neighbours:
            if neighbour.type == 'C':
                carbon_neighbours += 1
        if carbon_neighbours == 3:
            atom.chiral = None

    # add colouring to highlight tailoring reaction
    for atom in chain_intermediate.graph:
        if atom == atom_1:
            for bond in atom.bonds:
                for neighbour in bond.neighbours:
                    if neighbour == atom_2:
                        for coloured_atom in bond.neighbours:
                            coloured_atom.draw.colour = 'LIME'

    for atom in chain_intermediate.graph:
        if atom == atom_2:
            for bond in atom.bonds:
                for neighbour in bond.neighbours:
                    if neighbour == atom_1:
                        for coloured_atom in bond.neighbours:
                            coloured_atom.draw.colour = 'LIME'

    return chain_intermediate


def find_double_cc(structure):
    """
    Returns the double Bond between the two C's that needs to be reduced

    structure: PIKAChU Structure object
    """
    locations = structure.find_substructures(RECENT_DEHYDRATION.structure)
    bonds = []
    for match in locations:
        atom_1 = match.atoms[RECENT_DEHYDRATION.atom_1]
        atom_2 = match.atoms[RECENT_DEHYDRATION.atom_2]
        bond = structure.bond_lookup[atom_1][atom_2]
        bonds.append(bond)

    return bonds
