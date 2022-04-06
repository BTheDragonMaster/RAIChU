from class_domain import *
from pikachu.reactions.functional_groups import find_atoms, GroupDefiner


POLYKETIDE_S = GroupDefiner('Sulphur atom polyketide', 'SC(C)=O', 0)
NRP_C = GroupDefiner('C atom to attach to PCP domain', 'NCC(O)=O', 2)

def attach_to_domain_pk(polyketide, domain_type):
    """
    Attaches the sulphur atom in the input polyketide to a PKS domain and
    returns the attached structure as a PIKAChU Structure object

    domain_type: Str, domain type
    polyketide: PIKAChU Structure object, to-be attached structure
    """
    #Create domain
    next_atom_nr = polyketide.find_next_atom_nr()
    domain = make_domain(domain_type, next_atom_nr)
    domain.add_electron_shells()

    #Remove H atom from S in polyketide, to allow attachment to domain
    locations_sulphur = find_atoms(POLYKETIDE_S, polyketide)
    assert len(locations_sulphur) == 1
    s_atom = locations_sulphur[0]
    for neighbour in s_atom.neighbours:
        if neighbour.type == 'H':
            h_to_remove = neighbour
            break
    sh_bond = polyketide.bond_lookup[s_atom][h_to_remove]
    polyketide.break_bond(sh_bond)
    split = polyketide.split_disconnected_structures()
    one, two = split
    if len(one.graph) == 1:
        h_atom = one
        structure = two
    else:
        h_atom = two
        structure = one

    #Make new bond between S in polyketide and domain
    for i, neighbour in enumerate([s_atom]):
        next_bond_nr = structure.find_next_bond_nr()
        structure.make_bond(domain, neighbour, next_bond_nr)

    structure.set_atoms()

    return structure

def attach_to_domain_nrp(nrp, domain_type):
    """
    Attaches the input NRP to a PCP domain and returns the product as a
    PIKAChU Structure object

    domain_type: Str, domain type
    nrp: PIKAChU Structure object, to-be attached NRP
    """
    # Create domain
    next_atom_nr = nrp.find_next_atom_nr()
    domain = make_domain(domain_type, next_atom_nr)
    domain.add_electron_shells()

    # Remove OH group from carboxylic acid in NRP, to allow attachment
    # to domain
    locations_c_to_domain = find_atoms(NRP_C, nrp)
    assert len(locations_c_to_domain) == 1
    c_atom_to_domain = locations_c_to_domain[0]
    for neighbour in c_atom_to_domain.neighbours:
        if neighbour.type == 'O':
            if nrp.bond_lookup[c_atom_to_domain][neighbour].type == 'single':
                remove_o = neighbour
                bond_to_break = nrp.bond_lookup[c_atom_to_domain][neighbour]


    nrp.break_bond(bond_to_break)
    split = nrp.split_disconnected_structures()
    one, two = split
    if len(one.graph) == 2:
        hydroxyl = one
        structure = two
    else:
        hydroxyl = two
        structure = one

    # Add S atom to C in NRP
    structure.add_atom('S', [c_atom_to_domain])
    structure.set_atom_neighbours()
    for atom in structure.graph:
        if atom.type == 'S':
            for neighbour in atom.neighbours:
                if neighbour == c_atom_to_domain:
                    sulphur_to_pcp = atom

    # Make new bond between S in NRP and domain
    for i, neighbour in enumerate([sulphur_to_pcp]):
        next_bond_nr = structure.find_next_bond_nr()
        structure.make_bond(domain, neighbour, next_bond_nr)

    structure.set_atoms()

    for atom in structure.graph:
        if not hasattr(atom.annotations, 'in_central_chain'):
            for attribute in ATTRIBUTES:
                atom.annotations.add_annotation(attribute, False)

    return structure

