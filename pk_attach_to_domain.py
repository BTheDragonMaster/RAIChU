from class_domain import *
from pikachu.reactions.functional_groups import find_atoms, GroupDefiner

POLYKETIDE_S = GroupDefiner('Sulphur atom polyketide', 'SC(C)=O', 0)

def attach_to_domain(polyketide, domain_type):
    """
    Attaches the sulphur atom in the input polyketide to a PKS domain

    domain_type: Str, domain type
    """
    #Create domain
    next_atom_nr = polyketide.find_next_atom_nr()
    domain = make_domain(domain_type, next_atom_nr)
    domain.add_shell_layout()

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
        proton = one
        structure = two
    else:
        proton = two
        structure = one

    #Make new bond between S in polyketide and domain
    for i, neighbour in enumerate([s_atom]):
        next_bond_nr = structure.find_next_bond_nr()
        structure.make_bond(domain, neighbour, next_bond_nr)

    structure.set_connectivities()

    return structure

