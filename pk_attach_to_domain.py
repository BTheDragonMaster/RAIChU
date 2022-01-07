from pks_elongation_reactions import *
from class_domain import *



def attach_to_domain(structure, domain_type):
    """

    """
    #Create domain
    next_atom_nr = structure.find_next_atom_nr()
    domain = make_domain(domain_type, next_atom_nr)
    domain.add_shell_layout()

    #Remove H atom from S in polyketide, to allow attachment to domain
    #NEED TO CHANGE: USE FIND_SUBSTRUCTURE() INSTEAD OF TAKING A RANDOM S ATOM
    for atom in structure.graph:
        if atom.type == 'S':
            s_atom = atom
            for neighbour in s_atom.neighbours:
                if neighbour.type == 'H':
                    h_to_remove = neighbour
                    break
    sh_bond = structure.bond_lookup[s_atom][h_to_remove]
    structure.break_bond(sh_bond)
    split = structure.split_disconnected_structures()
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

