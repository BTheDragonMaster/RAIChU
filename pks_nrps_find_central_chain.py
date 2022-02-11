from pks_elongation_reactions import *
from pikachu.reactions.functional_groups import find_atoms, GroupDefiner

def find_central_chain_pks_nrps(pks_nrps_attached):
    """

    """
    # find atoms in the structure inside a cycle
    for atom in pks_nrps_attached.graph:
        if atom.in_ring(pks_nrps_attached):
            atom.inside_ring = True
            print(atom)
        else:
            atom.inside_ring = False

    #Identify starting point central chain attached NRP/PK
    for atom in pks_nrps_attached.graph:
        if atom.type == 'S':
            sulphur = atom
        if not hasattr(atom, 'in_central_chain'):
            atom.in_central_chain = False

    central_chain = [sulphur]
    visited = [sulphur]
    atom_central_chain = sulphur
    end_atom = False

    #Identify complete central chain from in_central_chain Atom attributes
    while not end_atom:
        for neighbour in atom_central_chain.neighbours:
            if neighbour.in_central_chain and neighbour not in visited:
                central_chain.append(neighbour)
                visited.append(neighbour)
                atom_central_chain = neighbour
            elif not neighbour.in_central_chain:
                neighbours = []
                for next_atom in atom_central_chain.neighbours:
                    neighbours.append(next_atom)
                if all(atom in visited for atom in neighbours):
                    end_atom = True
                else:
                    visited.append(neighbour)

    return central_chain