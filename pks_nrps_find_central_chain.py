from pks_elongation_reactions import *
from pikachu.reactions.functional_groups import find_atoms, GroupDefiner
from pikachu.smiles.smiles import *
from central_atoms_pk_starter import find_central_atoms_pk_starter

def find_central_chain_pks_nrps(pks_nrps_attached):
    """Identifies the central chain atoms in the input structure from the
    in_central_chain Atom attribute, and returns these Atom objects as a list

    pks_nrps_attached: PIKAChU Structure object, input (hybrid) PK/NRP structure
    """
    # If the structure is a PK starter unit, find the central carbon chain
    if not any(hasattr(atom, 'in_central_chain') for atom in pks_nrps_attached.graph) and\
            len(pks_nrps_attached.find_substructures(Smiles('C(=O)S').smiles_to_structure())) > 0:
        pks_nrps_attached = find_central_atoms_pk_starter(pks_nrps_attached)

    # Find atoms in the structure inside a cycle
    pks_nrps_attached.find_cycles()
    for atom in pks_nrps_attached.graph:
        if atom.in_ring(pks_nrps_attached):
            atom.inside_ring = True
        else:
            atom.inside_ring = False

    # Identify starting point central chain attached NRP/PK
    for atom in pks_nrps_attached.graph:
        if atom.type == 'S' and any(hasattr(neighbour, 'domain_type')
                                    for neighbour in atom.neighbours):
            sulphur = atom
        if not hasattr(atom, 'in_central_chain'):
            atom.in_central_chain = False


    central_chain = [sulphur]
    visited = [sulphur]
    atom_central_chain = sulphur
    end_atom = False

    # Identify complete central chain from in_central_chain Atom attributes
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
    print(central_chain)

    return central_chain

