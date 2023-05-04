
from pikachu.smiles.smiles import *
from pikachu.general import read_smiles
from raichu.central_chain_detection.label_central_chain import label_pk_central_chain


def reorder_central_chain(central_chain, drawer):
    stop_linearising = None
    rings = []
    current_ring = []
    current_ring_index = None

    for i, atom in enumerate(central_chain):
        if atom.inside_ring:
            ring_index = atom.get_ring_index(drawer.structure)
            if (current_ring and ring_index == current_ring_index) or not current_ring:
                current_ring.append(atom)
                current_ring_index = ring_index
            else:
                rings.append(current_ring)
                current_ring = [atom]
                current_ring_index = ring_index

    rings.append(current_ring)
    new_rings = []

    for ring in rings:

        if len(ring) > 3:
            alternative_path = drawer.find_shortest_path(ring[0], ring[-1], path_type='atom')
            if len(alternative_path) <= 3:
                new_backbone = []
                path_added = False
                for atom in central_chain:
                    if atom not in ring:
                        new_backbone.append(atom)
                    elif not path_added:
                        new_backbone += alternative_path
                        path_added = True
                central_chain = new_backbone
                new_rings.append(alternative_path)
            else:
                stop_linearising = ring[2]
                break
        else:
            new_rings.append(ring)

    rings = new_rings
    atom_to_ring = {}
    for ring in rings:
        for atom in ring:
            atom_to_ring[atom] = ring

    return central_chain, rings, atom_to_ring, stop_linearising


def find_central_chain(pks_nrps_attached):
    """Identifies the central chain atoms in the input structure from the
    in_central_chain Atom attribute, and returns these Atom objects as a list

    pks_nrps_attached: PIKAChU Structure object, input (hybrid) PK/NRP structure
    """
    # If the structure is a PK starter unit, find the central carbon chain

    if not any(atom.annotations.in_central_chain for atom in pks_nrps_attached.graph) and\
            len(pks_nrps_attached.find_substructures(read_smiles('C(=O)S'))) > 0:
        pks_nrps_attached = label_pk_central_chain(pks_nrps_attached)

    # Find atoms in the structure inside a cycle
    pks_nrps_attached.find_cycles()
    for atom in pks_nrps_attached.graph:
        if atom.in_ring(pks_nrps_attached):
            atom.inside_ring = True
        else:
            atom.inside_ring = False

    # Identify starting point central chain attached NRP/PK
    sulphur = None
    for atom in pks_nrps_attached.graph:
        if atom.type == 'S' and any(neighbour.annotations.domain_type
                                    for neighbour in atom.neighbours):
            sulphur = atom
        if not atom.annotations.in_central_chain:
            atom.annotations.in_central_chain = False
    assert sulphur

    central_chain = [sulphur]
    visited = [sulphur]
    atom_central_chain = sulphur
    end_atom = False

    # Identify complete central chain from in_central_chain Atom attributes
    while not end_atom:
        for neighbour in atom_central_chain.neighbours:
            if neighbour.annotations.in_central_chain and neighbour not in visited:
                central_chain.append(neighbour)
                visited.append(neighbour)
                atom_central_chain = neighbour
            elif not neighbour.annotations.in_central_chain:
                neighbours = []
                for next_atom in atom_central_chain.neighbours:
                    neighbours.append(next_atom)
                if all(atom in visited for atom in neighbours):
                    end_atom = True
                else:
                    visited.append(neighbour)

    return central_chain


def find_central_chain_ripp(ripp_attached):
    """Identifies the central chain atoms in the input structure from the
    in_central_chain Atom attribute, and returns these Atom objects as a list

    ripp_attached: PIKAChU Structure object, input (hybrid) ripp structure
    """
    # If the structure is a PK starter unit, find the central carbon chain

    if not any(atom.annotations.in_central_chain for atom in ripp_attached.graph) and\
            len(ripp_attached.find_substructures(read_smiles('C(=O)S'))) > 0:
        ripp_attached = label_pk_central_chain(ripp_attached)

    # Find atoms in the structure inside a cycle
    ripp_attached.find_cycles()
    for atom in ripp_attached.graph:
        if atom.in_ring(ripp_attached):
            atom.inside_ring = True
        else:
            atom.inside_ring = False

    # Identify starting point central chain
    nitrogen = None
    domains = [
        atom.annotations.domain_type for atom in ripp_attached.graph if atom.annotations.domain_type]
    if "Follower" in domains:
        for atom in ripp_attached.graph:
            if atom.annotations.domain_type:
                if atom.annotations.domain_type == "Follower":
                    nitrogen = atom.get_neighbour("N")

            if not atom.annotations.in_central_chain:
                atom.annotations.in_central_chain = False

    else:
        for atom in ripp_attached.graph:
            if atom.type == 'N' and any(neighbour.annotations.domain_type
                                        for neighbour in atom.neighbours):
                nitrogen = atom

            if not atom.annotations.in_central_chain:
                atom.annotations.in_central_chain = False

    assert nitrogen

    central_chain = [nitrogen]
    visited = [nitrogen]
    atom_central_chain = nitrogen
    end_atom = False

    # Identify complete central chain from in_central_chain Atom attributes
    while not end_atom:
        for neighbour in atom_central_chain.neighbours:
            if neighbour.annotations.in_central_chain and neighbour not in visited:
                central_chain.append(neighbour)
                visited.append(neighbour)
                atom_central_chain = neighbour
            elif not neighbour.annotations.in_central_chain:
                neighbours = []
                for next_atom in atom_central_chain.neighbours:
                    neighbours.append(next_atom)
                if all(atom in visited for atom in neighbours):
                    end_atom = True
                else:
                    visited.append(neighbour)

    return central_chain


def find_central_chain_not_attached(pks_nrps):
    """Identifies the central chain atoms in the input structure from the
    in_central_chain Atom attribute, and returns these Atom objects as a list

    pks_nrps: PIKAChU Structure object, input (hybrid) PK/NRP/RiPP structure
    """
    # Find atoms in the structure inside a cycle
    pks_nrps.find_cycles()
    for atom in pks_nrps.graph:
        if atom.in_ring(pks_nrps):
            atom.inside_ring = True
        else:
            atom.inside_ring = False

    # Identify starting point central chain attached NRP/PK
    carbon = None
    for atom in pks_nrps.graph:
        if atom.type == 'C' and atom.annotations.in_central_chain and [neighbour.type for neighbour in atom.neighbours].count("O") == 2:
            carbon = atom
        if not atom.annotations.in_central_chain:
            atom.annotations.in_central_chain = False
    assert carbon

    central_chain = [carbon]
    visited = [carbon]
    atom_central_chain = carbon
    end_atom = False

    # Identify complete central chain from in_central_chain Atom attributes
    while not end_atom:
        for neighbour in atom_central_chain.neighbours:
            if neighbour.annotations.in_central_chain and neighbour not in visited:
                central_chain.append(neighbour)
                visited.append(neighbour)
                atom_central_chain = neighbour
            elif not neighbour.annotations.in_central_chain:
                neighbours = []
                for next_atom in atom_central_chain.neighbours:
                    neighbours.append(next_atom)
                if all(atom in visited for atom in neighbours):
                    end_atom = True
                else:
                    visited.append(neighbour)

    return central_chain
