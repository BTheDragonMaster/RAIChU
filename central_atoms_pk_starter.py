from pikachu.reactions.functional_groups import GroupDefiner, find_atoms
from pikachu.general import read_smiles
from attributes import ATTRIBUTES
from attach_to_domain import attach_to_domain_pk

POLYKETIDE_S = GroupDefiner('Sulphur atom polyketide', 'SC(C)=O', 0)



def find_central_atoms_pk_starter(pk_starter_unit):
    """Finds the the atoms in the central chain of the polyketide starter unit,
    sets the in_central_chain Atom object attribute to True/False accordingly,
    after which the pk_starter_unit Structure object is returned.

    pk_starter_unit: PIKAChU Structure object of the PK starter unit
    """
    pk_starter_unit.find_cycles()
    central_chain = []
    visited = []

    # Check if structure is a polyketide, define sulphur attached to domain
    locations_sulphur = find_atoms(POLYKETIDE_S, pk_starter_unit)
    assert len(locations_sulphur) == 1
    sulphur_polyketide = locations_sulphur[0]
    central_chain.append(sulphur_polyketide)

    # Iterate over atom neighbours, check if the atom belongs to the central
    # carbon chain (i.e., not a methyl/ethyl/methoxy sidechain) and add to list
    for neighbour in sulphur_polyketide.neighbours:
        if neighbour.type == 'C':
            chain_carbon = neighbour
            end_carbon = False
            central_chain += [chain_carbon]
            visited.append(chain_carbon)
            while not end_carbon:
                for next_atom in chain_carbon.neighbours:
                    ethyl_branch = False
                    inside_cycle = False
                    methyl_group = False
                    next_atom_neighbour_types = []

                    # Build lists of neighbouring atom types
                    if next_atom.type == 'C' and next_atom not in visited:
                        c_neighbours = []
                        for neighbor in next_atom.neighbours:
                            next_atom_neighbour_types.append(neighbor.type)
                            if neighbor.type == 'C':
                                c_neighbours.append(neighbor)

                        # Confirm carbon doesn't belong to ethyl sidebranch
                        types = []
                        if len(c_neighbours) == 2 and 'O' not in \
                                next_atom_neighbour_types:
                            for atom in c_neighbours:
                                neighbours_atom = atom.neighbours
                                for neighbour in neighbours_atom:
                                    types.append(neighbour.type)
                            if types.count('H') == 4 and types.count('C') == 4:
                                ethyl_branch = True
                            for c_atom in c_neighbours:
                                if c_atom not in visited and c_atom != next_atom:
                                    ethyl_branch = False
                            visited.append(chain_carbon)

                        # Carbon in a terminal carboxylic acid group is the
                        # final carbon in the central chain
                        if next_atom_neighbour_types.count('O') == 2:
                            central_chain.append(next_atom)
                            end_carbon = True
                            visited.append(chain_carbon)

                        # Confirm carbon is not part of cycle
                        if len(c_neighbours) == 2:
                            if all(atom.in_ring(pk_starter_unit)
                                   for atom in c_neighbours):
                                visited.append(next_atom)
                                inside_cycle = True

                        # Case if carbon is part of a benzene ring
                        if len (c_neighbours) == 3 and \
                                len(next_atom_neighbour_types) == 3:
                            # central_chain.append(next_atom) #added 10/04/2021
                            visited.append(next_atom)
                            inside_cycle = True
                            end_carbon = True

                        # Case where the PK starter ends in two methyl branches
                        if next_atom_neighbour_types.count(
                                'C') == 3 and \
                                next_atom_neighbour_types.count('H') == 1:
                            nr_methyl_brances = 0
                            for c_neighbour in c_neighbours:
                                c_neighbour_neighbour_types = []
                                for c_neighbour_neighbour \
                                        in c_neighbour.neighbours:
                                    c_neighbour_neighbour_types.append(
                                        c_neighbour_neighbour.type)
                                if c_neighbour_neighbour_types.count('H') == 3:
                                    methyl_carbon = c_neighbour
                                    nr_methyl_brances += 1
                            if nr_methyl_brances == 2:
                                methyl_group = True
                                visited.append(next_atom)
                                central_chain.append(next_atom)
                                visited.append(methyl_carbon)
                                central_chain.append(methyl_carbon)
                                end_carbon = True

                        # Confirm carbon doesn't belong to methyl sidebranch
                        if next_atom_neighbour_types.count('H') == 3 or (
                                next_atom_neighbour_types.count(
                                        'H') == 2 and next_atom_neighbour_types.count(
                                '*') == 1):
                            methyl_group = True

                            # Terminal methylgroup is not a sidebranch:
                            count_c_neighbours_not_visited = 0
                            for neighbouring_atom in chain_carbon.neighbours:
                                if neighbouring_atom.type == 'C' and neighbouring_atom not in visited:
                                    count_c_neighbours_not_visited += 1
                            if count_c_neighbours_not_visited < 2:
                                methyl_group = False
                                visited.append(chain_carbon)
                                for neighbour in chain_carbon.neighbours:
                                    if neighbour.type == 'C' and neighbour not in visited:
                                        next_neighbours = []
                                        for next_neighbour in neighbour.neighbours:
                                            next_neighbours.append(
                                                next_neighbour.type)
                                        if next_neighbours.count('H') == 3 or (
                                                next_neighbours.count(
                                                        'H') == 2 and (
                                                next_neighbours.count(
                                                        '*')) == 1):
                                            central_chain.append(neighbour)
                                end_carbon = True

                            visited.append(chain_carbon)

                        # Identification of central chain carbon atoms through
                        # exclusion
                        if not ethyl_branch and not methyl_group and not \
                                end_carbon and not inside_cycle:
                            central_chain += [next_atom]
                            chain_carbon = next_atom
                            visited.append(chain_carbon)

    # Set in_central_chain Atom attribute
    for atom in pk_starter_unit.graph:
        if atom in central_chain:
            atom.annotations.in_central_chain = True
        else:
            atom.annotations.in_central_chain = False

    return pk_starter_unit

if __name__ == "__main__":
    starter_units_antismash = ['SC(=O)CC', 'SC(CC(O)=O)=O', 'SC(CC(O)=O)=O',
                               'SC(C(C(O)=O)CC)=O', 'SC(C(C(O)=O)OC)=O',
                               'SC(C*)=O',
                               'SC(C(C)CC)=O', 'SC(C1C(CCC1)C(=O)O)=O',
                               'SC(C)=O',
                               'SC(C1=CC=CC=C1)=O', 'SC(CC(C)C)=O',
                               'SC(C(C(=O)O)CC[Cl])=O']
    starter_unit_error = 'SC(C1=CC=CC=C1)=O'
    struct = read_smiles(starter_unit_error)
    struct.add_attributes(ATTRIBUTES, boolean=True)
    struct = attach_to_domain_pk(struct, 'ACP')
    cc = find_central_atoms_pk_starter(struct)
    cc.print_graph()
    for atom in cc.graph:
        print(atom, atom.annotations.in_central_chain)


