from pks_elongation_reactions import *
from pikachu.reactions.functional_groups import find_atoms, GroupDefiner

POLYKETIDE_S = GroupDefiner('Sulphur atom polyketide', 'SC(C)=O', 0)

def find_central_chain(polyketide):
    """
    Returns a list of the atoms that make up the central carbon chain in the
    polyketide, from end carbon to start carbon (the atom attached to the
    sulphur atom)

    polyketide: PIKAChU Structure object of a polyketide
    """
    central_chain = []
    visited = []

    #Check if structure is a polyketide, define sulphur attached to domain
    locations_sulphur = find_atoms(POLYKETIDE_S, polyketide)
    assert len(locations_sulphur) == 1
    sulphur_polyketide = locations_sulphur[0]

    #Iterate over atom neighbours, check if the atom belongs to the central
    #carbon chain (i.e., not a methyl/ethyl/methoxy sidechain) and add to list
    for neighbour in sulphur_polyketide.neighbours:
        if neighbour.type == 'C':
            chain_carbon = neighbour
            end_carbon = False
            central_chain += [chain_carbon]
            visited.append(chain_carbon)
            while not end_carbon:
                for next_atom in chain_carbon.neighbours:
                    ethyl_branch = False
                    methyl_group = False
                    next_atom_neighbour_types = []

                    #Build lists of neighbouring atom types
                    if next_atom.type == 'C' and next_atom not in visited:
                        c_neighbours = []
                        for neighbor in next_atom.neighbours:
                            next_atom_neighbour_types.append(neighbor.type)
                            if neighbor.type == 'C':
                                c_neighbours.append(neighbor)

                        #Confirm carbon doesn't belong to ethyl sidebranch
                        types = []
                        if len(c_neighbours) == 2 and 'O' not in \
                            next_atom_neighbour_types:
                            for atom in c_neighbours:
                                neighbours_atom = atom.neighbours
                                for neighbour in neighbours_atom:
                                    types.append(neighbour.type)
                            if types.count('H') == 4 and types.count('C') == 4:
                                ethyl_branch = True

                        #Carbon of terminal carboxylic acid group is the final
                        #carbon in the central chain
                        if next_atom_neighbour_types.count('O') == 2:
                            central_chain.append(next_atom)
                            end_carbon = True

                        #Confirm carbon doesn't belong to methyl sidebranch
                        if next_atom_neighbour_types.count('H') == 3:
                            methyl_group = True

                        #Identification of central chain carbon atoms through
                        #exclusion
                        if not ethyl_branch and not methyl_group and not \
                            end_carbon:
                            central_chain += [next_atom]
                            chain_carbon = next_atom
                            visited.append(chain_carbon)

    return list(reversed(central_chain))







if __name__ == "__main__":
    struct = Smiles('C(C(C(C(C(CC(=O)O)=O)CC)O)C)(S)=O').smiles_to_structure()
    struct_no_ethyl = Smiles('C(C(C(C(C(CC(=O)O)=O)C)O)C)(S)=O').smiles_to_structure()
    struct2 = Smiles('C(C(=CCC(CC(CC(=O)O)=O)O)C)(=O)S').smiles_to_structure()
    #print(struct.graph)
    #Drawer(struct)
    print(struct2.graph)

