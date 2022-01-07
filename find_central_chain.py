from pks_elongation_reactions import *
from pks_tailoring_reactions import *

def find_central_chain(polyketide):
    """
    Need to make it more specific: include some find_substructure searches to find beginning and end of polyketide

    """
    central_chain = []
    visited = []
    for atom in polyketide.graph:
        if atom.type == 'S':
            sulphur = atom
            visited.append(sulphur)
            for neighbour in sulphur.neighbours:
                if neighbour.type == 'C':
                    chain_carbon = neighbour
                    visited.append(chain_carbon)
                    flag = True
                    central_chain += [chain_carbon]
                    while flag:
                        visited.append(chain_carbon)
                        for next_atom in chain_carbon.neighbours:
                            if next_atom not in visited:
                                flag2 = True
                                next_atom_neighbour_types = []
                                if next_atom.type == 'C' and next_atom not in visited:
                                    for neighbor in next_atom.neighbours:
                                        if neighbor not in visited:
                                            next_atom_neighbour_types.append(neighbor.type)
                                            next_atom_neighbour_neighbour_types = []
                                            for next_atom_neighbour_neighbour in neighbor.neighbours:
                                                next_atom_neighbour_neighbour_types.append(next_atom_neighbour_neighbour.type)
                                            if next_atom_neighbour_neighbour_types.count('H') == 3:
                                                flag2 = False
                                        if next_atom_neighbour_types.count('H') == 1 and next_atom_neighbour_neighbour_types.count('C') == 3:
                                            flag2 = True
                                        if next_atom_neighbour_types.count('H') == 3:
                                            continue
                                        if next_atom_neighbour_types.count('O') == 2:
                                            central_chain.append(next_atom)
                                            flag = False
                                        else:
                                            if flag2:
                                                central_chain += [next_atom]
                                                visited.append(next_atom)
                                                chain_carbon = next_atom
    return list(reversed(central_chain))

def find_central_chain_improved(polyketide):
    """
    Need to make it more specific: include some find_substructure searches to
    find beginning and end of polyketide

    """
    central_chain = []
    visited = []
    for atom in polyketide.graph:
        if atom.type == 'S':
            sulphur = atom
            for neighbour in sulphur.neighbours:
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
                            if next_atom.type == 'C' and next_atom not in visited:
                                #build list of next_atom_neighbour_types
                                c_neighbours = []
                                for neighbor in next_atom.neighbours:
                                    next_atom_neighbour_types.append(neighbor.type)
                                    if neighbor.type == 'C':
                                        c_neighbours.append(neighbor)
                                types = []
                                if len(c_neighbours) == 2 and 'O' not in next_atom_neighbour_types:
                                    for atom in c_neighbours:
                                        neighbours_atom = atom.neighbours
                                        for neighbour in neighbours_atom:
                                            types.append(neighbour.type)
                                    if types.count('H') == 4 and types.count('C') == 4:
                                        ethyl_branch = True
                                if next_atom_neighbour_types.count('O') == 2:
                                    central_chain.append(next_atom)
                                    end_carbon = True
                                if next_atom_neighbour_types.count('H') == 3:
                                    methyl_group = True
                                if not ethyl_branch and not methyl_group and not end_carbon:
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
    print(find_central_chain_improved(struct2), 'CENTRAL CHAHN')
