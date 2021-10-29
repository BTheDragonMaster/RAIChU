from pks_modules_to_structure import *

def find_central_chain(polyketide):
    """
    Need to make it more specific: include some find_substructure searches to find beginning and end of polyketide

    """
    central_chain = []
    visited = []
    for atom in polyketide.graph:
        if atom.type == 'S':
            sulphur = atom
            for neighbour in sulphur.neighbours:
                if neighbour.type == 'C':
                    chain_carbon = neighbour
                    flag = True
                    central_chain += [chain_carbon]
                    while flag:
                        visited.append(chain_carbon)
                        for next_atom in chain_carbon.neighbours:
                            next_atom_neighbour_types = []
                            if next_atom.type == 'C' and next_atom not in visited:
                                for neighbour in next_atom.neighbours:
                                    next_atom_neighbour_types.append(neighbour.type)
                                if next_atom_neighbour_types.count('H') == 3:
                                    continue
                                if next_atom_neighbour_types.count('O') == 2:
                                    central_chain.append(next_atom)
                                    flag = False
                                else:
                                    central_chain += [next_atom]
                                    chain_carbon = next_atom
    return list(reversed(central_chain))






if __name__ == "__main__":
    start_unit = Smiles('SC(CC(=O)O)=O').smiles_to_structure()
    chain1 = add_methylmalonylunit(start_unit)
    chain2 = ketoreductase(chain1)
    chain3 = dehydratase(chain2)
    chain4 = add_malonylunit(chain3)
    Drawer(chain4)
    print(chain4.graph)
    print('central chain:', find_central_chain(chain4))