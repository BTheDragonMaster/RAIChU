from raichu_drawer import *
from pk_attach_to_domain import attach_to_domain_nrp

def find_central_chain_nrp(nrp):
    """
    NEED TO CLEAN THIS FUNCTION UP!!!!!!!!! WITH PROPER VARIABLE NAMES

    """
    visited = []
    found_the_sulphur = False
    for atom in nrp.graph:
        if atom.type == 'S':
            sulphur = atom
            for z in atom.neighbours:
                print(z)
                if hasattr(z, 'domain_type'):
                    pcp_domain = z
                    central_peptide_chain = [pcp_domain]
                    found_the_sulphur = True
            if found_the_sulphur:
                central_peptide_chain.append(sulphur)
                visited.append(sulphur)
                terminal_nitrogen = False
                for sulphur_neighbour in sulphur.neighbours:
                    if sulphur_neighbour.type =='C':
                        first_carbon = sulphur_neighbour
                        visited.append(first_carbon)
                        central_peptide_chain.append(first_carbon)
                        for u in first_carbon.neighbours:
                            if u.type == 'C':
                                second_carbon = u
                                current_atom = second_carbon
                while not terminal_nitrogen:
                    print(central_peptide_chain)
                    current_atom_neighbour_types = []
                    current_atom_neighbours = []
                    for current_atom_neighbour in current_atom.neighbours:
                        current_atom_neighbour_types.append(current_atom_neighbour.type)
                        current_atom_neighbours.append(current_atom_neighbour)
                    if current_atom.type == 'C':
                        if len(current_atom_neighbour_types) == 3 and\
                                current_atom_neighbour_types.count('O') == 1 and\
                                current_atom_neighbour_types.count('N') == 1 and\
                                current_atom_neighbour_types.count('C') == 1:
                            for x in current_atom_neighbours:
                                if x.type == 'N':
                                    if x in visited:
                                        visited.append(current_atom)
                                        central_peptide_chain.append(current_atom)
                                        for b in current_atom_neighbours:
                                            if b.type == 'C':
                                                if b not in visited:
                                                    current_atom = b
                        elif len(current_atom_neighbour_types) == 4 and\
                                current_atom_neighbour_types.count('H') == 1 and\
                                current_atom_neighbour_types.count('N') == 1:
                            for y in current_atom_neighbours:
                                if y.type == 'C':
                                    if y in visited:
                                        visited.append(current_atom)
                                        central_peptide_chain.append(current_atom)
                                        for c in current_atom_neighbours:
                                            if c.type == 'N':
                                                current_atom = c
                        elif len(current_atom_neighbour_types) == 4 and\
                                current_atom_neighbour_types.count('C') == 3 and\
                                current_atom_neighbour_types.count('N') == 1:
                            for y in current_atom_neighbours:
                                if y.type == 'C':
                                    if y in visited:
                                        visited.append(current_atom)
                                        central_peptide_chain.append(current_atom)
                                        for c in current_atom_neighbours:
                                            if c.type == 'N':
                                                current_atom = c
                    elif current_atom.type == 'N':
                        if current_atom_neighbour_types.count('H') == 2 and\
                                current_atom_neighbour_types.count('C') == 1:
                            for z in current_atom_neighbours:
                                if z.type == 'C':
                                    if z in visited:
                                        terminal_nitrogen = True
                                        visited.append(current_atom)
                                        central_peptide_chain.append(current_atom)
                        if current_atom_neighbour_types.count('H') == 1 and\
                                current_atom_neighbour_types.count('C') == 2:
                            visited.append(current_atom)
                            central_peptide_chain.append(current_atom)
                            for d in current_atom_neighbours:
                                if d.type == 'C':
                                    if d not in visited:
                                        current_atom = d


    return central_peptide_chain




if __name__ == "__main__":
    peptide = Smiles('C(C(NC(C(NC(C(=O)[O])C)=O)C(C)=O)=O)(N)CCC').smiles_to_structure()
    peptide = attach_to_domain_nrp(peptide, 'PCP')
    print(find_central_chain_nrp(peptide))