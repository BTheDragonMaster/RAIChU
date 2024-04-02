from pikachu.reactions.functional_groups import GroupDefiner, find_atoms, find_bonds
from raichu.reactions.general import initialise_atom_attributes

import timeout_decorator
from raichu.data.molecular_moieties import AMINO_FATTY_ACID, AMINO_ACID_BACKBONE, N_AMINO_ACID, C1_AMINO_ACID, \
    C2_AMINO_ACID, BETA_AMINO_ACID_BACKBONE, B_N_AMINO_ACID, B_C1_AMINO_ACID, B_C2_AMINO_ACID, B_C3_AMINO_ACID, \
    B_LEAVING_BOND, ACID_C1, OH_AMINO_ACID, N_AMINO_ACID_ATTACHED, C1_AMINO_ACID_ATTACHED, C2_AMINO_ACID_ATTACHED

POLYKETIDE_S = GroupDefiner('Sulphur atom polyketide', 'SC(C)=O', 0)
POLYKETIDE_S_INSERTED_O = GroupDefiner('Sulphur atom polyketide inserted O', 'SC(O)=O', 0)
POLYKETIDE_S_REDUCED_STARTER = GroupDefiner('Sulphur atom polyketide', 'SC(C)O', 0)


class CentralChain:
    def __init__(self, structure, chain_type):
        self.structure = structure
        if chain_type in ['ripp', 'nrps', 'pks']:
            self.chain_type = chain_type
        else:
            raise ValueError(f"Expected 'ripp', 'nrps', or 'pks' as central chain type; got {chain_type}")

    def label_central_chain(self):
        raise NotImplementedError


class PksCentralChain(CentralChain):
    def __init__(self, structure):
        super().__init__(structure, chain_type='pks')


class NrpsCentralChain(CentralChain):
    def __init__(self, structure):
        super().__init__(structure, chain_type='nrps')


class RippCentralChain(CentralChain):
    def __init__(self, structure):
        super().__init__(structure, chain_type='ripp')


@timeout_decorator.timeout(3)
def label_pk_central_chain(pk_starter_unit):
    """Finds the the atoms in the central chain of the polyketide starter unit,
    sets the in_central_chain Atom object attribute to True/False accordingly,
    after which the pk_starter_unit Structure object is returned.

    pk_starter_unit: PIKAChU Structure object of the PK starter unit
    """

    initialise_atom_attributes(pk_starter_unit)
    central_chain = []
    visited = []

    # Check if structure is a polyketide, define sulphur attached to domain
    locations_sulphur = find_atoms(POLYKETIDE_S, pk_starter_unit)
    if len(locations_sulphur) == 0:
        locations_sulphur = find_atoms(POLYKETIDE_S_INSERTED_O, pk_starter_unit)
        if len(locations_sulphur) == 0:
            locations_sulphur = find_atoms(POLYKETIDE_S_REDUCED_STARTER, pk_starter_unit)
    assert len(locations_sulphur) == 1
    
    sulphur_polyketide = locations_sulphur[0]
    central_chain.append(sulphur_polyketide)

    # Iterate over atom neighbours, check if the atom belongs to the central
    # carbon chain (i.e., not a methyl/ethyl/methoxy sidechain) and add to list
    for neighbour in sulphur_polyketide.neighbours:
        if neighbour.type == 'C':
            chain_carbon = neighbour
            end_carbon = False
            central_chain.append(chain_carbon)
            visited.append(chain_carbon)

            while not end_carbon:

                for next_atom in chain_carbon.neighbours:
                    ethyl_branch = False
                    inside_cycle = False
                    methyl_group = False

                    # Build lists of neighbouring atom types
                    if next_atom.type == 'C' and next_atom not in visited or \
                            (next_atom.type == 'O' and next_atom not in visited and
                             len(next_atom.get_neighbours('C')) == 2):

                        # Keep track of carbon neighbours

                        # Confirm carbon doesn't belong to ethyl sidebranch
                        types = []
                        c_neighbours = next_atom.get_neighbours('C')
                        if len(c_neighbours) == 2 and next_atom.get_neighbour('O') is None:

                            for carbon_neighbour in c_neighbours:
                                for branch_neighbour in carbon_neighbour.neighbours:
                                    types.append(branch_neighbour.type)

                            if types.count('H') == 4 and types.count('C') == 4:
                                ethyl_branch = True

                            for carbon_neighbour in c_neighbours:
                                if carbon_neighbour not in visited and carbon_neighbour != next_atom:
                                    ethyl_branch = False

                            visited.append(chain_carbon)

                        # Carbon in a terminal carboxylic acid group is the
                        # final carbon in the central chain
                        if len(next_atom.get_neighbours('O')) == 2:
                            for oxygen in next_atom.get_neighbours('O'):
                                if next_atom.get_bond(oxygen).type == 'single':
                                    central_chain.append(oxygen)
                                    visited.append(oxygen)
                            central_chain.append(next_atom)
                            end_carbon = True
                            visited.append(chain_carbon)

                        # Confirm carbon is not part of cycle
                        if len(c_neighbours) == 2:
                            if all(carbon.in_ring(pk_starter_unit) for carbon in c_neighbours):
                                visited.append(next_atom)
                                inside_cycle = True
                                end_carbon = True

                        # Case if carbon is part of a benzene ring
                        if len(c_neighbours) == 3 and len(next_atom.neighbours) == 3:
                            # central_chain.append(next_atom) #added 10/04/2021
                            visited.append(next_atom)
                            inside_cycle = True
                            end_carbon = True

                        # Case where the PK starter ends in two methyl branches
                        if len(next_atom.get_neighbours('C')) == 3 and len(next_atom.get_neighbours('H')) == 1:
                            nr_methyl_brances = 0
                            methyl_carbon = None
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
                                assert methyl_carbon is not None
                                methyl_group = True
                                visited.append(next_atom)
                                central_chain.append(next_atom)
                                visited.append(methyl_carbon)
                                central_chain.append(methyl_carbon)
                                end_carbon = True

                        # Confirm carbon doesn't belong to methyl sidebranch
                        if len(next_atom.get_neighbours('H')) == 3 or \
                                (len(next_atom.get_neighbours('H')) == 2 and len(next_atom.get_neighbours('*')) == 1):
                            methyl_group = True

                            # Terminal methyl group is not a side branch:
                            count_c_neighbours_not_visited = 0
                            for neighbouring_atom in chain_carbon.neighbours:
                                if neighbouring_atom.type == 'C' and neighbouring_atom not in visited:
                                    count_c_neighbours_not_visited += 1

                            if count_c_neighbours_not_visited < 2:
                                methyl_group = False
                                visited.append(chain_carbon)
                                for carbon_neighbour in chain_carbon.neighbours:
                                    if carbon_neighbour.type == 'C' and carbon_neighbour not in visited:
                                        next_neighbours = []
                                        for next_neighbour in carbon_neighbour.neighbours:
                                            next_neighbours.append(
                                                next_neighbour.type)
                                        if next_neighbours.count('H') == 3 or (
                                                next_neighbours.count(
                                                        'H') == 2 and (
                                                next_neighbours.count(
                                                        '*')) == 1):
                                            central_chain.append(carbon_neighbour)
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


def label_ripp_central_chain(peptide):
    initialise_atom_attributes(peptide)

    n_atoms_aa = find_atoms(N_AMINO_ACID_ATTACHED, peptide)
    c1_atoms_aa = find_atoms(C1_AMINO_ACID_ATTACHED, peptide)
    c2_atoms_aa = find_atoms(C2_AMINO_ACID_ATTACHED, peptide)
    oh_atoms_aa = find_atoms(OH_AMINO_ACID, peptide)

    for nitrogen in n_atoms_aa:
        nitrogen.annotations.in_central_chain = True
        nitrogen.annotations.n_atom_nmeth = True
    for carbon in c1_atoms_aa:
        carbon.annotations.in_central_chain = True
        carbon.annotations.chiral_c_ep = True
    for carbon in c2_atoms_aa:
        carbon.annotations.in_central_chain = True
    for oxygen in oh_atoms_aa:
        oxygen.annotations.in_central_chain = True


def label_nrp_central_chain(peptide, module_type='elongation'):
    initialise_atom_attributes(peptide)

    as_normal = False
    if module_type == 'starter':
        if peptide.find_substructures(AMINO_FATTY_ACID):
            acid = peptide
            label_acid_central_chain(acid)
        else:
            as_normal = True

    if as_normal or module_type == 'elongation':
        if peptide.find_substructures(AMINO_ACID_BACKBONE):
            n_atoms_aa = find_atoms(N_AMINO_ACID, peptide)
            c1_atoms_aa = find_atoms(C1_AMINO_ACID, peptide)
            c2_atoms_aa = find_atoms(C2_AMINO_ACID, peptide)

            assert len(n_atoms_aa) == 1
            assert len(c1_atoms_aa) == 1
            assert len(c2_atoms_aa) == 1

            n_atoms_aa[0].annotations.in_central_chain = True
            n_atoms_aa[0].annotations.n_atom_nmeth = True
            c1_atoms_aa[0].annotations.in_central_chain = True
            c1_atoms_aa[0].annotations.chiral_c_ep = True
            c2_atoms_aa[0].annotations.in_central_chain = True

        elif peptide.find_substructures(BETA_AMINO_ACID_BACKBONE):
            n_atoms_aa = find_atoms(B_N_AMINO_ACID, peptide)
            c1_atoms_aa = find_atoms(B_C1_AMINO_ACID, peptide)
            c2_atoms_aa = find_atoms(B_C2_AMINO_ACID, peptide)
            c3_atoms_aa = find_atoms(B_C3_AMINO_ACID, peptide)

            assert len(n_atoms_aa) == 1
            assert len(c1_atoms_aa) == 1
            assert len(c2_atoms_aa) == 1
            assert len(c3_atoms_aa) == 1

            if n_atoms_aa[0].in_ring(peptide) and c1_atoms_aa[0].in_ring(peptide) and c2_atoms_aa[0].in_ring(peptide):
                label_acid_central_chain(peptide)

            else:

                n_atoms_aa[0].annotations.in_central_chain = True
                n_atoms_aa[0].annotations.n_atom_nmeth = True
                c1_atoms_aa[0].annotations.in_central_chain = True
                c1_atoms_aa[0].annotations.chiral_c_ep = True
                c2_atoms_aa[0].annotations.in_central_chain = True
                c3_atoms_aa[0].annotations.in_central_chain = True

                bonds = find_bonds(B_LEAVING_BOND, peptide)
                assert len(bonds) == 1
                bond = bonds[0]
                if bond.atom_1.type == 'O':
                    bond.atom_1.annotations.leaving_oh_o = True

                    bond.atom_2.annotations.leaving_oh_h = True

                elif bond.atom_2.type == 'O':
                    bond.atom_2.annotations.leaving_oh_o = True

                    bond.atom_1.annotations.leaving_oh_h = True

        else:
            # If NRPS starter unit is a (fatty) acid instead of an amino acid
            acid = peptide
            label_acid_central_chain(acid)


def label_acid_central_chain(acid):
    c1_atom = None
    c2_atom = None
    oh_o_atom = None
    oh_h_atom = None
    backbone = []

    c1_atoms_aa = find_atoms(ACID_C1, acid)
    if c1_atoms_aa:
        c1_atom = c1_atoms_aa[0]

    assert c1_atom

    for neighbour in c1_atom.neighbours:
        neighbours = []

        if neighbour.type == 'C':
            for next_atom in neighbour.neighbours:
                neighbours.append(next_atom.type)
            if neighbours.count('O') == 2 and len(neighbours) == 3:
                c2_atom = neighbour

                # Find the leaving -OH group and annotate it

                for bond in c2_atom.bonds:
                    if bond.type == 'single' and (bond.has_neighbour('O')):
                        if bond.atom_1.type == 'O':
                            bond.atom_1.annotations.leaving_oh_o = True
                            oh_o_atom = bond.atom_1

                            bond.atom_2.annotations.leaving_oh_h = True
                            oh_h_atom = bond.atom_2

                        elif bond.atom_2.type == 'O':
                            bond.atom_2.annotations.leaving_oh_o = True
                            oh_o_atom = bond.atom_2

                            bond.atom_1.annotations.leaving_oh_h = True
                            oh_h_atom = bond.atom_1

    assert c2_atom
    assert oh_o_atom
    assert oh_h_atom

    c1_atom.annotations.in_central_chain = True
    c2_atom.annotations.in_central_chain = True
    c2_atom.annotations.c2_acid = True

    backbone.append(c1_atom)
    backbone.append(c2_atom)

    last_backbone_atom = c1_atom
    endpoint = False

    # Add atoms in chain with exactly 2 non-H neighbours to the central chain
    while not endpoint:

        neighbours = last_backbone_atom.get_non_hydrogen_neighbours()[:]
        for neighbour in last_backbone_atom.get_non_hydrogen_neighbours():
            if neighbour in backbone:
                neighbours.remove(neighbour)

        if not neighbours:
            endpoint = True

        else:
            # If the last backbone atom was the first ring member
            if last_backbone_atom.in_ring(acid):
                ring_members = last_backbone_atom.get_ring(acid)
                outgoing_pairs = []
                for atom in ring_members:
                    for neighbour in atom.get_non_hydrogen_neighbours():
                        if neighbour not in ring_members and neighbour not in backbone and not neighbour.in_ring(acid):
                            outgoing_pairs.append((atom, neighbour))

                if not outgoing_pairs:
                    endpoint = True

                else:
                    longest_index = -1
                    largest_subtree = 0
                    best_shortest_path = None
                    for i, pair in enumerate(outgoing_pairs):
                        ring_atom, outgoing_atom = pair
                        shortest_path = find_shortest_path_ring(ring_atom, last_backbone_atom, ring_members)
                        subtree_size = acid.get_subtree_size(outgoing_atom, set(ring_members))
                        if len(shortest_path) < 4 and subtree_size > largest_subtree:
                            longest_index = i
                            best_shortest_path = shortest_path

                    if longest_index >= 0:

                        backbone_ring_atom, next_backbone_atom = outgoing_pairs[longest_index]
                        for atom in best_shortest_path:
                            shortest_path_atom = acid.get_atom(atom)
                            if shortest_path_atom not in backbone:
                                shortest_path_atom.annotations.in_central_chain = True
                                backbone.append(shortest_path_atom)
                            next_backbone_atom = acid.get_atom(next_backbone_atom)
                            next_backbone_atom.annotations.in_central_chain = True
                            backbone.append(next_backbone_atom)
                            last_backbone_atom = next_backbone_atom
                    else:
                        endpoint = True

            else:
                longest_index = -1
                largest_subtree = 0

                for i, neighbour in enumerate(neighbours):
                    subtree_size = acid.get_subtree_size(neighbour, set(backbone))
                    if subtree_size > largest_subtree:
                        largest_subtree = subtree_size
                        longest_index = i

                if longest_index >= 0:
                    next_backbone_atom = acid.get_atom(neighbours[longest_index])
                else:
                    next_backbone_atom = acid.get_atom(neighbours[0])

                next_backbone_atom.annotations.in_central_chain = True
                backbone.append(next_backbone_atom)
                last_backbone_atom = next_backbone_atom


def find_shortest_path_ring(atom_1, atom_2, atom_list):
    max_distance = len(atom_list) // 2
    atom_1_idx = None
    atom_2_idx = None
    for i, atom in enumerate(atom_list):
        if atom == atom_1:
            atom_1_idx = i
        if atom == atom_2:
            atom_2_idx = i

    assert atom_1_idx is not None
    assert atom_2_idx is not None

    large_index = max([atom_1_idx, atom_2_idx])
    small_index = min([atom_1_idx, atom_2_idx])

    shortest_path = atom_list[small_index:large_index + 1]
    if len(shortest_path) - 1 > max_distance:
        shortest_path = atom_list[large_index:] + atom_list[:small_index + 1]

    return shortest_path
