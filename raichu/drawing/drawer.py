from collections import OrderedDict
from typing import Optional, List, Tuple, Dict, Set, Union

from pikachu.drawing.drawing import *
from pikachu.math_functions import *
from pikachu.chem.atom import Atom
from pikachu.chem.bond import Bond

from raichu.central_chain_detection.find_central_chain import find_central_chain, find_central_chain_not_attached, \
    find_central_chain_ripp, reorder_central_chain
from raichu.central_chain_detection.label_central_chain import label_ripp_central_chain


def get_ring(atom, structure):
    cycles = structure.cycles.all_cycles

    for i, cycle in enumerate(cycles):
        if atom in cycle:
            return cycle

    return None


class RaichuDrawer(Drawer):
    def __init__(self, structure, options=None, save_png=None, dont_show=False,
                 coords_only=True, dpi=100, save_svg=None, draw_Cs_in_pink=False, add_url=False, horizontal=False,
                 ripp=False, make_linear=True):
        self.dont_show = dont_show
        self.dpi = dpi
        self.draw_Cs_in_pink = draw_Cs_in_pink
        self.add_url = add_url
        self.horizontal = horizontal
        self.ripp = ripp
        self.make_linear = make_linear
        self.line_width = 2

        if options is None:
            self.options = Options()
        else:
            self.options = options
        if save_png is None:
            self.save_png = None
        else:
            # Check if filename is valid
            assert save_png.endswith('.png')
            self.save_png = save_png
        if save_svg is None:
            self.save_svg = None
        else:
            # Check if filename is valid
            assert save_svg.endswith('.svg')
            self.save_svg = save_svg
        super().__init__(structure, options=None, coords_only=coords_only)

    def center_on_carrier_domain(self):
        pass

    def find_clashing_atoms(self) -> List[Tuple[Atom, Atom]]:
        clashing_atoms = []
        for i, atom_1 in enumerate(self.drawn_atoms):
            for j in range(i + 1, len(self.drawn_atoms)):
                atom_2 = self.drawn_atoms[j]
                if not self.structure.bond_exists(atom_1, atom_2):
                    distance = Vector.subtract_vectors(atom_1.draw.position, atom_2.draw.position).get_squared_length()
                    if atom_1.has_neighbour('H') and atom_2.has_neighbour('H') and atom_1.type != 'C' and \
                            atom_2.type != 'C':
                        max_distance = self.options.bond_length_squared * 1.5
                    # elif atom_1.has_neighbour('H') or atom_2.has_neighbour('H'):
                    #     max_distance = self.options.bond_length_squared / 1.5
                    else:
                        max_distance = self.options.bond_length_squared / 2
                    if distance <= max_distance:
                        clashing_atoms.append((atom_1, atom_2))

        return clashing_atoms
    
    def find_out_of_bound_atoms(self, minimum_y) -> List[Atom]:
        clashing_atoms = []
        for i, atom in enumerate(self.drawn_atoms):
            if atom.draw.position.y < minimum_y and atom.type != 'I':
                clashing_atoms.append(atom)
        return clashing_atoms

    # TODO: add an option in pikachu Options object to set the min bond distance between which to consider fine-tuning
    def finetune_overlap_resolution(self, masked_bonds=None, highest_atom=None):

        if not masked_bonds:
            masked_bonds = set()

        else:
            masked_bonds = set(masked_bonds)

        # if self.total_overlap_score > self.options.overlap_sensitivity:
        clashing_atoms = self.find_clashing_atoms()

        best_bonds = []
        for atom_1, atom_2 in clashing_atoms:

            shortest_path = self.find_shortest_path(atom_1, atom_2)
            rotatable_bonds = []
            distances = []

            for i, bond in enumerate(shortest_path):
                distance_1 = i
                distance_2 = len(shortest_path) - i

                average_distance = len(shortest_path) / 2

                distance_metric = abs(
                    average_distance - distance_1) + abs(average_distance - distance_2)

                if self.bond_is_rotatable(bond) and bond not in masked_bonds:
                    rotatable_bonds.append(bond)
                    distances.append(distance_metric)

            best_bond = None
            optimal_distance = float('inf')
            for i, distance in enumerate(distances):
                if distance < optimal_distance:
                    # if highest_atom:
                    #     if atom_1.draw.position.y < highest_atom.draw.position.y and \
                    #             atom_2.draw.position.y < highest_atom.draw.position.y:
                    #         best_bond = rotatable_bonds[i]
                    #         optimal_distance = distance
                    # else:
                        best_bond = rotatable_bonds[i]
                        optimal_distance = distance

            if best_bond:
                best_bonds.append(best_bond)

        best_bonds = list(set(best_bonds))

        for best_bond in best_bonds:
            if self.total_overlap_score > self.options.overlap_sensitivity:

                subtree_size_1 = self.get_subgraph_size(
                    best_bond.atom_1, {best_bond.atom_2})
                subtree_size_2 = self.get_subgraph_size(
                    best_bond.atom_2, {best_bond.atom_1})

                if subtree_size_1 < subtree_size_2:
                    rotating_atom = best_bond.atom_1
                    parent_atom = best_bond.atom_2
                else:
                    rotating_atom = best_bond.atom_2
                    parent_atom = best_bond.atom_1

                overlap_score, _, _ = self.get_overlap_score()

                self.rotate_subtree(rotating_atom, parent_atom, math.radians(
                    15), parent_atom.draw.position)

                new_overlap_score, _, _ = self.get_overlap_score()

                scores = [new_overlap_score]

                for i in range(12):
                    self.rotate_subtree(rotating_atom, parent_atom, math.radians(
                        30), parent_atom.draw.position)
                    new_overlap_score, _, _ = self.get_overlap_score()
                    scores.append(new_overlap_score)

                assert len(scores) == 13
                scores = scores[:12]

                best_i = 0
                best_score = scores[0]

                for i, score in enumerate(scores):
                    if score < best_score:
                        best_score = score
                        best_i = i

                self.total_overlap_score = best_score

                self.rotate_subtree(rotating_atom, parent_atom, math.radians(
                    30 * best_i + 1), parent_atom.draw.position)
    
    def finetune_overlap_bubbles(self, masked_bonds=None, highest_atom=None):

        if not masked_bonds:
            masked_bonds = set()

        else:
            masked_bonds = set(masked_bonds)

        domain = None
        
        # Find minimum y

        for atom in self.structure.graph:
            if atom.type == 'I':
                if atom.annotations.domain_type in ["Leader", "Follower", "ACP", "PCP"]:
                    domain = atom
        # Add also padding for hydrogens drawn as a atom label but not as an explicit hydrogen (e.g. HN)
        assert domain is not None

        minimum_y = domain.draw.position.y + 5
        
        clashing_atoms = self.find_out_of_bound_atoms(minimum_y)

        best_bonds = []
        for atom_1 in clashing_atoms:

            shortest_path = self.find_shortest_path(atom_1, domain)
            rotatable_bonds = []
            distances = []

            for i, bond in enumerate(shortest_path):
                distance_1 = i
                distance_2 = len(shortest_path) - i

                average_distance = len(shortest_path) / 2

                distance_metric = abs(
                    average_distance - distance_1) + abs(average_distance - distance_2)

                if self.bond_is_rotatable(bond) and bond not in masked_bonds:
                    rotatable_bonds.append(bond)
                    distances.append(distance_metric)

            best_bond = None
            optimal_distance = float('inf')
            for i, distance in enumerate(distances):
                if distance < optimal_distance:
                    # if highest_atom:
                    #     if atom_1.draw.position.y < highest_atom.draw.position.y and \
                    #             atom_2.draw.position.y < highest_atom.draw.position.y:
                    #         best_bond = rotatable_bonds[i]
                    #         optimal_distance = distance
                    # else:
                    best_bond = rotatable_bonds[i]
                    optimal_distance = distance

            if best_bond:
                best_bonds.append(best_bond)

        best_bonds = list(set(best_bonds))

        for best_bond in best_bonds:
            if len(self.find_out_of_bound_atoms(minimum_y))>0:

                path1 = self.find_shortest_path(best_bond.atom_1, domain)
                path2 = self.find_shortest_path(best_bond.atom_2, domain)
                
                if len(path2) < len(path1):
                    rotating_atom = best_bond.atom_1
                    parent_atom = best_bond.atom_2
                else:
                    rotating_atom = best_bond.atom_2
                    parent_atom = best_bond.atom_1

                new_overlap_score, _, _ = self.get_overlap_score()

                scores = [new_overlap_score]
                overlappings = [len(self.find_out_of_bound_atoms(minimum_y))]

                self.rotate_subtree(rotating_atom, parent_atom, math.radians(
                    15), parent_atom.draw.position)

                for i in range(12):
                    self.rotate_subtree(rotating_atom, parent_atom, math.radians(
                        30), parent_atom.draw.position)
                    overlappings.append(len(self.find_out_of_bound_atoms(minimum_y)))
                    new_overlap_score, _, _ = self.get_overlap_score()
                    scores.append(new_overlap_score)

                assert len(scores) == 13
                scores = scores[:12]

                assert len(overlappings) == 13
                overlappings = overlappings[:12]

                best_i = 0
                best_overlapping = overlappings[0]
                best_score = self.total_overlap_score
                for i, overlapping in enumerate(overlappings):
                    if overlapping < best_overlapping and scores[i] < self.options.overlap_sensitivity + self.total_overlap_score:
                        best_overlapping = overlapping
                        best_i = i
                        best_score = scores[i]
                self.rotate_subtree(rotating_atom, parent_atom, math.radians(
                    30 * best_i + 1), parent_atom.draw.position)
                self.total_overlap_score = best_score
                

    def draw(self, coords_only: bool = False) -> None:

        if not self.options.draw_hydrogens:
            self.hide_hydrogens()

        self.get_atom_nr_to_atom()
        self.define_rings()

        if not self.multiple:
            self._process_structure()
            if self.make_linear:
                self.linearise()
            self.set_chiral_bonds()
            if not coords_only:
                self.draw_structure()
        else:
            self.restore_ring_information()

    # TODO: replace this entirely with PIKAChU's svg drawer
    def draw_structure(self):
        # Find the plotting dimensions of the molecule such that the canvas can be scaled to fit the molecule

        min_x = 100000000
        max_x = -100000000
        min_y = 100000000
        max_y = -100000000

        for atom in self.structure.graph:
            if atom.draw.positioned:
                if atom.draw.position.x < min_x:
                    min_x = atom.draw.position.x
                if atom.draw.position.y < min_y:
                    min_y = atom.draw.position.y
                if atom.draw.position.x > max_x:
                    max_x = atom.draw.position.x
                if atom.draw.position.y > max_y:
                    max_y = atom.draw.position.y

        height = max_y - min_y
        width = max_x - min_x

        fig, ax = plt.subplots(figsize=((width + 2 * self.options.padding) / 50.0,
                                        (height + 2 * self.options.padding) / 50.0), dpi=100)

        ax.set_aspect('equal', adjustable='box')
        ax.axis('off')

        ax.set_xlim([min_x - self.options.padding, max_x + self.options.padding])
        ax.set_ylim([min_y - self.options.padding, max_y + self.options.padding])
        plt.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)

        params = {'mathtext.default': 'regular'}
        plt.rcParams.update(params)

        ring_centers_x = []
        ring_centers_y = []

        for ring in self.rings:
            self.set_ring_center(ring)

            ring_centers_x.append(ring.center.x)
            ring_centers_y.append(ring.center.y)

        for bond_nr, bond in self.structure.bonds.items():
            if bond.atom_1.draw.positioned and bond.atom_2.draw.positioned:
                line = Line(bond.atom_1.draw.position, bond.atom_2.draw.position, bond.atom_1, bond.atom_2)
                midpoint = line.get_midpoint()

                if bond.type == 'single':
                    if bond in self.chiral_bonds:
                        orientation, chiral_center = self.chiral_bond_to_orientation[bond]
                        self._plot_chiral_bond(orientation, chiral_center, line, ax, midpoint)
                    else:
                        self._plot_halflines(line, ax, midpoint)
                elif bond.type == 'double':
                    if not self._is_terminal(bond.atom_1) and not self._is_terminal(bond.atom_2):
                        self._plot_halflines(line, ax, midpoint)

                        common_ring_numbers = self._get_common_rings(bond.atom_1, bond.atom_2)

                        if common_ring_numbers:
                            common_rings = []
                            for ring_nr in common_ring_numbers:
                                common_rings.append(self.get_ring(ring_nr))

                            common_rings.sort(key=lambda x: len(x.members))
                            common_ring = common_rings[0]
                            ring_centre = common_ring.center
                            second_line = line.double_line_towards_center(ring_centre, self.options.bond_spacing,
                                                                          self.options.double_bond_length)
                            second_line_midpoint = second_line.get_midpoint()
                            self._plot_halflines_double(second_line, ax, second_line_midpoint)

                        else:
                            bond_neighbours = bond.atom_1.drawn_neighbours + bond.atom_2.drawn_neighbours
                            if bond_neighbours:
                                vectors = [atom.draw.position for atom in bond_neighbours]
                                gravitational_point = Vector.get_average(vectors)
                                second_line = line.double_line_towards_center(gravitational_point,
                                                                              self.options.bond_spacing,
                                                                              self.options.double_bond_length)
                                second_line_midpoint = second_line.get_midpoint()
                                self._plot_halflines_double(second_line, ax, second_line_midpoint)
                            else:
                                print("Shouldn't happen!")
                    else:
                        if self._is_terminal(bond.atom_1) and self._is_terminal(bond.atom_2):
                            dummy_1 = Vector(bond.atom_1.draw.position.x + 1, bond.atom_1.draw.position.y + 1)
                            dummy_2 = Vector(bond.atom_1.draw.position.x - 1, bond.atom_1.draw.position.y - 1)
                            double_bond_line_1 = line.double_line_towards_center(dummy_1,
                                                                                 self.options.bond_spacing / 2.0,
                                                                                 self.options.double_bond_length)
                            double_bond_line_1_midpoint = double_bond_line_1.get_midpoint()
                            double_bond_line_2 = line.double_line_towards_center(dummy_2,
                                                                                 self.options.bond_spacing / 2.0,
                                                                                 self.options.double_bond_length)
                            double_bond_line_2_midpoint = double_bond_line_2.get_midpoint()

                            self._plot_halflines_double(double_bond_line_1, ax, double_bond_line_1_midpoint)
                            self._plot_halflines_double(double_bond_line_2, ax, double_bond_line_2_midpoint)

                        else:

                            if self._is_terminal(bond.atom_1):
                                terminal_atom = bond.atom_1
                                branched_atom = bond.atom_2
                            else:
                                terminal_atom = bond.atom_2
                                branched_atom = bond.atom_1

                            if len(branched_atom.drawn_neighbours) >= 3:
                                closest_two = self.get_sorted_distances_from_list(terminal_atom,
                                                                                  branched_atom.drawn_neighbours)
                                closest_atom_1 = closest_two[0][1]
                                closest_atom_2 = closest_two[1][1]

                                line = Line(terminal_atom.draw.position, branched_atom.draw.position, terminal_atom,
                                            branched_atom)

                                double_bond_line_1, double_bond_line_2 = line.get_perpendicular_lines(
                                    self.options.bond_spacing / 2.0)
                                terminal_atom_pos_1 = double_bond_line_1.get_atom_coords(terminal_atom)
                                terminal_atom_pos_2 = double_bond_line_2.get_atom_coords(terminal_atom)

                                closest_atom_to_pos_1 = terminal_atom_pos_1.get_closest_atom(closest_atom_1,
                                                                                             closest_atom_2)
                                closest_atom_to_pos_2 = terminal_atom_pos_2.get_closest_atom(closest_atom_1,
                                                                                             closest_atom_2)

                                bond_1_line = Line(branched_atom.draw.position, closest_atom_to_pos_1.draw.position,
                                                   branched_atom, closest_atom_to_pos_1)
                                bond_2_line = Line(branched_atom.draw.position, closest_atom_to_pos_2.draw.position,
                                                   branched_atom, closest_atom_to_pos_2)

                                double_bond_line_1_midpoint = double_bond_line_1.get_midpoint()
                                double_bond_line_2_midpoint = double_bond_line_2.get_midpoint()

                                intersection_1 = double_bond_line_1.find_intersection(bond_1_line)
                                intersection_2 = double_bond_line_2.find_intersection(bond_2_line)

                                if terminal_atom.draw.position.x > branched_atom.draw.position.x:
                                    # check for parallel lines
                                    if intersection_1 and intersection_1.x < 100000 and intersection_1.y < 100000:
                                        double_bond_line_1.point_1 = intersection_1
                                    if intersection_2 and intersection_2.x < 100000 and intersection_2.y < 100000:
                                        double_bond_line_2.point_1 = intersection_2

                                else:
                                    # check for parallel lines
                                    if intersection_1 and intersection_1.x < 100000 and intersection_1.y < 100000:
                                        double_bond_line_1.point_2 = intersection_1
                                    if intersection_2 and intersection_2.x < 100000 and intersection_2.y < 100000:
                                        double_bond_line_2.point_2 = intersection_2

                                self._plot_halflines(double_bond_line_1, ax, double_bond_line_1_midpoint)
                                self._plot_halflines(double_bond_line_2, ax, double_bond_line_2_midpoint)

                            else:
                                self._plot_halflines(line, ax, midpoint)

                                bond_neighbours = bond.atom_1.drawn_neighbours + bond.atom_2.drawn_neighbours
                                if bond_neighbours:
                                    vectors = [atom.draw.position for atom in bond_neighbours]
                                    gravitational_point = Vector.get_average(vectors)
                                    second_line = line.get_parallel_line(gravitational_point,
                                                                         self.options.bond_spacing)
                                    second_line_midpoint = second_line.get_midpoint()
                                    self._plot_halflines(second_line, ax, second_line_midpoint)
                                else:
                                    print("Shouldn't happen!")

                elif bond.type == 'triple':
                    self._plot_halflines(line, ax, midpoint)
                    line_1, line_2 = line.get_parallel_lines(self.options.bond_spacing)
                    line_1_midpoint = line_1.get_midpoint()
                    line_2_midpoint = line_2.get_midpoint()
                    self._plot_halflines(line_1, ax, line_1_midpoint)
                    self._plot_halflines(line_2, ax, line_2_midpoint)

        # If the starter unit contains an unknown moiety, number this R1
        for atom in self.structure.graph:
            if atom.type == '*' and not atom.annotations.unknown_index:
                atom.annotations.unknown_index = 1
        for atom in self.structure.graph:
            if ((atom.type != 'C' or self.draw_Cs_in_pink) and getattr(atom.annotations,"domain_type", None) not in ["Follower", "Leader"]) and atom.draw.positioned:
                text_h = ''
                text_h_pos = None

                if (atom.type != 'C' or self.draw_Cs_in_pink or atom.draw.draw_explicit) :
                    text = atom.type
                else:
                    text = ''

                if hasattr(atom.annotations, "domain_type") and atom.annotations.domain_type:
                    if atom.annotations.domain_type in ["Follower", "Leader"]:
                        text = ""
                    else:
                        text = atom.annotations.domain_type

                horizontal_alignment = 'center'

                orientation = self._get_hydrogen_text_orientation(atom)
                if orientation == 'H_above_atom':
                    text_h_pos = Vector(
                        atom.draw.position.x, atom.draw.position.y + 6)
                if orientation == 'H_below_atom':
                    text_h_pos = Vector(
                        atom.draw.position.x, atom.draw.position.y - 6)

                atom_draw_position = Vector(
                    atom.draw.position.x, atom.draw.position.y)

                if atom.type == '*':
                    neighbouring_c = atom.get_neighbour('C')
                    assert neighbouring_c
                    # In order to not let the number of the sidechain overlap
                    # with the bond, move Rx symbol along hor. axis depending on
                    # if group is added to the left or right of the chain
                    delta_x_r = self.get_delta_x_sidechain(
                        atom, neighbouring_c)
                    atom_draw_position.x += delta_x_r
                    text = fr'$R_{atom.annotations.unknown_index}$'

                if not atom.charge and (atom.type != 'C' or self.draw_Cs_in_pink or atom.draw.draw_explicit) and getattr(atom.annotations, "domain_type", None) not in ["Follower", "Leader"]:

                    if atom.draw.has_hydrogen:
                        hydrogen_count = 0
                        for neighbour in atom.neighbours:
                            if neighbour.type == 'H' and not neighbour.draw.is_drawn:
                                hydrogen_count += 1

                        if hydrogen_count and not (atom.type == 'C' and self.draw_Cs_in_pink or getattr(atom.annotations, "domain_type", None) in ["Follower", "Leader"]):

                            if hydrogen_count > 1:
                                if orientation == 'H_before_atom':
                                    text = r'$H_{hydrogens}{atom_type}$'.format(hydrogens=hydrogen_count,
                                                                                atom_type=atom.type)
                                    horizontal_alignment = 'right'
                                    atom_draw_position.x += 3
                                elif orientation == 'H_below_atom' or orientation == 'H_above_atom':
                                    text = atom.type
                                    text_h = r'$H_{hydrogens}$'.format(
                                        hydrogens=hydrogen_count)

                                else:
                                    text = r'${atom_type}H_{hydrogens}$'.format(hydrogens=hydrogen_count,
                                                                                atom_type=atom.type)
                                    horizontal_alignment = 'left'
                                    atom_draw_position.x -= 3
                            elif hydrogen_count == 1:
                                if orientation == 'H_before_atom':
                                    text = f'H{atom.type}'
                                    horizontal_alignment = 'right'
                                    atom_draw_position.x += 3
                                elif orientation == 'H_below_atom' or orientation == 'H_above_atom':
                                    text = atom.type
                                    text_h = 'H'
                                else:
                                    text = f'{atom.type}H'
                                    horizontal_alignment = 'left'
                                    atom_draw_position.x -= 3

                elif atom.charge:
                    if atom.charge > 0:
                        charge_symbol = '+'
                    else:
                        charge_symbol = '-'

                    hydrogen_count = 0
                    for neighbour in atom.neighbours:
                        if neighbour.type == 'H' and not neighbour.draw.is_drawn:
                            hydrogen_count += 1

                    if not hydrogen_count:

                        if abs(atom.charge) > 1:

                            text = r'${atom_type}^{charge}{charge_symbol}$'.format(charge=atom.charge,
                                                                                   atom_type=atom.type,
                                                                                   charge_symbol=charge_symbol)
                        elif abs(atom.charge) == 1:
                            text = r'${atom_type}^{charge_symbol}$'.format(atom_type=atom.type,
                                                                           charge_symbol=charge_symbol)

                        horizontal_alignment = 'left'
                        atom_draw_position.x -= 3
                    else:

                        if hydrogen_count > 1:
                            if orientation == 'H_before_atom':
                                if abs(atom.charge) > 1:
                                    text = r'$H_{hydrogens}{atom_type}^{charge}{charge_symbol}$'.format(
                                        hydrogens=hydrogen_count,
                                        atom_type=atom.type,
                                        charge=atom.charge,
                                        charge_symbol=charge_symbol)
                                elif abs(atom.charge) == 1:
                                    text = r'$H_{hydrogens}{atom_type}^{charge_symbol}$'.format(
                                        hydrogens=hydrogen_count,
                                        atom_type=atom.type,
                                        charge_symbol=charge_symbol)

                                horizontal_alignment = 'right'
                                atom_draw_position.x += 3
                            elif orientation == 'H_above_atom' or orientation == 'H_below_atom':
                                text_h = r'$H_{hydrogens}$'.format(
                                    hydrogens=hydrogen_count)
                                if abs(atom.charge) > 1:
                                    text = r'${atom_type}^{charge}{charge_symbol}$'.format(atom_type=atom.type,
                                                                                           charge=atom.charge,
                                                                                           charge_symbol=charge_symbol)
                                elif abs(atom.charge) == 1:
                                    text = r'${atom_type}^{charge_symbol}$'.format(atom_type=atom.type,
                                                                                   charge_symbol=charge_symbol)
                            else:
                                if abs(atom.charge) > 1:
                                    text = r'${atom_type}H_{hydrogens}^{charge}{charge_symbol}$'.format(
                                        hydrogens=hydrogen_count,
                                        atom_type=atom.type,
                                        charge=atom.charge,
                                        charge_symbol=charge_symbol)
                                elif abs(atom.charge) == 1:
                                    text = r'${atom_type}H_{hydrogens}^{charge_symbol}$'.format(
                                        hydrogens=hydrogen_count,
                                        atom_type=atom.type,
                                        charge_symbol=charge_symbol)

                                horizontal_alignment = 'left'
                                atom_draw_position.x -= 3
                        elif hydrogen_count == 1:
                            if orientation == 'H_before_atom':
                                if abs(atom.charge) > 1:

                                    text = r'$H{atom_type}^{charge}{charge_symbol}$'.format(atom_type=atom.type,
                                                                                            charge=atom.charge,
                                                                                            charge_symbol=charge_symbol)
                                elif abs(atom.charge) == 1:
                                    text = r'$H{atom_type}^{charge_symbol}$'.format(atom_type=atom.type,
                                                                                    charge_symbol=charge_symbol)
                                horizontal_alignment = 'right'
                                atom_draw_position.x += 3
                            elif orientation == 'H_above_atom' or orientation == 'H_below_atom':
                                text_h = 'H'
                                if abs(atom.charge) > 1:

                                    text = r'${atom_type}^{charge}{charge_symbol}$'.format(atom_type=atom.type,
                                                                                           charge=atom.charge,
                                                                                           charge_symbol=charge_symbol)
                                elif abs(atom.charge) == 1:
                                    text = r'${atom_type}^{charge_symbol}$'.format(atom_type=atom.type,
                                                                                   charge_symbol=charge_symbol)

                            else:
                                if abs(atom.charge) > 1:
                                    text = r'${atom_type}H^{charge}{charge_symbol}$'.format(atom_type=atom.type,
                                                                                            charge=atom.charge,
                                                                                            charge_symbol=charge_symbol)

                                elif abs(atom.charge) == 1:
                                    text = r'${atom_type}H^{charge_symbol}$'.format(atom_type=atom.type,
                                                                                    charge_symbol=charge_symbol)
                                horizontal_alignment = 'left'
                                atom_draw_position.x -= 3

                atom_color = atom.draw.colour
                if self.draw_Cs_in_pink:
                    if atom.type == 'C':
                        atom_color = "magenta"
                if self.add_url:
                    if text:
                        plt.text(atom_draw_position.x, atom_draw_position.y,
                                 text, url=str(atom),
                                 horizontalalignment=horizontal_alignment,
                                 verticalalignment='center',
                                 color=atom_color)
                    if text_h:
                        plt.text(text_h_pos.x, text_h_pos.y,
                                 text_h, url=str(atom),
                                 horizontalalignment='center',
                                 verticalalignment='center',
                                 color=atom_color)
                else:
                    if text:
                        plt.text(atom_draw_position.x, atom_draw_position.y,
                                 text,
                                 horizontalalignment=horizontal_alignment,
                                 verticalalignment='center',
                                 color=atom_color)
                    if text_h:
                        plt.text(text_h_pos.x, text_h_pos.y,
                                 text_h,
                                 horizontalalignment='center',
                                 verticalalignment='center',
                                 color=atom_color)

        # If a png filename is included in the initialization of the
        # Raichu_drawer object, don't show the structure, but do save it as a
        # png image to the provided filename
        if self.dont_show:
            return self

        elif self.save_png is None and not self.dont_show and self.save_svg is None:
            plt.show()

        else:
            if self.save_png:
                plt.savefig(self.save_png)
            elif self.save_svg:
                plt.savefig(self.save_svg)
            plt.clf()
            plt.close()

    def colour_modules(self):
        pass

    def place_top_atom(self, backbone_atoms, attachment_point):


        # Fix position of the S and C atom and the S-C angle
        top_atom = backbone_atoms[0]
        first_carbon = backbone_atoms[1]

        top_atom.draw.position.x = first_carbon.draw.position.x
        top_atom.draw.position.y = first_carbon.draw.position.y + 15

        angle = get_angle(first_carbon.draw.position,
                          top_atom.draw.position)

        angle_degrees = round(math.degrees(angle), 3)

        correct_angle_deg = -120
        delta_angle_deg = correct_angle_deg - angle_degrees
        delta_angle_rad = math.radians(delta_angle_deg)

        self.rotate_subtree(top_atom, first_carbon, delta_angle_rad,
                            first_carbon.draw.position)

        # Fix position domain straight above sulphur atom
        attachment_point.draw.position.x = top_atom.draw.position.x
        attachment_point.draw.position.y = top_atom.draw.position.y + 15


    @staticmethod
    def get_opposite_placement(placement):
        if placement == 'right':
            return 'left'
        elif placement == 'left':
            return 'right'
        else:
            raise ValueError("Only placements accepted are 'left' and 'right'.")

    def get_placements(self, backbone, atom_to_ring, stop_linearising):
        backbone_to_placement = OrderedDict()
        for i, atom in enumerate(backbone):
            if atom == stop_linearising:
                break

            if i == 0:
                backbone_to_placement[atom] = 'right'
            elif i == 1:
                backbone_to_placement[atom] = 'left'
            else:
                if atom not in backbone_to_placement:
                    previous_placement = backbone_to_placement.get(backbone[i - 1])
                    ring = atom_to_ring.get(atom)
                    if ring and len(ring) == 2:
                        if ring[0] in backbone_to_placement:
                            backbone_to_placement[ring[1]] = backbone_to_placement.get(ring[0])
                        else:
                            backbone_to_placement[atom] = self.get_opposite_placement(previous_placement)
                    elif hasattr(atom, "domain_type") and atom.annotations.domain_type == "Leader":
                        backbone_to_placement[atom] = previous_placement
                    else:
                        backbone_to_placement[atom] = self.get_opposite_placement(previous_placement)

        return backbone_to_placement

    def fix_rings(self, rings, backbone_atoms, backbone_to_placement):
        for ring in rings:
            if len(ring) == 1:
                atom_1 = ring[0]
                if atom_1 == backbone_atoms[-1]:
                    continue
                drawing_ring = self.get_ring(atom_1.draw.rings[0])
                center = Vector(0, 0)
                for atom in drawing_ring.members:
                    center.add(atom.draw.position)

                center.divide(len(drawing_ring.members))
                angle = get_angle(atom_1.draw.position,
                                  center)
                
                if backbone_to_placement[atom_1] == 'right':
                    desired_angle = 180
                else:
                    desired_angle = 0.0
                required_rotation_deg = desired_angle - angle
                required_rotation_rad = math.radians(required_rotation_deg)

                for atom in drawing_ring.members:
                    if atom != atom_1:
                        atom.draw.position.rotate_around_vector(required_rotation_rad, atom_1.draw.position)

            elif len(ring) == 2:
                atom_1, atom_2 = ring
                masked = {atom_1, atom_2}

                first_atom_cycle = None
                for next_atom in atom_1.neighbours:
                    if next_atom.type != 'H' and next_atom not in backbone_atoms:
                        first_atom_cycle = next_atom

                assert first_atom_cycle

                if backbone_to_placement[atom_1] == backbone_to_placement[atom_2] == 'right':
                    if atom_1.draw.position.x > first_atom_cycle.draw.position.x:
                        for atom in self.traverse_substructure(first_atom_cycle, masked):
                            delta_x = 2 * (atom_1.draw.position.x - atom.draw.position.x)
                            atom.draw.position.x += delta_x

                elif backbone_to_placement[atom_1] == backbone_to_placement[atom_2] == 'left':

                    if atom_1.draw.position.x < first_atom_cycle.draw.position.x:
                        for atom in self.traverse_substructure(first_atom_cycle, masked):
                            delta_x = 2 * (atom.draw.position.x - atom_1.draw.position.x)
                            atom.draw.position.x -= delta_x

    def get_overlap_score(self) -> Tuple[float, List[Tuple[float, Atom]], Dict[Atom, float]]:
        total = 0.0

        overlap_scores = {}
        for atom in self.drawn_atoms:
            overlap_scores[atom] = 0.0

        for i, atom_1 in enumerate(self.drawn_atoms):
            for j in range(i + 1, len(self.drawn_atoms)):
                atom_2 = self.drawn_atoms[j]
                distance = Vector.subtract_vectors(atom_1.draw.position, atom_2.draw.position).get_squared_length()
                if atom_1 not in atom_2.neighbours:
                    if distance < self.options.bond_length_squared + 1:

                        if atom_1.has_neighbour('H') and atom_2.has_neighbour('H') and atom_1.type != 'C' \
                                and atom_2.type != 'C':

                            weight = self.options.overlap_sensitivity + abs((self.options.bond_length - math.sqrt(distance))) / self.options.bond_length

                        else:
                            weight = abs((self.options.bond_length - math.sqrt(distance))) / self.options.bond_length
                        total += weight
                        overlap_scores[atom_1] += weight
                        overlap_scores[atom_2] += weight

        sorted_overlaps = []

        for atom in self.drawn_atoms:
            sorted_overlaps.append((overlap_scores[atom], atom))

        sorted_overlaps.sort(key=lambda x: x[0], reverse=True)

        return total, sorted_overlaps, overlap_scores

    def linearise(self):

        attachment_point = None
        domains = []

        for atom in self.structure.graph:
            if atom.type == 'I':
                if atom.annotations.domain_type in ["Leader", "Follower", "ACP", "PCP"]:
                    domains.append(atom)
                    # Ensures the follower peptide is set as attachment point if both leader and follower are present
                    if not attachment_point:
                        attachment_point = atom
                    elif attachment_point.annotations.domain_type == "Leader":
                        attachment_point = atom

        if attachment_point and attachment_point.annotations.domain_type == "Leader":
            horizontal_rotation = 'anticlockwise'
        else:
            horizontal_rotation = 'clockwise'

        if self.ripp:
            label_ripp_central_chain(self.structure)

        self.structure.refresh_structure()

        if attachment_point:
            if self.ripp:
                backbone_atoms = find_central_chain_ripp(self.structure)
            else:
                backbone_atoms = find_central_chain(self.structure)
        else:
            backbone_atoms = find_central_chain_not_attached(self.structure)

        if attachment_point:

            self.place_top_atom(backbone_atoms, attachment_point)

        self.structure.refresh_structure()
        backbone, full_rings, rings, atom_to_ring, stop_linearising = reorder_central_chain(backbone_atoms, self)
        if len(domains) == 2:
            for domain in domains:
                if domain.annotations.domain_type == 'Leader':
                    backbone.append(domain)

        backbone_to_placement = self.get_placements(backbone, atom_to_ring, stop_linearising)

        i = 0

        while i < len(backbone) - 1:

            atom_1 = backbone[i]
            atom_2 = backbone[i + 1]

            if atom_1 not in backbone_to_placement or atom_2 not in backbone_to_placement:
                break

            current_angle = get_angle(atom_1.draw.position, atom_2.draw.position)

            if backbone_to_placement[atom_1] == 'right' and backbone_to_placement[atom_2] == 'left':
                desired_angle = 60.0
            elif backbone_to_placement[atom_1] == 'left' and backbone_to_placement[atom_2] == 'right':
                desired_angle = 120.0
            elif backbone_to_placement[atom_1] == backbone_to_placement[atom_2]:
                desired_angle = 90.0
            else:
                raise ValueError("Only supported orientations are 'left' and 'right'.")

            required_rotation_deg = desired_angle - current_angle
            required_rotation_rad = math.radians(required_rotation_deg)

            if not (atom_2 in atom_to_ring and atom_1 in atom_to_ring and atom_to_ring[atom_2] == atom_to_ring[atom_1]):

                self.rotate_subtree(atom_2, atom_1, required_rotation_rad,
                                    atom_1.draw.position)
            else:
                if len(atom_to_ring[atom_2]) == 2:
                    self.rotate_subtree(atom_2, atom_1, required_rotation_rad,
                                        atom_1.draw.position)
                elif len(atom_to_ring[atom_2]) == 3:
                    if atom_2 == atom_to_ring[atom_2][1]:

                        self.rotate_subtree(atom_2, atom_1, required_rotation_rad,
                                            atom_1.draw.position)
                    elif atom_2 == atom_to_ring[atom_2][2]:

                        ring = atom_to_ring[atom_2]
                        ring_atoms = []
                        sidechain_atoms = []
                        for atom in atom_2.get_ring(self.structure):
                            if atom not in backbone:
                                ring_atoms.append(atom)

                        if abs(required_rotation_deg) > 90:
                            for atom in ring_atoms:
                                for neighbour in atom.neighbours:
                                    if neighbour.type != 'H':
                                        for sidechain_atom in self.traverse_substructure(neighbour,
                                                                                         set(ring_atoms + backbone)):
                                            if sidechain_atom not in sidechain_atoms and \
                                                    sidechain_atom not in ring_atoms and sidechain_atom not in backbone:
                                                sidechain_atoms.append(sidechain_atom)

                            for ring_atom in ring_atoms:
                                ring_atom.draw.position.mirror_about_line(ring[0].draw.position, ring[1].draw.position)
                            for sidechain_atom in sidechain_atoms:
                                sidechain_atom.draw.position.mirror_about_line(ring[0].draw.position, ring[1].draw.position)
                        self.rotate_subtree_ring(atom_2, set(ring_atoms + [atom_1]), required_rotation_rad,
                                                 atom_1.draw.position)

            i += 1

        atoms_in_rings = []
        for ring in full_rings:
            for atom in ring:
                atoms_in_rings.append(atom)
        atoms_in_rings = set(atoms_in_rings)

        self.fix_rings(rings, backbone, backbone_to_placement)

        self.position_sidechains(backbone, backbone_to_placement, atoms_in_rings)

        # self.resolve_primary_overlaps()
        self.total_overlap_score, sorted_overlap_scores, atom_to_scores = self.get_overlap_score()
        central_chain_bonds = set()

        for bond in self.structure.bonds.values():
            # if bond.atom_1.annotations.in_central_chain or \
            #         bond.atom_2.annotations.in_central_chain:
            if bond.atom_1 in backbone or bond.atom_2 in backbone:
                central_chain_bonds.add(bond)

        self.finetune_overlap_resolution(
            masked_bonds=central_chain_bonds, highest_atom=backbone[0])

        self.resolve_secondary_overlaps(sorted_overlap_scores)

        if self.horizontal:
            if horizontal_rotation == 'clockwise':
                self.rotate_structure(-1.5707)
            else:
                self.rotate_structure(1.5707)

    def resolve_overlaps(self, linear: bool = True, masked_bonds: Optional[Set["Bond"]] = None) -> None:
        """

        Parameters
        ----------
        linear: bool, if True, resolve overlaps for a linearised molecule; if False,
            resolve overlaps for a non-linearised molecule
        masked_bonds: set of bonds, bonds that are not allowed to be rotated

        """
        if not linear:

            self.resolve_primary_overlaps()
            self.total_overlap_score, sorted_overlap_scores, atom_to_scores = self.get_overlap_score()
            self.finetune_overlap_resolution()
            self.resolve_secondary_overlaps(sorted_overlap_scores)

        else:

            pass

    def find_high_atoms(self, atoms: List[Atom], max_y: float) -> List[Atom]:
        """
        Returns atoms that are drawn too high on the plot
        Parameters
        ----------
        atoms: list of [Atom, ->]
        max_y: float, maximum vertical y coordinate an atom is allowed to have

        Returns
        -------
        high_atoms: list of [Atom, ->], atoms with a y coordinate higher than max_y
        """
        high_atoms: List[Atom] = []
        for atom in atoms:
            if atom.draw.position.y > max_y:
                high_atoms.append(atom)
        return high_atoms

    def find_clashing_sidechain_atoms(self, sidechain_atoms):
        clashing_atoms: List[Tuple[Atom, Atom]] = []
        for i, atom_1 in enumerate(self.drawn_atoms):
            for j in range(i + 1, len(self.drawn_atoms)):
                atom_2 = self.drawn_atoms[j]
                if atom_1 in sidechain_atoms and not self.structure.bond_exists(atom_1, atom_2):
                    distance = Vector.subtract_vectors(atom_1.draw.position, atom_2.draw.position).get_squared_length()
                    if atom_1.has_neighbour('H') and atom_2.has_neighbour('H') and atom_1.type != 'C' and \
                            atom_2.type != 'C':
                        max_distance = self.options.bond_length_squared * 1.5
                    elif atom_1.has_neighbour('H') or atom_2.has_neighbour('H'):
                        max_distance = self.options.bond_length_squared
                    else:
                        max_distance = self.options.bond_length_squared * 0.8
                    if distance <= max_distance:
                        clashing_atoms.append((atom_1, atom_2))

        return clashing_atoms

    def resolve_sidechain_overlaps(self, backbone, backbone_to_placement, backbone_rings, central_chain_bonds):
        max_y = backbone[0].draw.position.y
        placed_atoms: List[Atom] = []
        backbone_to_sidechain: Dict[Atom, List[List[Atom]]] = {}
        backbone_to_rotatable_bonds: Dict[Atom, List[Bond]] = {}

        for backbone_atom in backbone:
            sidechains = self.get_sidechains(backbone, backbone_atom, backbone_rings)
            backbone_to_sidechain[backbone_atom] = sidechains
            beta_atoms: List[Atom] = self.get_beta_atoms(backbone, backbone_atom, backbone_rings)
            for sidechain in sidechains:
                beta_atom: Optional[Atom] = None

                for atom in sidechain:
                    if atom in beta_atoms:
                        beta_atom = atom
                        break

                assert beta_atom

                sidechain_bonds = self.get_sidechain_bonds(sidechain)
                masked_bond = self.structure.bond_lookup[beta_atom][backbone_atom]
                rotatable_bonds: List[Bond] = []

                for bond in sidechain_bonds:
                    if self.bond_is_rotatable(bond) and bond != masked_bond:
                        rotatable_bonds.append(bond)
                backbone_to_rotatable_bonds[backbone_atom] = rotatable_bonds


            if len(sidechains) == 2:
                pass
            elif len(sidechains) == 1:
                sidechain = sidechains[0]
                atoms_to_adjust: List[Atom] = self.find_high_atoms(sidechain, max_y)
                clashing_atoms = self.find_clashing_sidechain_atoms(sidechain)
                if clashing_atoms or atoms_to_adjust:
                    beta_atoms: List[Atom] = self.get_beta_atoms(backbone, backbone_atom, backbone_rings)
                    assert len(beta_atoms) == 1
                    beta_atom = beta_atoms[0]

                    sidechain_bonds = self.get_sidechain_bonds(sidechain)
                    rotatable_bonds: List[Bond] = []

                    masked_bond = self.structure.bond_lookup[beta_atom][backbone_atom]
                    for bond in sidechain_bonds:
                        if self.bond_is_rotatable(bond) and bond != masked_bond and \
                                (bond.atom_1 in beta_atoms or bond.atom_2 in beta_atoms):
                            rotatable_bonds.append(bond)

                    if not rotatable_bonds:
                        for bond in sidechain_bonds:
                            if self.bond_is_rotatable(bond) and bond != masked_bond:
                                rotatable_bonds.append(bond)

    def resolve_secondary_overlaps(self, sorted_scores: List[Tuple[float, Atom]]) -> None:
        for score, atom in sorted_scores:
            if score > self.options.overlap_sensitivity:
                if len(atom.drawn_neighbours) <= 1:
                    if atom.drawn_neighbours and atom.drawn_neighbours[0]._adjacent_to_stereobond():
                        continue

                    closest_atom = self.get_closest_atom(atom)

                    drawn_neighbours = closest_atom.drawn_neighbours

                    if len(drawn_neighbours) <= 1:
                        if not closest_atom.draw.previous_position:
                            closest_position = drawn_neighbours[0].draw.position
                        else:
                            closest_position = closest_atom.draw.previous_position

                    else:
                        if not closest_atom.draw.previous_position:
                            closest_position = drawn_neighbours[0].draw.position
                        else:
                            closest_position = closest_atom.draw.position

                    if not atom.draw.previous_position:
                        atom_previous_position = atom.drawn_neighbours[0].draw.position
                    else:
                        atom_previous_position = atom.draw.previous_position
                    #
                    # atom.draw.position.rotate_away_from_vector(closest_position, atom_previous_position,
                    #                                            math.radians(20))

    def rotate_subtree_ring(self, root, masked, angle, center):
        for atom in self.traverse_substructure(root, masked):
            atom.draw.position.rotate_around_vector(angle, center)
            for anchored_ring in atom.draw.anchored_rings:
                if anchored_ring.center:
                    anchored_ring.center.rotate_around_vector(angle, center)

    def get_beta_atoms(self, backbone: List[Atom], backbone_atom: Atom, backbone_rings: List[Atom]) -> List[Atom]:
        """
        Returns atoms directly adjacent to the backbone atom that is not in the same ring as a backbone atom

        Parameters
        ----------
        backbone: list of [Atom, ->] with each atom a backbone atom
        backbone_atom: Atom instance, backbone atom for which the sidechains are determined
        backbone_rings: list of [Atom, ->] with each atom a member of a ring that overlaps with the backbone

        Returns
        -------
        beta_atoms: list of [Atom, ->], representing the atoms directly adjacent to the backbone atom
        """
        beta_atoms: List[Atom] = []
        for neighbour in backbone_atom.neighbours:
            if neighbour.type != 'H' and not neighbour.annotations.domain_type:
                if neighbour not in backbone and neighbour not in backbone_rings:
                    beta_atoms.append(neighbour)

        return beta_atoms



    def get_sidechains(self, backbone: List[Atom], backbone_atom: Atom, backbone_rings: List[Atom]) -> List[List[Atom]]:
        """
        Returns the non-ring sidechains of a backbone atom

        Parameters
        ----------
        backbone: list of [Atom, ->] with each atom a backbone atom
        backbone_atom: Atom instance, backbone atom for which the sidechains are determined
        backbone_rings: list of [Atom, ->] with each atom a member of a ring that overlaps with the backbone

        Returns
        -------
        sidechains: list of [[Atom, ->] ->], representing the sidechains coming out of the backbone atom

        """
        sidechains: List[List[Atom]] = []
        beta_atoms = self.get_beta_atoms(backbone, backbone_atom, backbone_rings)
        for beta_atom in beta_atoms:
            sidechain: List[Atom] = []

            for atom in self.traverse_substructure(beta_atom, {backbone_atom}):
                sidechain.append(atom)
            sidechains.append(sidechain)
        return sidechains

    def get_sidechain_bonds(self, sidechain: List[Atom]) -> List[Bond]:
        """
        Returns a list of bonds that exist between sidechain atoms

        Parameters
        ----------
        sidechain: list of [Atom, ->]

        Returns
        -------
        sidechain_bonds: list of [Bond, ->]

        """
        bonds: Set[Bond] = set()

        for atom in sidechain:
            for neighbour in atom.neighbours:
                if neighbour in sidechain:
                    bonds.add(self.structure.bond_lookup[atom][neighbour])

        sidechain_bonds: List[Bond] = list(bonds)

        return sidechain_bonds

    def position_sidechains(self, backbone, backbone_to_placement, atoms_in_rings):
        backbone_atom_before = None
        backbone_to_neighbours = {}

        for backbone_atom in backbone:
            neighbours = []
            for neighbour in backbone_atom.neighbours:
                if neighbour.type != 'H' and not neighbour.annotations.domain_type:
                    if neighbour not in backbone and neighbour not in atoms_in_rings \
                            and neighbour not in neighbours:
                        neighbours.append(neighbour)
            backbone_to_neighbours[backbone_atom] = neighbours

        for backbone_atom, neighbours in backbone_to_neighbours.items():
            placement = backbone_to_placement[backbone_atom]
            if len(neighbours) == 0:
                pass
            elif len(neighbours) == 1:
                neighbour = neighbours[0]

                current_angle = get_angle(backbone_atom.draw.position, neighbour.draw.position)

                if placement == 'right':
                    desired_angle = 180.0
                elif placement == 'left':
                    desired_angle = 0.0
                else:
                    raise ValueError("Placement must be 'left' or 'right'.")
                # if ring in central chain but last member or tertiary carbon
                if backbone_atom in atoms_in_rings and backbone_atom_before in atoms_in_rings and backbone_atom_before:
                    desired_angle += 180

                required_rotation_deg = desired_angle - current_angle
                required_rotation_rad = math.radians(required_rotation_deg)

                self.rotate_subtree(neighbour, backbone_atom, required_rotation_rad, backbone_atom.draw.position)

            elif len(neighbours) == 2:
                neighbour_1, neighbour_2 = neighbours

                current_angle_1 = get_angle(backbone_atom.draw.position, neighbour_1.draw.position)
                current_angle_2 = get_angle(backbone_atom.draw.position, neighbour_2.draw.position)

                if placement == 'right':
                    desired_angle_1 = 140.0
                    desired_angle_2 = 220.0
                elif placement == 'left':
                    desired_angle_1 = 40.0
                    desired_angle_2 = -40.0
                else:
                    raise ValueError("Placement must be 'left' or 'right'.")

                required_rotation_deg_1 = desired_angle_1 - current_angle_1
                required_rotation_rad_1 = math.radians(required_rotation_deg_1)

                required_rotation_deg_2 = desired_angle_2 - current_angle_2
                required_rotation_rad_2 = math.radians(required_rotation_deg_2)

                self.rotate_subtree(neighbour_1, backbone_atom, required_rotation_rad_1, backbone_atom.draw.position)
                self.rotate_subtree(neighbour_2, backbone_atom, required_rotation_rad_2, backbone_atom.draw.position)
            else:
                raise ValueError("Backbone atoms can only have two non-backbone neighbours at most")
            backbone_atom_before = backbone_atom


    def plot_halflines_s_domain(self, line, ax, midpoint):
        truncated_line = line.get_truncated_line(
            self.options.short_bond_length)
        self.plot_line_dashed(truncated_line, ax, color='#a6a6a6')

    def plot_line_dashed(self, line, ax, color='grey'):
        with matplotlib.rc_context({'path.sketch': (5, 10, 1)}):
            ax.plot([line.point_1.x, line.point_2.x],
                    [line.point_1.y, line.point_2.y], color=color,
                    linewidth=self.line_width/1.5)

    def get_delta_x_sidechain(self, sidechain_atom, neighbouring_carbon):
        """Returns the change tot the x-position of the sidegroup symbol in
        the structure drawing as float
        sidechain_atom: Atom object of the unknown sidechain group
        neighoburin_carbon: Atom object of the carbon in the central chain
        bound to the unknown sidechain group
        """
        sidechainatom_x = sidechain_atom.draw.position.x
        carbon_x = neighbouring_carbon.draw.position.x
        if (carbon_x - sidechainatom_x) < -10:
            delta_x = 1.5
        elif (carbon_x - sidechainatom_x) > 10:
            delta_x = -1.5
        else:
            delta_x = 0.0

        return delta_x

    def draw_svg(self, annotation: Union[None, str] = None, numbered_atoms: List = None) -> str:

        self.set_annotation_for_grouping(annotation)

        ring_centers_x = []
        ring_centers_y = []

        for ring in self.rings:
            self.set_ring_center(ring)

            ring_centers_x.append(ring.center.x)
            ring_centers_y.append(ring.center.y)

        for bond_nr, bond in self.structure.bonds.items():
            if bond.atom_1.draw.positioned and bond.atom_2.draw.positioned:
                line = Line(bond.atom_1.draw.position, bond.atom_2.draw.position, bond.atom_1, bond.atom_2)
                midpoint = line.get_midpoint()

                if bond.type == 'single':
                    if bond in self.chiral_bonds:

                        orientation, chiral_center = self.chiral_bond_to_orientation[bond]
                        self._draw_chiral_bond(orientation, chiral_center, line, midpoint)
                    else:
                        self._draw_halflines(line, midpoint)
                elif bond.type == 'double':
                    if not self._is_terminal(bond.atom_1) and not self._is_terminal(bond.atom_2):
                        self._draw_halflines(line, midpoint)

                        common_ring_numbers = self._get_common_rings(bond.atom_1, bond.atom_2)

                        if common_ring_numbers:
                            common_rings = []
                            for ring_nr in common_ring_numbers:
                                common_rings.append(self.get_ring(ring_nr))

                            common_rings.sort(key=lambda x: len(x.members))
                            common_ring = common_rings[0]
                            ring_centre = common_ring.center
                            second_line = line.double_line_towards_center(ring_centre, self.options.bond_spacing,
                                                                          self.options.double_bond_length)
                            second_line_midpoint = second_line.get_midpoint()
                            self._draw_halflines_double(second_line, second_line_midpoint)

                        else:
                            bond_neighbours = bond.atom_1.drawn_neighbours + bond.atom_2.drawn_neighbours
                            if bond_neighbours:
                                vectors = [atom.draw.position for atom in bond_neighbours]
                                gravitational_point = Vector.get_average(vectors)
                                second_line = line.double_line_towards_center(gravitational_point,
                                                                              self.options.bond_spacing,
                                                                              self.options.double_bond_length)
                                second_line_midpoint = second_line.get_midpoint()
                                self._draw_halflines_double(second_line, second_line_midpoint)
                            else:
                                print("Shouldn't happen!")
                    else:
                        if self._is_terminal(bond.atom_1) and self._is_terminal(bond.atom_2):
                            dummy_1 = Vector(bond.atom_1.draw.position.x + 1, bond.atom_1.draw.position.y + 1)
                            dummy_2 = Vector(bond.atom_1.draw.position.x - 1, bond.atom_1.draw.position.y - 1)
                            double_bond_line_1 = line.double_line_towards_center(dummy_1,
                                                                                 self.options.bond_spacing / 2.0,
                                                                                 self.options.double_bond_length)
                            double_bond_line_1_midpoint = double_bond_line_1.get_midpoint()
                            double_bond_line_2 = line.double_line_towards_center(dummy_2,
                                                                                 self.options.bond_spacing / 2.0,
                                                                                 self.options.double_bond_length)
                            double_bond_line_2_midpoint = double_bond_line_2.get_midpoint()

                            self._draw_halflines_double(double_bond_line_1, double_bond_line_1_midpoint)
                            self._draw_halflines_double(double_bond_line_2, double_bond_line_2_midpoint)

                        else:

                            if self._is_terminal(bond.atom_1):
                                terminal_atom = bond.atom_1
                                branched_atom = bond.atom_2
                            else:
                                terminal_atom = bond.atom_2
                                branched_atom = bond.atom_1

                            if len(branched_atom.drawn_neighbours) >= 3:
                                closest_two = self.get_sorted_distances_from_list(terminal_atom,
                                                                                  branched_atom.drawn_neighbours)
                                closest_atom_1 = closest_two[0][1]
                                closest_atom_2 = closest_two[1][1]

                                line = Line(terminal_atom.draw.position, branched_atom.draw.position, terminal_atom,
                                            branched_atom)

                                double_bond_line_1, double_bond_line_2 = line.get_perpendicular_lines(
                                    self.options.bond_spacing / 2.0)
                                terminal_atom_pos_1 = double_bond_line_1.get_atom_coords(terminal_atom)
                                terminal_atom_pos_2 = double_bond_line_2.get_atom_coords(terminal_atom)

                                closest_atom_to_pos_1 = terminal_atom_pos_1.get_closest_atom(closest_atom_1,
                                                                                             closest_atom_2)
                                closest_atom_to_pos_2 = terminal_atom_pos_2.get_closest_atom(closest_atom_1,
                                                                                             closest_atom_2)

                                bond_1_line = Line(branched_atom.draw.position, closest_atom_to_pos_1.draw.position,
                                                   branched_atom, closest_atom_to_pos_1)
                                bond_2_line = Line(branched_atom.draw.position, closest_atom_to_pos_2.draw.position,
                                                   branched_atom, closest_atom_to_pos_2)

                                double_bond_line_1_midpoint = double_bond_line_1.get_midpoint()
                                double_bond_line_2_midpoint = double_bond_line_2.get_midpoint()

                                intersection_1 = double_bond_line_1.find_intersection(bond_1_line)
                                intersection_2 = double_bond_line_2.find_intersection(bond_2_line)

                                if terminal_atom.draw.position.x > branched_atom.draw.position.x:
                                    # check for parallel lines
                                    if intersection_1 and intersection_1.x < 100000 and intersection_1.y < 100000:
                                        double_bond_line_1.point_1 = intersection_1
                                    if intersection_2 and intersection_2.x < 100000 and intersection_2.y < 100000:
                                        double_bond_line_2.point_1 = intersection_2

                                else:
                                    # check for parallel lines
                                    if intersection_1 and intersection_1.x < 100000 and intersection_1.y < 100000:
                                        double_bond_line_1.point_2 = intersection_1
                                    if intersection_2 and intersection_2.x < 100000 and intersection_2.y < 100000:
                                        double_bond_line_2.point_2 = intersection_2

                                self._draw_halflines(double_bond_line_1, double_bond_line_1_midpoint)
                                self._draw_halflines(double_bond_line_2, double_bond_line_2_midpoint)

                            else:
                                self._draw_halflines(line, midpoint)

                                bond_neighbours = bond.atom_1.drawn_neighbours + bond.atom_2.drawn_neighbours
                                if bond_neighbours:
                                    vectors = [atom.draw.position for atom in bond_neighbours]
                                    gravitational_point = Vector.get_average(vectors)
                                    second_line = line.get_parallel_line(gravitational_point,
                                                                         self.options.bond_spacing)
                                    second_line_midpoint = second_line.get_midpoint()
                                    self._draw_halflines(second_line, second_line_midpoint)
                                else:
                                    print("Shouldn't happen!")

                elif bond.type == 'triple':
                    self._draw_halflines(line, midpoint)
                    line_1, line_2 = line.get_parallel_lines(self.options.bond_spacing)
                    line_1_midpoint = line_1.get_midpoint()
                    line_2_midpoint = line_2.get_midpoint()
                    self._draw_halflines(line_1, line_1_midpoint)
                    self._draw_halflines(line_2, line_2_midpoint)

        for atom in self.structure.graph:
            if atom.draw.positioned:
                svg_text = ''
                svg_h_text = ''
                svg_charge_text = ''
                svg_h_count_text = ''

                if atom.type != 'C' or atom.draw.draw_explicit or atom.charge:
                    if atom.type == 'C' and not atom.charge:
                        svg_text = self._draw_text('.', atom.draw.position.x, atom.draw.position.y - 2,
                                                   color=atom.draw.colour)
                    else:
                        svg_text = self._draw_text(atom.type, atom.draw.position.x, atom.draw.position.y,
                                                   color=atom.draw.colour)
                if hasattr(atom, "domain_type") and atom.annotations.domain_type in ["Leader", "Follower", "ACP", "PCP"]:
                    svg_text = ''

                        # TODO: Make this possible in svg writing
                        # text = self.set_r_group_indices_subscript(atom.type)

                orientation = self._get_hydrogen_text_orientation(atom)

                # Swap up-down orientation due to swapped svg coordinate system
                if orientation == 'H_below_atom':
                    orientation = 'H_above_atom'
                elif orientation == 'H_above_atom':
                    orientation = 'H_below_atom'

                if atom.type != 'C' or atom.draw.draw_explicit or atom.charge:

                    hydrogen_count = 0

                    for neighbour in atom.neighbours:
                        if neighbour.type == 'H' and not neighbour.draw.is_drawn:
                            hydrogen_count += 1

                    h_x = atom.draw.position.x
                    h_y = atom.draw.position.y
                    h_subscript_x = atom.draw.position.x
                    h_subscript_y = atom.draw.position.y + 3
                    charge_x = atom.draw.position.x + 6
                    charge_y = atom.draw.position.y - 3

                    if orientation == 'H_below_atom':
                        h_y = atom.draw.position.y + 7
                        h_subscript_x += 2
                        h_subscript_y += 10
                    elif orientation == 'H_above_atom':
                        h_y = atom.draw.position.y - 7
                        h_subscript_x += 2
                        h_subscript_y -= 4
                    elif orientation == 'H_before_atom':
                        if hydrogen_count > 1:
                            h_x -= 10
                            h_subscript_x -= 5
                        elif hydrogen_count == 1:
                            h_x -= 6
                    else:
                        h_x += len(atom.type) * 6
                        h_subscript_x += len(atom.type) * 6 + 5

                    if abs(atom.charge):
                        if hydrogen_count and orientation == 'H_after_atom':
                            if hydrogen_count > 1:
                                charge_x += len(atom.type) * 6 + 7
                            else:
                                charge_x += len(atom.type) * 6 + 2

                    h_pos = Vector(h_x, h_y)
                    h_subscript_pos = Vector(h_subscript_x, h_subscript_y)
                    charge_pos = Vector(charge_x, charge_y)

                    if hydrogen_count:
                        svg_h_text = self._draw_text('H', h_pos.x, h_pos.y,
                                                     color=atom.draw.colour)
                        if hydrogen_count > 1:
                            svg_h_count_text = self._draw_text(hydrogen_count, h_subscript_pos.x, h_subscript_pos.y,
                                                               font_size=self.options.svg_font_size_small,
                                                               color=atom.draw.colour)
                    if atom.charge:
                        if atom.charge > 0:
                            charge_symbol = '+'
                        else:
                            charge_symbol = '-'

                        if abs(atom.charge) > 1:
                            charge_text = f"{abs(atom.charge)}{charge_symbol}"
                        else:
                            charge_text = charge_symbol

                        svg_charge_text = self._draw_text(charge_text, charge_pos.x, charge_pos.y,
                                                          font_size=self.options.svg_font_size_small,
                                                          color=atom.draw.colour)

                if svg_text or svg_h_text or svg_charge_text or svg_h_count_text:
                    if self.structure_id:
                        text_group = f'<g id="atom_{atom}_{self.structure_id}_text">\n'
                    else:
                        text_group = f'<g id="atom_{atom}_text">\n'
                    if svg_text:
                        text_group += svg_text
                        text_group += '\n'
                    if svg_h_text:
                        text_group += svg_h_text
                        text_group += '\n'
                    if svg_charge_text:
                        text_group += svg_charge_text
                        text_group += '\n'
                    if svg_h_count_text:
                        text_group += svg_h_count_text
                        text_group += '\n'
                    text_group += '</g>'
                    self._add_svg_element(text_group, atom)

        if numbered_atoms:
            self.show_atom_numbers(numbered_atoms)

        svg = self._assemble_svg()
        return svg


def get_angle(vector1, vector2):
    difference = Vector.subtract_vectors(vector1, vector2)
    difference_angle = difference.angle()
    return math.degrees(difference_angle)
