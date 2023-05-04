from pikachu.drawing.drawing import *
from pikachu.math_functions import *

from raichu.central_chain_detection.find_central_chain import find_central_chain, find_central_chain_not_attached, \
    find_central_chain_ripp, reorder_central_chain
from raichu.central_chain_detection.label_central_chain import label_nrp_central_chain


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

    def export_to_json(self):
        json_structure = {"atoms": [],
                          "bonds": [],
                          "modules": []}

        ring_centers_x = []
        ring_centers_y = []

        for ring in self.rings:
            self.set_ring_center(ring)

            ring_centers_x.append(ring.center.x)
            ring_centers_y.append(ring.center.y)

        for bond_nr, bond in self.structure.bonds.items():
            if bond.atom_1.draw.positioned and bond.atom_2.draw.positioned:
                json_bond = {"id": bond.nr,
                             "atom_1": bond.atom_1.nr,
                             "atom_2": bond.atom_2.nr,
                             "type": bond.type,
                             "chiral": False,
                             "lines": []}

                if bond in self.chiral_bonds:

                    orientation, chiral_center = self.chiral_bond_to_orientation[bond]
                    json_bond["chiral"] = True
                    json_bond["wedge_orientation"] = orientation
                    json_bond["wedge_origin"] = chiral_center.nr
                elif (bond.atom_1.type == 'S' and bond.atom_2.annotations.domain_type or
                      (bond.atom_2.type == 'S' and bond.atom_1.annotations.domain_type)):
                    json_bond["type"] = "squiggle"

        # If the starter unit contains an unknown moiety, number this R1
        for atom in self.structure.graph:
            if atom.type == '*' and not atom.annotations.unknown_index:
                atom.annotations.unknown_index = 1
        for atom in self.structure.graph:
            if atom.draw.positioned:
                json_atom = {"id": atom.nr,
                             "type": atom.type,
                             "text": '',
                             "atom_x": atom.draw.position.x,
                             "atom_y": atom.draw.position.y}

                if atom.type != 'C' and atom.draw.positioned:
                    text_h = ''
                    text_h_pos = None

                    if atom.type != 'C' or atom.draw.draw_explicit:
                        text = atom.type
                    else:
                        text = ''

                    if atom.annotations.domain_type:
                        text = atom.annotations.domain_type

                    horizontal_alignment = 'center'

                    orientation = self.get_hydrogen_text_orientation(atom)
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

                    if not atom.charge and (atom.type != 'C' or atom.draw.draw_explicit or self.draw_Cs_in_pink ):

                        if atom.draw.has_hydrogen:
                            hydrogen_count = 0
                            for neighbour in atom.neighbours:
                                if neighbour.type == 'H' and not neighbour.draw.is_drawn:
                                    hydrogen_count += 1

                            if hydrogen_count:

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

    def find_clashing_atoms(self) -> List[Tuple[Atom, Atom]]:
        clashing_atoms = []
        for i, atom_1 in enumerate(self.drawn_atoms):
            for j in range(i + 1, len(self.drawn_atoms)):
                atom_2 = self.drawn_atoms[j]
                if not self.structure.bond_exists(atom_1, atom_2):
                    distance = Vector.subtract_vectors(atom_1.draw.position, atom_2.draw.position).get_squared_length()
                    # if distance < self.options.bond_length_squared:
                    if distance < self.options.bond_length / 2:
                        clashing_atoms.append((atom_1, atom_2))

        return clashing_atoms

    def finetune_overlap_resolution(self, masked_bonds=None, highest_atom=None):

        if not masked_bonds:
            masked_bonds = set()

        else:
            masked_bonds = set(masked_bonds)

        if self.total_overlap_score > self.options.overlap_sensitivity:
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
                        if highest_atom:
                            if not atom_1.draw.position.y > highest_atom.draw.position.y and not atom_2.draw.position.y > highest_atom.draw.position.y:
                                best_bond = rotatable_bonds[i]
                                optimal_distance = distance
                        else:
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

                    scores = [overlap_score]

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

    def draw(self, coords_only: bool = False) -> None:

        if not self.options.draw_hydrogens:
            self.hide_hydrogens()

        self.get_atom_nr_to_atom()
        self.define_rings()

        if not self.multiple:
            self.process_structure()
            if self.make_linear:
                self.linearise()
            self.set_chiral_bonds()
            if not coords_only:
                self.draw_structure()
        else:
            self.restore_ring_information()

    def draw_structure(self):
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
        self.line_width = 2

        fig, ax = plt.subplots(figsize=((width + 2 * self.options.padding) /
                                        50.0, (height + 2 * self.options.padding) / 50.0), dpi=self.dpi)

        ax.set_aspect('equal', adjustable='box')
        ax.axis('off')

        ax.set_xlim(
            [min_x - self.options.padding, max_x + self.options.padding])
        ax.set_ylim(
            [min_y - self.options.padding, max_y + self.options.padding])
        plt.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0,
                            hspace=0)

        params = {'mathtext.default': 'regular', }
        plt.rcParams.update(params)

        ring_centers_x = []
        ring_centers_y = []

        for ring in self.rings:
            self.set_ring_center(ring)

            ring_centers_x.append(ring.center.x)
            ring_centers_y.append(ring.center.y)

        for bond_nr, bond in self.structure.bonds.items():
            if bond.atom_1.draw.positioned and bond.atom_2.draw.positioned:
                line = Line(bond.atom_1.draw.position,
                            bond.atom_2.draw.position, bond.atom_1,
                            bond.atom_2)
                midpoint = line.get_midpoint()
                truncated_line = line.get_truncated_line(
                    self.options.short_bond_length)
                if bond.type == 'single':
                    if bond in self.chiral_bonds:
                        orientation, chiral_center = self.chiral_bond_to_orientation[bond]
                        self.plot_chiral_bond(orientation, chiral_center, line,
                                              ax, midpoint)
                    else:
                        if hasattr(bond.atom_1.annotations, "domain_type") and hasattr(bond.atom_2.annotations, "domain_type"):
                            if (bond.atom_1.type == 'S' and
                                    bond.atom_2.annotations.domain_type or
                                    (bond.atom_2.type == 'S' and
                                     bond.atom_1.annotations.domain_type)):
                                self.plot_halflines_s_domain(
                                    line, ax, midpoint)
                            else:
                                self.plot_halflines(line, ax, midpoint)
                        else:
                            self.plot_halflines(line, ax, midpoint)
                elif bond.type == 'double':
                    if not self.is_terminal(
                            bond.atom_1) and not self.is_terminal(bond.atom_2):
                        self.plot_halflines(line, ax, midpoint)

                        common_ring_numbers = self.get_common_rings(
                            bond.atom_1, bond.atom_2)

                        if common_ring_numbers:
                            common_rings = []
                            for ring_nr in common_ring_numbers:
                                common_rings.append(self.get_ring(ring_nr))

                            common_rings.sort(key=lambda x: len(x.members))
                            common_ring = common_rings[0]
                            ring_centre = common_ring.center
                            second_line = line.double_line_towards_center(
                                ring_centre, self.options.bond_spacing,
                                self.options.double_bond_length)
                            second_line_midpoint = second_line.get_midpoint()
                            self.plot_halflines_double(second_line, ax,
                                                       second_line_midpoint)

                        else:
                            bond_neighbours = bond.atom_1.drawn_neighbours +\
                                bond.atom_2.drawn_neighbours
                            if bond_neighbours:
                                vectors = [atom.draw.position for atom in
                                           bond_neighbours]
                                gravitational_point = Vector.get_average(
                                    vectors)
                                second_line = line.double_line_towards_center(
                                    gravitational_point,
                                    self.options.bond_spacing,
                                    self.options.double_bond_length)
                                second_line_midpoint = second_line.get_midpoint()
                                self.plot_halflines_double(second_line, ax,
                                                           second_line_midpoint)
                            else:
                                print("Shouldn't happen!")
                    else:
                        if self.is_terminal(bond.atom_1) and self.is_terminal(
                                bond.atom_2):
                            dummy_1 = Vector(bond.atom_1.draw.position.x + 1,
                                             bond.atom_1.draw.position.y + 1)
                            dummy_2 = Vector(bond.atom_1.draw.position.x - 1,
                                             bond.atom_1.draw.position.y - 1)
                            double_bond_line_1 = line.double_line_towards_center(
                                dummy_1,
                                self.options.bond_spacing / 2.0,
                                self.options.double_bond_length)
                            double_bond_line_1_midpoint = double_bond_line_1.get_midpoint()
                            double_bond_line_2 = line.double_line_towards_center(
                                dummy_2,
                                self.options.bond_spacing / 2.0,
                                self.options.double_bond_length)
                            double_bond_line_2_midpoint = double_bond_line_2.get_midpoint()

                            self.plot_halflines_double(double_bond_line_1, ax,
                                                       double_bond_line_1_midpoint)
                            self.plot_halflines_double(double_bond_line_2, ax,
                                                       double_bond_line_2_midpoint)

                        else:

                            if self.is_terminal(bond.atom_1):
                                terminal_atom = bond.atom_1
                                branched_atom = bond.atom_2
                            else:
                                terminal_atom = bond.atom_2
                                branched_atom = bond.atom_1

                            if len(branched_atom.drawn_neighbours) >= 3:
                                closest_two = self.get_sorted_distances_from_list(
                                    terminal_atom,
                                    branched_atom.drawn_neighbours)
                                closest_atom_1 = closest_two[0][1]
                                closest_atom_2 = closest_two[1][1]

                                line = Line(terminal_atom.draw.position,
                                            branched_atom.draw.position,
                                            terminal_atom, branched_atom)

                                bond_1_line = Line(branched_atom.draw.position,
                                                   closest_atom_1.draw.position,
                                                   branched_atom,
                                                   closest_atom_1)
                                bond_2_line = Line(branched_atom.draw.position,
                                                   closest_atom_2.draw.position,
                                                   branched_atom,
                                                   closest_atom_1)

                                double_bond_line_1 = line.double_line_towards_center(
                                    closest_atom_1.draw.position,
                                    self.options.bond_spacing / 2.0,
                                    self.options.double_bond_length)
                                double_bond_line_2 = line.double_line_towards_center(
                                    closest_atom_2.draw.position,
                                    self.options.bond_spacing / 2.0,
                                    self.options.double_bond_length)

                                double_bond_line_1_midpoint = \
                                    double_bond_line_1.get_midpoint()
                                double_bond_line_2_midpoint = \
                                    double_bond_line_2.get_midpoint()

                                intersection_1 = double_bond_line_1.find_intersection(
                                    bond_1_line)
                                intersection_2 = double_bond_line_2.find_intersection(
                                    bond_2_line)

                                if terminal_atom.draw.position.x >\
                                        branched_atom.draw.position.x:
                                    double_bond_line_1.point_1 = intersection_1
                                    double_bond_line_2.point_1 = intersection_2
                                else:
                                    double_bond_line_1.point_2 = intersection_1
                                    double_bond_line_2.point_2 = intersection_2

                                self.plot_halflines(double_bond_line_1, ax,
                                                    double_bond_line_1_midpoint)
                                self.plot_halflines(double_bond_line_2, ax,
                                                    double_bond_line_2_midpoint)

                            else:
                                pass

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

                orientation = self.get_hydrogen_text_orientation(atom)
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

    def place_top_atom(self, backbone_atoms):
        attachment_point = None

        for atom in self.structure.graph:
            if atom.annotations.domain_type:
                attachment_point = atom

        assert attachment_point

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
        backbone_to_placement = {}
        for i, atom in enumerate(backbone):
            if atom == stop_linearising:
                break

            if i == 0:
                backbone_to_placement[atom] = 'right'
            elif i == 1:
                backbone_to_placement[atom] = 'left'
            else:
                if atom not in backbone_to_placement:
                    previous_placement = backbone_to_placement[backbone[i - 1]]
                    ring = atom_to_ring.get(atom)
                    if ring and len(ring) == 2:
                        if ring[0] in backbone_to_placement:
                            backbone_to_placement[ring[1]] = backbone_to_placement[ring[0]]
                        else:
                            backbone_to_placement[atom] = self.get_opposite_placement(previous_placement)
                    else:
                        backbone_to_placement[atom] = self.get_opposite_placement(previous_placement)

        return backbone_to_placement

    def fix_rings(self, rings, backbone_atoms, backbone_to_placement):
        for ring in rings:
            if len(ring) == 1:
                atom_1 = ring[0]
                drawing_ring = self.get_ring(atom_1.draw.rings[0])
                angle = get_angle(atom_1.draw.position,
                                  drawing_ring.center)

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

    def linearise(self):

        attached_to_domain = False

        for atom in self.structure.graph:
            if atom.type == 'I':
                if atom.annotations.domain_type in ["Leader", "Follower", "ACP", "PCP"]:
                    attached_to_domain = True

        if self.ripp:
            label_nrp_central_chain(self.structure, is_ripp=True)

        if attached_to_domain:
            if self.ripp:
                backbone_atoms = find_central_chain_ripp(self.structure)
            else:
                backbone_atoms = find_central_chain(self.structure)
        else:
            backbone_atoms = find_central_chain_not_attached(self.structure)

        if attached_to_domain:

            self.place_top_atom(backbone_atoms)

        self.structure.refresh_structure()

        backbone, rings, atom_to_ring, stop_linearising = reorder_central_chain(backbone_atoms, self)
        backbone_to_placement = self.get_placements(backbone, atom_to_ring, stop_linearising)

        i = 0
        print(backbone_atoms)
        print(backbone)

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
                        print(i, ring)
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

                            print(sidechain_atoms)

                            for ring_atom in ring_atoms:
                                ring_atom.draw.position.mirror_about_line(ring[0].draw.position, ring[1].draw.position)
                            for sidechain_atom in sidechain_atoms:
                                sidechain_atom.draw.position.mirror_about_line(ring[0].draw.position, ring[1].draw.position)
                        self.rotate_subtree_ring(atom_2, set(ring_atoms + [atom_1]), required_rotation_rad,
                                                 atom_1.draw.position)

            i += 1

        self.fix_rings(rings, backbone, backbone_to_placement)
        self.position_sidechains(backbone, backbone_to_placement)

        self.resolve_primary_overlaps()
        self.total_overlap_score, sorted_overlap_scores, atom_to_scores = self.get_overlap_score()
        central_chain_bonds = set()

        for bond in self.structure.bonds.values():
            if bond.atom_1.annotations.in_central_chain or \
                    bond.atom_2.annotations.in_central_chain:
                central_chain_bonds.add(bond)

        self.finetune_overlap_resolution(
            masked_bonds=central_chain_bonds, highest_atom=backbone[0])

        self.resolve_secondary_overlaps(sorted_overlap_scores)

        if self.horizontal:
            self.rotate_structure(-1.5707)

    def resolve_secondary_overlaps(self, sorted_scores: List[Tuple[float, Atom]]) -> None:
        for score, atom in sorted_scores:
            if score > self.options.overlap_sensitivity:
                if len(atom.drawn_neighbours) <= 1:
                    if atom.drawn_neighbours and atom.drawn_neighbours[0].adjacent_to_stereobond():
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

    def position_sidechains(self, backbone, backbone_to_placement):

        sidechain_to_placement = {}
        sidechain_to_neighbour = {}

        for backbone_atom in backbone:
            for neighbour in backbone_atom.neighbours:
                if neighbour.type != 'H' and not neighbour.annotations.domain_type:
                    if neighbour not in backbone and not neighbour.inside_ring and neighbour not in sidechain_to_placement:

                        sidechain_to_placement[neighbour] = backbone_to_placement[backbone_atom]
                        sidechain_to_neighbour[neighbour] = backbone_atom

        for atom, placement in sidechain_to_placement.items():
            backbone_atom = sidechain_to_neighbour[atom]
            current_angle = get_angle(backbone_atom.draw.position, atom.draw.position)
            if placement == 'right':
                desired_angle = 180.0
            elif placement == 'left':
                desired_angle = 0.0
            else:
                raise ValueError("Placement must be 'left' or 'right'.")

            required_rotation_deg = desired_angle - current_angle
            required_rotation_rad = math.radians(required_rotation_deg)

            self.rotate_subtree(atom, backbone_atom, required_rotation_rad, backbone_atom.draw.position)

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


def get_angle(vector1, vector2):
    difference = Vector.subtract_vectors(vector1, vector2)
    difference_angle = difference.angle()
    return math.degrees(difference_angle)
