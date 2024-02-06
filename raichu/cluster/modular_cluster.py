from typing import List
import os

from pikachu.drawing.drawing import Drawer
from pikachu.reactions.functional_groups import find_atoms

from raichu.reactions.chain_release import cyclic_release, find_o_betapropriolactone
from raichu.drawing.drawer import RaichuDrawer
from raichu.module import _Module
from raichu.drawing.bubbles import draw_bubbles
from raichu.substrate import PKSSubstrate
from raichu.data.trans_at import TRANSATOR_CLADE_TO_TAILORING_REACTIONS, TRANSATOR_CLADE_TO_STARTER_SUBSTRATE, \
    TRANSATOR_CLADE_TO_ELONGATING
from raichu.domain.domain import TailoringDomain
from raichu.data.molecular_moieties import O_OH, N_NH, KETO_GROUP_C
from raichu.cluster.base_cluster import Cluster
from raichu.representations import TailoringRepresentation, MacrocyclizationRepresentation


class ModularCluster(Cluster):
    def __init__(self, modules: List[_Module], tailoring_representations: List[TailoringRepresentation] = None) -> None:
        super().__init__(tailoring_representations)
        self.modules = modules
        self.tailoring_enzymes = tailoring_representations
        self.chain_intermediate = None

        self.tailored_product = None
        self.handle_transat()
        self.cyclisation_type = None
        self.modular = True

    def handle_transat(self):
        for i, module in enumerate(self.modules):
            if module.type.name == "PKS" and module.subtype.name == "PKS_TRANS":
                substrate = PKSSubstrate("MALONYL_COA")
                if i < len(self.modules) - 1:
                    j = i + 1
                    next_module = self.modules[j]
                    if next_module.is_broken:
                        while next_module.is_broken:
                            next_module = self.modules[j]
                            j += 1
                    if next_module.type.name == "PKS":
                        if module.is_starter_module and next_module.subtype.name == "PKS_TRANS":
                            if next_module.synthesis_domain:
                                if next_module.synthesis_domain.subtype and next_module.synthesis_domain.subtype.name != "UNKNOWN":
                                    substrate_name = TRANSATOR_CLADE_TO_STARTER_SUBSTRATE.get(next_module.synthesis_domain.subtype.name)
                                    if substrate_name is not None:
                                        substrate = PKSSubstrate(substrate_name)
                                    else:
                                        substrate = PKSSubstrate("ACETYL_COA")
                                else:
                                    substrate = PKSSubstrate("ACETYL_COA")
                            else:
                                substrate = PKSSubstrate("ACETYL_COA")
                                    # TODO: Transfer names of substrates to a dictionary in raichu.substrate
                        elif module.is_starter_module:
                            substrate = PKSSubstrate("ACETYL_COA")
                        if next_module.synthesis_domain:
                            if not module.is_termination_module and next_module.subtype.name == "PKS_TRANS" and next_module.synthesis_domain.subtype and next_module.synthesis_domain.subtype.name != "UNKNOWN" and next_module.synthesis_domain.subtype.name != "MISCELLANEOUS" and next_module.synthesis_domain.subtype.name != "NON_ELONGATING":
                                # Ignore tailoring domains in the module itself
                                module.tailoring_domains = []
                                for dummy_domain_type, dummy_domain_subtype in TRANSATOR_CLADE_TO_TAILORING_REACTIONS[next_module.synthesis_domain.subtype.name]:
                                    self.modules[i].tailoring_domains.append(TailoringDomain(dummy_domain_type,
                                                                                             dummy_domain_subtype))
                if module.synthesis_domain:
                    if module.synthesis_domain.subtype and module.synthesis_domain.subtype.name != "UNKNOWN" and module.synthesis_domain.subtype.name != "MISCELLANEOUS":
                        self.modules[i].synthesis_domain.is_elongating = TRANSATOR_CLADE_TO_ELONGATING[module.synthesis_domain.subtype.name]
                self.modules[i].recognition_domain.substrate = substrate

    def compute_structures(self, compute_cyclic_products=True):
        for module in self.modules:
            structure = module.run_module(self.chain_intermediate)
            self.modular_intermediates.append(structure.deepcopy())
            self.chain_intermediate = structure
            if module.is_termination_module:
                self.linear_product = module.release_chain(structure)
                self.chain_intermediate = self.linear_product
                if module.termination_domain.type.name in ["TE", 'DUMMY_TE']:
                    self.cyclisation_type = 'condensative'
                elif module.termination_domain.type.name in ["TD", "DUMMY_TD"]:
                    self.cyclisation_type = 'oxidative'
                else:
                    raise ValueError(f"RAIChU does not recognise termination domain type {module.termination_domain.type.name}.")
                break
        else:
            raise ValueError("Cluster must contain at least one termination module.")

        if compute_cyclic_products:
            self.cyclise_all()

    def cyclise(self, atom):

        cyclised_product = cyclic_release(self.chain_intermediate, atom)
        self.cyclic_product = cyclised_product
        return cyclised_product

    def cyclise_all(self):

        self.macrocyclisation_representations = []
        self.cyclised_products = []

        o_not_to_use = find_o_betapropriolactone(self.chain_intermediate, release_type=self.cyclisation_type)

        candidate_target_atoms = []
        terminal_o = None
        terminal_c = None

        for atom in find_atoms(N_NH, self.chain_intermediate):
            if len(atom.get_neighbours('H')) == 2:
                candidate_target_atoms.append(str(atom))

        if self.cyclisation_type == 'condensative':
            for atom in find_atoms(O_OH, self.chain_intermediate):

                if atom.has_neighbour('H') and atom != o_not_to_use and not getattr(atom.annotations, 'terminal_o'):
                    candidate_target_atoms.append(str(atom))
                elif getattr(atom.annotations, 'terminal_o'):
                    terminal_o = str(atom)

            assert terminal_o
            for atom in candidate_target_atoms:
                macrocyclisation_representation = MacrocyclizationRepresentation(atom, terminal_o, 'condensative')
                self.macrocyclisation_representations.append(macrocyclisation_representation)

        elif self.cyclisation_type == 'oxidative':
            for atom in find_atoms(O_OH, self.chain_intermediate):
                if atom.has_neighbour('H') and atom != o_not_to_use:
                    candidate_target_atoms.append(str(atom))

            for atom in find_atoms(KETO_GROUP_C, self.chain_intermediate):
                if atom.has_neighbour('H') and getattr(atom.annotations, 'terminal_c'):
                    terminal_c = str(atom)

            assert terminal_c

            for atom in candidate_target_atoms:
                macrocyclisation_representation = MacrocyclizationRepresentation(atom, terminal_c, 'oxidative')
                self.macrocyclisation_representations.append(macrocyclisation_representation)

        else:
            raise ValueError(f"Macrocyclisation must be oxidative or condensative. '{self.cyclisation_type}' given.")

        self.do_macrocyclization(sequential=False)

    def draw_all_products(self, out_dir):
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)

        last_i = 0

        for i, product in enumerate(self.cyclic_products):
            out_file = os.path.join(out_dir, f'product_{i}.svg')
            drawing = Drawer(product)
            drawing.write_svg(out_file)
            last_i = i

        drawing = Drawer(self.linear_product)
        out_file = os.path.join(out_dir, f'product_{last_i + 1}.svg')
        drawing.write_svg(out_file)

    # TODO: Use PIKAChU's SVG drawer
    # def draw_spaghettis(self):
    #     spaghetti_svgs = []
    #     for structure in self.modular_intermediates:
    #         drawing = RaichuDrawer(structure, dont_show=True)
    #         drawing.draw_structure()
    #         svg_string = drawing.save_svg_string()
    #         spaghetti_svgs.append(svg_string)
    #
    #     linear_drawing = Drawer(self.linear_product)
    #     linear_svg = linear_drawing.save_svg_string()
    #
    #     return spaghetti_svgs + [linear_svg]

    def get_spaghettis(self, whitespace=30):
        drawings = []
        widths = []

        for i, structure in enumerate(self.modular_intermediates):

            drawing = RaichuDrawer(structure, dont_show=True)

            drawing.flip_y_axis()
            drawing.finetune_overlap_bubbles()
            drawing.move_to_positive_coords()
            drawing.convert_to_int()

            carrier_domain_pos = None

            for atom in drawing.structure.graph:
                if atom.annotations.domain_type:
                    carrier_domain_pos = atom.draw.position
                    atom.draw.positioned = False

            assert carrier_domain_pos

            min_x = 100000000
            max_x = -100000000

            for atom in drawing.structure.graph:
                if atom.draw.positioned:
                    if atom.draw.position.x < min_x:
                        min_x = atom.draw.position.x
                    if atom.draw.position.x > max_x:
                        max_x = atom.draw.position.x

            width = (carrier_domain_pos.x - min_x + 0.5 * whitespace,
                     max_x - carrier_domain_pos.x + 0.5 * whitespace)
            widths.append(width)
            drawings.append(drawing)

        return drawings, widths

    def draw_cluster(self, as_string=True, out_file=None, colour_by_module=True):
        drawings, widths = self.get_spaghettis()
        bubble_svg, bubble_positions, last_domain_coord = draw_bubbles(self, widths)
        min_x = 100000000
        max_x = -100000000
        min_y = 100000000
        max_y = -100000000

        svg_strings = []
        squiggly_svgs = []
        padding = None
        svg_style = None

        for i, drawing in enumerate(drawings):
            drawing.set_structure_id(f"s{i}")
            drawing.set_annotation_for_grouping('module_nr')
            if colour_by_module:
                drawing.colour_by_annotation()
            svg_style = drawing.svg_style

            padding = drawing.options.padding

            carrier_domain_pos = None

            for atom in drawing.structure.graph:
                if atom.annotations.domain_type:
                    carrier_domain_pos = atom.draw.position
                    atom.draw.positioned = False

            assert carrier_domain_pos
            bubble_x, bubble_y = bubble_positions[i]
            x_translation = bubble_x - carrier_domain_pos.x
            y_translation = bubble_y - carrier_domain_pos.y

            drawing.move_structure(x_translation, y_translation + 15)
            svg = drawing.draw_svg(annotation='module_nr')
            svg_strings.append(svg)

            sulphur_pos = None
            carrier_pos = None

            for atom in drawing.structure.graph:
                if atom.draw.positioned:
                    if atom.draw.position.x < min_x:
                        min_x = atom.draw.position.x
                    if atom.draw.position.y < min_y:
                        min_y = atom.draw.position.y
                    if atom.draw.position.x > max_x:
                        max_x = atom.draw.position.x
                    if atom.draw.position.y > max_y:
                        max_y = atom.draw.position.y
                if atom.annotations.domain_type:
                    sulphur_pos = atom.get_neighbour('S').draw.position
                    carrier_pos = atom.draw.position

            squiggly_svg = f'<path d="M {sulphur_pos.x} {sulphur_pos.y - 5} Q {sulphur_pos.x - 5} {sulphur_pos.y - (sulphur_pos.y - 5 - carrier_pos.y)/2}, {carrier_pos.x} {sulphur_pos.y - 5 - (sulphur_pos.y - 5 - carrier_pos.y)/2} T {carrier_pos.x} {carrier_pos.y}" stroke="grey" fill="white"/>'
            squiggly_svgs.append(squiggly_svg)

        assert padding is not None

        x1 = 0
        x2 = max([max_x + padding, last_domain_coord + padding])
        y1 = 0
        y2 = max_y + padding

        width = x2
        height = y2

        svg_string = f"""<svg width="{width}" height="{height}" viewBox="{x1} {y1} {x2} {y2}" xmlns="http://www.w3.org/2000/svg">\n"""
        if svg_style:
            svg_string += f"{svg_style}\n"
        svg_string += bubble_svg
        for string in svg_strings:
            svg_string += string
        for squiggly_svg in squiggly_svgs:
            svg_string += squiggly_svg
        svg_string += "</svg>"

        if as_string:
            return svg_string
        else:
            if out_file is None:
                raise ValueError("Must provide output svg directory if 'as_string' is set to False.")
            else:
                with open(out_file, 'w') as svg_out:
                    svg_out.write(svg_string)
