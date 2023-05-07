from pikachu.drawing.drawing import Drawer

from raichu.drawing.drawer import RaichuDrawer
from raichu.drawing.pathway import draw_arrow_and_text
from raichu.data.attributes import ENZYME_TO_NAME
from raichu.tailoring_enzymes import TailoringEnzyme
from raichu.reactions.general_tailoring_reactions import cyclisation, oxidative_bond_formation


class Cluster:
    def __init__(self, tailoring_representations=None, macrocyclisation_representations=None) -> None:
        self.tailoring_representations = tailoring_representations
        self.macrocyclisation_representations = macrocyclisation_representations

        self.chain_intermediate = None
        self.tailored_product = None

        self.linear_product = None
        self.cyclic_product = None
        self.cyclic_products = []

        self.modular_intermediates = []
        self.cyclic_intermediates = []
        self.tailoring_intermediates = []

    def initialize_macrocyclization(self):
        initialized_macrocyclization_atoms = []
        if self.macrocyclisation_representations:
            for cyclization in self.macrocyclisation_representations:
                atom_1 = None
                atom_2 = None
                for atom in self.chain_intermediate.atoms.values():
                    if str(atom) == cyclization.atom_1:
                        atom_1 = atom
                    elif str(atom) == cyclization.atom_2:
                        atom_2 = atom
                atoms = [atom_1, atom_2]
                # atoms = [atom for atom in self.chain_intermediate.atoms.values()
                #          if str(atom) in [cyclization.atom_1, cyclization.atom_2]]
                initialized_macrocyclization_atoms.append((atoms, cyclization.type))
        return initialized_macrocyclization_atoms

    def do_macrocyclization(self, sequential=True):
        initialized_macrocyclization_atoms = self.initialize_macrocyclization()

        for index, macrocyclization_atoms in enumerate(initialized_macrocyclization_atoms):
            atoms, cyclisation_type = macrocyclization_atoms
            if [str(atom) for atom in atoms] != [self.macrocyclisation_representations[index].atom_1,
                                                 self.macrocyclisation_representations[index].atom_2]:
                raise ValueError(
                    f'Not all atoms {str(initialized_macrocyclization_atoms)} \
                    for macrocyclisation exist in the structure.')

            if sequential:
                structure = self.chain_intermediate
            else:
                structure = self.chain_intermediate.deepcopy()

            atom1 = structure.get_atom(atoms[0])
            atom2 = structure.get_atom(atoms[1])

            if cyclisation_type == 'oxidative':
                structure = oxidative_bond_formation(atom1, atom2, structure)

            elif cyclisation_type == 'condensative':
                structure = cyclisation(structure, atom1, atom2)
            else:
                raise ValueError("Cyclisation type must be oxidative or condensative.")

            if sequential:
                self.chain_intermediate = structure
                self.cyclic_intermediates.append(structure.deepcopy())
            else:
                self.cyclic_products.append(structure.deepcopy())

            self.cyclic_product = structure

    def initialize_modification_sites(self, modification_sites):
        modification_sites_initialized = []
        for atom_representations_for_reaction in modification_sites:
            atom_indices_and_types = map(
                lambda atom: [int(''.join(filter(str.isdigit, atom))), atom], atom_representations_for_reaction)
            atoms_in_structure = list(map(
                str, self.chain_intermediate.atoms.values()))
            atoms_for_reaction_initialized = [self.chain_intermediate.atoms[atom[0]]
                                              for atom in atom_indices_and_types if atom[1] in atoms_in_structure]
            atoms_for_reaction_initialized = list(
                filter(lambda atom: atom is not None, atoms_for_reaction_initialized))
            modification_sites_initialized += [
                atoms_for_reaction_initialized]
        return modification_sites_initialized

    def do_tailoring(self):
        if self.tailoring_representations:
            self.tailoring_intermediates.append(self.chain_intermediate.deepcopy())

            for tailoring_enzyme_representation in self.tailoring_representations:
                modification_sites = self.initialize_modification_sites(
                    tailoring_enzyme_representation.modification_sites)

                if [[str(atom) for atom in atoms_for_reaction] for atoms_for_reaction in modification_sites] != tailoring_enzyme_representation.modification_sites:
                    raise ValueError(
                        f'Not all atoms {tailoring_enzyme_representation.modification_sites} for \
                        {tailoring_enzyme_representation.type} exist in the structure.')

                tailoring_enzyme = TailoringEnzyme(tailoring_enzyme_representation.gene_name,
                                                   tailoring_enzyme_representation.type, modification_sites,
                                                   tailoring_enzyme_representation.substrate)
                self.tailored_product = tailoring_enzyme.do_tailoring(self.chain_intermediate)
                self.tailoring_intermediates.append(self.tailored_product.deepcopy())
                self.chain_intermediate = self.tailored_product

    def get_drawings_tailoring(self, whitespace=30):

        drawings = []
        widths = []
        centre_points = []
        max_height = -1000000000

        for i, structure in enumerate(self.tailoring_intermediates):

            drawing = Drawer(structure)

            drawing.flip_y_axis()
            drawing.move_to_positive_coords()
            drawing.convert_to_int()

            min_x = 100000000
            max_x = -100000000
            min_y = 100000000
            max_y = -100000000

            for atom in drawing.structure.graph:
                if atom.draw.positioned:
                    if atom.draw.position.x < min_x:
                        min_x = atom.draw.position.x
                    if atom.draw.position.x > max_x:
                        max_x = atom.draw.position.x
                    if atom.draw.position.y < min_y:
                        min_y = atom.draw.position.y
                    if atom.draw.position.y > max_y:
                        max_y = atom.draw.position.y

            width = max_x - min_x + whitespace
            height = max_y - min_y + whitespace
            if height > max_height:
                max_height = height

            widths.append(width)
            centre_points.append(((max_x + min_x) / 2, (max_y + min_y) / 2))
            drawings.append(drawing)

        print(len(drawings))

        return drawings, widths, max_height, centre_points

    def draw_tailoring(self, arrow_size=20, mode='gene_name', as_string=False, out_file=None):

        assert self.tailoring_representations
        assert mode in ['gene_name', 'reaction_name']
        drawings, widths, max_height, centre_points = self.get_drawings_tailoring()

        target_x = 0
        target_y = max_height / 2
        arrow_svgs = []
        structure_svgs = []
        padding = None
        svg_style = None

        for i, drawing in enumerate(drawings):
            enzyme = self.tailoring_representations[i]
            current_x, current_y = centre_points[i]
            width = widths[i]
            translation_x = (target_x + 0.5 * width) - current_x
            translation_y = target_y - current_y

            drawing.set_structure_id(f"s{i}")
            svg_style = drawing.svg_style

            target_x += width
            if i != len(drawings) - 1:

                arrow_start = target_x
                arrow_end = target_x + arrow_size
                arrow_y = target_y
                if mode == 'gene_name':
                    text = enzyme.gene_name
                elif mode == 'reaction_name':
                    text = ENZYME_TO_NAME[enzyme.type]
                else:
                    raise ValueError("Drawing mode must be 'gene_name' or 'reaction_name'.")

                arrow_svg = draw_arrow_and_text(arrow_start, arrow_end, arrow_y, i, text)
                arrow_svgs.append(arrow_svg)

            drawing.move_structure(translation_x, translation_y)
            structure_svg = drawing.draw_svg()
            structure_svgs.append(structure_svg)
            padding = drawing.options.padding

        assert svg_style and padding

        x1 = 0
        x2 = target_x + padding
        y1 = 0
        y2 = max_height + padding

        width = x2
        height = y2

        svg_string = f"""<svg width="{width}" height="{height}" viewBox="{x1} {y1} {x2} {y2}" xmlns="http://www.w3.org/2000/svg">\n"""
        if svg_style:
            svg_string += f"{svg_style}\n"
        svg_string += "</svg>"

        svg_string += f"""<g id="arrows">\n"""
        for arrow in arrow_svgs:
            svg_string += arrow
        svg_string += "</g>\n"

        for string in structure_svgs:
            svg_string += string

        if as_string:
            return svg_string
        else:
            if out_file is None:
                raise ValueError("Must provide output svg path if 'as_string' is set to False.")
            else:
                with open(out_file, 'w') as svg_out:
                    svg_out.write(svg_string)

    def draw_product(self, as_string=False, out_file=None):
        assert self.chain_intermediate
        drawing = RaichuDrawer(self.chain_intermediate, dont_show=True, add_url=True, draw_Cs_in_pink=False,
                               make_linear=False)
        drawing.draw_structure()
        svg_string = drawing.save_svg_string()

        if as_string:
            return svg_string
        else:
            if out_file is None:
                raise ValueError("Must provide output svg directory if 'as_string' is set to False.")
            else:
                with open(out_file, 'w') as svg_out:
                    svg_out.write(svg_string)
