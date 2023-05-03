import math

from pikachu.reactions.functional_groups import find_bonds
from pikachu.general import read_smiles
from pikachu.math_functions import Vector


from raichu.data.molecular_moieties import PEPTIDE_BOND
from raichu.data.attributes import AMINOACID_ONE_LETTER_TO_NAME, AMINOACID_ONE_LETTER_TO_SMILES
from raichu.reactions.general_tailoring_reactions import proteolytic_cleavage, cyclisation
from raichu.tailoring_enzymes import TailoringEnzyme
from raichu.drawing.drawer import RaichuDrawer
from raichu.drawing.colours import AMINO_ACID_FILL_COLOURS, AMINO_ACID_OUTLINE_COLOURS
from raichu.attach_to_domain import attach_to_follower_ripp, attach_to_leader_ripp
from raichu.central_chain_detection.label_central_chain import label_nrp_central_chain

def make_circle(x_coord, y_coord, size, amino_acid):
    """Easy function to draw circle for the domain visualization. Returns
    matplotlib.patches.circle object

    x_coord: int, x-coordinate of the center of the circle to be drawn
    amino_acid : amino acid to be visualized
    size: radius of circle
    """
    colour = AMINO_ACID_FILL_COLOURS[amino_acid]
    outline_colour = AMINO_ACID_OUTLINE_COLOURS[amino_acid]
    circle = f"""<circle r="{size}" cx="{int(x_coord)}" cy="{int(y_coord)}" stroke="{outline_colour}"
        fill="{colour}" stroke-width="1"/>"""
    return circle


class RiPP_Cluster:
    def __init__(self, gene_name_precursor: str, full_amino_acid_sequence: str, amino_acid_sequence_for_structure_prediction: str, cleavage_sites: list() = None, macrocyclisations: list() = None, tailoring_enzymes_representation=None) -> None:
        self.gene_name = gene_name_precursor
        self.full_amino_acid_sequence = full_amino_acid_sequence.upper()
        self.amino_acid_seqence__for_structure_prediction = amino_acid_sequence_for_structure_prediction.upper()
        self.cleavage_sites = cleavage_sites
        self.cleavage_bonds = []
        self.macrocyclisations = macrocyclisations
        self.tailoring_enzymes_representation = tailoring_enzymes_representation
        self.chain_intermediate = None
        self.linear_product = None
        self.tailored_product = None
        self.cyclised_product = None
        self.final_product = None
        self.initialized_macrocyclization_atoms = []

    def make_peptide(self):
        smiles_peptide_chain = ""
        for amino_acid in self.amino_acid_seqence__for_structure_prediction:
            if amino_acid in AMINOACID_ONE_LETTER_TO_SMILES:
                substrate = AMINOACID_ONE_LETTER_TO_SMILES[amino_acid]
            else:
                raise ValueError(f"Unknown amino acid: {amino_acid}")
            smiles_peptide_chain += str(substrate)
        smiles_peptide_chain += "O"
        self.linear_product = read_smiles(smiles_peptide_chain)
        label_nrp_central_chain(
            self.linear_product, module_type='elongation')
        self.chain_intermediate = self.linear_product

    def initialize_cleavage_sites_on_structure(self) -> list:
        for cleavage_site in self.cleavage_sites:
            amino_acid_cleavage = cleavage_site.position_amino_acid
            number_cleavage = cleavage_site.position_index
            if self.amino_acid_seqence__for_structure_prediction[number_cleavage-1] == amino_acid_cleavage:
                peptide_bonds = find_bonds(PEPTIDE_BOND, self.linear_product)
                peptide_bonds = sorted(peptide_bonds, key=lambda bond: bond.nr)
                cleavage_bond = peptide_bonds[number_cleavage]
                self.cleavage_bonds += [[cleavage_bond,
                    cleavage_site.structure_to_keep]]
            else:
                raise ValueError(
                    f"No {AMINOACID_ONE_LETTER_TO_NAME[amino_acid_cleavage]} in position {number_cleavage} for cleavage.")

    def do_proteolytic_claevage(self):
        self.initialize_cleavage_sites_on_structure()
        for bond, structure_to_keep in self.cleavage_bonds:
            self.chain_intermediate = proteolytic_cleavage(
                bond, self.chain_intermediate, structure_to_keep=structure_to_keep)
        self.final_product = self.chain_intermediate

    def initialize_macrocyclization_on_structure(self) -> list(list()):
        if self.macrocyclisations:
            for cyclization in self.macrocyclisations:
                atoms = [atom for atom in self.chain_intermediate.atoms.values() if str(
                    atom) in [cyclization.atom1, cyclization.atom2]]
                self.initialized_macrocyclization_atoms += [atoms]

    def do_macrocyclization(self):
        self.initialize_macrocyclization_on_structure()
        for index, macrocyclization_atoms in enumerate(self.initialized_macrocyclization_atoms):
            if [str(atom) for atom in macrocyclization_atoms] != [self.macrocyclisations[index].atom1, self.macrocyclisations[index].atom2]:
                raise ValueError(
                    f'Not all atoms {str(self.initialized_macrocyclization_atoms)} for macrocyclisation exist in the structure.')
            atom1 = self.chain_intermediate.get_atom(macrocyclization_atoms[0])
            atom2 = self.chain_intermediate.get_atom(macrocyclization_atoms[1])
            self.chain_intermediate = cyclisation(
                self.chain_intermediate, atom1, atom2)
        self.cyclised_product = self.chain_intermediate

    def initialize_modification_sites_on_structure(self, modification_sites):
        modification_sites_initialized = []
        for atoms_for_reaction in modification_sites:
            atoms_for_reaction_with_numbers = map(
                lambda atom: [int(''.join(filter(str.isdigit, atom))), atom], atoms_for_reaction)
            atoms_in_structure = list(map(
                str, self.chain_intermediate.atoms.values()))
            atoms_for_reaction_initialized = [self.chain_intermediate.atoms[atom[0]]
                                                for atom in atoms_for_reaction_with_numbers if atom[1] in atoms_in_structure]
            atoms_for_reaction_initialized = list(
                filter(lambda atom: atom is not None, atoms_for_reaction_initialized))
            modification_sites_initialized += [
                atoms_for_reaction_initialized]
        return modification_sites_initialized

    def do_tailoring(self):
        if self.tailoring_enzymes_representation:
            for tailoring_enzyme_representation in self.tailoring_enzymes_representation:
                modification_sites = self.initialize_modification_sites_on_structure(
                    tailoring_enzyme_representation.modification_sites)
                if [[str(atom) for atom in atoms_for_reaction]for atoms_for_reaction in modification_sites] != tailoring_enzyme_representation.modification_sites:
                    raise ValueError(
                        f'Not all atoms {tailoring_enzyme_representation.modification_sites} for {tailoring_enzyme_representation.type} exist in the structure.')
                tailoring_enzyme = TailoringEnzyme(
                    tailoring_enzyme_representation.gene_name, tailoring_enzyme_representation.type, modification_sites, tailoring_enzyme_representation.substrate)
                self.tailored_product = tailoring_enzyme.do_tailoring(
                    self.chain_intermediate)
                self.chain_intermediate = self.tailored_product

    def draw_product(self, as_string=True, out_file=None, draw_straightened = True):
        assert self.chain_intermediate
        drawing = RaichuDrawer(self.chain_intermediate, dont_show=True, add_url=True,
                               draw_Cs_in_pink=False, draw_straightened=draw_straightened, horizontal=False)
        drawing.draw_structure()
        svg_string = drawing.save_svg_string()
        if as_string:
            return svg_string
        else:
            if out_file is None:
                raise ValueError(
                    "Must provide output svg directory if 'as_string' is set to False.")
            else:
                with open(out_file, 'w') as svg_out:
                    svg_out.write(svg_string)

    def draw_precursor(self, fold=10, size=14, as_string=True, out_file=None, amino_acid_sequence=None, leader=True, x_translation=0, y_translation=0, min_padding_y =0):
        if amino_acid_sequence == None:
            amino_acid_sequence = self.full_amino_acid_sequence
        # set begin of chain so that last amino acid is forward
        circles = []
        texts = []
        text_colour = "#000000"
        if leader:
            current_y = min_padding_y+size
            current_x = min([fold,len(amino_acid_sequence)])*size*2+size*0.5
            forward = False
            # reverse sequence to go backwards
            amino_acid_sequence = amino_acid_sequence[::-1]
            for index, amino_acid in enumerate(amino_acid_sequence):
                # TODO: create dictionary for amino acid colors
                circles.append(make_circle(current_x, current_y, size, amino_acid))
                text = f"""<text x="{current_x}" y="{current_y}" fill="{text_colour}" text-anchor="middle" font-family="verdana" font-size = "{10}">\
                    <tspan y="{current_y}" dy="0.35em">{amino_acid}</tspan></text>"""
                texts.append(text)
                step = size*2
                if index % fold == 0 and index != 0:
                    forward = not forward
                if (index % fold > fold-2 or index % fold < 1) and index > 2:
                    step = 2**(1/2) * step/2
                    current_y += step
                if forward:
                    current_x += step
                else:
                    current_x -= step

        if not leader:
            current_y = size+y_translation
            current_x = size+x_translation
            forward = True
            for index, amino_acid in enumerate(amino_acid_sequence):
                circles.append(make_circle(current_x, current_y, size, amino_acid))
                text = f"""<text x="{current_x}" y="{current_y}" fill="{text_colour}" text-anchor="middle" font-family="verdana" font-size = "{10}">\
                    <tspan y="{current_y}" dy="0.35em">{amino_acid}</tspan></text>"""
                texts.append(text)
                step = size*2
                if index % fold == 0 and index != 0:
                    forward = not forward
                if (index % fold > fold-2 or index % fold < 1) and index > 2:
                    step = 2**(1/2) * step/2
                    current_y += step
                if forward:
                    current_x += step
                else:
                    current_x -= step

        svg = ''
        svg += f"""<g id="domain_circles_{leader}">\n"""
        for i, circle in enumerate(circles):
            text = texts[i]
            svg += f"""<g id="domain_bubble_{i}_{leader}">\n"""
            svg += f"{circle}\n"
            svg += f"{text}\n"
            svg += "</g>\n"
        svg += "</g>\n"
        if as_string:
            return svg
        else:
            svg_string = f"""<svg width="{fold*size*2+size*2}" height="{math.ceil(len(amino_acid_sequence)/fold)*size*(2+2**0.5)+size}" viewBox="{0} {0} {fold*size*2+size} {math.ceil(len(amino_acid_sequence)/fold)*size*3+size}" xmlns="http://www.w3.org/2000/svg">\n"""
            svg_string += svg
            svg_string += "</svg>"
            if out_file is None:
                raise ValueError(
                    "Must provide output svg directory if 'as_string' is set to False.")
            else:
                with open(out_file, 'w') as svg_out:
                    svg_out.write(svg_string)

    def draw_precursor_with_modified_product(self, fold=10, size=14, as_string=True, out_file=None, add_url=True, draw_Cs_in_pink=False):
        leader_pos = Vector(0, 0)
        follower_pos = Vector(0, 0)
        amino_acid_sequence_without_core= self.full_amino_acid_sequence.split(
            self.amino_acid_seqence__for_structure_prediction)
        if len(amino_acid_sequence_without_core)!=2:
            raise ValueError("Core peptide not in precursor.")
        [amino_acid_sequence_leader,
        amino_acid_sequence_follower] = amino_acid_sequence_without_core
        structure = self.chain_intermediate.deepcopy()
        if len(amino_acid_sequence_follower) > 0:
            structure = attach_to_follower_ripp(
                structure)
        if len(amino_acid_sequence_leader) > 0:
            structure = attach_to_leader_ripp(
                structure)
        drawing = RaichuDrawer(structure, dont_show=True, add_url=add_url,
                                draw_straightened=True, horizontal=True, draw_Cs_in_pink=draw_Cs_in_pink)
        drawing.flip_y_axis()
        drawing.move_to_positive_coords()
        drawing.convert_to_int()
        drawing.draw_structure()
        svg_style = drawing.svg_style
        x_translation_leader = 0
        y_translation_leader = 0
        max_x = 0
        max_y = 0
        for atom in drawing.structure.graph:
            if atom.draw.position.x > max_x:
                max_x = atom.draw.position.x
            if atom.draw.position.y > max_y:
                max_y = atom.draw.position.y
            if atom.annotations.domain_type:
                if atom.annotations.domain_type == "Follower":
                    follower_pos = atom.draw.position
                if atom.annotations.domain_type == "Leader":
                    leader_pos = atom.draw.position
        
        svg_bubbles_leader = self.draw_precursor(fold=fold, size = size, as_string=True, amino_acid_sequence=amino_acid_sequence_leader, leader=True, min_padding_y= max(leader_pos.y, follower_pos.y))
        position_x_first_bubble_leader = min([fold,len(amino_acid_sequence_leader)])*size*2+size*0.5
        position_y_first_bubble_leader =   max(leader_pos.y, follower_pos.y)+size
        x_translation_leader = position_x_first_bubble_leader - leader_pos.x + 10 - 0.2
        if leader_pos.x != 0:
            y_translation_leader = position_y_first_bubble_leader - leader_pos.y
        else:
            y_translation_leader = position_y_first_bubble_leader - follower_pos.y
        drawing.move_structure(x_translation_leader, y_translation_leader)
        svg = drawing.draw_svg()
        x_translation = position_x_first_bubble_leader - leader_pos.x + 10 - 0.2 + max_x + x_translation_leader
        y_translation = -position_y_first_bubble_leader
        if follower_pos.x != 0:
            y_translation = follower_pos.y -size
            x_translation = follower_pos.x
            svg_bubbles_follower = self.draw_precursor(
                fold=fold, size=size, as_string=True, amino_acid_sequence=amino_acid_sequence_follower, leader=False, x_translation=x_translation, y_translation=y_translation)
        else:
            svg_bubbles_follower = ""


        x1 = 0
        x2 = x_translation + 2 * size * min(len(amino_acid_sequence_follower), fold) + 2*size
        y1 = 0
        y2 = max(max_y + y_translation_leader + 10, (math.ceil(len(amino_acid_sequence_leader)/fold) * size * (1+2**0.5) + max_y + 10),
                 y_translation + 2*size + math.ceil(len(amino_acid_sequence_follower)/fold) * size * (1+2**0.5))
        svg_string = f"""<svg width="{x2}" height="{y2}" viewBox="{x1} {y1} {x2} {y2}" xmlns="http://www.w3.org/2000/svg">\n"""
        if svg_style:
            svg_string += f"{svg_style}\n"
        svg_string += svg_bubbles_leader
        if follower_pos.x != 0:
            svg_string += svg_bubbles_follower
        svg_string += svg
        # # add connections
        # if len(amino_acid_sequence_leader) > 0:
        #     svg_string += f"""<line x1="{position_x_first_bubble_leader+size}" y1="{position_y_first_bubble_leader}" x2="{position_x_first_bubble_leader+size+5}" y2="{position_y_first_bubble_leader}" stroke = "black" />"""
        # if follower_pos:
        #     svg_string += f"""<line x1="{x_translation -3}" y1="{y_translation +size}" x2="{x_translation}" y2="{y_translation +size}" stroke = "black" />"""
        svg_string += "</svg>"

        if as_string:
            return svg_string
        else:
            if out_file is None:
                raise ValueError(
                    "Must provide output svg directory if 'as_string' is set to False.")
            else:
                with open(out_file, 'w') as svg_out:
                    svg_out.write(svg_string)
