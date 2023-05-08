import math
from typing import List
import os

from pikachu.reactions.functional_groups import find_bonds
from pikachu.general import read_smiles
from pikachu.math_functions import Vector

from raichu.data.molecular_moieties import PEPTIDE_BOND
from raichu.data.attributes import AMINOACID_ONE_LETTER_TO_NAME, AMINOACID_ONE_LETTER_TO_SMILES
from raichu.reactions.general_tailoring_reactions import proteolytic_cleavage
from raichu.drawing.drawer import RaichuDrawer
from raichu.attach_to_domain import attach_to_follower_ripp, attach_to_leader_ripp
from raichu.central_chain_detection.label_central_chain import label_nrp_central_chain
from raichu.representations import CleavageSiteRepresentation, MacrocyclizationRepresentation, TailoringRepresentation
from raichu.cluster.base_cluster import Cluster
from raichu.drawing.ripp_drawing import make_circle


class RiPPCluster(Cluster):
    def __init__(self, gene_name_precursor: str, full_peptide: str,
                 core_peptide: str, cleavage_sites: List = None,
                 macrocyclisations: List = None, tailoring_representations=None) -> None:

        super().__init__(tailoring_representations, macrocyclisations)
        self.gene_name_precursor = gene_name_precursor
        self.full_peptide = full_peptide.upper()
        self.core_peptide = core_peptide.upper()
        self.cleavage_sites = cleavage_sites

        self.cleavage_bonds = []

        self.linear_product = None
        self.final_product = None

        self.initialized_macrocyclization_atoms = []
        self.ripp = True

    def write_cluster(self, out_dir):
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)
        cluster_out = os.path.join(out_dir, 'cluster.txt')
        with open(cluster_out, 'w') as out:
            out.write(f"gene_name\t{self.gene_name_precursor}\n")
            out.write(f"full_peptide\t{self.full_peptide}\n")
            out.write(f"core_peptide\t{self.core_peptide}\n")
            cleavage_strings = []
            for cleavage_site in self.cleavage_sites:
                cleavage_strings.append(f"{cleavage_site.amino_acid}|{cleavage_site.amino_acid_index}|{cleavage_site.peptide_fragment_to_keep}")
            cyclic_strings = []
            for cyclic in self.macrocyclisation_representations:
                cyclic_strings.append(f"{cyclic.atom_1}|{cyclic.atom_2}")
            out.write(f"cleavage_sites\t{':'.join(cleavage_strings)}")
            out.write(f"macrocyclisations\t{':'.join(cyclic_strings)}")

        if self.tailoring_representations:
            tailoring_out = os.path.join(out_dir, 'tailoring.txt')
            with open(tailoring_out, 'w') as tailoring:
                tailoring.write('gene_name\ttype\tsubstrate\tmodification_sites\n')

                for enzyme in self.tailoring_representations:

                    if enzyme.modification_sites:
                        site_reprs = []
                        for modification_site in enzyme.modification_sites:

                            site_repr = '|'.join(list(map(str, modification_site)))
                            site_reprs.append(site_repr)

                        site_str = ':'.join(site_reprs)
                    else:
                        site_str = str(None)

                    tailoring.write(f"{enzyme.gene_name}\t{enzyme.type}\t{enzyme.substrate}\t{site_str}\n")

    @classmethod
    def from_file(cls, in_dir):
        in_cluster = os.path.join(in_dir, 'cluster.txt')
        in_tailoring = os.path.join(in_dir, 'tailoring.txt')
        gene_name = None
        full_peptide = None
        core_peptide = None
        cleavage_sites = []
        macrocyclisations = []

        with open(in_cluster, 'r') as cluster:
            for line in cluster:
                line = line.strip()
                if line:
                    category, value = line.split('\t')
                    if category == 'gene_name':
                        gene_name = value
                    elif category == 'full_peptide':
                        full_peptide = value
                    elif category == 'core_peptide':
                        core_peptide = value
                    elif category == 'cleavage_sites':
                        for cleavage_string in value.split(':'):
                            aa, aa_index, fragment_to_keep = cleavage_string.split('|')
                            cleavage_sites.append(CleavageSiteRepresentation(aa, int(aa_index), fragment_to_keep))
                    elif category == 'macrocyclisations':
                        for cyclic_string in value.split(':'):
                            atom_1, atom_2 = cyclic_string.split('|')
                            macrocyclisations.append(MacrocyclizationRepresentation(atom_1, atom_2))

        if os.path.exists(in_tailoring):
            tailoring_enzymes = []
            with open(in_tailoring, 'r') as tailoring:
                tailoring.readline()
                for line in tailoring:
                    line_info = line.split('\t')
                    line_info_cleaned = []
                    for i, entry in enumerate(line_info):
                        if line_info[i] == "None":
                            line_info_cleaned.append(None)
                        else:
                            line_info_cleaned.append(line_info[i])

                    gene_name, enzyme_type, substrate, site_str = line_info_cleaned

                    sites = []

                    if site_str:
                        site = site_str.split(':')
                        sites.append(site.split('|'))
                    tailoring_representation = TailoringRepresentation(gene_name, enzyme_type, sites, substrate)
                    tailoring_enzymes.append(tailoring_representation)

        else:
            tailoring_enzymes = None

        if not macrocyclisations:
            macrocyclisations = None
        if not cleavage_sites:
            cleavage_sites = None

        return cls(gene_name, full_peptide, core_peptide, cleavage_sites, macrocyclisations, tailoring_enzymes)

    def make_peptide(self):
        smiles_peptide_chain = ""
        for amino_acid in self.core_peptide:
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

    def initialize_cleavage_sites(self) -> None:
        if self.cleavage_sites:
            for cleavage_site in self.cleavage_sites:
                amino_acid_cleavage = cleavage_site.amino_acid
                number_cleavage = cleavage_site.amino_acid_index
                if self.core_peptide[number_cleavage - 1] == amino_acid_cleavage:
                    peptide_bonds = find_bonds(PEPTIDE_BOND, self.linear_product)
                    peptide_bonds = sorted(peptide_bonds, key=lambda bond: bond.nr)
                    cleavage_bond = peptide_bonds[number_cleavage]
                    self.cleavage_bonds += [[cleavage_bond,
                                             cleavage_site.peptide_fragment_to_keep]]
                else:
                    raise ValueError(
                        f"No {AMINOACID_ONE_LETTER_TO_NAME[amino_acid_cleavage]} in position {number_cleavage} for cleavage.")

    def do_proteolytic_cleavage(self):
        self.initialize_cleavage_sites()
        self.cleaved_intermediates.append(self.chain_intermediate.deepcopy())
        for bond, structure_to_keep in self.cleavage_bonds:
            self.chain_intermediate = proteolytic_cleavage(
                bond, self.chain_intermediate, structure_to_keep=structure_to_keep)
            self.cleaved_intermediates.append(self.chain_intermediate.deepcopy())
        self.final_product = self.chain_intermediate

    def draw_precursor(self, fold=10, size=14, as_string=True, out_file=None, amino_acid_sequence=None, leader=True,
                       x_translation=0, y_translation=0, min_padding_y=0):
        if amino_acid_sequence is None:
            amino_acid_sequence = self.full_peptide
        # set begin of chain so that last amino acid is forward
        circles = []
        texts = []
        text_colour = "#000000"
        if leader:

            current_y = min_padding_y + size
            current_x = min([fold, len(amino_acid_sequence)]) * size * 2 + size * 0.5

            forward = False
            # reverse sequence to go backwards
            amino_acid_sequence = amino_acid_sequence[::-1]
            for index, amino_acid in enumerate(amino_acid_sequence):
                # TODO: create dictionary for amino acid colors
                circles.append(make_circle(current_x, current_y, size, amino_acid))
                text = f"""<text x="{current_x}" y="{current_y}" fill="{text_colour}" text-anchor="middle" font-family="verdana" font-size = "{10}">\
                    <tspan y="{current_y}" dy="0.35em">{amino_acid}</tspan></text>"""
                texts.append(text)
                step = size * 2
                if index % fold == 0 and index != 0:
                    forward = not forward
                if (index % fold > fold - 2 or index % fold < 1) and index > 2:
                    step = 2 ** (1 / 2) * step / 2
                    current_y += step
                if forward:
                    current_x += step
                else:
                    current_x -= step

        if not leader:
            current_y = size + y_translation
            current_x = size + x_translation
            forward = True
            for index, amino_acid in enumerate(amino_acid_sequence):
                circles.append(make_circle(current_x, current_y, size, amino_acid))
                text = f"""<text x="{current_x}" y="{current_y}" fill="{text_colour}" text-anchor="middle" font-family="verdana" font-size = "{10}">\
                    <tspan y="{current_y}" dy="0.35em">{amino_acid}</tspan></text>"""
                texts.append(text)
                step = size * 2
                if index % fold == 0 and index != 0:
                    forward = not forward
                if (index % fold > fold - 2 or index % fold < 1) and index > 2:
                    step = 2 ** (1 / 2) * step / 2
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
            svg_string = f"""<svg width="{fold * size * 2 + size * 2}" height="{math.ceil(len(amino_acid_sequence) / fold) * size * (2 + 2 ** 0.5) + size}" viewBox="{0} {0} {fold * size * 2 + size} {math.ceil(len(amino_acid_sequence) / fold) * size * 3 + size}" xmlns="http://www.w3.org/2000/svg">\n"""
            svg_string += svg
            svg_string += "</svg>"
            if out_file is None:
                raise ValueError(
                    "Must provide output svg directory if 'as_string' is set to False.")
            else:
                with open(out_file, 'w') as svg_out:
                    svg_out.write(svg_string)

    def draw_cluster(self, fold=10, size=14, as_string=True, out_file=None, add_url=True,
                     draw_cs_in_pink=False):
        leader_pos = Vector(0, 0)
        follower_pos = Vector(0, 0)

        amino_acid_sequence_without_core = self.full_peptide.split(
            self.core_peptide)
        if len(amino_acid_sequence_without_core) != 2:
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

        if not self.cyclic_product:

            drawing = RaichuDrawer(structure, dont_show=True, add_url=add_url, ripp=True, horizontal=True,
                                   draw_Cs_in_pink=draw_cs_in_pink)
        else:
            drawing = RaichuDrawer(structure, dont_show=True, add_url=add_url, ripp=True,
                                   make_linear=False, horizontal=True, draw_Cs_in_pink=draw_cs_in_pink)

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

        svg_bubbles_leader = self.draw_precursor(fold=fold, size=size, as_string=True,
                                                 amino_acid_sequence=amino_acid_sequence_leader, leader=True,
                                                 min_padding_y=max(leader_pos.y, follower_pos.y))
        position_x_first_bubble_leader = min([fold, len(amino_acid_sequence_leader)]) * size * 2 + size * 0.5
        position_y_first_bubble_leader = max(leader_pos.y, follower_pos.y) + size

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

            y_translation = follower_pos.y - size

            x_translation = follower_pos.x
            svg_bubbles_follower = self.draw_precursor(
                fold=fold, size=size, as_string=True, amino_acid_sequence=amino_acid_sequence_follower, leader=False,
                x_translation=x_translation, y_translation=y_translation)
        else:
            svg_bubbles_follower = ""

        x1 = 0
        x2 = x_translation + 2 * size * min(len(amino_acid_sequence_follower), fold) + 2 * size
        y1 = 0

        y2 = max(max_y + y_translation_leader + 10,
                 (math.ceil(len(amino_acid_sequence_leader) / fold) * size * (1 + 2 ** 0.5) + max_y + 10),
                 y_translation + 2 * size + math.ceil(len(amino_acid_sequence_follower) / fold) * size * (1 + 2 ** 0.5))
        # print(x1, x2, y1, y2)

        svg_string = f"""<svg width="{x2}" height="{y2}" viewBox="{x1} {y1} {x2} {y2}" xmlns="http://www.w3.org/2000/svg">\n"""
        if svg_style:
            svg_string += f"{svg_style}\n"
        svg_string += svg_bubbles_leader
        if follower_pos.x != 0:
            svg_string += svg_bubbles_follower
        svg_string += svg
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