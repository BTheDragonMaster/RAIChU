from pikachu.general import read_smiles

from raichu.substrate import AminoAcidSubstrate
from raichu.tailoring_enzymes import TailoringEnzyme
from raichu.drawing.drawer import RaichuDrawer

class Alkaloid_Cluster:
    def __init__(self, start_amino_acid: str, tailoring_enzymes_representation = None) -> None:
        self.start_amino_acid = start_amino_acid
        self.tailoring_enzymes_representation = tailoring_enzymes_representation
        self.chain_intermediate = None
        self.linear_product = None
        self.tailored_product = None
        self.final_product = None

    def make_scaffold(self):
        substrate = AminoAcidSubstrate(self.start_amino_acid)
        self.linear_product = read_smiles(substrate.smiles)
        self.chain_intermediate = self.linear_product

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

    def draw_product(self, as_string=True, out_file=None):
            assert self.chain_intermediate
            for atom in self.chain_intermediate.graph:
                if atom.chiral:
                    atom.chiral = False
            drawing = RaichuDrawer(self.chain_intermediate, dont_show=True, add_url=True, draw_Cs_in_pink=False, draw_straightened=False)
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

