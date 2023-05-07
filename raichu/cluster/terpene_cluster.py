from typing import List

from pikachu.general import read_smiles

from raichu.substrate import TerpeneCyclaseSubstrate
from raichu.representations import MacrocyclizationRepresentation, TailoringRepresentation
from raichu.reactions.general_tailoring_reactions import dephosphorylation, oxidative_bond_formation
from raichu.drawing.drawer import RaichuDrawer
from raichu.cluster.base_cluster import Cluster


class TerpeneCluster(Cluster):
    def __init__(self, cyclase_name: str, precursor: str, cyclase_type: str = None,
                 macrocyclisations: List[MacrocyclizationRepresentation] = None,
                 tailoring_representations: List[TailoringRepresentation] = None) -> None:
        super().__init__(tailoring_representations)
        self.cyclase_name = cyclase_name
        self.cyclase_type = cyclase_type
        self.precursor = precursor
        self.macrocyclisations = macrocyclisations
        self.tailoring_enzymes_representation = tailoring_representations
        self.chain_intermediate = None
        self.tailored_product = None
        self.cyclised_product = None
        self.final_product = None

    def create_precursor(self) -> None:
        substrate = TerpeneCyclaseSubstrate(self.precursor)
        self.chain_intermediate = read_smiles(substrate.smiles)

    def do_macrocyclization(self):
        initialized_macrocyclization_atoms = self.initialize_macrocyclization()

        if self.cyclase_type != "Class_2":
            self.chain_intermediate = dephosphorylation(self.chain_intermediate)

        for macrocyclization_atoms in initialized_macrocyclization_atoms:
            atom1 = self.chain_intermediate.get_atom(macrocyclization_atoms[0])
            atom2 = self.chain_intermediate.get_atom(macrocyclization_atoms[1])
            self.chain_intermediate = oxidative_bond_formation(atom1, atom2, self.chain_intermediate)

        self.cyclised_product = self.chain_intermediate

    def draw_product(self, as_string=True, out_file=None):
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
