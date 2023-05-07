from typing import List

from pikachu.general import read_smiles

from raichu.substrate import TerpeneCyclaseSubstrate
from raichu.representations import MacrocyclizationRepresentation, TailoringRepresentation
from raichu.reactions.general_tailoring_reactions import dephosphorylation, oxidative_bond_formation
from raichu.cluster.base_cluster import Cluster


class TerpeneCluster(Cluster):
    def __init__(self, cyclase_name: str, precursor: str, cyclase_type: str = None,
                 macrocyclisations: List[MacrocyclizationRepresentation] = None,
                 tailoring_representations: List[TailoringRepresentation] = None) -> None:
        super().__init__(tailoring_representations, macrocyclisations)
        self.cyclase_name = cyclase_name
        self.cyclase_type = cyclase_type
        self.precursor = precursor

        self.chain_intermediate = None
        self.tailored_product = None
        self.final_product = None

    def create_precursor(self) -> None:
        substrate = TerpeneCyclaseSubstrate(self.precursor)
        self.chain_intermediate = read_smiles(substrate.smiles)

    def do_macrocyclization(self, sequential=True):
        initialized_macrocyclization_atoms = self.initialize_macrocyclization()
        self.cyclic_intermediates.append(self.chain_intermediate.deepcopy())

        if self.cyclase_type != "Class_2":
            self.chain_intermediate = dephosphorylation(self.chain_intermediate)

        for macrocyclization_atoms, cyclisation_type in initialized_macrocyclization_atoms:
            atom1 = self.chain_intermediate.get_atom(macrocyclization_atoms[0])
            atom2 = self.chain_intermediate.get_atom(macrocyclization_atoms[1])
            self.chain_intermediate = oxidative_bond_formation(atom1, atom2, self.chain_intermediate)
            self.cyclic_intermediates.append(self.chain_intermediate.deepcopy())

        self.cyclic_product = self.chain_intermediate

