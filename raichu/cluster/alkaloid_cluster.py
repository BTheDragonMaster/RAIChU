from pikachu.general import read_smiles

from raichu.substrate import AminoAcidSubstrate
from raichu.tailoring_enzymes import TailoringEnzyme
from raichu.drawing.drawer import RaichuDrawer
from raichu.cluster.base_cluster import Cluster


class AlkaloidCluster(Cluster):
    def __init__(self, start_amino_acid: str, tailoring_representations=None) -> None:
        super().__init__(tailoring_representations)
        self.start_amino_acid = start_amino_acid
        self.linear_product = None
        self.final_product = None

    def make_scaffold(self):
        substrate = AminoAcidSubstrate(self.start_amino_acid)
        self.linear_product = read_smiles(substrate.smiles)
        self.chain_intermediate = self.linear_product


