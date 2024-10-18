import pytest
import os
from pikachu.general import structure_to_smiles
from raichu.representations import (
    MacrocyclizationRepresentation,
    TailoringRepresentation,
    IsomerizationRepresentation,
    MethylShiftRepresentation,
)
from raichu.cluster.terpene_cluster import TerpeneCluster


def test_terpene():
    cluster = TerpeneCluster(
        "TPS21",
        "GERANYLGERANYL_PYROPHOSPHATE",
        macrocyclisations=[
            MacrocyclizationRepresentation("C_15", "C_11"),
            MacrocyclizationRepresentation("C_10", "C_12"),
            MacrocyclizationRepresentation("C_9", "C_14"),
        ],
        cyclase_type="Class_1",
        methyl_shifts=[MethylShiftRepresentation(["C_16", "C_17"])],
        double_bond_isomerisations=[
            IsomerizationRepresentation(["C_24", "C_25", "C_25", "C_27"])
        ],
        tailoring_representations=[
            TailoringRepresentation("CYP71BN1", "EPOXIDASE", [["C_19", "C_20"]]),
            TailoringRepresentation("CYP71BN1", "HYDROXYLASE", [["C_19"]]),
            TailoringRepresentation("CYP71BN1", "REDUCTIVE_LYASE", [["C_19", "O_5"]]),
        ],
    )
    cluster.create_precursor()
    cluster.do_macrocyclization()
    cluster.do_tailoring()

    # Check if final structure is correct
    assert (
        structure_to_smiles(cluster.chain_intermediate)
        == r"CC(CC(O4)C4(CCCC(=C)C)C)C3C1CC2C3(C)C2C1"
    )


if __name__ == "__main__":
    test_terpene()
