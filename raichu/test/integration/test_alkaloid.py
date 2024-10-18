import pytest
import os
from pikachu.general import structure_to_smiles
from raichu.representations import (
    MacrocyclizationRepresentation,
    TailoringRepresentation,
    IsomerizationRepresentation,
    MethylShiftRepresentation,
)
from raichu.cluster.alkaloid_cluster import AlkaloidCluster
from raichu.run_raichu import get_tailoring_sites, get_tailoring_sites_atom_names


def test_alkaloid():
    cluster = AlkaloidCluster(
        start_amino_acid="tryptophan",
        tailoring_representations=[
            TailoringRepresentation("bufA", "DECARBOXYLASE", [["C_12"]]),
            TailoringRepresentation("bufA", "METHYLTRANSFERASE", [["N_15"], ["N_15"]]),
            TailoringRepresentation("bufA", "HYDROXYLASE", [["C_0"]]),
        ],
    )
    cluster.make_scaffold()
    cluster.do_tailoring()

    # Check if final structure is correct
    assert (
        structure_to_smiles(cluster.chain_intermediate) == r"Oc1ccc([nH]cc2CCN(C)C)c2c1"
    )


if __name__ == "__main__":
    test_alkaloid()
