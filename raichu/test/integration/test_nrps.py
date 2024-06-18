import pytest
import os
from pikachu.general import structure_to_smiles
from raichu.run_raichu import build_cluster
from raichu.representations import (
    ClusterRepresentation,
    ModuleRepresentation,
    DomainRepresentation,
    TailoringRepresentation,
)


def test_NRPS():
    current_dir = os.path.dirname(os.path.abspath(__file__))
    # Construct the path to the test data
    cluster_dir = os.path.join(current_dir, "../test_data/nrps/hormaomycin/")
    cluster_representation = ClusterRepresentation.from_file(cluster_dir)
    cluster = build_cluster(cluster_representation, strict=False)
    cluster.compute_structures(compute_cyclic_products=False)
    cluster.do_tailoring()

    # Check if all tayloring reactions are performed
    for module in cluster.modules:
        for domain in module.tailoring_domains:
            if domain.supertype == "TAILORING":
                if domain.active == True:
                    assert domain.used == True
    # Check if final structure is correct
    assert (
        structure_to_smiles(cluster.chain_intermediate)
        == "OC(=O)[C@H](N(C(=O)[C@@H](NC(=O)[C@@H](NC(=O)[C@H](NC(=O)[C@@H](NC(=O)[C@H](NC(=O)[C@@H](NC(c6ccc(Cl)[nH]6)=O)C[C@@H]5C[C@@H]5[N+]([O-])=O)[C@H](O)C)[C@H](c4ccccc4)C)C[C@@H]2C[C@@H]2[N+]([O-])=O)[C@H](c3ccccc3)C)[C@H](C)CC)C1)C[C@H]1/C=C\C"
    )


if __name__ == "__main__":
    test_NRPS()
