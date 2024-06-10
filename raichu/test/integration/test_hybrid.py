import pytest
import os
from pikachu.general import structure_to_smiles
from raichu.run_raichu import build_cluster
from raichu.representations import ClusterRepresentation


def test_hybrid_from_mibig():
    cluster_dir = os.path.realpath("../test_data/hybrid/microcystin/")
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
        == "OC(=O)[C@@H](NC(=O)C([*])NC(=O)C([*])NC(=O)[C@H](NC(=O)[C@@H](N(C)C(=O)C([*])NC(=O)CC(=O)C=CCC(=O)CC(=O)C([*])N)CO)C)CCCN=C(N)N"
    )


if __name__ == "__main__":
    test_hybrid_from_mibig()
