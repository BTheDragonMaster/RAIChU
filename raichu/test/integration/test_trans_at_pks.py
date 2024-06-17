import pytest
import os
from pikachu.general import structure_to_smiles
from raichu.run_raichu import build_cluster
from raichu.representations import ClusterRepresentation


def test_extensive_KS_subtypes():
    cluster_dir = os.path.realpath("../test_data/trans_at_pks/")
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
        == "OC(=O)C[C@H](O)CC(=C)C[C@@H](O1)CC(C[*])CC1C(C)C=C[C@H](O2)CC(OC)C2C(=O)C=C[C@H](O)C"
    )


if __name__ == "__main__":
    test_extensive_KS_subtypes()
