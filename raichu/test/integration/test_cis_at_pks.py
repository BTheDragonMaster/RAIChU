import pytest
import os
from pikachu.general import structure_to_smiles
from raichu.run_raichu import build_cluster
from raichu.representations import ClusterRepresentation


def test_cis_at_PKS():
    current_dir = os.path.dirname(os.path.abspath(__file__))
    cluster_dir = os.path.join(current_dir, "../test_data/cis_at_pks/erythromycin/")
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
        == "OC([C@H]([C@H]([C@H]([C@@H](C(C[C@H](C([C@H]([C@@H]([C@@H]([C@@H](CC)O)C)O)C)=O)C)C)O)C)O)C)=O"
    )


if __name__ == "__main__":
    test_cis_at_PKS()
