import pytest
import os

from raichu.validation.drawing_validation.validate_central_chain import check_central_chains, \
    check_central_chains_from_file
from raichu.validation.create_random_clusters import generate_modular_cluster, AA_STARTER_CHOICES

# Remove substrate that causes known central chain error
AA_STARTER_CHOICES.remove("dehydrophenylalanine")


# NRPS


def test_central_chain_hormaomycin():
    current_dir = os.path.dirname(os.path.abspath(__file__))
    # Construct the path to the test data
    cluster_dir = os.path.join(current_dir, "../test_data/nrps/hormaomycin/")
    central_chain_error, faulty_spaghettis, total_backbone_errors = check_central_chains_from_file(cluster_dir)

    assert not central_chain_error
    assert faulty_spaghettis == 0
    assert total_backbone_errors == 0


# trans-AT PKS

def test_central_chain_trans_at():
    current_dir = os.path.dirname(os.path.abspath(__file__))
    # Construct the path to the test data
    cluster_dir = os.path.join(current_dir, "../test_data/trans_at_pks/")
    central_chain_error, faulty_spaghettis, total_backbone_errors = check_central_chains_from_file(cluster_dir)

    assert not central_chain_error
    assert faulty_spaghettis == 0
    assert total_backbone_errors == 0


# cis-AT PKS

def test_central_chain_coelimycin():
    current_dir = os.path.dirname(os.path.abspath(__file__))
    # Construct the path to the test data
    cluster_dir = os.path.join(current_dir, "../test_data/cis_at_pks/coelimycin/")
    central_chain_error, faulty_spaghettis, total_backbone_errors = check_central_chains_from_file(cluster_dir)

    assert not central_chain_error
    assert faulty_spaghettis == 0
    assert total_backbone_errors == 0


def test_central_chain_erythromycin():
    current_dir = os.path.dirname(os.path.abspath(__file__))
    # Construct the path to the test data
    cluster_dir = os.path.join(current_dir, "../test_data/cis_at_pks/erythromycin/")
    central_chain_error, faulty_spaghettis, total_backbone_errors = check_central_chains_from_file(cluster_dir)

    assert not central_chain_error
    assert faulty_spaghettis == 0
    assert total_backbone_errors == 0


# Hybrid

def test_central_chain_microcystin():
    current_dir = os.path.dirname(os.path.abspath(__file__))
    # Construct the path to the test data
    cluster_dir = os.path.join(current_dir, "../test_data/hybrid/microcystin/")
    central_chain_error, faulty_spaghettis, total_backbone_errors = check_central_chains_from_file(cluster_dir)

    assert not central_chain_error
    assert faulty_spaghettis == 0
    assert total_backbone_errors == 0


# Random

def test_random():
    for i in range(10):
        cluster = generate_modular_cluster(3)
        central_chain_error, faulty_spaghettis, total_backbone_errors = check_central_chains(cluster)

        assert not central_chain_error
        assert faulty_spaghettis == 0
        assert total_backbone_errors == 0


if __name__ == "__main__":
    test_central_chain_hormaomycin()
    test_central_chain_coelimycin()
    test_central_chain_erythromycin()
    test_central_chain_trans_at()
    test_central_chain_microcystin()
    test_random()
