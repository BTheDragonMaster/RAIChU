import pytest
import os
from pikachu.general import structure_to_smiles
from raichu.representations import (
    TailoringRepresentation,
)
from raichu.cluster.ripp_cluster import RiPPCluster
from raichu.run_raichu import get_tailoring_sites


def _helper_build_substrate(tailoring_representations=None) -> RiPPCluster:
    cluster = RiPPCluster(
        "truE",
        "",
        "AKHDSTNCT",
        tailoring_representations=tailoring_representations,
    )
    cluster.make_peptide()
    cluster.do_tailoring()
    return cluster


def test_double_bond_isomerase():
    cluster = _helper_build_substrate(
        [
            TailoringRepresentation("truF", "THREONINE_SERINE_DEHYDRATASE", [["O_48"]]),
            TailoringRepresentation(
                "truF", "DOUBLE_BOND_ISOMERASE", [["C_46", "C_44", "C_46", "C_49"]]
            ),
        ]
    )
    assert (
        structure_to_smiles(cluster.chain_intermediate)
        == r"N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(NC(C(N[C@H](C(N[C@H](C(N[C@H](C(O)=O)[C@@H](C)O)=O)CS)=O)CC(N)=O)=O)C=C)=O)CO)=O)CC(O)=O)=O)Cc1cnc[nH]1)=O)CCCCN)=O)C"
    )


def test_threonin_serin_dehydratase():
    cluster = _helper_build_substrate(
        [
            TailoringRepresentation("truF", "THREONINE_SERINE_DEHYDRATASE", [["O_48"]]),
        ]
    )
    assert (
        structure_to_smiles(cluster.chain_intermediate)
        == r"N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(NC(C(N[C@H](C(N[C@H](C(N[C@H](C(O)=O)[C@@H](C)O)=O)CS)=O)CC(N)=O)=O)=CC)=O)CO)=O)CC(O)=O)=O)Cc1cnc[nH]1)=O)CCCCN)=O)C"
    )


def test_prenyltransferase():
    cluster = _helper_build_substrate(
        [
            TailoringRepresentation(
                "truF", "PRENYLTRANSFERASE", [["C_10"]], "SQUALENE"
            ),
        ]
    )

    assert (
        structure_to_smiles(cluster.chain_intermediate)
        == r"N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(O)=O)[C@@H](C)O)=O)CS)=O)CC(N)=O)=O)[C@@H](C)O)=O)CO)=O)CC(O)=O)=O)Cc1cnc[nH]1)=O)CC(CC(C)=CCC/C(/C)=C/CC/C(/C)=C/CC/C=C(\C)/CC/C=C(\C)/CCC=C(C)C)CCN)=O)C"
    )


def test_aminotransferase():
    cluster = _helper_build_substrate(
        [
            TailoringRepresentation("truF", "AMINOTRANSFERASE", [["C_34"]]),
        ]
    )

    assert (
        structure_to_smiles(cluster.chain_intermediate)
        == r"N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N)N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(O)=O)[C@@H](C)O)=O)CS)=O)CC(N)=O)=O)[C@@H](C)O)=O)CO)CC(O)=O)=O)Cc1cnc[nH]1)=O)CCCCN)=O)C"
    )


def test_hydroxylase():
    cluster = _helper_build_substrate(
        [
            TailoringRepresentation("truF", "HYDROXYLASE", [["C_12"]]),
        ]
    )

    assert (
        structure_to_smiles(cluster.chain_intermediate)
        == r"N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(O)=O)[C@@H](C)O)=O)CS)=O)CC(N)=O)=O)[C@@H](C)O)=O)CO)=O)CC(O)=O)=O)Cc1cnc[nH]1)=O)CCCC(O)N)=O)C"
    )


def test_halogenase():
    cluster = _helper_build_substrate(
        [
            TailoringRepresentation("truF", "HALOGENASE", [["C_12"]], "F"),
        ]
    )

    assert (
        structure_to_smiles(cluster.chain_intermediate)
        == r"N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(O)=O)[C@@H](C)O)=O)CS)=O)CC(N)=O)=O)[C@@H](C)O)=O)CO)=O)CC(O)=O)=O)Cc1cnc[nH]1)=O)CCCC(F)N)=O)C"
    )


def test_acetyltransferase():
    cluster = _helper_build_substrate(
        [
            TailoringRepresentation("truF", "ACETYLTRANSFERASE", [["C_12"]]),
        ]
    )
    assert (
        structure_to_smiles(cluster.chain_intermediate)
        == r"N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(O)=O)[C@@H](C)O)=O)CS)=O)CC(N)=O)=O)[C@@H](C)O)=O)CO)=O)CC(O)=O)=O)Cc1cnc[nH]1)=O)CCCC(C(=O)C)N)=O)C"
    )


def test_keto_reduction():
    cluster = _helper_build_substrate(
        [
            TailoringRepresentation("truF", "KETO_REDUCTION", [["C_4"]]),
        ]
    )
    assert (
        structure_to_smiles(cluster.chain_intermediate)
        == r"N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(O)=O)[C@@H](C)O)=O)CS)=O)CC(N)=O)=O)[C@@H](C)O)=O)CO)=O)CC(O)=O)=O)Cc1cnc[nH]1)=O)CCCCN)O)C"
    )


def test_acyltransferase():
    cluster = _helper_build_substrate(
        [
            TailoringRepresentation(
                "truF", "ACYLTRANSFERASE", [["C_10"]], "OLEIC_ACID"
            ),
        ]
    )
    assert (
        structure_to_smiles(cluster.chain_intermediate)
        == r"N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(O)=O)[C@@H](C)O)=O)CS)=O)CC(N)=O)=O)[C@@H](C)O)=O)CO)=O)CC(O)=O)=O)Cc1cnc[nH]1)=O)CC(CCCCCCCCC=CCCCCCCCC(O)=O)CCN)=O)C"
    )


def test_double_bond_reductase():
    cluster = _helper_build_substrate(
        [
            TailoringRepresentation("truF", "THREONINE_SERINE_DEHYDRATASE", [["O_48"]]),
            TailoringRepresentation(
                "truF", "DOUBLE_BOND_REDUCTASE", [["C_44", "C_46"]]
            ),
        ]
    )
    assert (
        structure_to_smiles(cluster.chain_intermediate)
        == r"N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(NC(C(N[C@H](C(N[C@H](C(N[C@H](C(O)=O)[C@@H](C)O)=O)CS)=O)CC(N)=O)=O)CC)=O)CO)=O)CC(O)=O)=O)Cc1cnc[nH]1)=O)CCCCN)=O)C"
    )


def test_dehydrogenase():
    cluster = _helper_build_substrate(
        [
            TailoringRepresentation("truF", "DEHYDROGENASE", [["C_11", "C_12"]]),
        ]
    )
    assert (
        structure_to_smiles(cluster.chain_intermediate)
        == r"N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(O)=O)[C@@H](C)O)=O)CS)=O)CC(N)=O)=O)[C@@H](C)O)=O)CO)=O)CC(O)=O)=O)Cc1cnc[nH]1)=O)CCC=CN)=O)C"
    )


def test_dehydratase():
    cluster = _helper_build_substrate(
        [
            TailoringRepresentation("truF", "DEHYDRATASE", [["C_37", "C_39"]]),
        ]
    )

    assert (
        structure_to_smiles(cluster.chain_intermediate)
        == r"N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(NC(C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(O)=O)[C@@H](C)O)=O)CS)=O)CC(N)=O)=O)[C@@H](C)O)=O)=C)=O)CC(O)=O)=O)Cc1cnc[nH]1)=O)CCCCN)=O)C"
    )


def test_peptidase():
    cluster = _helper_build_substrate(
        [
            TailoringRepresentation("truF", "PEPTIDASE", [["N_43", "C_41"]]),
        ]
    )

    assert (
        structure_to_smiles(cluster.chain_intermediate)
        == r"SC[C@@H](C(N[C@H](C(O)=O)[C@@H](C)O)=O)NC(=O)[C@H](CC(N)=O)NC(=O)[C@H]([C@@H](C)O)N"
    )


def test_protease():
    cluster = _helper_build_substrate(
        [
            TailoringRepresentation("truF", "PROTEASE", [["N_52", "C_50"]]),
        ]
    )

    assert (
        structure_to_smiles(cluster.chain_intermediate)
        == r"SC[C@@H](C(N[C@H](C(O)=O)[C@@H](C)O)=O)NC(=O)[C@H](CC(N)=O)N"
    )


def test_alcohol_dehydrogenase():
    cluster = _helper_build_substrate(
        [
            TailoringRepresentation("truF", "ALCOHOL_DEHYDROGENASE", [["O_40"]]),
        ]
    )

    assert (
        structure_to_smiles(cluster.chain_intermediate)
        == r"N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(O)=O)[C@@H](C)O)=O)CS)=O)CC(N)=O)=O)[C@@H](C)O)=O)C=O)=O)CC(O)=O)=O)Cc1cnc[nH]1)=O)CCCCN)=O)C"
    )


def test_epoxidase():
    cluster = _helper_build_substrate(
        [
            TailoringRepresentation("truF", "THREONINE_SERINE_DEHYDRATASE", [["O_48"]]),
            TailoringRepresentation("truF", "EPOXIDASE", [["C_46", "C_44"]]),
        ]
    )

    assert (
        structure_to_smiles(cluster.chain_intermediate)
        == r"N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(NC(O2)(C(N[C@H](C(N[C@H](C(N[C@H](C(O)=O)[C@@H](C)O)=O)CS)=O)CC(N)=O)=O)C2C)=O)CO)=O)CC(O)=O)=O)Cc1cnc[nH]1)=O)CCCCN)=O)C"
    )


def test_methyl_mutase():
    cluster = _helper_build_substrate(
        [
            TailoringRepresentation("truF", "METHYL_MUTASE", [["C_3", "C_64"]]),
        ]
    )

    assert (
        structure_to_smiles(cluster.chain_intermediate)
        == r"NCC(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(O)=O)[C@@H](C)O)=O)C(C)S)=O)CC(N)=O)=O)[C@@H](C)O)=O)CO)=O)CC(O)=O)=O)Cc1cnc[nH]1)=O)CCCCN)=O"
    )


def test_monoamine_oxidase():
    cluster = _helper_build_substrate(
        [
            TailoringRepresentation("truF", "MONOAMINE_OXIDASE", [["N_13"]]),
        ]
    )
    assert (
        structure_to_smiles(cluster.chain_intermediate)
        == r"N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(O)=O)[C@@H](C)O)=O)CS)=O)CC(N)=O)=O)[C@@H](C)O)=O)CO)=O)CC(O)=O)=O)Cc1cnc[nH]1)=O)CCCC=O)=O)C"
    )


if __name__ == "__main__":
    test_acyltransferase()
