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

def _helper_build_substrate_with_arginin(tailoring_representations=None) -> RiPPCluster:
    cluster = RiPPCluster(
        "truE",
        "",
        "AKRSTRCT",
        tailoring_representations=tailoring_representations,
    )
    cluster.make_peptide()
    cluster.do_tailoring()
    return cluster

def _helper_build_substrate_with_aspartate(tailoring_representations=None) -> RiPPCluster:
    cluster = RiPPCluster(
        "truE",
        "",
        "AKDSTRCT",
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
        == "N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(NC(C(N[C@H](C(N[C@H](C(N[C@H](C(O)=O)[C@@H](C)O)=O)CS)=O)CC(N)=O)=O)C=C)=O)CO)=O)CC(O)=O)=O)Cc1cnc[nH]1)=O)CCCCN)=O)C"
    )


def test_threonin_serin_dehydratase():
    cluster = _helper_build_substrate(
        [
            TailoringRepresentation("truF", "THREONINE_SERINE_DEHYDRATASE", [["O_48"]]),
        ]
    )
    assert (
        structure_to_smiles(cluster.chain_intermediate)
        == "N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(NC(C(N[C@H](C(N[C@H](C(N[C@H](C(O)=O)[C@@H](C)O)=O)CS)=O)CC(N)=O)=O)=CC)=O)CO)=O)CC(O)=O)=O)Cc1cnc[nH]1)=O)CCCCN)=O)C"
    )

def test_prenyltransferase():
    cluster = _helper_build_substrate(
        [
            TailoringRepresentation("truF", "PRENYLTRANSFERASE", [["C_10"]], "SQUALENE" ),
        ]
    )
    assert (
        structure_to_smiles(cluster.chain_intermediate)
        == "N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(O)=O)[C@@H](C)O)=O)CS)=O)CC(N)=O)=O)[C@@H](C)O)=O)CO)=O)CC(O)=O)=O)Cc1cnc[nH]1)=O)CC(CC(C)=CCC/C(/C)=C/CC/C(/C)=C/CC/C=C(\C)/CC/C=C(\C)/CCC=C(C)C)CCN)=O)C"
    )


def test_aminotransferase():
    cluster = _helper_build_substrate(
        [
            TailoringRepresentation("truF", "AMINOTRANSFERASE", [["C_34"]]),
        ]
    )
    assert (
        structure_to_smiles(cluster.chain_intermediate)
        == "N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N)N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(O)=O)[C@@H](C)O)=O)CS)=O)CC(N)=O)=O)[C@@H](C)O)=O)CO)CC(O)=O)=O)Cc1cnc[nH]1)=O)CCCCN)=O)C"
    )

def test_hydroxylase():
    cluster = _helper_build_substrate(
        [
            TailoringRepresentation("truF", "HYDROXYLASE", [["C_12"]]),
        ]
    )
    assert (
        structure_to_smiles(cluster.chain_intermediate)
        == "N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(O)=O)[C@@H](C)O)=O)CS)=O)CC(N)=O)=O)[C@@H](C)O)=O)CO)=O)CC(O)=O)=O)Cc1cnc[nH]1)=O)CCCC(O)N)=O)C"
    )

def test_halogenase():
    cluster = _helper_build_substrate(
        [
            TailoringRepresentation("truF", "HALOGENASE", [["C_12"]], "F" ),
        ]
        
    )
    assert (
        structure_to_smiles(cluster.chain_intermediate)
        == "N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(O)=O)[C@@H](C)O)=O)CS)=O)CC(N)=O)=O)[C@@H](C)O)=O)CO)=O)CC(O)=O)=O)Cc1cnc[nH]1)=O)CCCC(F)N)=O)C"
    )

def test_acetyltransferase():
    cluster = _helper_build_substrate(
        [
            TailoringRepresentation("truF", "ACETYLTRANSFERASE", [["C_12"]] ),
        ]
        
    )
    assert (
        structure_to_smiles(cluster.chain_intermediate)
        == "N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(O)=O)[C@@H](C)O)=O)CS)=O)CC(N)=O)=O)[C@@H](C)O)=O)CO)=O)CC(O)=O)=O)Cc1cnc[nH]1)=O)CCCC(C(=O)C)N)=O)C"
    )

def test_keto_reduction():
    cluster = _helper_build_substrate(
        [
            TailoringRepresentation("truF", "KETO_REDUCTION", [["C_4"]] ),
        ]
        
    )
    assert (
        structure_to_smiles(cluster.chain_intermediate)
        == "N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(O)=O)[C@@H](C)O)=O)CS)=O)CC(N)=O)=O)[C@@H](C)O)=O)CO)=O)CC(O)=O)=O)Cc1cnc[nH]1)=O)CCCCN)O)C"
    )
def test_acyltransferase():
    cluster = _helper_build_substrate(
        [
            TailoringRepresentation("truF", "ACYLTRANSFERASE", [["C_10"]], "PALAMITIC_ACID" ),
        ]
        
    )
    assert (
        structure_to_smiles(cluster.chain_intermediate)
        == "N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(O)=O)[C@@H](C)O)=O)CS)=O)CC(N)=O)=O)[C@@H](C)O)=O)CO)=O)CC(O)=O)=O)Cc1cnc[nH]1)=O)CC(CCCCCCCCCCCCCCCC(=O)[O-])CCN)=O)C"
    )

def test_double_bond_reductase():
    cluster = _helper_build_substrate(
        [
            TailoringRepresentation("truF", "THREONINE_SERINE_DEHYDRATASE", [["O_48"]]),
            TailoringRepresentation("truF", "DOUBLE_BOND_REDUCTASE", [[ "C_44", "C_46"]] ),
        ]

    )
    assert (
        structure_to_smiles(cluster.chain_intermediate)
        == "N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(NC(C(N[C@H](C(N[C@H](C(N[C@H](C(O)=O)[C@@H](C)O)=O)CS)=O)CC(N)=O)=O)CC)=O)CO)=O)CC(O)=O)=O)Cc1cnc[nH]1)=O)CCCCN)=O)C"
    )
def test_dehydrogenase():
    cluster = _helper_build_substrate(
        [
            TailoringRepresentation("truF", "DEHYDROGENASE", [[ "C_11", "C_12"]] ),
        ]
        
    )
    assert (
        structure_to_smiles(cluster.chain_intermediate)
        == "N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(O)=O)[C@@H](C)O)=O)CS)=O)CC(N)=O)=O)[C@@H](C)O)=O)CO)=O)CC(O)=O)=O)Cc1cnc[nH]1)=O)CCC=CN)=O)C"
    )
def test_dehydratase():
    cluster = _helper_build_substrate(
        [
            TailoringRepresentation("truF", "DEHYDRATASE", [[ "C_37", "C_39"]] ),
        ]   
    )
    assert (
        structure_to_smiles(cluster.chain_intermediate)
        == "N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(NC(C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(O)=O)[C@@H](C)O)=O)CS)=O)CC(N)=O)=O)[C@@H](C)O)=O)=C)=O)CC(O)=O)=O)Cc1cnc[nH]1)=O)CCCCN)=O)C"
    )
def test_peptidase():
    cluster = _helper_build_substrate(
        [
            TailoringRepresentation("truF", "PEPTIDASE", [[ "N_43", "C_41"]]),
        ]
    )
    assert (
        structure_to_smiles(cluster.chain_intermediate)
        == "SC[C@@H](C(N[C@H](C(O)=O)[C@@H](C)O)=O)NC(=O)[C@H](CC(N)=O)NC(=O)[C@H]([C@@H](C)O)N"
    )
def test_protease():
    cluster = _helper_build_substrate(
        [
            TailoringRepresentation("truF", "PROTEASE", [[ "N_52", "C_50"]]),
        ]
        
    )
    assert (
        structure_to_smiles(cluster.chain_intermediate)
        == "SC[C@@H](C(N[C@H](C(O)=O)[C@@H](C)O)=O)NC(=O)[C@H](CC(N)=O)N"
    )

def test_alcohol_dehydrogenase():
    cluster = _helper_build_substrate(
        [
            TailoringRepresentation("truF", "ALCOHOL_DEHYDROGENASE", [[ "O_40"]]),
        ]
    )
    assert (
        structure_to_smiles(cluster.chain_intermediate)
        == "N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(O)=O)[C@@H](C)O)=O)CS)=O)CC(N)=O)=O)[C@@H](C)O)=O)C=O)=O)CC(O)=O)=O)Cc1cnc[nH]1)=O)CCCCN)=O)C"
    )
def test_epoxidase():
    cluster = _helper_build_substrate(
        [
            TailoringRepresentation("truF", "THREONINE_SERINE_DEHYDRATASE", [["O_48"]]),
            TailoringRepresentation("truF", "EPOXIDASE", [[ "C_46" , "C_44"]]),
        ]    
    )
    assert (
        structure_to_smiles(cluster.chain_intermediate)
        == "N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(NC(O2)(C(N[C@H](C(N[C@H](C(N[C@H](C(O)=O)[C@@H](C)O)=O)CS)=O)CC(N)=O)=O)C2C)=O)CO)=O)CC(O)=O)=O)Cc1cnc[nH]1)=O)CCCCN)=O)C"
    )

def test_methyl_mutase():
    cluster = _helper_build_substrate(
        [
            TailoringRepresentation("truF", "METHYL_MUTASE", [[ "C_3", "C_64"]]),
        ]  
    )
    assert (
        structure_to_smiles(cluster.chain_intermediate)
        == "NCC(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(O)=O)[C@@H](C)O)=O)C(C)S)=O)CC(N)=O)=O)[C@@H](C)O)=O)CO)=O)CC(O)=O)=O)Cc1cnc[nH]1)=O)CCCCN)=O"
    )

def test_monoamine_oxidase():
    cluster = _helper_build_substrate(
        [
            TailoringRepresentation("truF", "MONOAMINE_OXIDASE", [[ "N_13"]]),
        ]
        
    )
    assert (
        structure_to_smiles(cluster.chain_intermediate)
        == "N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(O)=O)[C@@H](C)O)=O)CS)=O)CC(N)=O)=O)[C@@H](C)O)=O)CO)=O)CC(O)=O)=O)Cc1cnc[nH]1)=O)CCCC=O)=O)C"
    )

def test_decarboxylase():
    cluster = _helper_build_substrate(
        [
            TailoringRepresentation("truF", "DECARBOXYLASE", [[ "C_75"]]),
        ]  
    )
    assert (
        structure_to_smiles(cluster.chain_intermediate)
        == "N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(NC[C@@H](C)O)=O)CS)=O)CC(N)=O)=O)[C@@H](C)O)=O)CO)=O)CC(O)=O)=O)Cc1cnc[nH]1)=O)CCCCN)=O)C"
    )


def test_splicease():
    cluster = _helper_build_substrate(
        [
            TailoringRepresentation("truF", "SPLICEASE", [[ "C_19", "N_43"]]),
        ]   
    )
    assert (
        structure_to_smiles(cluster.chain_intermediate)
        == "SC[C@@H](C(N[C@H](C(O)=O)[C@@H](C)O)=O)NC(=O)[C@H](CC(N)=O)NC(=O)[C@H]([C@@H](C)O)NCc1cnc[nH]1"

    )
def test_arginase():
    cluster = _helper_build_substrate_with_arginin(
        [
           TailoringRepresentation("truF", "ARGINASE", [[ "N_50"]]),
        ]
        
    )
    assert (
        structure_to_smiles(cluster.chain_intermediate)
        == "N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(O)=O)[C@@H](C)O)=O)CS)=O)CCCN)=O)[C@@H](C)O)=O)CO)=O)CCCNC(N)=N)=O)CCCCN)=O)C"

    )

def test_oxidative_bond_synthase():
    cluster = _helper_build_substrate_with_arginin(
        [
           TailoringRepresentation("truF", "OXIDATIVE_BOND_SYNTHASE", [[ "O_72", "C_41"]]),
        ]
    )
    assert (
        structure_to_smiles(cluster.chain_intermediate)
        == "N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C1=O)[C@@H](COC(=O)[C@H]([C@@H](C)O)NC(=O)[C@H](CS)NC(=O)[C@H](CCCNC(N)=N)N1)O)=O)CO)=O)CCCNC(N)=N)=O)CCCCN)=O)C"

    )
def test_methyltransferase():
    cluster = _helper_build_substrate_with_arginin(
        [
           TailoringRepresentation("truF", "METHYLTRANSFERASE", [["C_9"]]),
        ]
    )
    assert (
        structure_to_smiles(cluster.chain_intermediate)
        == "N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(O)=O)[C@@H](C)O)=O)CS)=O)CCCNC(N)=N)=O)[C@@H](C)O)=O)CO)=O)CCCNC(N)=N)=O)C(C)CCCN)=O)C"

    )
def test_c_methyltransferase():
    cluster = _helper_build_substrate_with_arginin(
        [
           TailoringRepresentation("truF", "C_METHYLTRANSFERASE", [["C_49"]]),
        ]   
    )
    assert (
        structure_to_smiles(cluster.chain_intermediate)
        == "N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(O)=O)[C@@H](C)O)=O)CS)=O)CCC(C)NC(N)=N)=O)[C@@H](C)O)=O)CO)=O)CCCNC(N)=N)=O)CCCCN)=O)C"

    )
def test_n_methyltransferase():
    cluster = _helper_build_substrate_with_arginin(
        [
           TailoringRepresentation("truF", "N_METHYLTRANSFERASE", [["N_50"]]),
        ]   
    )
    assert (
        structure_to_smiles(cluster.chain_intermediate)
        == "N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(O)=O)[C@@H](C)O)=O)CS)=O)CCCN(C)C(N)=N)=O)[C@@H](C)O)=O)CO)=O)CCCNC(N)=N)=O)CCCCN)=O)C"

    )
def test_o_methyltransferase():
    cluster = _helper_build_substrate_with_arginin(
        [
           TailoringRepresentation("truF", "O_METHYLTRANSFERASE", [["O_32"]]),
        ]  
    )
    assert (
        structure_to_smiles(cluster.chain_intermediate)
        == "N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(O)=O)[C@@H](C)O)=O)CS)=O)CCCNC(N)=N)=O)[C@@H](C)O)=O)COC)=O)CCCNC(N)=N)=O)CCCCN)=O)C"

    )

def test_amino_acid_epimerase():
    cluster = _helper_build_substrate_with_arginin(
        [
           TailoringRepresentation("truF", "AMINO_ACID_EPIMERASE", [["C_64"]]),
        ]   
    )
    assert (
        structure_to_smiles(cluster.chain_intermediate)
        == "N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@@H](C(O)=O)[C@@H](C)O)=O)CS)=O)CCCNC(N)=N)=O)[C@@H](C)O)=O)CO)=O)CCCNC(N)=N)=O)CCCCN)=O)C"

    )
def test_omega_ester():
    cluster = _helper_build_substrate_with_aspartate(
        [
           TailoringRepresentation("truF", "OMEGA_ESTER", [["O_22","O_29"]]),
        ]   
    )
    assert (
        structure_to_smiles(cluster.chain_intermediate)
        == "N[C@H](C(N[C@H](C(N[C@H](C1=O)CC(OC[C@@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(O)=O)[C@@H](C)O)=O)CS)=O)CCCNC(N)=N)=O)[C@@H](C)O)=O)N1)=O)=O)CCCCN)=O)C"

    )
def test_omega_thioester():
    cluster = _helper_build_substrate_with_aspartate(
        [
           TailoringRepresentation("truF", "OMEGA_THIOESTER", [["O_22","S_57"]]),
        ]  
    )
    assert (
        structure_to_smiles(cluster.chain_intermediate)
        == "N[C@H](C(N[C@H](C(N[C@H](C1=O)CC(SC[C@@H](C(N[C@H](C(O)=O)[C@@H](C)O)=O)NC(=O)[C@H](CCCNC(N)=N)NC(=O)[C@H]([C@@H](C)O)NC(=O)[C@H](CO)N1)=O)=O)CCCCN)=O)C"

    )

def test_omega_amide():
    cluster = _helper_build_substrate_with_aspartate(
        [
           TailoringRepresentation("truF", "OMEGA_AMIDE", [["O_22","N_13"]]),
        ]   
    )
    assert (
        structure_to_smiles(cluster.chain_intermediate)
        == "N[C@H](C(N[C@H]1CCCCNC(=O)C[C@@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(O)=O)[C@@H](C)O)=O)CS)=O)CCCNC(N)=N)=O)[C@@H](C)O)=O)CO)=O)NC1=O)=O)C"

    )
def test_hydrolase():
    cluster = _helper_build_substrate_with_aspartate(
        [
           TailoringRepresentation("truF", "OMEGA_ESTER", [["O_22","O_29"]]),
           TailoringRepresentation("truF", "HYDROLASE", [["C_28","O_29"]]),
        ]  
    )
    assert (
        structure_to_smiles(cluster.chain_intermediate)
        == "OC[C@@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(O)=O)[C@@H](C)O)=O)CS)=O)CCCNC(N)=N)=O)[C@@H](C)O)=O)NC(=O)[C@H](CC(O)=O)NC(=O)[C@H](CCCCN)NC(=O)[C@H](C)N"

    )

def test_macrolactamidination():
    cluster = _helper_build_substrate_with_aspartate(
        [
           TailoringRepresentation("truF", "MACROLACTAMIDINATION", [["C_58", "N_13"]]),
        ]  
    )
    assert (
        structure_to_smiles(cluster.chain_intermediate)
        == "N[C@H](C(N[C@H]1CCCCNC(=N[C@H](C(O)=O)[C@@H](C)O)[C@H](CS)NC(=O)[C@H](CCCNC(N)=N)NC(=O)[C@H]([C@@H](C)O)NC(=O)[C@H](CO)NC(=O)[C@H](CC(O)=O)NC1=O)=O)C"

    )


def test_lanthipeptide_cyclase():
    cluster = _helper_build_substrate_with_aspartate(
        [
            TailoringRepresentation("truF", "THREONINE_SERINE_DEHYDRATASE", [["O_37"]]),
            TailoringRepresentation("truF", "LANTHIPEPTIDE_CYCLASE", [["S_57", "C_35"]]),
        ]   
    )
    assert (
        structure_to_smiles(cluster.chain_intermediate)
        == "N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(NC(C1=O)C(SC[C@@H](C(N[C@H](C(O)=O)[C@@H](C)O)=O)NC(=O)[C@H](CCCNC(N)=N)N1)C)=O)CO)=O)CC(O)=O)=O)CCCCN)=O)C"

    )
def test_lanthionine_synthetase():
    cluster = _helper_build_substrate_with_aspartate(
        [
            TailoringRepresentation("truF", "LANTHIONINE_SYNTHETASE", [["S_57", "C_28"]]),
        ]   
    )
    assert (
        structure_to_smiles(cluster.chain_intermediate)
        == "N[C@H](C(N[C@H](C(N[C@H](C(N[C@H]1CSC[C@@H](C(N[C@H](C(O)=O)[C@@H](C)O)=O)NC(=O)[C@H](CCCNC(N)=N)NC(=O)[C@H]([C@@H](C)O)NC1=O)=O)CC(O)=O)=O)CCCCN)=O)C"

    )
def test_thioamidation():
    cluster = _helper_build_substrate_with_aspartate(
        [
            TailoringRepresentation("truF", "THIOAMIDATION", [["C_23"]]),
        ] 
    )
    assert (
        structure_to_smiles(cluster.chain_intermediate)
        == "N[C@H](C(N[C@H](C(N[C@H](C(=S)N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(O)=O)[C@@H](C)O)=O)CS)=O)CCCNC(N)=N)=O)[C@@H](C)O)=O)CO)CC(O)=O)=O)CCCCN)=O)C"

    )
def test_reductive_lyase():
    cluster = _helper_build_substrate_with_aspartate(
        [
            TailoringRepresentation("truF", "REDUCTIVE_LYASE", [["C_26", "C_30"]]),
        ] 
    )
    assert (
        structure_to_smiles(cluster.chain_intermediate)
        == "N[C@H](C(N[C@H](C(N[C@H](C(NCCO)=O)CC(O)=O)=O)CCCCN)=O)C"
    )






def test_macrolactam_synthetase():
    cluster = _helper_build_substrate_with_aspartate(
        [
           TailoringRepresentation("truF", "MACROLACTAM_SYNTHETASE", [["O_22"]]),
        ]   
    )
    get_tailoring_sites(
    cluster.chain_intermediate,
    enzyme_name="MACROLACTAM_SYNTHETASE",
    out_file="test.svg",
    )
    #print(structure_to_smiles(cluster.chain_intermediate))
    assert (
        structure_to_smiles(cluster.chain_intermediate)
        == "N[C@H](C(N[C@H](C(N[C@H](C1=O)CC(NC(=N)NCCC[C@@H](C(N[C@H](C(N[C@H](C(O)=O)[C@@H](C)O)=O)CS)=O)NC(=O)[C@H]([C@@H](C)O)NC(=O)[C@H](CO)N1)=O)=O)CCCCN)=O)C"

    )

def test_thiopeptide_cyclase():
    cluster = _helper_build_substrate_with_aspartate(
        [
            TailoringRepresentation("truF", "THREONINE_SERINE_DEHYDRATASE", [["O_37"]]),
            TailoringRepresentation("truF", "THIOPEPTIDE_CYCLASE", [["C_35"]]),
        ]
        
    )
    get_tailoring_sites(
        cluster.chain_intermediate,
        enzyme_name="THIOPEPTIDE_CYCLASE",
        out_file="test.svg",
    )
    #print(structure_to_smiles(cluster.chain_intermediate))
    assert (
        structure_to_smiles(cluster.chain_intermediate)
        == "N[C@H](C(N[C@H](C(N[C@H](C(N[C@H]1CSC[C@@H](C(N[C@H](C(O)=O)[C@@H](C)O)=O)NC(=O)[C@H](CCCNC(N)=N)NC(=O)[C@H]([C@@H](C)O)NC1=O)=O)CC(O)=O)=O)CCCCN)=O)C"

    )
def test_cyclodehydrase():
    cluster = _helper_build_substrate_with_aspartate(
        #[
            #TailoringRepresentation("truF", "THREONINE_SERINE_DEHYDRATASE", [["O_37"]]),
            #TailoringRepresentation("truF", "CYCLODEHYDRASE", [["C_35"]]),
        #]
        
    )
    get_tailoring_sites(
        cluster.chain_intermediate,
        enzyme_name="CYCLODEHYDRASE",
        out_file="test.svg",
    )
    #print(structure_to_smiles(cluster.chain_intermediate))
    assert (
        structure_to_smiles(cluster.chain_intermediate)
        == "N[C@H](C(N[C@H](C(N[C@H](C(N[C@H]1CSC[C@@H](C(N[C@H](C(O)=O)[C@@H](C)O)=O)NC(=O)[C@H](CCCNC(N)=N)NC(=O)[C@H]([C@@H](C)O)NC1=O)=O)CC(O)=O)=O)CCCCN)=O)C"

    )
if __name__ == "__main__":
    test_thiopeptide_cyclase()


