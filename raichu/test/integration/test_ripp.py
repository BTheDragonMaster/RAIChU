from pikachu.general import structure_to_smiles
from raichu.representations import (
    MacrocyclizationRepresentation,
    TailoringRepresentation,
)
from raichu.run_raichu import get_tailoring_sites
from raichu.cluster.ripp_cluster import RiPPCluster


def test_ripp():
    cyanobactin_cluster_trunkamide = RiPPCluster(
        "truE",
        "MNKKNILPQLGQPVIRLTAGQLSSQLAELSEEALGGVDASTSIAPFCSYDGVDASSYDGVDASSYDD",
        "TSIAPFC",
        macrocyclisations=[
            MacrocyclizationRepresentation("N_0", "O_59", "condensative")
        ],
        tailoring_representations=[
            TailoringRepresentation("truD", "CYCLODEHYDRASE", [["S_56"]]),
            TailoringRepresentation(
                "truF", "PRENYLTRANSFERASE", [["O_13"], ["O_5"]], "3_METHYL_1_BUTENYL"
            ),
        ],
    )

    sacti_peptide_cluster_thurincin = RiPPCluster(
        "thnA",
        "METPVVQPRDWTCWSCLVCAACSVELLNLVTAATGASTAS",
        "DWTCWSCLVCAACSVELLNLVTAATGASTAS",
        tailoring_representations=[
            TailoringRepresentation(
                "thnB",
                "OXIDATIVE_BOND_SYNTHASE",
                [
                    ["S_109", "C_222"],
                    ["S_66", "C_203"],
                    ["S_90", "C_182"],
                    ["S_37", "C_156"],
                ],
            )
        ],
    )
    thiopeptide_cluster_thiomuracin = RiPPCluster(
        "tpdA",
        "MDLSDLPMDVFELADDGVAVESLTAGHGMTEVGASCNCFCYICCSCSSA",
        "SCNCFCYICCSCSS",
        tailoring_representations=[
            TailoringRepresentation(
                "tpdD",
                "THREONINE_SERINE_DEHYDRATASE",
                [["O_104"], ["O_90"], ["O_4"], ["O_111"]],
            ),
            TailoringRepresentation(
                "tpdE",
                "CYCLODEHYDRASE",
                [["S_97"], ["S_11"], ["S_76"], ["S_46"], ["S_83"], ["S_27"]],
            ),
            TailoringRepresentation("tpdF", "THIOPEPTIDE_CYCLASE", [["C_3", "C_89"]]),
            TailoringRepresentation("tpdF", "METHYLTRANSFERASE", [["C_26"]]),
            TailoringRepresentation("tpdF", "HYDROXYLASE", [["C_33"]]),
            TailoringRepresentation(
                "tpdF", "DEHYDROGENASE", [["C_67", "C_68"], ["C_3", "C_89"]]
            ),
            TailoringRepresentation("tpdF", "EPOXIDASE", [["C_67", "C_68"]]),
        ],
    )

    thiopeptide_cluster_thiomuracin.make_peptide()
    thiopeptide_cluster_thiomuracin.do_macrocyclization()
    thiopeptide_cluster_thiomuracin.do_tailoring()
    # Check if final structure is correct
    assert (
        structure_to_smiles(thiopeptide_cluster_thiomuracin.chain_intermediate)
        == r"NC12C=CC(C%10=N[C@H](C(NC(C(NC(C(O)=O)=C)=O)=C)=O)CS%10)=NC2(O)[C@H](CS4)N=C4[C@H](CS5)N=C5[C@H]([C@@H](C)C6CO6)NC(=O)[C@H](Cc(cc%11)ccc%11O)NC(=O)[C@H](CS7)N=C7[C@H](C(O)c8ccccc8)NC(=O)[C@H](C(C)S9)N=C9[C@H](CC(N)=O)NC(=O)[C@H]3N=C1SC3"
    )

    sacti_peptide_cluster_thurincin.make_peptide()
    sacti_peptide_cluster_thurincin.do_macrocyclization()
    sacti_peptide_cluster_thurincin.do_tailoring()

    # Check if final structure is correct
    assert (
        structure_to_smiles(sacti_peptide_cluster_thurincin.chain_intermediate)
        == r"N[C@H](C(N[C@H](C(N[C@H](C(N[C@H]3CSC(C4=O)(CC(N)=O)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCC(O)=O)NC(=O)[C@H](C(C)C)NC(=O)[C@H](CO)NC(=O)[C@H](CSC(C(N[C@H](C(N[C@H](C(N[C@H](C(O)=O)CO)=O)C)=O)[C@@H](C)O)=O)(CO)NC(=O)[C@H](C)NC(=O)CNC(=O)C6([C@@H](C)O)NC(=O)[C@H](C)NC(=O)[C@H](C)NC(=O)C5([C@@H](C)O)NC(=O)[C@H](C(C)C)NC(=O)[C@H](CC(C)C)N4)NC(=O)[C@H](C)NC(=O)[C@H](C)NC(=O)[C@H](CS5)NC(=O)[C@H](C(C)C)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CS6)NC(=O)[C@H](CO)NC(=O)[C@H](Cc7c[nH]c(cccc8)c78)NC3=O)=O)[C@@H](C)O)=O)Cc1c[nH]c(cccc2)c12)=O)CC(O)=O"
    )

    cyanobactin_cluster_trunkamide.make_peptide()
    cyanobactin_cluster_trunkamide.do_macrocyclization()
    cyanobactin_cluster_trunkamide.do_tailoring()

    # Check if final structure is correct
    assert (
        structure_to_smiles(cyanobactin_cluster_trunkamide.chain_intermediate)
        == r"C=CC(OC[C@@H]1NC(=O)[C@H]([C@@H](C)OC(C=C)(C)C)NC(=O)[C@H](CS2)N=C2[C@H](Cc3ccccc3)NC(=O)[C@H](CCC4)N4C(=O)[C@H](C)NC(=O)[C@H]([C@@H](C)CC)NC1=O)(C)C"
    )


if __name__ == "__main__":
    test_ripp()
