from raichu.cluster.ripp_cluster import RiPPCluster
from raichu.representations import MacrocyclizationRepresentation, TailoringRepresentation, CleavageSiteRepresentation
from raichu.run_raichu import draw_ripp_structure

cyanobactin_cluster_trunkamide = RiPPCluster("truE",
                                             "MNKKNILPQLGQPVIRLTAGQLSSQLAELSEEALGGVDASTSIAPFCSYDGVDASSYDGVDASSYDD",
                                             "TSIAPFC",
                                             macrocyclisations=[MacrocyclizationRepresentation("N_0", "O_59",
                                                                                               "condensative")],
                                             tailoring_representations=[TailoringRepresentation("truD",
                                                                                                "CYCLODEHYDRATION",
                                                                                                [["S_56"]]),
                                                                        TailoringRepresentation("truF",
                                                                                                "PRENYLTRANSFERASE",
                                                                                                [['O_13'], ['O_5']],
                                                                                                "3_METHYL_1_BUTENYL")])

draw_ripp_structure(cyanobactin_cluster_trunkamide, "cyanobactin")

sacti_peptide_cluster_thurincin = RiPPCluster("thnA", "METPVVQPRDWTCWSCLVCAACSVELLNLVTAATGASTAS",
                                              "DWTCWSCLVCAACSVELLNLVTAATGASTAS",
                                              tailoring_representations=[TailoringRepresentation(
                                                  "thnB", "OXIDATIVE_BOND_SYNTHASE",
                                                  [['S_109', "C_222"], ['S_66', "C_203"], ['S_90', "C_182"],
                                                   ['S_37', "C_156"]])]
                                                  )

# draw_ripp_structure(sacti_peptide_cluster_thurincin, "thurincin")


thiopeptide_cluster_thiomuracin = RiPPCluster("tpdA", "MDLSDLPMDVFELADDGVAVESLTAGHGMTEVGASCNCFCYICCSCSSA",
                                              "SCNCFCYICCSCSS",
                                              tailoring_representations=[
                                                  TailoringRepresentation("tpdD", "THREONINE_SERINE_DEHYDRATASE",
                                                                          [['O_104'], ['O_90'], ['O_4'], ['O_111']]),
                                                  TailoringRepresentation("tpdE", "CYCLODEHYDRATION",
                                                                          [['S_97'], ['S_11'], ['S_76'], ['S_46'],
                                                                           ['S_83'], ['S_27']]),
                                                  TailoringRepresentation("tpdF", "THIOPEPTIDE_CYCLASE",
                                                                          [["C_3", "C_89"]]),
                                                  TailoringRepresentation("tpdF", "METHYLTRANSFERASE", [["C_26"]]),
                                                  TailoringRepresentation("tpdF", "HYDROXYLASE", [["C_33"]]),
                                                  TailoringRepresentation(
                                                      "tpdF", "DEHYDROGENASE",
                                                      [["C_67", "C_68"], ["C_3", "C_89"]]),
                                                  TailoringRepresentation("tpdF", "EPOXIDASE", [["C_67", "C_68"]]),
                                                  TailoringRepresentation("tpdF", "MONOAMINE_OXIDASE", [["N_0"]]),
                                                  TailoringRepresentation(
                                                      "tpdF", "DEHYDRATASE", [["C_1", "C_84"]])

                                              ])
# draw_ripp_structure(thiopeptide_cluster_thiomuracin, "thiomuracin")

ripp_cluster = RiPPCluster("", "mkaekslkayawyiwyaha", "slkayawyiwy",
                           cleavage_sites=[CleavageSiteRepresentation("Y", 5, "end")])


draw_ripp_structure(ripp_cluster, "cleavage_example_ripp")
