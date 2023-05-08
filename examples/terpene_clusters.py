from raichu.representations import MacrocyclizationRepresentation, TailoringRepresentation
from raichu.cluster.terpene_cluster import TerpeneCluster

terpene_lycosantalonol = TerpeneCluster("BcBOT2", "GERANYLGERANYL_PYROPHOSPHATE",
                                        macrocyclisations=[MacrocyclizationRepresentation(
                                            "C_15", "C_11"), MacrocyclizationRepresentation(
                                            "C_10", "C_12"), MacrocyclizationRepresentation(
                                            "C_9", "C_14")],
                                        cyclase_type="Class_1",
                                        tailoring_representations=[TailoringRepresentation("cyclase", "DOUBLE_BOND_REDUCTION",
                                                                    [["C_20", "C_19"], ["C_15", "C_14"], ["C_11", "C_10"]]),
                                                                   TailoringRepresentation("hydroxylase", "HYDROXYLATION", [["C_19"], ["C_20"]]),
                                                                   TailoringRepresentation("dehydrogenase", "ALCOHOL_DEHYDROGENASE", [["O_71"]])])
terpene_lycosantalonol.create_precursor()
terpene_lycosantalonol.do_tailoring()
terpene_lycosantalonol.do_macrocyclization()
terpene_lycosantalonol.draw_pathway(out_file="lycosantalonol_pathway.svg")
terpene_lycosantalonol.draw_pathway(out_file="lycosantalonol_pathway_summary.svg", summary=True)
terpene_lycosantalonol.draw_product(out_file="lycosantalonol.svg")