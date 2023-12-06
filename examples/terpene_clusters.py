from raichu.representations import MacrocyclizationRepresentation, TailoringRepresentation
from raichu.cluster.terpene_cluster import TerpeneCluster
from raichu.run_raichu import get_tailoring_sites_atom_names

terpene_lycosantalonol = TerpeneCluster("TPS21", "GERANYLGERANYL_PYROPHOSPHATE",
                                        macrocyclisations=[MacrocyclizationRepresentation(
                                            "C_15", "C_11"), MacrocyclizationRepresentation(
                                            "C_10", "C_12"), MacrocyclizationRepresentation(
                                            "C_9", "C_14")],
                                        cyclase_type="Class_1",
                                        tailoring_representations=[
                                            TailoringRepresentation("CYP71BN1", "EPOXIDASE",[["C_19","C_20"]]),
                                            TailoringRepresentation("CYP71BN1", "HYDROXYLASE",[["C_19"]]),
                                            TailoringRepresentation("CYP71BN1", "REDUCTIVE_BOND_BREAKAGE",[["C_19", "O_5"]]),
                                            TailoringRepresentation("CYP71BN1", "ALCOHOL_DEHYDROGENASE", [["O_4"]])
                                            ])
terpene_lycosantalonol.create_precursor()
terpene_lycosantalonol.do_macrocyclization()
terpene_lycosantalonol.do_tailoring()
terpene_lycosantalonol.draw_pathway(order=("cyclisation", "tailoring"),out_file= "lycosantalonol_summary.svg", summary=True)
terpene_lycosantalonol.draw_pathway(order=("cyclisation", "tailoring"),out_file= "lycosantalonol_pathway.svg")
