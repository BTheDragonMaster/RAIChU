from raichu.representations import MacrocyclizationRepresentation, TailoringRepresentation
from raichu.cluster.terpene_cluster import TerpeneCluster

terpene_lycosantalonol = TerpeneCluster("BcBOT2", "GERANYLGERANYL_PYROPHOSPHATE",
                                        macrocyclisations=[MacrocyclizationRepresentation(
                                            "C_15", "C_11"), MacrocyclizationRepresentation(
                                            "C_10", "C_12"), MacrocyclizationRepresentation(
                                            "C_9", "C_14")],
                                        terpene_cyclase_type="Class_1",
                                        tailoring_enzymes_representation=[
                                            TailoringRepresentation("hydroxylase", "HYDROXYLATION",[["C_19"],["C_20"]]),
                                            TailoringRepresentation("keto", "ALCOHOL_DEHYDROGENASE", [["O_71"]])])
terpene_lycosantalonol.create_precursor()

terpene_lycosantalonol.do_macrocyclization()
terpene_lycosantalonol.do_tailoring()
terpene_lycosantalonol.draw_product(as_string=False, out_file= "lycosantalonol.svg")
