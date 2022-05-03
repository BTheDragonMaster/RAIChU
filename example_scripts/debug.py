from raichu.visualize_cluster import *
from pikachu.general import draw_structure, structure_to_smiles, draw_smiles, svg_from_structure


clust = [["Module 1", "starter_module_pks", "SC(=O)CC"],
         ["Module 2", "elongation_module_pks", "methylmalonylcoa", ["KR_A1"]],
         ["Module 3", "elongation_module_nrps", "Threonine", []],
         ["Module 4", "elongation_module_pks", "malonylcoa", ["KR", "DH", "ER"]]]

erythromycin_bgc = [['module 0', 'starter_module_pks', 'SC(=O)CC'],
                    ['module 1', 'elongation_module_pks', 'methylmalonylcoa', ['KR_B2']],
                    ['module 2', 'elongation_module_pks', 'methylmalonylcoa', ['KR_A1']],
                    ['module 3', 'elongation_module_pks', 'methylmalonylcoa', ['KR_C2']],
                    ['module 4', 'elongation_module_pks', 'methylmalonylcoa', ['KR', 'DH','ER']],
                    ['module 5', 'elongation_module_pks', 'methylmalonylcoa', ['KR_A1']],
                    ['module 6', 'terminator_module_pks', 'methylmalonylcoa', ['KR_A1']]]


# draw_cluster(clust, interactive=True)

product = cluster_to_structure(erythromycin_bgc, attach_to_acp=True)
drawer = RaichuDrawer(product, save_svg="erythromycin_linear.svg")
drawer.draw_structure()

#
# linear_product = thioesterase_linear_product(product)

all_products = thioesterase_all_products(product, out_folder="erythromycin")
#
# RaichuDrawer(lin_te)
# draw_structure(lin_te)