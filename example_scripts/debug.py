from raichu.visualize_cluster import *
from pikachu.general import draw_structure, structure_to_smiles, draw_smiles, svg_from_structure
import os


clust = [["Module 1", "starter_module_pks", "SC(=O)CC"],
         ["Module 2", "elongation_module_pks", "methylmalonylcoa", ["KR"]],
         ["Module 3", "elongation_module_nrps", "Threonine", ['E', 'nMT']],
         ["Module 4", "elongation_module_pks", "malonylcoa", ["KR", "DH", "ER"]]]

clust = [["Module 1", "starter_module_pks", "SC(=O)CC"],
         ["Module 2", "elongation_module_pks", "malonylcoa", ["KR", "DH", "ER"]],
         ["Module 3", "elongation_module_nrps", "Threonine", ['E', 'nMT']],]

erythromycin_bgc = [['module 0', 'starter_module_pks', 'SC(=O)CC'],
                    ['module 1', 'elongation_module_pks', 'methylmalonylcoa', ['KR_B2']],
                    ['module 2', 'elongation_module_pks', 'methylmalonylcoa', ['KR_A1']],
                    ['module 3', 'elongation_module_pks', 'methylmalonylcoa', ['KR_C2']],
                    ['module 4', 'elongation_module_pks', 'methylmalonylcoa', ['KR', 'DH', 'ER']],
                    ['module 5', 'elongation_module_pks', 'methylmalonylcoa', ['KR_A1']],
                    ['module 6', 'terminator_module_pks', 'methylmalonylcoa', ['KR_A1']]]


# draw_cluster(clust, interactive=True)

products = cluster_to_structure(clust, visualization_mechanism=True, attach_to_acp=True, draw_mechanism_per_module=True)
print(products)
product_nr = 0
for i, product_group in enumerate(products):
    for product in product_group:
        product_nr += 1
        product.dont_show = False
        for atom in product.structure.graph:
            atom.draw.colour = 'black'

        product.save_svg = os.path.join(os.getcwd(), f"figure_1_product_{product_nr}.svg")
        product.draw_structure()


all_products = thioesterase_all_products(products[-1][0].structure, out_folder="products")
#
# RaichuDrawer(lin_te)
# draw_structure(lin_te)