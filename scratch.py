from pks_modules_to_structure import *
from pks_thioesterase_reactions import *
from pikachu.general import structure_to_smiles
from visualize_pks_cluster import draw_pks_cluster


erythromycin_cluster = [['module_1', 'starter_module', 'SC(=O)CC'],
                       ['module_2', 'elongation_module', 'methylmalonylcoa', ['KR_B2']],
                       ['module_3', 'elongation_module', 'methylmalonylcoa', ['KR_A1']],
                       ['module_4', 'elongation_module', 'methylmalonylcoa', ['KR_C2']],
                       ['module_5', 'elongation_module', 'methylmalonylcoa', ['KR', 'DH', 'ER']],
                       ['module_6', 'elongation_module', 'methylmalonylcoa', ['KR_A1']],
                       ['module_7', 'terminator_module', 'methylmalonylcoa', ['KR_A1']]]
product = pks_cluster_to_structure(erythromycin_cluster)
# thioesterase_all_products(product)

#draw_pks_cluster(erythromycin_cluster, interactive=True)

weird_structure_cluster =  [['pks module 1 [5148:6510]', 'starter_module', 'SC(CC(O)=O)=O'],
                           ['pks module 2 [6591:11160]', 'elongation_module', 'methylmalonylcoa', ['KR_B1']],
                           ['pks module 3 [11277:15915]', 'elongation_module', 'malonylcoa', ['KR_A1']],
                           ['pks module 4 [16029:19218]', 'elongation_module', 'methylmalonylcoa', []],
                           ['pks module 5 [19606:24247]', 'elongation_module', 'methylmalonylcoa', ['KR_A2']],
                           ['pks module 6 [24325:29008]', 'elongation_module', 'ethylmalonylcoa', ['KR_B1']],
                           ['pks module 7 [29086:34684]', 'elongation_module', 'malonylcoa', ['KR', 'DH']]]

#weird_product = pks_cluster_to_structure(weird_structure_cluster)


cluster =                  [['pks module 1 [5148:6510]', 'starter_module', 'SC(CC(O)=O)=O'],
                           ['pks module 2 [6591:11160]', 'elongation_module', 'methylmalonylcoa', ['KR_A2']],
                           ['pks module 3 [6591:11160]', 'elongation_module', 'ethylmalonylcoa', ['KR_B1']],
                           ['pks module 4 [11277:15915]', 'elongation_module', 'malonylcoa', ['KR' , 'DH']]]
#pks_cluster_to_structure(cluster)

wrong_cluster =            [['pks module 1', 'starter_module', 'SC(CC(O)=O)=O'],
                           ['pks module 2', 'elongation_module', 'methylmalonylcoa', ['KR_B1']],
                           ['pks module 3', 'elongation_module', 'malonylcoa', ['KR_A1']],
                           ['pks module 4', 'elongation_module', 'methylmalonylcoa', []],
                           ['pks module 5', 'elongation_module', 'methylmalonylcoa', ['KR_A2']],
                           ['pks module 6', 'elongation_module', 'ethylmalonylcoa', ['KR_B1']],
                           ['pks module 7', 'elongation_module', 'malonylcoa', ['KR', 'DH']]]

#pks_cluster_to_structure(wrong_cluster, attach_to_acp=True)





