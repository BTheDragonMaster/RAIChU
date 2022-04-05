"""
This script contains examples of the general functionalities of RAIChU
"""
from visualize_cluster import *


if __name__ == "__main__":
    # Simple elongation reaction
    input_polyketide_intermediate = read_smiles('SC(=O)CC')
    Drawer(input_polyketide_intermediate)
    product = add_malonylunit(input_polyketide_intermediate)
    # Or: add_methylmalonylunit, add_ethylmalonylunit, add_methoxymalonylunit
    Drawer(product)

    # Simple tailoring reactions
    product2 = ketoreductase(product)
    Drawer(product2)
    product3 = dehydratase(product2)
    Drawer(product3)
    product4 = enoylreductase(product3)
    Drawer(product4)


    #Calculate structure directly from input PKS cluster
    erythromycin_cluster = [['module_1', 'starter_module', 'SC(=O)CC'],
                           ['module_2', 'elongation_module', 'methylmalonylcoa', ['KR_B2']],
                           ['module_3', 'elongation_module', 'methylmalonylcoa', ['KR_A1']],
                           ['module_4', 'elongation_module', 'methylmalonylcoa', ['KR_C2']],
                           ['module_5', 'elongation_module', 'methylmalonylcoa', ['KR', 'DH', 'ER']],
                           ['module_6', 'elongation_module', 'methylmalonylcoa', ['KR_A1']],
                           ['module_7', 'terminator_module', 'methylmalonylcoa', ['KR_A1']]]
    final_product = cluster_to_structure(erythromycin_cluster, attach_to_acp=True)

    #Visualise PKS cluster
    draw_pks_cluster(erythromycin_cluster)

    #Visualise PKS cluster (interactive mode)
    draw_pks_cluster(erythromycin_cluster, interactive=True)

    #Visualise NRPS cluster
    nrps_cluster = [['NRPS module 1', 'starter_module_nrps', 'd-threonine'],
                    ['NRPS module 2', 'elongation_module_nrps', 'valine'],
                    ['NRPS module 3', 'elongation_module_nrps', 'serine'],
                    ['NRPS module 4', 'elongation_module_nrps', 'cysteine'],
                    ['NRPS module 5', 'elongation_module_nrps', 'glutamicacid'],
                    ['NRPS module 6', 'elongation_module_nrps', 'alanine'],
                    ['NRPS module 7', 'terminator_module_nrps', 'valine']]
    draw_pks_cluster(nrps_cluster)

    #Visualise hybrid PKS/NRPS cluster
    nrps_cluster = [['PKS module 1', 'starter_module', 'SC(=O)CC'],
                    ['NRPS module 2', 'elongation_module_nrps', 'valine'],
                    ['NRPS module 3', 'elongation_module_nrps', 'serine'],
                    ['PKS module 4', 'elongation_module', 'malonylcoa', ['KR','DH']],
                    ['NRPS module 5', 'elongation_module_nrps', 'glutamicacid'],
                    ['PKS module 6', 'elongation_module', 'methoxymalonylacp', ['KR','DH','ER']],
                    ['NRPS module 7', 'terminator_module_nrps', 'valine']]
    draw_pks_cluster(nrps_cluster)