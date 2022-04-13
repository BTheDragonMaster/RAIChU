from visualize_cluster import *
from nrps_tailoring_reactions import nrps_epimerization, nrps_methylation
from pikachu.general import draw_structure


clust = [['module 1', 'starter_module_nrps', 'threonine'],
         ['module 2', 'elongation_module_nrps', 'glycine'],
         ['module 3', 'elongation_module_pks', 'malonylcoa', ['KR']],
         ['module 4', 'elongation_module_pks', 'pk', []]]
#
# clust = [['module 1', 'starter_module_pks', 'SC(=O)CC'],
#          ['module 2', 'elongation_module_pks', 'ethylmalonylcoa', ['KR']],
#          ['module 3', 'elongation_module_pks', 'methylmalonylcoa', ['KR']],
#          ['module 4', 'elongation_module_pks', 'methylmalonylcoa', ['KR']],
#          ['module 5', 'elongation_module_pks', 'pk', []]]

clust = [['module 1', 'starter_module_pks', 'SC(C1=CC=CC=C1)=O'],
         ['module 2', 'elongation_module_pks', 'methylmalonylcoa', []],
         ['module 3', 'elongation_module_pks', 'methylmalonylcoa', ['KR', 'DH', 'ER']],
         ['module 4', 'elongation_module_nrps', 'glutamine', ['E']],
         ['module 5', 'elongation_module_pks', 'malonylcoa', ['KR', 'DH']],
         ['module 6', 'elongation_module_nrps', 'valine', ['E', 'nMT']],
         ['module 7', 'terminator_module_nrps', 'serine', []]]

clust2 = [['module 1', 'starter_module_nrps', 'glycine'],
          ['module 2', 'terminator_module_nrps', 'serine', ['E']],
          ['module 3', 'terminator_module_nrps', 'asparagine', ['nMT', 'E']],
          ['module 4', 'terminator_module_nrps', 'valine', []]]

draw_cluster(clust, interactive=True)
# draw_cluster(clust2, interactive=True)
# product = cluster_to_structure(clust2)
# for atom in product.graph:
#     print(atom, atom.annotations.chiral_c_ep)

# ep_product = nrps_epimerization(product)
# RaichuDrawer(ep_product)
# draw_cluster(clust)
#
# # asparagine = read_smiles('CC(C)[C@@H](C(=O)O)N')
# asparagine = read_smiles(('C([C@@H](C(=O)O)N)C(=O)N'))
# asparagine.add_attributes(ATTRIBUTES, boolean=True)
# product = condensation_nrps(asparagine, intermediate)
#
#
# product = attach_to_domain_nrp(product, 'PCP')
#
# RaichuDrawer(product)