from visualize_cluster import *

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
         ['module 2', 'elongation_module_pks', 'methylmalonylcoa', ['KR', 'DH', 'ER']],
         ['module 3', 'elongation_module_nrps', 'glutamine'],
         ['module 4', 'elongation_module_pks', 'ethylmalonylcoa', ['KR']],
         ['module 5', 'elongation_module_nrps', 'glutamine'],
         ['module 6', 'elongation_module_pks', 'malonylcoa', ['KR', 'DH']],
         ['module 7', 'elongation_module_nrps', 'glycine'],
         ['module 8', 'terminator_module_nrps', 'serine']]

intermediate = cluster_to_structure(clust, attach_to_acp=True)

# asparagine = read_smiles('CC(C)[C@@H](C(=O)O)N')
asparagine = read_smiles(('C([C@@H](C(=O)O)N)C(=O)N'))
asparagine.add_attributes(ATTRIBUTES, boolean=True)
product = condensation_nrps(asparagine, intermediate)


product = attach_to_domain_nrp(product, 'PCP')

RaichuDrawer(product)