from visualize_cluster import *

clust = [['module 1', 'starter_module_nrps', 'threonine'],
         ['module 2', 'elongation_module_nrps', 'glycine'],
         ['module 3', 'elongation_module_pks', 'malonylcoa', ['KR']],
         ['module 4', 'elongation_module_pks', 'pk', []]]

intermediate = cluster_to_structure(clust, attach_to_acp=True)
#asparagine = read_smiles('CC(C)[C@@H](C(=O)O)N')+
asparagine = read_smiles(('C([C@@H](C(=O)O)N)C(=O)N'))
asparagine.add_attributes(ATTRIBUTES, boolean=True)
product = condensation_nrps(asparagine, intermediate)
print('product graph', product.graph)
print('product bond lookup', product.bond_lookup)

#attach_to_domain_nrp(product, 'PCP')
RaichuDrawer(product)