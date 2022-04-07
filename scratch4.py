from visualize_cluster import *

clust = [['module 1', 'starter_module_pks', 'SC(C1C(CCC1)C(=O)O)=O'],
         ['module 2', 'elongation_module_pks', 'methoxymalonylacp', ['KR', 'DH', 'ER']]]

#draw_cluster(clust)

start = read_smiles('SC(C1C(CCC1)C(=O)O)=O')
start.add_attributes(ATTRIBUTES, boolean=True)
start = attach_to_domain_pk(start, 'ACP')
start = add_methoxymalonylunit(start)
start = ketoreductase(start)


start = dehydratase(start)
print(start.graph)
RaichuDrawer(start)
# RaichuDrawer(start)
start = enoylreductase(start)
RaichuDrawer(start)