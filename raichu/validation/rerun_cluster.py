from raichu.representations import ClusterRepresentation
from raichu.run_raichu import draw_cluster, write_linear_product
from raichu.general import draw_cluster_from_modular_cluster_representation, draw_products_from_modular_cluster_representation
from sys import argv
import os

cluster = ClusterRepresentation.from_file(argv[1])


out_path = argv[1]
if not os.path.exists(out_path):
    os.mkdir(out_path)

out_pathway = os.path.join(out_path, "pathway.svg")
out_smiles = os.path.join(out_path, "smiles.txt")

out_file = os.path.join(argv[1], "cluster.svg")
draw_cluster(cluster, out_file, colour_by_module=False)

draw_cluster_from_modular_cluster_representation(cluster, out_pathway)
draw_products_from_modular_cluster_representation(cluster, out_path)
write_linear_product(cluster, out_smiles)
