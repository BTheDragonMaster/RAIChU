from raichu.representations import ClusterRepresentation
from raichu.run_raichu import draw_cluster
from raichu.general import draw_cluster_from_modular_cluster_representation, draw_products_from_modular_cluster_representation
from sys import argv
import os

cluster = ClusterRepresentation.from_file(argv[1])


out_path = argv[1]
if not os.path.exists(out_path):
    os.mkdir(out_path)

out_pathway = os.path.join(out_path, "pathway.svg")

draw_cluster_from_modular_cluster_representation(cluster, out_pathway)
draw_products_from_modular_cluster_representation(cluster, out_path)

out_file = os.path.join(argv[1], "cluster.svg")

draw_cluster(cluster, out_file)