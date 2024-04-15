from sys import argv
import os

from raichu.antismash import load_antismash_gbk
from raichu.general import draw_cluster_from_modular_cluster_representation, draw_products_from_modular_cluster_representation

cluster = load_antismash_gbk(argv[1])
cluster_path = argv[1].split('.gbk')[0]
if not os.path.exists(cluster_path):
    os.mkdir(cluster_path)
pathway_out = os.path.join(cluster_path, "pathway.svg")
cluster_out = os.path.join(cluster_path, "raichu_cluster")

cluster.write_cluster(cluster_out)
draw_cluster_from_modular_cluster_representation(cluster, pathway_out)
draw_products_from_modular_cluster_representation(cluster, cluster_path)