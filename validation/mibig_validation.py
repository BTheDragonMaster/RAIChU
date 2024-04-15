import os
from sys import argv

from raichu.antismash import load_antismash_gbk
from raichu.general import draw_cluster_from_modular_cluster_representation, \
    draw_products_from_modular_cluster_representation

input_dir = argv[1]
out_dir = argv[2]
if not os.path.exists(out_dir):
    os.mkdir(out_dir)
for genbank_file in os.listdir(input_dir):
    if genbank_file.endswith('.gbk'):
        path = os.path.join(input_dir, genbank_file)
        cluster = load_antismash_gbk(path)
        out_name = genbank_file.split('.gbk')[0]
        print(f"Working on {out_name}...")

        out_bgc = os.path.join(out_dir, out_name)
        if not os.path.exists(out_bgc):
            os.mkdir(out_bgc)

        out_path = os.path.join(out_bgc, "pathway.svg")
        cluster_out = os.path.join(out_bgc, "raichu_cluster")
        if not os.path.exists(cluster_out):
            os.mkdir(cluster_out)
        cluster.write_cluster(cluster_out)
        draw_cluster_from_modular_cluster_representation(cluster, out_path)
        draw_products_from_modular_cluster_representation(cluster, out_bgc)


