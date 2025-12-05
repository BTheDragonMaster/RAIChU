from raichu.representations import ClusterRepresentation
from raichu.run_raichu import draw_cluster, write_linear_product
from raichu.general import draw_cluster_from_modular_cluster_representation, \
    draw_products_from_modular_cluster_representation
from sys import argv
import os

colour_dictionary = {"module_01": '#000000',
                     "module_02": '#D5BFDE',
                     "module_03": "#A075B1",
                     "module_04": "#1475BC",
                     "module_05": "#662483",
                     "module_06": "#818386",
                     "module_07": "#DEDCF0",
                     "module_08": "#1475BC",
                     "module_09": "#744696",
                     "module_10": "#1475BC",
                     "module_11": "#818386",
                     "module_12": "#0E5B88",
                     "module_13": "#4AA0D9",
                     "module_14": "#595E67"}




cluster = ClusterRepresentation.from_file(argv[1])


out_path = argv[1]
if not os.path.exists(out_path):
    os.mkdir(out_path)

out_pathway = os.path.join(out_path, "pathway.svg")
out_smiles = os.path.join(out_path, "smiles.txt")

out_file = os.path.join(argv[1], "cluster.svg")
draw_cluster(cluster, out_file, colour_by_module=True, colour_dictionary=colour_dictionary)

draw_cluster_from_modular_cluster_representation(cluster, out_pathway)
draw_products_from_modular_cluster_representation(cluster, out_path)
write_linear_product(cluster, out_smiles)
