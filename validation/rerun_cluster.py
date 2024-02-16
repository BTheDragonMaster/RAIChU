from raichu.representations import ClusterRepresentation
from raichu.run_raichu import draw_cluster
from sys import argv
import os

cluster = ClusterRepresentation.from_file(argv[1])

out_file = os.path.join(argv[1], "cluster.svg")

draw_cluster(cluster, out_file)