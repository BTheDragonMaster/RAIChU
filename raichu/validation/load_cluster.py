from raichu.representations import ClusterRepresentation
from raichu.run_raichu import draw_cluster

from sys import argv

cluster_repr = ClusterRepresentation.from_file(argv[1])
draw_cluster(cluster_repr)