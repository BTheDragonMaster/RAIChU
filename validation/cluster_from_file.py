from raichu.representations import ClusterRepresentation
from raichu.run_raichu import build_cluster
from sys import argv
import os

file_location = argv[1]
file_name = os.path.basename(file_location)


cluster_representation = ClusterRepresentation.from_file(file_location)
cluster = build_cluster(cluster_representation, strict=True)

out_path = os.path.join(file_location, "pathway.svg")


cluster.compute_structures(compute_cyclic_products=False)
cluster.do_tailoring()
cluster.draw_cluster(as_string=False, out_file=out_path, colour_by_module=True)
cluster.cyclise_all()
cluster.draw_all_products(file_location)
