from raichu.representations import ClusterRepresentation
from raichu.run_raichu import build_cluster
from pikachu.general import structure_to_smiles
from sys import argv
import os

file_location = argv[1]
file_name = os.path.basename(file_location)


cluster_representation = ClusterRepresentation.from_file(file_location)
cluster = build_cluster(cluster_representation, strict=True)

out_path = os.path.join(file_location, "pathway.svg")
out_intermediates = os.path.join(file_location, "intermediates.txt")


cluster.compute_structures(compute_cyclic_products=False)
with open(out_intermediates, 'w') as out_i:
    out_i.write('module_nr\tsmiles\n')
    for i, intermediate in enumerate(cluster.modular_intermediates):
        smiles = structure_to_smiles(intermediate)
        out_i.write(f"{i}\t{smiles}\n")

cluster.do_tailoring()
cluster.draw_cluster(as_string=False, out_file=out_path, colour_by_module=True)
cluster.cyclise_all()
cluster.draw_all_products(file_location)
