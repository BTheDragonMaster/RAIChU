import os
from sys import argv

from pikachu.math_functions import Vector

from raichu.representations import ClusterRepresentation
from raichu.general import build_cluster
from raichu.central_chain_detection.find_central_chain import reorder_central_chain, find_central_chain


def check_central_chains_from_file(cluster_folder):
    cluster_repr = ClusterRepresentation.from_file(cluster_folder)
    return check_central_chains(cluster_repr)


def check_central_chains(cluster_repr):

    cluster = build_cluster(cluster_repr, strict=False)
    cluster.compute_structures(compute_cyclic_products=False)
    cluster.do_tailoring()
    drawings, widths = cluster.get_spaghettis()
    faulty_spaghettis = 0
    total_backbone_errors = 0

    for drawing in drawings:
        carrier_domain_pos = None
        central_chain = find_central_chain(drawing.structure)
        backbone, _, _, _, last_atom_to_check = reorder_central_chain(central_chain, drawing)

        for atom in drawing.structure.graph:
            if atom.annotations.domain_type:
                carrier_domain_pos = atom.draw.position
                break

        assert carrier_domain_pos

        # According to the 30-60-90 triangle rule, atoms of the central chain should either have the same coordinate
        # as the carrier domain, or be half a bond length positioned to the right

        alternative_position = Vector.subtract_vectors(carrier_domain_pos, Vector(drawing.options.bond_length / 2, 0))
        max_deviation = float(drawing.options.bond_length) / 10

        errors_in_backbone = 0

        for atom in backbone:

            if atom == last_atom_to_check:
                break
            if alternative_position.x - max_deviation < atom.draw.position.x < alternative_position.x + max_deviation or \
                    carrier_domain_pos.x - max_deviation < atom.draw.position.x < carrier_domain_pos.x + max_deviation:
                continue
            else:
                errors_in_backbone += 1

        if errors_in_backbone:
            faulty_spaghettis += 1
            total_backbone_errors += errors_in_backbone

    central_chain_error = False

    if faulty_spaghettis:
        central_chain_error = True

    return central_chain_error, faulty_spaghettis, total_backbone_errors


def check_clusters(folder):

    faulty_clusters = 0

    for i, cluster_name in enumerate(os.listdir(folder)):
        if i % 10 == 0 and i != 0:
            print(f"Checked {i} clusters, of which {faulty_clusters} have an incorrectly drawn central chain.")
        cluster_folder = os.path.join(folder, cluster_name)
        if os.path.isdir(cluster_folder) and 'cluster' in cluster_name:
            central_chain_error, faulty_spaghettis, total_backbone_errors = check_central_chains_from_file(cluster_folder)
            if central_chain_error:
                faulty_clusters += 1
                print(f"Central chain error in {cluster_name}.")
                print(f"Faulty spaghettis: {faulty_spaghettis}")
                print(f"Total backbone errors: {total_backbone_errors}")
                print('\n')



    return faulty_clusters


if __name__ == "__main__":
    print(check_clusters(argv[1]))
