from joblib import dump
import os
from pathlib import Path

from raichu.run_raichu import build_cluster, ClusterRepresentation, ModuleRepresentation, DomainRepresentation
import raichu.test.test_data

# TEST_FOLDER = os.path.dirname(raichu.test.test_data.__file__)
TEST_FOLDER = Path(__file__).resolve().parent

TEST_DATA_FOLDER = TEST_FOLDER / "test_data"


MARFORMYCIN_CLUSTER = ClusterRepresentation([ModuleRepresentation("NRPS", None, "valine",
                                                                  [DomainRepresentation("MfnC", 'A', None, None,
                                                                                        True,
                                                                                        True),
                                                                   DomainRepresentation("MfnC", 'PCP', None, None,
                                                                                        True,
                                                                                        True)
                                                                   ]),
                                             ModuleRepresentation("NRPS", None, "threonine",
                                                                  [DomainRepresentation("MfnC", 'C', None, None,
                                                                                        True,
                                                                                        True),
                                                                   DomainRepresentation("MfnC", 'A', None, None,
                                                                                        True,
                                                                                        True),
                                                                   DomainRepresentation("MfnC", 'PCP', None, None,
                                                                                        True,
                                                                                        True)
                                                                   ]),
                                             ModuleRepresentation("NRPS", None, "tyrosine",
                                                                  [DomainRepresentation("MfnC", 'C', None, None,
                                                                                        True,
                                                                                        True),
                                                                   DomainRepresentation("MfnC", 'A', None, None,
                                                                                        True,
                                                                                        True),
                                                                   DomainRepresentation("MfnC", 'PCP', None, None,
                                                                                        True,
                                                                                        True),

                                                                   DomainRepresentation("MfnC", 'E', None, None,
                                                                                        True,
                                                                                        True),
                                                                   ]),
                                             ModuleRepresentation("NRPS", None, "isoleucine",
                                                                  [DomainRepresentation("MfnD", 'C', None, None,
                                                                                        True,
                                                                                        True),
                                                                   DomainRepresentation("MfnD", 'A', None, None,
                                                                                        True,
                                                                                        True),
                                                                   DomainRepresentation("MfnD", 'PCP', None, None,
                                                                                        True,
                                                                                        True),
                                                                   DomainRepresentation("MfnD", 'E', None, None,
                                                                                        True,
                                                                                        True)
                                                                   ]),
                                             ModuleRepresentation("NRPS", None, "piperazic acid",
                                                                  [DomainRepresentation("MfnE", 'C', None, None,
                                                                                        True,
                                                                                        True),
                                                                   DomainRepresentation("MfnE", 'A', None, None,
                                                                                        True,
                                                                                        True),
                                                                   DomainRepresentation("MfnE", 'PCP', None, None,
                                                                                        True,
                                                                                        True)
                                                                   ]),
                                             ModuleRepresentation("NRPS", None, "leucine",
                                                                  [DomainRepresentation("MfnE", 'A', None, None,
                                                                                        True,
                                                                                        True),
                                                                   DomainRepresentation("MfnE", 'PCP', None, None,
                                                                                        True, False),
                                                                   DomainRepresentation("MfnE", 'C', None, None,
                                                                                        True,
                                                                                        True),
                                                                   DomainRepresentation("MfnE", 'PCP', None, None,
                                                                                        True, True)
                                                                   ]),
                                             ModuleRepresentation("NRPS", None, "valine",
                                                                  [DomainRepresentation("MfnE", 'C', None, None,
                                                                                        True,
                                                                                        True),
                                                                   DomainRepresentation("MfnE", 'A', None, None,
                                                                                        True,
                                                                                        True),
                                                                   DomainRepresentation("MfnE", 'nMT', None, None,
                                                                                        True,
                                                                                        True),
                                                                   DomainRepresentation("MfnE", 'PCP', None, None,
                                                                                        True,
                                                                                        True),
                                                                   DomainRepresentation("MfnE", 'TE', None, None,
                                                                                        True,
                                                                                        True)
                                                                   ])
                                             ]
                                            )


def generate_test_nrps(cluster_representation, cluster_name):
    nrps_folder = os.path.join(TEST_DATA_FOLDER, "nrps")
    if not os.path.exists(nrps_folder):
        os.mkdir(nrps_folder)

    out_folder = os.path.join(nrps_folder, cluster_name)
    if not os.path.exists(out_folder):
        os.mkdir(out_folder)

    cluster = build_cluster(cluster_representation, strict=False)
    cluster.compute_structures(compute_cyclic_products=False)

    for i, intermediate in enumerate(cluster.modular_intermediates):
        dump(intermediate, os.path.join(out_folder, f"module {i + 1}.structure"))


if __name__ == "__main__":
    generate_test_nrps(MARFORMYCIN_CLUSTER, "marformycin")