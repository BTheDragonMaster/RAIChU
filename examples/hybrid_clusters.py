from raichu.run_raichu import ClusterRepresentation, ModuleRepresentation, DomainRepresentation, draw_cluster

# epoxomicin

epoxomicin_cluster = ClusterRepresentation([ModuleRepresentation("NRPS", None, "acetic acid",
                                                                 [DomainRepresentation("EpxD", 'A',
                                                                                       None, None, True, True),
                                                                  DomainRepresentation("EpxD", 'PCP',
                                                                                       None, None, True, True)
                                                                  ]),
                                            ModuleRepresentation("NRPS", None, "isoleucine",
                                                                 [DomainRepresentation("EpxD", 'C',
                                                                                       None, None, True, True),
                                                                  DomainRepresentation("EpxD", 'A',
                                                                                       None, None, True, True),
                                                                  DomainRepresentation("EpxD", 'nMT',
                                                                                       None, None, True, True),
                                                                  DomainRepresentation("EpxD", 'PCP',
                                                                                       None, None, True, True)
                                                                  ]),
                                            ModuleRepresentation("NRPS", None, "isoleucine",
                                                                 [DomainRepresentation("EpxD", 'C',
                                                                                       None, None, True, True),
                                                                  DomainRepresentation("EpxD", 'A',
                                                                                       None, None, True, True),
                                                                  DomainRepresentation("EpxD", 'PCP',
                                                                                       None, None, True, True)
                                                                  ]),
                                            ModuleRepresentation("NRPS", None, "threonine",
                                                                 [DomainRepresentation("EpxD", 'C',
                                                                                       None, None, True, True),
                                                                  DomainRepresentation("EpxD", 'A',
                                                                                       None, None, True, True),
                                                                  DomainRepresentation("EpxD", 'PCP',
                                                                                       None, None, True, True)
                                                                  ]),
                                            ModuleRepresentation("NRPS", None, "leucine",
                                                                 [DomainRepresentation("EpxD", 'C',
                                                                                       None, None, True, True),
                                                                  DomainRepresentation("EpxD", 'A',
                                                                                       None, None, True, True),
                                                                  DomainRepresentation("EpxD", 'PCP',
                                                                                       None, None, True, True)
                                                                  ]),
                                            ModuleRepresentation("PKS", "PKS_CIS", "MALONYL_COA",
                                                                 [DomainRepresentation("EpxE", 'KS',
                                                                                       None, None, True, True),
                                                                  DomainRepresentation("EpxE", 'AT',
                                                                                       None, None, True, True),
                                                                  DomainRepresentation("EpxE", 'UNKNOWN',
                                                                                       None, "cMT", True, True),
                                                                  DomainRepresentation("EpxE", 'ACP',
                                                                                       None, None, True, True),
                                                                  DomainRepresentation("EpxE", 'TE',
                                                                                       None, None, True, True)
                                                                  ])
                                            ]
                                           )

draw_cluster(epoxomicin_cluster, "epoxomicin.svg")

