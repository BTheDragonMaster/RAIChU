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

# paper example figure 3

cluster = ClusterRepresentation([ModuleRepresentation('PKS', 'PKS_CIS', "PROPIONYL_COA",
                                                      [DomainRepresentation('gene 1', 'AT'),
                                                       DomainRepresentation('gene 1', 'ACP')]),
                                 ModuleRepresentation('PKS', 'PKS_CIS', "METHYLMALONYL_COA",
                                                      [DomainRepresentation('gene 1', 'KS'),
                                                       DomainRepresentation('gene 1', 'AT'),
                                                       DomainRepresentation('gene 1', 'ACP')]),
                                 ModuleRepresentation("NRPS", None, "valine",
                                                      [DomainRepresentation("gene 1", 'C'),
                                                       DomainRepresentation("gene 1", 'A'),
                                                       DomainRepresentation("gene 1", 'PCP')]),
                                 ModuleRepresentation('PKS', 'PKS_CIS', "MALONYL_COA",
                                                      [DomainRepresentation('gene 1', 'KS'),
                                                       DomainRepresentation('gene 1', 'AT'),
                                                       DomainRepresentation('gene 1', 'DH'),
                                                       DomainRepresentation('gene 1', 'ER'),
                                                       DomainRepresentation('gene 1', 'KR', 'A1'),
                                                       DomainRepresentation('gene 1', 'ACP'),
                                                       DomainRepresentation('gene 1', 'TE')]),
                                 ])

draw_cluster(cluster, 'paper_figure_3.svg')


