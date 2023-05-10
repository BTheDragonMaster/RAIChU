from raichu.run_raichu import ClusterRepresentation, ModuleRepresentation, DomainRepresentation, draw_cluster, \
    draw_products

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


cluster = ClusterRepresentation([ModuleRepresentation('PKS', 'PKS_CIS', "PROPIONYL_COA",
                                                      [DomainRepresentation('gene 1', 'AT'),
                                                       DomainRepresentation('gene 1', 'ACP')]),
                                 ModuleRepresentation('PKS', 'PKS_CIS', "METHYLMALONYL_COA",
                                                      [DomainRepresentation('gene 1', 'KS'),
                                                       DomainRepresentation('gene 1', 'AT'),
                                                       DomainRepresentation('gene 1', 'KR', 'A1'),
                                                       DomainRepresentation('gene 1', 'ACP')]),
                                 ModuleRepresentation("NRPS", None, "threonine",
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

draw_cluster(cluster, 'paper_figure_4.svg')
draw_products(cluster, 'paper_figure_4')

leinamycin_cluster = ClusterRepresentation([ModuleRepresentation("NRPS", None, "alanine",
                                                                 [DomainRepresentation("LnmQ", 'A'),
                                                                  DomainRepresentation("LnmP", 'PCP')
                                                                  ]),
                                            ModuleRepresentation("NRPS", None, "cysteine",
                                                                 [DomainRepresentation("LnmI", 'C'),
                                                                  DomainRepresentation("LnmI", 'CYC'),
                                                                  DomainRepresentation("LnmI", 'A'),
                                                                  DomainRepresentation("LnmI", 'OX'),
                                                                  DomainRepresentation("LnmI", 'PCP')
                                                                  ]),
                                            ModuleRepresentation("PKS", 'PKS_TRANS', "MALONYL_COA",
                                                                 [DomainRepresentation("LnmI", 'KS', 'ST'),
                                                                  DomainRepresentation("LnmI", 'ACP')
                                                                  ]),
                                            ModuleRepresentation("PKS", 'PKS_TRANS', "MALONYL_COA",
                                                                 [DomainRepresentation("LnmI", 'KS', 'DB'),
                                                                  DomainRepresentation("LnmI", 'UNKNOWN', name="Do",
                                                                                       used=False),
                                                                  DomainRepresentation("LnmI", 'KR'),
                                                                  DomainRepresentation("LnmI", 'ACP')
                                                                  ]),
                                            ModuleRepresentation("PKS", 'PKS_TRANS', "MALONYL_COA",
                                                                 [DomainRepresentation("LnmI", 'KS', 'DB'),
                                                                  DomainRepresentation("LnmJ", 'DH'),
                                                                  DomainRepresentation("LnmJ", 'ACP'),
                                                                  DomainRepresentation("LnmJ", 'KR'),
                                                                  ]),
                                            ModuleRepresentation("PKS", 'PKS_TRANS', "MALONYL_COA",
                                                                 [DomainRepresentation("LnmJ", 'KS', 'KETO'),
                                                                  DomainRepresentation("LnmJ", 'UNKNOWN', name="Do",
                                                                                       used=False),
                                                                  DomainRepresentation("LnmJ", 'ACP'),
                                                                  ]),
                                            ModuleRepresentation("PKS", 'PKS_TRANS', "MALONYL_COA",
                                                                 [DomainRepresentation("LnmJ", 'KS', 'ALPHAME_EDB'),
                                                                  DomainRepresentation("LnmJ", 'UNKNOWN', name="Do",
                                                                                       used=False),
                                                                  DomainRepresentation("LnmJ", 'DH'),
                                                                  DomainRepresentation("LnmJ", 'KR'),
                                                                  DomainRepresentation("LnmJ", 'ACP', used=False),
                                                                  DomainRepresentation("LnmJ", 'UNKNOWN', name="cMT",
                                                                                       used=False),
                                                                  DomainRepresentation("LnmJ", 'ACP')
                                                                  ]),
                                            ModuleRepresentation("PKS", 'PKS_TRANS', "MALONYL_COA",
                                                                 [DomainRepresentation("LnmJ", 'KS', 'RED'),
                                                                  DomainRepresentation("LnmJ", 'UNKNOWN', name="Do",
                                                                                       used=False),
                                                                  DomainRepresentation("LnmJ", 'DH'),
                                                                  DomainRepresentation("LnmJ", 'KR'),
                                                                  DomainRepresentation("LnmJ", 'ACP')
                                                                  ]),
                                            ModuleRepresentation("PKS", 'PKS_TRANS', "MALONYL_COA",
                                                                 [DomainRepresentation("LnmJ", 'KS', 'RED'),
                                                                  DomainRepresentation("LnmJ", 'DH',
                                                                                       active=False, used=False),
                                                                  DomainRepresentation("LnmJ", 'ACP', used=False),
                                                                  DomainRepresentation("LnmJ", 'ACP'),
                                                                  DomainRepresentation("LnmJ", 'UNKNOWN',
                                                                                       used=False),
                                                                  DomainRepresentation("LnmJ", 'UNKNOWN', name="SH",
                                                                                       used=False),
                                                                  DomainRepresentation("LnmN", 'TE')

                                                                  ])
                                            ])

draw_cluster(leinamycin_cluster, 'leinamycin.svg')

cmc_thuggacin_cluster = ClusterRepresentation([ModuleRepresentation('PKS', 'PKS_CIS', "ACETYL_COA",
                                                                    [DomainRepresentation('TugA', 'AT'),
                                                                     DomainRepresentation('TugA', 'ACP')]
                                                                    ),
                                               ModuleRepresentation('PKS', 'PKS_CIS', "METHYLMALONYL_COA",
                                                                    [DomainRepresentation('TugA', 'KS'),
                                                                     DomainRepresentation('TugA', 'AT'),
                                                                     DomainRepresentation('TugA', 'DH'),
                                                                     DomainRepresentation('TugA', 'KR'),
                                                                     DomainRepresentation('TugA', 'ACP')]
                                                                    ),
                                               ModuleRepresentation('PKS', 'PKS_CIS', "METHYLMALONYL_COA",
                                                                    [DomainRepresentation('TugA', 'KS'),
                                                                     DomainRepresentation('TugA', 'AT'),
                                                                     DomainRepresentation('TugA', 'DH'),
                                                                     DomainRepresentation('TugA', 'KR', 'B1'),
                                                                     DomainRepresentation('TugA', 'ACP')]
                                                                    ),
                                               ModuleRepresentation('PKS', 'PKS_CIS', "METHYLMALONYL_COA",
                                                                    [DomainRepresentation('TugA', 'KS'),
                                                                     DomainRepresentation('TugA', 'AT'),
                                                                     DomainRepresentation('TugA', 'DH'),
                                                                     DomainRepresentation('TugA', 'ER'),
                                                                     DomainRepresentation('TugA', 'KR', 'B1'),
                                                                     DomainRepresentation('TugA', 'ACP'),
                                                                     DomainRepresentation('TugA', 'UNKNOWN', name="Do",
                                                                                          used=False)]
                                                                    ),
                                               ModuleRepresentation('PKS', 'PKS_CIS', "MALONYL_COA",
                                                                    [DomainRepresentation('TugB', 'KS'),
                                                                     DomainRepresentation('TugB', 'AT'),
                                                                     DomainRepresentation('TugB', 'KR', 'B1'),
                                                                     DomainRepresentation('TugB', 'ACP')]
                                                                    ),
                                               ModuleRepresentation('PKS', 'PKS_CIS', "MALONYL_COA",
                                                                    [DomainRepresentation('TugB', 'KS'),
                                                                     DomainRepresentation('TugB', 'AT'),
                                                                     DomainRepresentation('TugB', 'KR', 'B1'),
                                                                     DomainRepresentation('TugB', 'ACP'),
                                                                     DomainRepresentation('TugB', 'UNKNOWN', name="Do",
                                                                                          used=False)]
                                                                    ),
                                               ModuleRepresentation('PKS', 'PKS_CIS', "MALONYL_COA",
                                                                    [DomainRepresentation('TugC', 'UNKNOWN', name="Do",
                                                                                          used=False),
                                                                     DomainRepresentation('TugC', 'KS'),
                                                                     DomainRepresentation('TugC', 'AT'),
                                                                     DomainRepresentation('TugC', 'DH'),
                                                                     DomainRepresentation('TugC', 'KR', 'A1'),
                                                                     DomainRepresentation('TugC', 'ACP')]
                                                                    ),
                                               ModuleRepresentation('PKS', 'PKS_CIS', "MALONYL_COA",
                                                                    [DomainRepresentation('TugC', 'KS'),
                                                                     DomainRepresentation('TugC', 'AT'),
                                                                     DomainRepresentation('TugC', 'DH'),
                                                                     DomainRepresentation('TugC', 'KR', 'B1'),
                                                                     DomainRepresentation('TugC', 'ACP'),
                                                                     DomainRepresentation('TugC', 'UNKNOWN', name="Do",
                                                                                          used=False)]
                                                                    ),
                                               ModuleRepresentation('PKS', 'PKS_CIS', "MALONYL_COA",
                                                                    [DomainRepresentation('TugD', 'KS'),
                                                                     DomainRepresentation('TugD', 'AT'),
                                                                     DomainRepresentation('TugD', 'KR', 'B1'),
                                                                     DomainRepresentation('TugD', 'ACP')]
                                                                    ),
                                               ModuleRepresentation('PKS', 'PKS_CIS', "METHYLMALONYL_COA",
                                                                    [DomainRepresentation('TugD', 'KS'),
                                                                     DomainRepresentation('TugD', 'AT'),
                                                                     DomainRepresentation('TugD', 'KR', 'A1'),
                                                                     DomainRepresentation('TugD', 'ACP')]
                                                                    ),
                                               ModuleRepresentation("NRPS", None, "cysteine",
                                                                    [DomainRepresentation("TugD", 'C'),
                                                                     DomainRepresentation("TugD", 'CYC'),
                                                                     DomainRepresentation("TugD", 'A'),
                                                                     DomainRepresentation("TugD", 'OX'),
                                                                     DomainRepresentation("TugD", 'PCP')]
                                                                    ),
                                               ModuleRepresentation('PKS', 'PKS_CIS', "METHYLMALONYL_COA",
                                                                    [DomainRepresentation('TugD', 'KS'),
                                                                     DomainRepresentation('TugD', 'AT'),
                                                                     DomainRepresentation('TugD', 'DH'),
                                                                     DomainRepresentation('TugD', 'KR', 'B1'),
                                                                     DomainRepresentation('TugD', 'ACP')]
                                                                    )])

draw_cluster(cmc_thuggacin_cluster, 'cmc_thuggacin.svg')
draw_products(cmc_thuggacin_cluster, 'cmc_thuggacin')