from raichu.run_raichu import ClusterRepresentation, ModuleRepresentation, DomainRepresentation, draw_cluster, \
    draw_product, draw_products

# trans-AT

# cis-AT

# Phenylnannolone A

phenylnannolone_cluster = ClusterRepresentation([ModuleRepresentation("NRPS", None, "cinnamic acid",
                                                                      [DomainRepresentation("Phn2", 'A',
                                                                                            None, None, True, True),
                                                                       DomainRepresentation("Phn2", 'PCP',
                                                                                            None, None, True, True)]),
                                                 ModuleRepresentation("PKS", "PKS_CIS", "ETHYLMALONYL_COA",
                                                                      [DomainRepresentation("Phn2", 'KS',
                                                                                            None, None, True, True),
                                                                       DomainRepresentation("Phn2", 'AT',
                                                                                            None, None, True, True),
                                                                       DomainRepresentation("Phn2", 'DH',
                                                                                            None, None, True, True),
                                                                       DomainRepresentation("Phn2", 'KR',
                                                                                            None, None, True, True),
                                                                       DomainRepresentation("Phn2", 'ACP',
                                                                                            None, None, True, True)
                                                                       ]),
                                                 ModuleRepresentation("PKS", "PKS_CIS", "MALONYL_COA",
                                                                      [DomainRepresentation("Phn2", 'KS',
                                                                                            None, None, True, True),
                                                                       DomainRepresentation("Phn2", 'AT',
                                                                                            None, None, True, True),
                                                                       DomainRepresentation("Phn2", 'DH',
                                                                                            None, None, True, True),
                                                                       DomainRepresentation("Phn2", 'KR',
                                                                                            None, None, True, True),
                                                                       DomainRepresentation("Phn2", 'ACP',
                                                                                            None, None, True, True)]),
                                                 ModuleRepresentation("PKS", "PKS_CIS", "MALONYL_COA",
                                                                      [DomainRepresentation("Phn2", 'KS',
                                                                                            None, None, True, True),
                                                                       DomainRepresentation("Phn2", 'AT',
                                                                                            None, None, True, True),
                                                                       DomainRepresentation("Phn2", 'ACP',
                                                                                            None, None, True, True)]),
                                                 ModuleRepresentation("PKS", "PKS_CIS", "MALONYL_COA",
                                                                      [DomainRepresentation("Phn2", 'KS',
                                                                                            None, None, True, True),
                                                                       DomainRepresentation("Phn2", 'AT',
                                                                                            None, None, True, True),
                                                                       DomainRepresentation("Phn2", 'DH',
                                                                                            None, None, True, True),
                                                                       DomainRepresentation("Phn2", 'KR',
                                                                                            None, None, True, True),
                                                                       DomainRepresentation("Phn2", 'ACP',
                                                                                            None, None, True, True),
                                                                       DomainRepresentation("Phn2", 'TE',
                                                                                            None, None, True, True)
                                                                       ])

                                                 ]
                                                )
# erythromycin

erythromycin_cluster = ClusterRepresentation([ModuleRepresentation('PKS', 'PKS_CIS', 'PROPIONYL_COA',
                                                                   [DomainRepresentation("eryAI", 'AT'),
                                                                    DomainRepresentation("eryAI", 'ACP')
                                                                    ]),
                                              ModuleRepresentation('PKS', 'PKS_CIS', 'METHYLMALONYL_COA',
                                                                   [DomainRepresentation("eryAI", 'KS'),
                                                                    DomainRepresentation("eryAI", 'AT'),
                                                                    DomainRepresentation("eryAI", 'DH',
                                                                                         active=False, used=False),
                                                                    DomainRepresentation("eryAI", 'KR', 'B2'),
                                                                    DomainRepresentation("eryAI", 'ACP')
                                                                    ]),
                                              ModuleRepresentation('PKS', 'PKS_CIS', 'METHYLMALONYL_COA',
                                                                   [DomainRepresentation("eryAI", 'KS'),
                                                                    DomainRepresentation("eryAI", 'AT'),
                                                                    DomainRepresentation("eryAI", 'KR', 'A1'),
                                                                    DomainRepresentation("eryAI", 'ACP')
                                                                    ]),
                                              ModuleRepresentation('PKS', 'PKS_CIS', 'METHYLMALONYL_COA',
                                                                   [DomainRepresentation("eryAII", 'KS'),
                                                                    DomainRepresentation("eryAII", 'AT'),
                                                                    DomainRepresentation("eryAII", 'KR', 'C2'),
                                                                    DomainRepresentation("eryAII", 'ACP')
                                                                    ]),
                                              ModuleRepresentation('PKS', 'PKS_CIS', 'METHYLMALONYL_COA',
                                                                   [DomainRepresentation("eryAII", 'KS'),
                                                                    DomainRepresentation("eryAII", 'AT'),
                                                                    DomainRepresentation("eryAII", 'DH'),
                                                                    DomainRepresentation("eryAII", 'ER'),
                                                                    DomainRepresentation("eryAII", 'KR'),
                                                                    DomainRepresentation("eryAII", 'ACP')
                                                                    ]),
                                              ModuleRepresentation('PKS', 'PKS_CIS', 'METHYLMALONYL_COA',
                                                                   [DomainRepresentation("eryAIII", 'KS'),
                                                                    DomainRepresentation("eryAIII", 'AT'),
                                                                    DomainRepresentation("eryAIII", 'KR', 'A1'),
                                                                    DomainRepresentation("eryAIII", 'ACP')
                                                                    ]),
                                              ModuleRepresentation('PKS', 'PKS_CIS', 'METHYLMALONYL_COA',
                                                                   [DomainRepresentation("eryAIII", 'KS'),
                                                                    DomainRepresentation("eryAIII", 'AT'),
                                                                    DomainRepresentation("eryAIII", 'KR', 'A1'),
                                                                    DomainRepresentation("eryAIII", 'ACP'),
                                                                    DomainRepresentation("eryAIII", 'TE')
                                                                    ]),
                                              ])

draw_cluster(phenylnannolone_cluster, "phenylnannolone.svg")
draw_cluster(erythromycin_cluster, "erythromycin.svg")
draw_product(erythromycin_cluster, "erythromycin_product.svg")
draw_products(erythromycin_cluster, "erythromycin")
erythromycin_cluster.write_cluster("erythromycin")

# iterative