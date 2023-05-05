from raichu.run_raichu import ClusterRepresentation, ModuleRepresentation, DomainRepresentation, draw_cluster

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

draw_cluster(phenylnannolone_cluster, "phenylnannolone.svg")

# iterative