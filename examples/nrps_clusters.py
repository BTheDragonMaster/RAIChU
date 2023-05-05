from raichu.run_raichu import ClusterRepresentation, ModuleRepresentation, DomainRepresentation, draw_cluster

# marformycin A

marformycin_cluster = ClusterRepresentation([ModuleRepresentation("NRPS", None, "valine",
                                                               [DomainRepresentation("MfnC", 'A', None, None, True,
                                                                                     True),
                                                                DomainRepresentation("MfnC", 'PCP', None, None, True,
                                                                                     True)
                                                                ]),
                                             ModuleRepresentation("NRPS", None, "threonine",
                                                                  [DomainRepresentation("MfnC", 'C', None, None, True,
                                                                                        True),
                                                                   DomainRepresentation("MfnC", 'A', None, None, True,
                                                                                        True),
                                                                   DomainRepresentation("MfnC", 'PCP', None, None, True,
                                                                                        True)
                                                                   ]),
                                             ModuleRepresentation("NRPS", None, "tyrosine",
                                                               [DomainRepresentation("MfnC", 'C', None, None, True,
                                                                                     True),
                                                                DomainRepresentation("MfnC", 'A', None, None, True,
                                                                                     True),
                                                                DomainRepresentation("MfnC", 'PCP', None, None, True,
                                                                                     True),

                                                                DomainRepresentation("MfnC", 'E', None, None, True,
                                                                                     True),
                                                                ]),
                                             ModuleRepresentation("NRPS", None, "isoleucine",
                                                                  [DomainRepresentation("MfnD", 'C', None, None, True,
                                                                                        True),
                                                                   DomainRepresentation("MfnD", 'A', None, None, True,
                                                                                        True),
                                                                   DomainRepresentation("MfnD", 'PCP', None, None, True,
                                                                                        True),
                                                                   DomainRepresentation("MfnD", 'E', None, None, True,
                                                                                        True)
                                                                   ]),
                                             ModuleRepresentation("NRPS", None, "piperazic acid",
                                                                  [DomainRepresentation("MfnE", 'C', None, None, True,
                                                                                        True),
                                                                   DomainRepresentation("MfnE", 'A', None, None, True,
                                                                                        True),
                                                                   DomainRepresentation("MfnE", 'PCP', None, None, True,
                                                                                        True)
                                                                   ]),
                                             ModuleRepresentation("NRPS", None, "leucine",
                                                                  [DomainRepresentation("MfnE", 'A', None, None, True,
                                                                                        True),
                                                                   DomainRepresentation("MfnE", 'PCP', None, None,
                                                                                        True, False),
                                                                   DomainRepresentation("MfnE", 'C', None, None, True,
                                                                                        True),
                                                                   DomainRepresentation("MfnE", 'PCP', None, None,
                                                                                        True, True)
                                                                   ]),
                                             ModuleRepresentation("NRPS", None, "valine",
                                                                  [DomainRepresentation("MfnE", 'C', None, None, True,
                                                                                        True),
                                                                   DomainRepresentation("MfnE", 'A', None, None, True,
                                                                                        True),
                                                                   DomainRepresentation("MfnE", 'nMT', None, None, True,
                                                                                        True),
                                                                   DomainRepresentation("MfnE", 'PCP', None, None, True,
                                                                                        True),
                                                                   DomainRepresentation("MfnE", 'TE', None, None, True,
                                                                                        True)
                                                                   ])
                                             ]
                                            )

draw_cluster(marformycin_cluster, "marformycin_A.svg")
