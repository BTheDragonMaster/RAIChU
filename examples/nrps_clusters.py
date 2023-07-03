from raichu.run_raichu import ClusterRepresentation, ModuleRepresentation, DomainRepresentation, draw_cluster, \
    draw_products, build_cluster, TailoringRepresentation, get_tailoring_sites
from pikachu.general import svg_from_structure

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

#draw_cluster(marformycin_cluster, "marformycin_A.svg")

hormaomycin_cluster = ClusterRepresentation([ModuleRepresentation('NRPS', None, "5‐chloropyrrole‐2‐carboxylic acid",
                                                                  [DomainRepresentation('HrmK', 'A'),
                                                                   DomainRepresentation('HrmL', 'PCP')]),
                                             ModuleRepresentation('NRPS', None, "3-[(1r,2r)-2-nitrocyclopropyl]-l-alanine",
                                                                  [DomainRepresentation('HrmO', 'C'),
                                                                   DomainRepresentation('HrmO', 'A'),
                                                                   DomainRepresentation('HrmO', 'PCP')]),
                                             ModuleRepresentation('NRPS', None, "threonine",
                                                                  [DomainRepresentation('HrmO', 'C'),
                                                                   DomainRepresentation('HrmO', 'A'),
                                                                   DomainRepresentation('HrmO', 'PCP'),
                                                                   DomainRepresentation('HrmO', 'E')]),
                                             ModuleRepresentation('NRPS', None, "R-beta-methylphenylalanine",
                                                                  [DomainRepresentation('HrmO', 'C'),
                                                                   DomainRepresentation('HrmO', 'A'),
                                                                   DomainRepresentation('HrmO', 'PCP')]),
                                             ModuleRepresentation('NRPS', None, "3-[(1r,2r)-2-nitrocyclopropyl]-l-alanine",
                                                                  [DomainRepresentation('HrmO', 'C'),
                                                                   DomainRepresentation('HrmO', 'A'),
                                                                   DomainRepresentation('HrmO', 'PCP'),
                                                                   DomainRepresentation('HrmO', 'E')]),
                                             ModuleRepresentation('NRPS', None, "R-beta-methylphenylalanine",
                                                                  [DomainRepresentation('HrmP', 'C'),
                                                                   DomainRepresentation('HrmP', 'A'),
                                                                   DomainRepresentation('HrmP', 'PCP')]),
                                             ModuleRepresentation('NRPS', None, "isoleucine",
                                                                  [DomainRepresentation('HrmP', 'C'),
                                                                   DomainRepresentation('HrmP', 'A'),
                                                                   DomainRepresentation('HrmP', 'PCP')]),
                                             ModuleRepresentation('NRPS', None, "4S-propenylproline",
                                                                  [DomainRepresentation('HrmP', 'C'),
                                                                   DomainRepresentation('HrmP', 'A'),
                                                                   DomainRepresentation('HrmP', 'PCP')])],
                                            tailoring_enzymes=[
                                                TailoringRepresentation(
                                                    "unknown", "HYDROXYLATION", [["N_135"]])
                                             ])
cluster = build_cluster(hormaomycin_cluster, strict=False)
cluster.compute_structures(compute_cyclic_products=False)
print(get_tailoring_sites(cluster.chain_intermediate,
                          enzyme_name='N_METHYLTRANSFERASE', out_file="hormaomycin_n.svg"))
draw_cluster(hormaomycin_cluster, "hormaomycin.svg")
print("drawn cluster")
draw_products(hormaomycin_cluster, "hormaomycin")
print("Done hormaomycin")

figure_1_example = ClusterRepresentation([ModuleRepresentation('NRPS', None, "threonine",
                                                               [DomainRepresentation('gene 1', 'A'),
                                                                DomainRepresentation('gene 1', 'PCP')]),
                                          ModuleRepresentation('NRPS', None, "glycine",
                                                               [DomainRepresentation('gene 1', 'C'),
                                                                DomainRepresentation('gene 1', 'A'),
                                                                DomainRepresentation('gene 1', 'PCP',),
                                                                DomainRepresentation('gene 1', 'TE')])])

draw_cluster(figure_1_example, "figure_1_example_TE.svg")
draw_products(figure_1_example, "figure_1_example_TE")

figure_1_example = ClusterRepresentation([ModuleRepresentation('NRPS', None, "glycine",
                                                               [DomainRepresentation('gene 1', 'A'),
                                                                DomainRepresentation('gene 1', 'PCP')]),
                                          ModuleRepresentation('NRPS', None, "threonine",
                                                               [DomainRepresentation('gene 1', 'C'),
                                                                DomainRepresentation('gene 1', 'CYC'),
                                                                DomainRepresentation('gene 1', 'A'),
                                                                DomainRepresentation('gene 1', 'OX'),
                                                                DomainRepresentation('gene 1', 'PCP',),
                                                                DomainRepresentation('gene 1', 'TE')])])

draw_cluster(figure_1_example, "figure_1_example_cycox.svg")
draw_products(figure_1_example, "figure_1_example_cycox")

figure_1_example = ClusterRepresentation([ModuleRepresentation('NRPS', None, "glycine",
                                                               [DomainRepresentation('gene 1', 'A'),
                                                                DomainRepresentation('gene 1', 'PCP')]),
                                          ModuleRepresentation('NRPS', None, "threonine",
                                                               [DomainRepresentation('gene 1', 'C'),
                                                                DomainRepresentation('gene 1', 'CYC'),
                                                                DomainRepresentation('gene 1', 'A'),
                                                                DomainRepresentation('gene 1', 'PCP',),
                                                                DomainRepresentation('gene 1', 'TE')])])

draw_cluster(figure_1_example, "figure_1_example_cyc.svg")
draw_products(figure_1_example, "figure_1_example_cyc")

figure_1_example = ClusterRepresentation([ModuleRepresentation('NRPS', None, "glycine",
                                                               [DomainRepresentation('gene 1', 'A'),
                                                                DomainRepresentation('gene 1', 'PCP')]),
                                          ModuleRepresentation('NRPS', None, "threonine",
                                                               [DomainRepresentation('gene 1', 'C'),
                                                                DomainRepresentation('gene 1', 'A'),
                                                                DomainRepresentation('gene 1', 'PCP',),
                                                                DomainRepresentation('gene 1', 'TE')])])

draw_cluster(figure_1_example, "figure_1_example_glythr.svg")
draw_products(figure_1_example, "figure_1_example_glythr")

figure_1_example = ClusterRepresentation([ModuleRepresentation('NRPS', None, "threonine",
                                                               [DomainRepresentation('gene 1', 'A'),
                                                                DomainRepresentation('gene 1', 'PCP')]),
                                          ModuleRepresentation('NRPS', None, "glycine",
                                                               [DomainRepresentation('gene 1', 'C'),
                                                                DomainRepresentation('gene 1', 'A'),
                                                                DomainRepresentation('gene 1', 'PCP',),
                                                                DomainRepresentation('gene 1', 'TD')])])

draw_cluster(figure_1_example, "figure_1_example_TD.svg")
draw_products(figure_1_example, "figure_1_example_TD")
