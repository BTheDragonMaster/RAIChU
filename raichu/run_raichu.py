import os

from raichu.cluster.modular_cluster import ModularCluster
from raichu.cluster.ripp_cluster import RiPPCluster
from raichu.cluster.terpene_cluster import TerpeneCluster
from raichu.domain.domain import (
    TailoringDomain,
    CarrierDomain,
    SynthesisDomain,
    RecognitionDomain,
    TerminationDomain,
    UnknownDomain,
    Domain,
)
from raichu.module import (
    PKSModuleSubtype,
    NRPSModule,
    LinearPKSModule,
    IterativePKSModule,
    TransATPKSModule,
    ModuleType,
)
from raichu.domain.domain_types import (
    TailoringDomainType,
    TerminationDomainType,
    CarrierDomainType,
    SynthesisDomainType,
    RecognitionDomainType,
)
from raichu.cluster.alkaloid_cluster import AlkaloidCluster
from raichu.tailoring_enzymes import TailoringEnzyme, TailoringEnzymeType
from raichu.representations import (
    ClusterRepresentation,
    DomainRepresentation,
    ModuleRepresentation,
    TailoringRepresentation,
    CleavageSiteRepresentation,
    MacrocyclizationRepresentation,
)

DOMAIN_TO_SUPERTYPE = {}
for domain_name in TailoringDomainType.__members__:
    DOMAIN_TO_SUPERTYPE[domain_name] = TailoringDomain
for domain_name in CarrierDomainType.__members__:
    DOMAIN_TO_SUPERTYPE[domain_name] = CarrierDomain
for domain_name in SynthesisDomainType.__members__:
    DOMAIN_TO_SUPERTYPE[domain_name] = SynthesisDomain
for domain_name in RecognitionDomainType.__members__:
    DOMAIN_TO_SUPERTYPE[domain_name] = RecognitionDomain
for domain_name in TerminationDomainType.__members__:
    DOMAIN_TO_SUPERTYPE[domain_name] = TerminationDomain
DOMAIN_TO_SUPERTYPE["UNKNOWN"] = UnknownDomain


def make_domain(
    domain_repr: DomainRepresentation, substrate: str, strict: bool = True
) -> Domain:
    domain_class = DOMAIN_TO_SUPERTYPE.get(domain_repr.type)
    if not domain_repr.name:
        domain_repr.name = domain_repr.type
    if domain_class:
        if domain_class == RecognitionDomain:
            domain = domain_class(
                domain_repr.type,
                substrate,
                domain_subtype=domain_repr.subtype,
                active=domain_repr.active,
                used=domain_repr.used,
            )
        elif domain_class == UnknownDomain:
            domain = UnknownDomain(domain_repr.name)
        else:
            domain = domain_class(
                domain_repr.type,
                domain_subtype=domain_repr.subtype,
                active=domain_repr.active,
                used=domain_repr.used,
            )
    elif strict:
        raise ValueError(f"Unrecognised domain type: {domain_repr.type}")
    else:
        domain = UnknownDomain(domain_repr.name)

    return domain


def build_cluster(
    cluster_repr: ClusterRepresentation, strict: bool = True
) -> ModularCluster:

    genes = set()

    modules = []
    previous_domain = None
    for i, module_repr in enumerate(cluster_repr.modules):
        if i == 0:
            starter = True
        else:
            starter = False

        if i == len(cluster_repr.modules) - 1:
            terminator = True
        else:
            terminator = False

        domains = []
        for domain_repr in module_repr.domains:
            domain = make_domain(domain_repr, module_repr.substrate, strict=strict)

            if domain_repr.gene_name is not None:
                if previous_domain:
                    previous_gene = previous_domain.gene
                    if (
                        previous_gene != domain_repr.gene_name
                        and domain_repr.gene_name in genes
                    ):
                        raise ValueError(
                            f"Gene name '{previous_gene}' already assigned to upstream domain(s)."
                        )

                domain.set_gene(domain_repr.gene_name)
                genes.add(domain_repr.gene_name)
            domains.append(domain)
            previous_domain = domain
        module_type = ModuleType.from_string(module_repr.type)
        if module_type.name == "PKS":
            if module_repr.subtype is not None:
                module_subtype = PKSModuleSubtype.from_string(module_repr.subtype)

                if module_subtype.name == "PKS_CIS":
                    module = LinearPKSModule(
                        i, domains, starter=starter, terminator=terminator
                    )
                elif module_subtype.name == "PKS_TRANS":
                    module = TransATPKSModule(
                        i, domains, starter=starter, terminator=terminator
                    )
                elif module_subtype.name == "PKS_ITER":
                    module = IterativePKSModule(
                        i,
                        domains,
                        starter=starter,
                        terminator=terminator,
                        iterations=module_repr.iterations,
                    )
                else:
                    raise ValueError(
                        f"Unrecognised PKS module subtype: {module_subtype}."
                    )
            else:
                raise ValueError("PKS module subtype must be specified.")

        elif module_type.name == "NRPS":
            if module_repr.subtype is None:
                module = NRPSModule(i, domains, starter=starter, terminator=terminator)
            else:
                raise ValueError(
                    "NRPS module subtypes are currently not supported. Please pass None."
                )
        else:
            raise ValueError(f"Unrecognised module type: {module_repr.type}")

        if (
            not any([domain.type in ["CP", "ACP", "PCP"] for domain in module.domains])
            and strict == True
        ):
            raise ValueError(
                f"Module {module} does not contain a carrier protein domain."
            )
        modules.append(module)
    cluster = ModularCluster(modules, cluster_repr.tailoring_enzymes)

    return cluster


def draw_cluster(
    cluster_repr: ClusterRepresentation, out_file=None, colour_by_module=True
) -> None:
    cluster = build_cluster(cluster_repr, strict=False)
    cluster.compute_structures(compute_cyclic_products=False)
    cluster.do_tailoring()

    if out_file:
        return cluster.draw_cluster(
            as_string=False, out_file=out_file, colour_by_module=colour_by_module
        )
    else:
        return cluster.draw_cluster()


def draw_products(cluster_repr: ClusterRepresentation, out_dir) -> None:
    cluster = build_cluster(cluster_repr)
    cluster.compute_structures(compute_cyclic_products=False)
    cluster.do_tailoring()
    cluster.cyclise_all()
    cluster.draw_all_products(out_dir)


def draw_product(cluster_repr: ClusterRepresentation, out_file) -> None:
    cluster = build_cluster(cluster_repr)
    cluster.compute_structures(compute_cyclic_products=True)
    cluster.do_tailoring()
    cluster.draw_product(as_string=False, out_file=out_file)


def draw_ripp_structure(ripp_cluster: RiPPCluster, out_folder: str) -> None:
    if not os.path.exists(out_folder):
        os.mkdir(out_folder)

    ripp_cluster.make_peptide()
    ripp_cluster.draw_product(
        as_string=False, out_file=os.path.join(out_folder, "peptide_test_ripp.svg")
    )
    ripp_cluster.draw_cluster(
        as_string=False, out_file=os.path.join(out_folder, "ripp_inline.svg")
    )

    order = []

    ripp_cluster.do_tailoring()
    if ripp_cluster.tailoring_representations:
        ripp_cluster.draw_tailoring(
            out_file=os.path.join(out_folder, "ripp_tailoring_steps.svg")
        )
        ripp_cluster.draw_product(
            as_string=False, out_file=os.path.join(out_folder, "ripp_tailoring.svg")
        )
        order.append("tailoring")
    if ripp_cluster.macrocyclisation_representations:
        ripp_cluster.do_macrocyclization()
        ripp_cluster.draw_product(
            as_string=False,
            out_file=os.path.join(out_folder, "ripp_macrocyclisation.svg"),
        )
        order.append("cyclisation")

    if ripp_cluster.cleavage_sites:
        print("here")
        ripp_cluster.do_proteolytic_cleavage()
        ripp_cluster.draw_product(
            as_string=False, out_file=os.path.join(out_folder, "cleavage_test_ripp.svg")
        )
        order.append("cleavage")

    ripp_cluster.draw_pathway(
        out_file=os.path.join(out_folder, "ripp_pathway.svg"), order=tuple(order)
    )


def draw_terpene_structure(terpene_cluster: TerpeneCluster) -> None:
    terpene_cluster.create_precursor()
    terpene_cluster.draw_product(as_string=False, out_file="precursor_test_terpene.svg")
    terpene_cluster.do_macrocyclization()
    terpene_cluster.draw_product(
        as_string=False, out_file="macroyclisation_test_terpene.svg"
    )
    terpene_cluster.do_tailoring()
    terpene_cluster.draw_product(as_string=False, out_file="tailoring_test_terpene.svg")
    print(get_tailoring_sites(terpene_cluster.chain_intermediate))


def draw_alkaloid_structure(alkaloid_cluster: AlkaloidCluster) -> None:
    alkaloid_cluster.make_scaffold()
    alkaloid_cluster.draw_product(
        as_string=False, out_file="precursor_test_alkaloid.svg"
    )
    alkaloid_cluster.do_tailoring()
    alkaloid_cluster.draw_product(
        as_string=False, out_file="tailoring_test_alkaloid.svg"
    )


def get_spaghettis(cluster_repr: ClusterRepresentation) -> list[str]:

    cluster = build_cluster(cluster_repr)
    cluster.compute_structures(compute_cyclic_products=False)
    cluster.do_tailoring()
    cluster.draw_cluster()
    spaghettis = cluster.get_spaghettis()

    return spaghettis


def get_tailoring_sites(structure, enzyme_name="all", out_file=None):
    tailoring_sites = {}
    if enzyme_name == "all":
        for enzyme_type in TailoringEnzymeType:
            tailoring_enzyme = TailoringEnzyme("gene", enzyme_type.name)
            tailoring_sites[enzyme_type.name] = tailoring_enzyme.get_possible_sites(
                structure
            )
    else:
        for enzyme_type in TailoringEnzymeType:
            if enzyme_type.name == enzyme_name:
                tailoring_enzyme = TailoringEnzyme("gene", enzyme_type.name)
                tailoring_sites[enzyme_type.name] = tailoring_enzyme.get_possible_sites(
                    structure, out_file=out_file
                )

    return tailoring_sites


def get_tailoring_sites_atom_names(structure):
    tailoring_sites = {}
    for enzyme_type in TailoringEnzymeType:
        tailoring_enzyme = TailoringEnzyme("gene", enzyme_type.name)
        tailoring_sites[enzyme_type.name] = [
            str(atom) if type(atom) != list else [str(subatom) for subatom in atom]
            for atom in tailoring_enzyme.get_possible_sites(structure)
        ]
    return tailoring_sites


def read_in_cluster(input_dir):
    cluster_representation = ClusterRepresentation.from_file(input_dir)
    print(cluster_representation)
    cluster = build_cluster(cluster_representation, strict=False)
    return cluster, cluster_representation


if __name__ == "__main__":
    # atropopeptide

    ripp_cluster = RiPPCluster(
        "best_ripp(tryptorubin)_encoding_gene",
        "mkaekslkayawyiwyaha",
        "slkayawyiwy",
        cleavage_sites=[CleavageSiteRepresentation("Y", 10, "end")],
        tailoring_representations=[
            TailoringRepresentation(
                "p450", "DOUBLE_BOND_REDUCTASE", [["C_139", "C_138"]]
            ),
            TailoringRepresentation(
                "p450",
                "OXIDATIVE_BOND_SYNTHASE",
                [["C_139", "N_134"], ["C_120", "N_102"], ["C_138", "C_107"]],
            ),
        ],
    )

    lanthipeptide_type_I_cluster_catenulipeptin = RiPPCluster(
        "Caci_4240",
        "MTEEMTLLDLQGMEQTETDSWGGSGHGGGGDSGLSVTGCNGHSGISLLCDL",
        "GHGGGGDSGLSVTGCNGHSGISLLCDL",
        tailoring_representations=[
            TailoringRepresentation(
                "Caci_4239", "THREONINE_SERINE_DEHYDRATASE", [["O_64"], ["O_144"]]
            ),
            TailoringRepresentation(
                "Caci_4239",
                "LANTHIPEPTIDE_CYCLASE",
                [
                    ["C_43", "C_63"],
                    ["S_92", "C_63"],
                    ["C_122", "C_143"],
                    ["S_169", "C_143"],
                ],
            ),
        ],
    )
    lanthipeptide_type_III_cluster_catenulipeptin = RiPPCluster(
        "Caci_4240",
        "MTEEMTLLDLQGMEQTETDSWGGSGHGGGGDSGLSVTGCNGHSGISLLCDL",
        "GHGGGGDSGLSVTGCNGHSGISLLCDL",
        tailoring_representations=[
            TailoringRepresentation(
                "Caci_4239",
                "LANTHIONINE_SYNTHETASE",
                [
                    ["C_43", "C_63"],
                    ["S_92", "C_63"],
                    ["C_122", "C_143"],
                    ["S_169", "C_143"],
                ],
            )
        ],
    )
    proteusins_cluster_polytheonamide_a = RiPPCluster(
        "poyA",
        "MADSDNTPTSRKDFETAIIAKAWKDPEYLRRLRSNPREVLQEELEALHPGAQLPDDLGISIHEEDENHVHLVMPRHPQNVSDQTLTDDDLDQAAGGTGIGVVVAVVAGAVANTGAGVNQVAGGNINVVGNINVNANVSVNMNQTT",
        "TGIGVVVAVVAGAVANTGAGVNQVAGGNINVVGNINVNANVSVNMNQTT",
        tailoring_representations=[
            TailoringRepresentation(
                "rSAM epimerase",
                "AMINO_ACID_EPIMERASE",
                [
                    ["C_269"],
                    ["C_144"],
                    ["C_284"],
                    ["C_163"],
                    ["C_36"],
                    ["C_301"],
                    ["C_52"],
                    ["C_90"],
                    ["C_185"],
                    ["C_316"],
                    ["C_66"],
                    ["C_204"],
                    ["C_334"],
                    ["C_84"],
                    ["C_221"],
                    ["C_353"],
                    ["C_104"],
                    ["C_233"],
                    ["C_252"],
                ],
            )
        ],
    )
    sliceotide_cluster = RiPPCluster(
        "plpA",
        "NVSVNMNQTTRVSVNMNQTNVSVNVSVNMNQTNVSVNVSVNMNQTNVSVNVSVN",
        "NVSVNMNQTTR",
        tailoring_representations=[
            TailoringRepresentation("rSAM epimerase", "SPLICEASE", [["C_30", "C_22"]]),
            TailoringRepresentation("arginase", "ARGINASE", [["N_93"]]),
        ],
    )
    thiopeptide_cluster_thiomuracin = RiPPCluster(
        "tpdA",
        "MDLSDLPMDVFELADDGVAVESLTAGHGMTEVGASCNCFCYICCSCSSA",
        "SCNCFCYICCSCSS",
        tailoring_representations=[
            TailoringRepresentation(
                "tpdD",
                "THREONINE_SERINE_DEHYDRATASE",
                [["O_104"], ["O_90"], ["O_4"], ["O_111"]],
            ),
            TailoringRepresentation(
                "tpdE",
                "CYCLODEHYDRASE",
                [["S_97"], ["S_11"], ["S_76"], ["S_46"], ["S_83"], ["S_27"]],
            ),
            TailoringRepresentation("tpdF", "THIOPEPTIDE_CYCLASE", [["C_3", "C_89"]]),
            TailoringRepresentation("tpdF", "METHYLTRANSFERASE", [["C_26"]]),
            TailoringRepresentation("tpdF", "HYDROXYLASE", [["C_33"]]),
            TailoringRepresentation(
                "tpdF", "DEHYDROGENASE", [["C_67", "C_68"], ["C_3", "C_89"]]
            ),
            TailoringRepresentation("tpdF", "EPOXIDASE", [["C_67", "C_68"]]),
            TailoringRepresentation("tpdF", "MONOAMINE_OXIDASE", [["N_0"]]),
            TailoringRepresentation("tpdF", "DEHYDRATASE", [["C_1", "C_84"]]),
        ],
    )
    cyanobactin_cluster_trunkamide = RiPPCluster(
        "truE",
        "MNKKNILPQLGQPVIRLTAGQLSSQLAELSEEALGGVDASTSIAPFCSYDGVDASSYDGVDASSYDD",
        "TSIAPFC",
        macrocyclisations=[MacrocyclizationRepresentation("N_0", "O_59")],
        tailoring_representations=[
            TailoringRepresentation("truD", "CYCLODEHYDRASE", [["S_56"]]),
            TailoringRepresentation(
                "truF", "PRENYLTRANSFERASE", [["O_13"], ["O_5"]], "3_METHYL_1_BUTENYL"
            ),
        ],
    )

    lasso_peptide_cluster = RiPPCluster(
        "A1S42_RS12075",
        "MKYCKPTFESIATFKKDTKGLWTGKFRDIFGGRAIVRIRIEF",
        "MKYCKPTFESIATFKKDTKGLWTGKFRDIFGGRAIVRIRIEF",
        tailoring_representations=[
            TailoringRepresentation("lasB", "PROTEASE", [["N_180", "C_178"]]),
            TailoringRepresentation("lasC", "MACROLACTAM_SYNTHETASE", [["O_261"]]),
        ],
    )
    sacti_peptide_cluster_thurincin = RiPPCluster(
        "thnA",
        "METPVVQPRDWTCWSCLVCAACSVELLNLVTAATGASTAS",
        "DWTCWSCLVCAACSVELLNLVTAATGASTAS",
        tailoring_representations=[
            TailoringRepresentation(
                "thnB",
                "OXIDATIVE_BOND_SYNTHASE",
                [
                    ["S_109", "C_222"],
                    ["S_66", "C_203"],
                    ["S_90", "C_182"],
                    ["S_37", "C_156"],
                ],
            )
        ],
    )

    terpene_cluster = TerpeneCluster(
        "limonene_synthase",
        "GERANYL_PYROPHOSPHATE",
        macrocyclisations=[MacrocyclizationRepresentation("C_13", "C_8")],
        cyclase_type="Class_1",
        tailoring_representations=[
            TailoringRepresentation(
                "pseudo_isomerase",
                "DOUBLE_BOND_SHIFT",
                [["C_13", "C_14", "C_14", "C_15"]],
            )
        ],
    )

    alkaloid_cluster = AlkaloidCluster(
        "phenylalanine",
        tailoring_representations=[
            TailoringRepresentation("pseudo_decarboxylase", "DECARBOXYLASE", [["C_9"]]),
            TailoringRepresentation(
                "pseudo_hydroxylase", "PRENYLTRANSFERASE", [["C_7"]], "DIMETHYLALLYL"
            ),
            TailoringRepresentation(
                "pseudo_decarboxylase", "HALOGENASE", [["C_10"]], "Cl"
            ),
            TailoringRepresentation("pseudo_hydroxylase", "HYDROXYLASE", [["C_6"]]),
            TailoringRepresentation(
                "methyltransferase", "METHYLTRANSFERASE", [["N_12"], ["C_7"], ["O_25"]]
            ),
        ],
    )

    cluster_repr = ClusterRepresentation(
        [
            ModuleRepresentation(
                "PKS",
                "PKS_CIS",
                "ACETYL_COA",
                [
                    DomainRepresentation("Gene 1", "AT", None, None, True, True),
                    DomainRepresentation("Gene 1", "ACP", None, None, True, True),
                ],
            ),
            ModuleRepresentation(
                "PKS",
                "PKS_ITER",
                "METHYLMALONYL_COA",
                [
                    DomainRepresentation("Gene 1", "KS", None, None, True, True),
                    DomainRepresentation("Gene 1", "AT", None, None, True, True),
                    DomainRepresentation("Gene 1", "AT", None, None, True, False),
                    DomainRepresentation("Gene 1", "DH", None, None, True, True),
                    DomainRepresentation("Gene 1", "ER", None, None, True, True),
                    DomainRepresentation("Gene 1", "ACP", None, None, True, True),
                ],
                5,
            ),
            ModuleRepresentation(
                "PKS",
                "PKS_CIS",
                "METHYLMALONYL_COA",
                [
                    DomainRepresentation("Gene 1", "KS", None, None, True, True),
                    DomainRepresentation("Gene 1", "AT", None, None, True, True),
                    DomainRepresentation("Gene 1", "TE", None, None, True, True),
                ],
            ),
        ]
    )

    trans_at_ks_cluster_repr = ClusterRepresentation(
        [
            ModuleRepresentation(
                "PKS",
                "PKS_TRANS",
                "ACETYL_COA",
                [
                    DomainRepresentation("Gene 1", "AT", None, None, True, True),
                    DomainRepresentation("Gene 1", "ACP", None, None, True, True),
                ],
            ),
            ModuleRepresentation(
                "PKS",
                "PKS_TRANS",
                "METHYLMALONYL_COA",
                [
                    DomainRepresentation("Gene 1", "KS", "BETA_D_OH", None, True, True),
                    DomainRepresentation("Gene 1", "DH", None, None, True, True),
                    DomainRepresentation("Gene 1", "ER", None, None, True, True),
                    DomainRepresentation("Gene 1", "ACP", None, None, True, True),
                ],
            ),
            ModuleRepresentation(
                "PKS",
                "PKS_TRANS",
                "METHYLMALONYL_COA",
                [
                    DomainRepresentation("Gene 1", "KS", "BETA_D_OH", None, True, True),
                    DomainRepresentation("Gene 1", "DH", None, None, True, True),
                    DomainRepresentation("Gene 1", "ER", None, None, True, True),
                    DomainRepresentation("Gene 1", "ACP", None, None, True, True),
                ],
            ),
            ModuleRepresentation(
                "PKS",
                "PKS_TRANS",
                "METHYLMALONYL_COA",
                [
                    DomainRepresentation("Gene 1", "KS", "BETA_D_OH", None, True, True),
                    DomainRepresentation("Gene 1", "DH", None, None, True, True),
                    DomainRepresentation("Gene 1", "ER", None, None, True, True),
                    DomainRepresentation("Gene 1", "ACP", None, None, True, True),
                ],
            ),
            ModuleRepresentation(
                "PKS",
                "PKS_TRANS",
                "METHYLMALONYL_COA",
                [
                    DomainRepresentation("Gene 1", "KS", "PYR", None, True, True),
                    DomainRepresentation("Gene 1", "AT", None, None, True, True),
                    DomainRepresentation("Gene 1", "AT", None, None, True, False),
                    DomainRepresentation("Gene 1", "DH", None, None, True, True),
                    DomainRepresentation("Gene 1", "ER", None, None, True, True),
                    DomainRepresentation("Gene 1", "ACP", None, None, True, True),
                ],
            ),
            ModuleRepresentation(
                "PKS",
                "PKS_CIS",
                "METHYLMALONYL_COA",
                [
                    DomainRepresentation("Gene 1", "KS", None, None, True, True),
                    DomainRepresentation("Gene 1", "AT", None, None, True, True),
                    DomainRepresentation("Gene 1", "TE", None, None, True, True),
                ],
            ),
        ]
    )
    # draw_cluster(trans_at_ks_cluster_repr, outfile = "iterative_pks.svg")
    # draw_ripp_structure(ripp_cluster)

    cluster_nrps = ClusterRepresentation(
        [
            ModuleRepresentation(
                "NRPS",
                None,
                "valine",
                [
                    DomainRepresentation("Gene 1", "A", None, None, True, True),
                    DomainRepresentation("Gene 1", "PCP", None, None, True, True),
                ],
            ),
            ModuleRepresentation(
                "NRPS",
                None,
                "proline",
                [
                    DomainRepresentation("Gene 1", "C", None, None, True, True),
                    DomainRepresentation("Gene 1", "A", None, None, True, True),
                    DomainRepresentation("Gene 1", "PCP", None, None, True, True),
                ],
            ),
            ModuleRepresentation(
                "NRPS",
                None,
                "threonine",
                [
                    DomainRepresentation("Gene 1", "C", None, None, True, True),
                    DomainRepresentation("Gene 1", "CYC", None, None, True, True),
                    DomainRepresentation("Gene 1", "A", None, None, True, True),
                    DomainRepresentation("Gene 1", "OX", None, None, True, True),
                    DomainRepresentation("Gene 1", "PCP", None, None, True, True),
                ],
            ),
            ModuleRepresentation(
                "NRPS",
                None,
                "proline",
                [
                    DomainRepresentation("Gene 1", "C", None, None, True, True),
                    DomainRepresentation("Gene 1", "A", None, None, True, True),
                    DomainRepresentation("Gene 1", "PCP", None, None, True, True),
                    DomainRepresentation("Gene 1", "TE", None, None, True, True),
                ],
            ),
            ModuleRepresentation(
                "NRPS",
                None,
                "valine",
                [
                    DomainRepresentation("Gene 1", "C", None, None, True, True),
                    DomainRepresentation("Gene 1", "A", None, None, True, True),
                    DomainRepresentation("Gene 1", "PCP", None, None, True, True),
                ],
            ),
            ModuleRepresentation(
                "NRPS",
                None,
                "valine",
                [
                    DomainRepresentation("Gene 1", "C", None, None, True, True),
                    DomainRepresentation("Gene 1", "A", None, None, True, True),
                    DomainRepresentation("Gene 1", "PCP", None, None, True, True),
                ],
            ),
        ]
    )

    nrps_cluster = ClusterRepresentation(
        [
            ModuleRepresentation(
                "NRPS",
                None,
                "serine",
                [
                    DomainRepresentation("Gene 1", "C", None, None, True, True),
                    DomainRepresentation("Gene 1", "A", None, None, True, True),
                    DomainRepresentation("Gene 1", "PCP", None, None, True, True),
                ],
            ),
            ModuleRepresentation(
                "NRPS",
                None,
                "threonine",
                [
                    DomainRepresentation("Gene 1", "C", None, None, True, True),
                    DomainRepresentation("Gene 1", "A", None, None, True, True),
                    DomainRepresentation("Gene 1", "PCP", None, None, True, True),
                ],
            ),
            ModuleRepresentation(
                "NRPS",
                None,
                "tryptophan",
                [
                    DomainRepresentation("Gene 1", "C", None, None, True, True),
                    DomainRepresentation("Gene 1", "A", None, None, True, True),
                    DomainRepresentation("Gene 1", "PCP", None, None, True, True),
                ],
            ),
            ModuleRepresentation(
                "NRPS",
                None,
                "aspartic acid",
                [
                    DomainRepresentation("Gene 1", "C", None, None, True, True),
                    DomainRepresentation("Gene 1", "A", None, None, True, True),
                    DomainRepresentation("Gene 1", "PCP", None, None, True, True),
                ],
            ),
            ModuleRepresentation(
                "NRPS",
                None,
                "aspartic acid",
                [
                    DomainRepresentation("Gene 1", "C", None, None, True, True),
                    DomainRepresentation("Gene 1", "A", None, None, True, True),
                    DomainRepresentation("Gene 1", "PCP", None, None, True, True),
                ],
            ),
            ModuleRepresentation(
                "NRPS",
                None,
                "4-hydroxyphenylglycine",
                [
                    DomainRepresentation("Gene 1", "C", None, None, True, True),
                    DomainRepresentation("Gene 1", "A", None, None, True, True),
                    DomainRepresentation("Gene 1", "PCP", None, None, True, True),
                ],
            ),
            ModuleRepresentation(
                "NRPS",
                None,
                "aspartic acid",
                [
                    DomainRepresentation("Gene 1", "C", None, None, True, True),
                    DomainRepresentation("Gene 1", "A", None, None, True, True),
                    DomainRepresentation("Gene 1", "PCP", None, None, True, True),
                ],
            ),
            ModuleRepresentation(
                "NRPS",
                None,
                "alanine",
                [
                    DomainRepresentation("Gene 1", "C", None, None, True, True),
                    DomainRepresentation("Gene 1", "A", None, None, True, True),
                    DomainRepresentation("Gene 1", "PCP", None, None, True, True),
                ],
            ),
            ModuleRepresentation(
                "NRPS",
                None,
                "asparagine",
                [
                    DomainRepresentation("Gene 1", "C", None, None, True, True),
                    DomainRepresentation("Gene 1", "A", None, None, True, True),
                    DomainRepresentation("Gene 1", "PCP", None, None, True, True),
                ],
            ),
            ModuleRepresentation(
                "NRPS",
                None,
                "glutamic acid",
                [
                    DomainRepresentation("Gene 1", "C", None, None, True, True),
                    DomainRepresentation("Gene 1", "A", None, None, True, True),
                    DomainRepresentation("Gene 1", "PCP", None, None, True, True),
                ],
            ),
            ModuleRepresentation(
                "NRPS",
                None,
                "tryptophan",
                [
                    DomainRepresentation("Gene 1", "C", None, None, True, True),
                    DomainRepresentation("Gene 1", "A", None, None, True, True),
                    DomainRepresentation("Gene 1", "PCP", None, None, True, True),
                ],
            ),
        ]
    )
    draw_cluster(nrps_cluster, out_file="antismash_test.svg")
    # draw_ripp_structure(ripp_cluster)
    # ripp_cluster.draw_precursor(as_string= False, out_file= "bubbles.svg")

    # cyanobactin_cluster_trunkamide.make_peptide()
    # # print(get_tailoring_sites_atom_names(
    # #     cyanobactin_cluster_trunkamide.chain_intermediate))
    # cyanobactin_cluster_trunkamide.draw_product(
    #     as_string=False, out_file="peptide_test_cyanobactin_peptide.svg")
    # cyanobactin_cluster_trunkamide.do_tailoring()
    # cyanobactin_cluster_trunkamide.draw_product(
    #     as_string=False, out_file="tailored_test_cyanobactin_peptide.svg", draw_straightened=False)
    # cyanobactin_cluster_trunkamide.do_macrocyclization()
    # cyanobactin_cluster_trunkamide.draw_product(
    #     as_string=False, out_file="final_peptide_test_cyanobactin_peptide.svg", draw_straightened=False)
    # lanthipeptide_type_I_cluster_catenulipeptin.make_peptide()
    # # print(get_tailoring_sites_atom_names(
    # #     lanthipeptide_type_III_cluster_catenulipeptin.chain_intermediate))
    # lanthipeptide_type_I_cluster_catenulipeptin.draw_product(
    #     as_string=False, out_file="peptide_test_lanthipeptide_peptide.svg")
    # lanthipeptide_type_I_cluster_catenulipeptin.do_tailoring()
    # lanthipeptide_type_I_cluster_catenulipeptin.draw_product(
    #     as_string=False, out_file="tailored_test_lanthipeptide_peptide.svg", draw_straightened=False)

    # thiopeptide_cluster_thiomuracin.make_peptide()
    #
    # thiopeptide_cluster_thiomuracin.draw_product(
    #     as_string=False, out_file="peptide_test_thiopeptide_peptide.svg")
    # thiopeptide_cluster_thiomuracin.do_tailoring()
    # thiopeptide_cluster_thiomuracin.draw_product(
    #     as_string=False, out_file="tailored_test_thiopeptide_peptide.svg")
    # get_tailoring_sites(thiopeptide_cluster_thiomuracin.chain_intermediate,
    #                     enzyme_name="DEHYDROGENASE", out_file="methylation_options_thiopeptide.svg")
    # # print(get_tailoring_sites_atom_names(
    # #     thiopeptide_cluster_thiomuracin.chain_intermediate))

    # get_tailoring_sites(thiopeptide_cluster_thiomuracin.chain_intermediate, enzyme_name="CYCLODEHYDRASE",
    #                     out_file="cyclodehydration.svg")
    # thiopeptide_cluster_thiomuracin.draw_precursor_with_modified_product(as_string=False, out_file="bubbles.svg")

    # print(get_tailoring_sites_atom_names(
    #     thiopeptide_cluster_thiomuracin.chain_intermediate))
    # sacti_peptide_cluster_thurincin.make_peptide()
    # # # print(get_tailoring_sites_atom_names(
    # # # sancti_peptide_cluster_thurincin.chain_intermediate))
    # # # sancti_peptide_cluster_thurincin.draw_product(
    # # #     as_string=False, out_file="peptide_test_sancti_peptide.svg")
    # sacti_peptide_cluster_thurincin.do_tailoring()
    # drawer = Drawer(sacti_peptide_cluster_thurincin.chain_intermediate)
    # drawer.show_molecule()

    # # sancti_peptide_cluster_thurincin.draw_product(
    # #     as_string=False, out_file="tailored_test_sancti_peptide.svg", draw_straightened=False)
    # # print(get_tailoring_sites_atom_names(
    # #     sancti_peptide_cluster_thurincin.chain_intermediate))
    # # proteusins_cluster_polytheonamide_a.make_peptide()
    # # proteusins_cluster_polytheonamide_a.draw_product(
    # #     as_string=False, out_file="peptide_test_proteusin_peptide.svg")
    # # proteusins_cluster_polytheonamide_a.do_tailoring()
    # # proteusins_cluster_polytheonamide_a.draw_product(
    # #     as_string=False, out_file="tailored_test_proteusin_peptide.svg", draw_straightened=False)
    # print(get_tailoring_sites_atom_names(
    #     thiopeptide_cluster_thiomuracin.chain_intermediate))

    # sancti_peptide_cluster_thurincin.make_peptide()
    # sancti_peptide_cluster_thurincin.draw_precursor_with_modified_product(as_string=False, out_file="bubbles.svg")
    # print(get_tailoring_sites_atom_names(
    # sancti_peptide_cluster_thurincin.chain_intermediate))
    # sancti_peptide_cluster_thurincin.draw_product(
    #     as_string=False, out_file="peptide_test_sancti_peptide.svg")
    # sancti_peptide_cluster_thurincin.do_tailoring()
    # sancti_peptide_cluster_thurincin.draw_product(
    #     as_string=False, out_file="tailored_test_sancti_peptide.svg", draw_straightened=False)
    # print(get_tailoring_sites_atom_names(
    #     sancti_peptide_cluster_thurincin.chain_intermediate))
    # proteusins_cluster_polytheonamide_a.make_peptide()
    # proteusins_cluster_polytheonamide_a.draw_product(
    #     as_string=False, out_file="peptide_test_proteusin_peptide.svg")
    # proteusins_cluster_polytheonamide_a.do_tailoring()
    # proteusins_cluster_polytheonamide_a.draw_product(
    #     as_string=False, out_file="tailored_test_proteusin_peptide.svg", draw_straightened=False)
    # sliceotide_cluster.make_peptide()
    # sliceotide_cluster.draw_product(as_string=False, out_file = "ripp-test.svg")
    # # print(get_tailoring_sites_atom_names(
    # #     sliceotide_cluster.chain_intermediate))
    # # sliceotide_cluster.draw_product(
    # #     as_string=False, out_file="peptide_test_sliceotide_cluster.svg")
    # sliceotide_cluster.draw_precursor_with_modified_product(as_string=False, out_file="bubbles.svg")
    # sliceotide_cluster.do_tailoring()
    # sliceotide_cluster.draw_product(
    #     as_string=False, out_file="tailored_test_sliceotide_cluster.svg", draw_straightened=False)
    # lasso_peptide_cluster.do_tailoring()
    ## draw_terpene_structure(terpene_cluster)
    # draw_alkaloid_structure(alkaloid_cluster)
    # cyanobactin_cluster_trunkamide.make_peptide()
    # #cyanobactin_cluster_trunkamide.draw_precursor_with_modified_product(as_string=False, out_file= "trunkamide_peptide.svg")
    # cyanobactin_cluster_trunkamide.do_tailoring()
    # cyanobactin_cluster_trunkamide.do_macrocyclization()
    # drawer = Drawer(cyanobactin_cluster_trunkamide.chain_intermediate)
    # # drawer.show_molecule()

    # cyanobactin_cluster_trunkamide.draw_product(as_string=False, out_file= "trunkamide.svg")
    # draw_terpene_structure(terpene_cluster)
    # draw_alkaloid_structure(alkaloid_cluster)

    # draw_cluster(trans_at_ks_cluster_repr, outfile="iterative_pks.svg")
