from typing import List, Union

from pikachu.drawing.drawing import *
from raichu.cluster import Cluster
from raichu.ripp import RiPP_Cluster
from raichu.terpene import Terpene_Cluster
from raichu.domain.domain import TailoringDomain, CarrierDomain, SynthesisDomain, RecognitionDomain, \
    TerminationDomain, UnknownDomain, Domain
from raichu.module import PKSModuleSubtype, NRPSModule, LinearPKSModule, IterativePKSModule, TransATPKSModule,\
    ModuleType
from raichu.domain.domain_types import TailoringDomainType, TerminationDomainType, CarrierDomainType, \
    SynthesisDomainType, RecognitionDomainType
from dataclasses import dataclass
from raichu.alkaloid import Alkaloid_Cluster
from raichu.tailoring_enzymes import TailoringEnzyme, TailoringEnzymeType

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

@dataclass
class MacrocyclizationRepresentation:
    atom1: str
    atom2: str


@dataclass
class CleavageSiteRepresentation:
    position_amino_acid: str
    position_index: int
    structure_to_keep: str


@dataclass
class TailoringRepresentation:
    gene_name: str
    type: str
    # Some tailoring reactions involve more than one atom
    modification_sites: List[List[str]]
    substrate: Union[str, None] = None


@dataclass
class DomainRepresentation:
    gene_name: Union[str, None]
    type: str
    subtype: Union[str, None]
    name: Union[str, None]
    active: bool
    used: bool


@dataclass
class ModuleRepresentation:
    type: str
    subtype: Union[str, None]
    substrate: str
    domains: List[DomainRepresentation]
    iterations: int = 1


@dataclass
class ClusterRepresentation:
    modules: List[ModuleRepresentation]
    tailoring_enzymes: Union[List[TailoringRepresentation], None] = None


def make_domain(domain_repr: DomainRepresentation, substrate: str, strict: bool = True) -> Domain:
    domain_class = DOMAIN_TO_SUPERTYPE.get(domain_repr.type)
    if not domain_repr.name:
        domain_repr.name = domain_repr.type
    if domain_class:
        if domain_class == RecognitionDomain:
            domain = domain_class(domain_repr.type, substrate, domain_subtype=domain_repr.subtype,
                                  active=domain_repr.active,
                                  used=domain_repr.used)
        elif domain_class == UnknownDomain:
            domain = UnknownDomain(domain_repr.name)
        else:
            domain = domain_class(domain_repr.type, domain_subtype=domain_repr.subtype, active=domain_repr.active,
                                  used=domain_repr.used)
    elif strict:
        raise ValueError(f"Unrecognised domain type: {domain_repr.type}")
    else:
        domain = UnknownDomain(domain_repr.name)

    return domain


def build_cluster(cluster_repr: ClusterRepresentation, strict: bool = True) -> Cluster:

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
            domain = make_domain(domain_repr, module_repr.substrate,
                                 strict=strict)

            if domain_repr.gene_name is not None:
                if previous_domain:
                    previous_gene = previous_domain.gene
                    if previous_gene != domain_repr.gene_name and domain_repr.gene_name in genes:
                        raise ValueError(
                            f"Gene name '{previous_gene}' already assigned to upstream domain(s).")

                domain.set_gene(domain_repr.gene_name)
                genes.add(domain_repr.gene_name)
            domains.append(domain)
            previous_domain = domain
        module_type = ModuleType.from_string(module_repr.type)
        if module_type.name == 'PKS':
            if module_repr.subtype is not None:
                module_subtype = PKSModuleSubtype.from_string(
                    module_repr.subtype)

                if module_subtype.name == 'PKS_CIS':
                    module = LinearPKSModule(
                        i, domains, starter=starter, terminator=terminator)
                elif module_subtype.name == 'PKS_TRANS':
                    module = TransATPKSModule(
                        i, domains, starter=starter, terminator=terminator)
                elif module_subtype.name == 'PKS_ITER':
                    module = IterativePKSModule(
                        i, domains, starter=starter, terminator=terminator, iterations=module_repr.iterations)
                else:
                    raise ValueError(
                        f"Unrecognised PKS module subtype: {module_subtype}.")
            else:
                raise ValueError("PKS module subtype must be specified.")

        elif module_type.name == 'NRPS':
            if module_repr.subtype is None:
                module = NRPSModule(
                    i, domains, starter=starter, terminator=terminator)
            else:
                raise ValueError(
                    "NRPS module subtypes are currently not supported. Please pass None.")
        else:
            raise ValueError(f"Unrecognised module type: {module_repr.type}")

        modules.append(module)
    cluster = Cluster(modules, cluster_repr.tailoring_enzymes)

    return cluster


def draw_cluster(cluster_repr: ClusterRepresentation, outfile=None) -> None:
    cluster = build_cluster(cluster_repr)
    cluster.compute_structures(compute_cyclic_products=False)
    cluster.do_tailoring()
    cluster.draw_product(as_string=False, out_file="tailoring_test.svg")
    if outfile:
        return cluster.draw_cluster(as_string=False, out_file=outfile)
    else:
        cluster.draw_cluster()


def draw_ripp_structure(ripp_cluster: RiPP_Cluster) -> None:
    ripp_cluster.make_peptide()
    ripp_cluster.draw_product(
        as_string=False, out_file="peptide_test_ripp.svg")
    ripp_cluster.do_tailoring()
    ripp_cluster.draw_product(
        as_string=False, out_file="tailoring_test_ripp.svg")
    ripp_cluster.do_macrocyclization()
    ripp_cluster.draw_product(
        as_string=False, out_file="macrocyclisation_test_ripp.svg")
    ripp_cluster.do_proteolytic_claevage()
    ripp_cluster.draw_product(
        as_string=False, out_file="cleavage_test_ripp.svg")


def draw_terpene_structure(terpene_cluster: Terpene_Cluster) -> None:
    terpene_cluster.create_precursor()
    terpene_cluster.draw_product(
        as_string=False, out_file="precursor_test_terpene.svg")
    terpene_cluster.do_macrocyclization()
    terpene_cluster.draw_product(
        as_string=False, out_file="macroyclisation_test_terpene.svg")
    terpene_cluster.do_tailoring()
    terpene_cluster.draw_product(
        as_string=False, out_file="tailoring_test_terpene.svg")
    print(get_tailoring_sites(terpene_cluster.chain_intermediate))


def draw_alkaloid_structure(alkaloid_cluster: Alkaloid_Cluster) -> None:
    alkaloid_cluster.make_scaffold()
    alkaloid_cluster.draw_product(
        as_string=False, out_file="precursor_test_alkaloid.svg")
    alkaloid_cluster.do_tailoring()
    alkaloid_cluster.draw_product(
        as_string=False, out_file="tailoring_test_alkaloid.svg")



def get_spaghettis(cluster_repr: ClusterRepresentation) -> List[str]:

    cluster = build_cluster(cluster_repr)
    cluster.compute_structures(compute_cyclic_products=False)
    cluster.do_tailoring()
    cluster.draw_cluster()
    spaghettis = cluster.draw_spaghettis()

    return spaghettis


def get_tailoring_sites(structure):
    tailoring_sites = {}
    for enzyme_type in TailoringEnzymeType:
        tailoring_enzyme = TailoringEnzyme("gene", enzyme_type.name)
        tailoring_sites[enzyme_type.name] = tailoring_enzyme.get_possible_sites(
            structure)
    return tailoring_sites


def get_tailoring_sites_atom_names(structure):
    tailoring_sites = {}
    for enzyme_type in TailoringEnzymeType:
        tailoring_enzyme = TailoringEnzyme("gene", enzyme_type.name)
        tailoring_sites[enzyme_type.name] = [str(atom) if type(atom) != list else [str(subatom) for subatom in atom] for atom in tailoring_enzyme.get_possible_sites(
            structure)]
    return tailoring_sites


if __name__ == "__main__":
    # atropopeptide
    ripp_cluster = RiPP_Cluster("best_ripp(tryptorubin)_encoding_gene", "mkaekslkayawyiwy", "mkaekslkayawyiwy",
                                cleavage_sites=[CleavageSiteRepresentation("Y", 10, "follower")],
                                tailoring_enzymes_representation=[TailoringRepresentation("p450", "REDUCTASE_DOUBLE_BOND_REDUCTION", [["C_139", "C_138"]]), TailoringRepresentation("p450", "P450_OXIDATIVE_BOND_FORMATION", [["C_139", "N_134"], ["C_120", "N_102"], ["C_138", "C_107"]])])
    lanthipeptide_type_I_cluster_catenulipeptin = RiPP_Cluster("Caci_4240", "MTEEMTLLDLQGMEQTETDSWGGSGHGGGGDSGLSVTGCNGHSGISLLCDL", "GHGGGGDSGLSVTGCNGHSGISLLCDL",
                                                               tailoring_enzymes_representation=[TailoringRepresentation(
                                                                   "Caci_4239", "THREONINE_SERINE_DEHYDRATASE", [['O_64'], ['O_144']]), TailoringRepresentation(
                                                                   "Caci_4239", "LANTHIPEPTIDE_CYCLASE", [["C_43", "C_63"], ["S_92", "C_63"], ["C_122", "C_143"], ["S_169", "C_143"]])])
    lanthipeptide_type_III_cluster_catenulipeptin = RiPP_Cluster("Caci_4240", "MTEEMTLLDLQGMEQTETDSWGGSGHGGGGDSGLSVTGCNGHSGISLLCDL", "GHGGGGDSGLSVTGCNGHSGISLLCDL",
                                                                tailoring_enzymes_representation=[TailoringRepresentation(
                                                                    "Caci_4239", "LANTHIONINE_SYNTHETASE", [["C_43", "C_63"], ["S_92", "C_63"], ["C_122", "C_143"], ["S_169", "C_143"]])])
    proteusins_cluster_polytheonamide_a = RiPP_Cluster(
        "poyA", "MADSDNTPTSRKDFETAIIAKAWKDPEYLRRLRSNPREVLQEELEALHPGAQLPDDLGISIHEEDENHVHLVMPRHPQNVSDQTLTDDDLDQAAGGTGIGVVVAVVAGAVANTGAGVNQVAGGNINVVGNINVNANVSVNMNQTT",
        "TGIGVVVAVVAGAVANTGAGVNQVAGGNINVVGNINVNANVSVNMNQTT", tailoring_enzymes_representation=[TailoringRepresentation("rSAM epimerase", "AMINO_ACID_EPIMERASE", [['C_269'], ['C_144'], ['C_284'], ['C_163'], ['C_36'], ['C_301'], ['C_52'], ['C_90'],  ['C_185'],  ['C_316'], ['C_66'], ['C_204'], ['C_334'], ['C_84'],  ['C_221'], ['C_353'], ['C_104'], ['C_233'], ['C_252']])])
    sliceotide_cluster = RiPP_Cluster("plpA", "HGGVGFNVSVNMNQTTRFGVGHHS", "NVSVNMNQTTR", tailoring_enzymes_representation=[TailoringRepresentation(
        "rSAM epimerase", "SPLICEASE", [["C_30", "C_22"]]), TailoringRepresentation("arginase", "ARGINASE", [['N_93']])])
    thiopeptide_cluster_thiomuracin = RiPP_Cluster(
        "tpdA", "MDLSDLPMDVFELADDGVAVESLTAGHGMTEVGASCNCFCYICCSCSSA", "SCNCFCYICCSCSS",
        tailoring_enzymes_representation=[TailoringRepresentation(
            "tpdD", "THREONINE_SERINE_DEHYDRATASE", [['O_104'], ['O_90'], ['O_4'], ['O_111']]), TailoringRepresentation(
            "tpdE", "CYCLODEHYDRATION", [['S_97'], ['S_11'], ['S_76'], ['S_46'], ['S_83'], ['S_27']]),
            TailoringRepresentation("tpdF", "THIOPEPTIDE_CYCLASE",[["C_3", "C_89"]],)
            ]
        )
    cyanobactin_cluster_trunkamide = RiPP_Cluster("truE", "MNKKNILPQLGQPVIRLTAGQLSSQLAELSEEALGGVDASTSIAPFCSYDGVDASTSIAPFCSYDGVDASTSIAPFCSYDD", "TSIAPFC",
                                                  macrocyclisations=[MacrocyclizationRepresentation("N_0","O_59")],
                                                  tailoring_enzymes_representation=[TailoringRepresentation("truD", "CYCLODEHYDRATION", [["S_56"]]),
                                                                                    TailoringRepresentation("truF", "PRENYLTRANSFERASE", [['O_13'], ['O_5']], "3_METHYL_1_BUTENYL")]
                                )

    lasso_peptide_cluster = RiPP_Cluster("A1S42_RS12075", "MKYCKPTFESIATFKKDTKGLWTGKFRDIFGGRAIVRIRIEF", "MKYCKPTFESIATFKKDTKGLWTGKFRDIFGGRAIVRIRIEF",
                                         tailoring_enzymes_representation=[TailoringRepresentation("lasB", "PROTEASE", [["N_180", "C_178"]]), TailoringRepresentation("lasC", "MACROLACTAM_SYNTHETASE", [["O_261"]])]
                                )
    sancti_peptide_cluster_thurincin = RiPP_Cluster("thnA", "METPVVQPRDWTCWSCLVCAACSVELLNLVTAATGASTAS", "DWTCWSCLVCAACSVELLNLVTAATGASTAS",
                                                    tailoring_enzymes_representation=[TailoringRepresentation(
                                                        "thnB", "OXIDATIVE_BOND_FORMATION", [['S_109', "C_222"], ['S_66', "C_203"], ['S_90', "C_182"], ['S_37', "C_156"]])]
                                    )

    terpene_cluster = Terpene_Cluster("limonene_synthase", "GERANYL_PYROPHOSPHATE",
                                      macrocyclisations=[MacrocyclizationRepresentation("C_13", "C_8")],
                                      terpene_cyclase_type="Class_1",
                                      tailoring_enzymes_representation=[
                                          TailoringRepresentation("pseudo_isomerase", "ISOMERASE_DOUBLE_BOND_SHIFT",
                                                                  [["C_13", "C_14", "C_14", "C_15"]]),
                                          TailoringRepresentation("prenyltransferase", "PRENYLTRANSFERASE",
                                                                  [["C_16"]], "DIMETHYLALLYL")])

    alkaloid_cluster = Alkaloid_Cluster("phenylalanine",
                                        tailoring_enzymes_representation=[TailoringRepresentation("pseudo_decarboxylase", "DECARBOXYLASE", [["C_9"]]),
                                                                          TailoringRepresentation(
                                                                              "pseudo_hydroxylase", "PRENYLTRANSFERASE", [["C_7"]], "DIMETHYLALLYL"),
                                                                          TailoringRepresentation(
                                                                              "pseudo_decarboxylase", "HALOGENASE", [["C_10"]], "Cl"),
                                                                          TailoringRepresentation(
                                                                              "pseudo_hydroxylase", "HYDROXYLATION", [["C_6"]]),
                                                                          TailoringRepresentation("methyltransferase", "METHYLTRANSFERASE", [["N_12"], ["C_7"], ["O_25"]]),

                                                                          ])

    cluster_repr = ClusterRepresentation([ModuleRepresentation("PKS", "PKS_CIS", "ACETYL_COA",
                                                               [DomainRepresentation("Gene 1", 'AT', None, None, True,
                                                                                     True),
                                                                DomainRepresentation("Gene 1", 'ACP', None, None, True,
                                                                                     True)
                                                                ]),
                                          ModuleRepresentation("PKS", "PKS_ITER", "METHYLMALONYL_COA",
                                                               [DomainRepresentation("Gene 1", 'KS',
                                                                                     None, None, True,
                                                                                     True),
                                                                DomainRepresentation("Gene 1", 'AT', None, None, True,
                                                                                     True),
                                                                DomainRepresentation("Gene 1", 'AT', None, None, True,
                                                                                     False),
                                                                DomainRepresentation("Gene 1", 'DH', None, None, True,
                                                                                     True),
                                                                DomainRepresentation("Gene 1", 'ER', None, None, True,
                                                                                     True),
                                                                DomainRepresentation("Gene 1", 'ACP', None, None, True,
                                                                                     True)
                                                                ], 5),
                                         ModuleRepresentation("PKS", "PKS_CIS", "METHYLMALONYL_COA",
                                                              [DomainRepresentation("Gene 1", 'KS',
                                                                                    None, None, True,
                                                                                    True),
                                                               DomainRepresentation("Gene 1", 'AT', None, None, True,
                                                                                    True),
                                                               DomainRepresentation("Gene 1", 'TE',
                                                                                    None, None, True,
                                                                                    True)

                                                               ])]
                                         )
    # draw_cluster(cluster_repr, outfile = "iterative_pks.svg")
    # draw_ripp_structure(ripp_cluster)
    #ripp_cluster.draw_precursor(as_string= False, out_file= "bubbles.svg")

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
    # print(get_tailoring_sites_atom_names(
    #     thiopeptide_cluster_thiomuracin.chain_intermediate))
    # thiopeptide_cluster_thiomuracin.draw_product(
    #     as_string=False, out_file="peptide_test_thiopeptide_peptide.svg")
    # thiopeptide_cluster_thiomuracin.do_tailoring()
    # thiopeptide_cluster_thiomuracin.draw_product(
    #     as_string=False, out_file="tailored_test_thiopeptide_peptide.svg", draw_straightened=False)
    # print(get_tailoring_sites_atom_names(
    #     thiopeptide_cluster_thiomuracin.chain_intermediate))

    # sancti_peptide_cluster_thurincin.make_peptide()
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
    # print(get_tailoring_sites_atom_names(
    #     sliceotide_cluster.chain_intermediate))
    # sliceotide_cluster.draw_product(
    #     as_string=False, out_file="peptide_test_sliceotide_cluster.svg")
    # sliceotide_cluster.draw_precursor_with_modified_product(as_string=False, out_file="bubbles.svg")
    # sliceotide_cluster.do_tailoring()
    # sliceotide_cluster.draw_product(
    #     as_string=False, out_file="tailored_test_sliceotide_cluster.svg", draw_straightened=False)
    #lasso_peptide_cluster.do_tailoring()
    # draw_terpene_structure(terpene_cluster)
    draw_alkaloid_structure(alkaloid_cluster)
