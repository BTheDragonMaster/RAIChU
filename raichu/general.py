from raichu.cluster.modular_cluster import ModularCluster
from raichu.cluster.ripp_cluster import RiPPCluster
from raichu.cluster.terpene_cluster import TerpeneCluster
from raichu.run_raichu import make_domain, build_cluster
from raichu.domain.domain import TailoringDomain, CarrierDomain, SynthesisDomain, RecognitionDomain, \
    TerminationDomain, UnknownDomain, Domain
from raichu.module import PKSModuleSubtype, NRPSModule, LinearPKSModule, IterativePKSModule, TransATPKSModule,\
    ModuleType
from raichu.domain.domain_types import TailoringDomainType, TerminationDomainType, CarrierDomainType, \
    SynthesisDomainType, RecognitionDomainType
from raichu.cluster.alkaloid_cluster import AlkaloidCluster
from raichu.tailoring_enzymes import TailoringEnzyme, TailoringEnzymeType
from raichu.representations import *
from raichu.antismash import load_antismash_gbk, parse_antismash_to_cluster_file

def draw_cluster_from_antiSMASH_output(antismash_file: str, out_file: str) -> None:
    draw_cluster_from_modular_cluster_representation(
        load_antismash_gbk(antismash_file),
        out_file=out_file,
    )


def draw_cluster_from_modular_cluster_representation(cluster_repr: ClusterRepresentation, out_file=None, colour_by_module=True) -> None:
    cluster = build_cluster(cluster_repr, strict= False)
    cluster.compute_structures(compute_cyclic_products=False)
    cluster.do_tailoring()

    if out_file:
        return cluster.draw_cluster(as_string=False, out_file=out_file, colour_by_module=colour_by_module)
    else:
        return cluster.draw_cluster()


def draw_products_from_modular_cluster_representation(cluster_repr: ClusterRepresentation, out_dir) -> None:
    cluster = build_cluster(cluster_repr)
    cluster.compute_structures(compute_cyclic_products=False)
    cluster.do_tailoring()
    cluster.cyclise_all()
    cluster.draw_all_products(out_dir)


def draw_product_from_modular_cluster_representation(cluster_repr: ClusterRepresentation, out_file) -> None:
    cluster = build_cluster(cluster_repr)
    cluster.compute_structures(compute_cyclic_products=True)
    cluster.do_tailoring()
    cluster.draw_product(as_string=False, out_file=out_file)


def draw_ripp_structure_from_ripp_cluster(ripp_cluster: RiPPCluster, out_folder: str) -> None:
    if not os.path.exists(out_folder):
        os.mkdir(out_folder)

    ripp_cluster.make_peptide()
    ripp_cluster.draw_product(
        as_string=False, out_file=os.path.join(out_folder, "peptide_test_ripp.svg"))
    ripp_cluster.draw_cluster(as_string=False, out_file=os.path.join(out_folder, 'ripp_inline.svg'))

    order = []

    ripp_cluster.do_tailoring()
    if ripp_cluster.tailoring_representations:
        ripp_cluster.draw_tailoring(out_file=os.path.join(out_folder, "ripp_tailoring_steps.svg"))
        ripp_cluster.draw_product(
            as_string=False, out_file=os.path.join(out_folder, "ripp_tailoring.svg"))
        order.append("tailoring")
    if ripp_cluster.macrocyclisation_representations:
        ripp_cluster.do_macrocyclization()
        ripp_cluster.draw_product(
            as_string=False, out_file=os.path.join(out_folder, "ripp_macrocyclisation.svg"))
        order.append('cyclisation')

    if ripp_cluster.cleavage_sites:
        ripp_cluster.do_proteolytic_cleavage()
        ripp_cluster.draw_product(
            as_string=False, out_file=os.path.join(out_folder, "cleavage_test_ripp.svg"))
        order.append('cleavage')

    ripp_cluster.draw_pathway(out_file=os.path.join(out_folder, "ripp_pathway.svg"), order=tuple(order))


def draw_terpene_structure_from_terpene_cluster(terpene_cluster: TerpeneCluster, out_folder: str) -> None:
    terpene_cluster.create_precursor()
    terpene_cluster.draw_product(
        as_string=False, out_file=os.path.join(out_folder, "precursor_test_terpene.svg"))
    terpene_cluster.do_macrocyclization()
    terpene_cluster.draw_product(
        as_string=False, out_file=os.path.join(out_folder, "macroyclisation_test_terpene.svg"))
    terpene_cluster.do_tailoring()
    terpene_cluster.draw_product(
        as_string=False, out_file=os.path.join(out_folder, "tailoring_test_terpene.svg"))


def draw_alkaloid_structure_from_alkaloid_cluster(alkaloid_cluster: AlkaloidCluster, out_folder: str) -> None:
    alkaloid_cluster.make_scaffold()
    alkaloid_cluster.draw_product(
        as_string=False, out_file=os.path.join(out_folder, "precursor_test_alkaloid.svg"))
    alkaloid_cluster.do_tailoring()
    alkaloid_cluster.draw_product(
        as_string=False, out_file=os.path.join(out_folder, "tailoring_test_alkaloid.svg"))


def get_tailoring_sites(structure, enzyme_name='all', out_file=None):
    tailoring_sites = {}
    if enzyme_name == 'all':
        for enzyme_type in TailoringEnzymeType:
            tailoring_enzyme = TailoringEnzyme("gene", enzyme_type.name)
            tailoring_sites[enzyme_type.name] = tailoring_enzyme.get_possible_sites(
                structure)
    else:
        for enzyme_type in TailoringEnzymeType:
            if enzyme_type.name == enzyme_name:
                tailoring_enzyme = TailoringEnzyme("gene", enzyme_type.name)
                tailoring_sites[enzyme_type.name] = tailoring_enzyme.get_possible_sites(structure, out_file=out_file)
    return tailoring_sites


def get_tailoring_sites_atom_names(structure):
    tailoring_sites = {}
    for enzyme_type in TailoringEnzymeType:
        tailoring_enzyme = TailoringEnzyme("gene", enzyme_type.name)
        tailoring_sites[enzyme_type.name] = [str(atom) if type(atom) != list else [str(subatom) for subatom in atom] for atom in tailoring_enzyme.get_possible_sites(
            structure)]
    return tailoring_sites
