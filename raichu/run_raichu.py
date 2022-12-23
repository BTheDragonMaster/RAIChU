from typing import List, Union

from raichu.cluster import Cluster
from raichu.ripp import RiPP_Cluster
from raichu.domain.domain import TailoringDomain, CarrierDomain, SynthesisDomain, RecognitionDomain, \
    TerminationDomain, UnknownDomain, Domain
from raichu.module import PKSModuleSubtype, NRPSModule, LinearPKSModule, IterativePKSModule, TransATPKSModule,\
    ModuleType
from raichu.domain.domain_types import TailoringDomainType, TerminationDomainType, CarrierDomainType, \
    SynthesisDomainType, RecognitionDomainType
from dataclasses import dataclass
from raichu.reactions.chain_release import find_all_o_n_atoms_for_cyclization
from raichu.reactions.general_tailoring_reactions import find_atoms_for_tailoring


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
    atoms: List[List[str]] # Some tailoring reactions involve more than one atom
    
    
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
                        raise ValueError(f"Gene name '{previous_gene}' already assigned to upstream domain(s).")

                domain.set_gene(domain_repr.gene_name)
                genes.add(domain_repr.gene_name)
            domains.append(domain)
            previous_domain = domain
        module_type = ModuleType.from_string(module_repr.type)
        if module_type.name == 'PKS':
            if module_repr.subtype is not None:
                module_subtype = PKSModuleSubtype.from_string(module_repr.subtype)

                if module_subtype.name == 'PKS_CIS':
                    module = LinearPKSModule(i, domains, starter=starter, terminator=terminator)
                elif module_subtype.name == 'PKS_TRANS':
                    module = TransATPKSModule(i, domains, starter=starter, terminator=terminator)
                elif module_subtype.name == 'PKS_ITER':
                    module = IterativePKSModule(i, domains, starter=starter, terminator=terminator)
                else:
                    raise ValueError(f"Unrecognised PKS module subtype: {module_subtype}.")
            else:
                raise ValueError("PKS module subtype must be specified.")

        elif module_type.name == 'NRPS':
            if module_repr.subtype is None:
                module = NRPSModule(i, domains, starter=starter, terminator=terminator)
            else:
                raise ValueError("NRPS module subtypes are currently not supported. Please pass None.")
        else:
            raise ValueError(f"Unrecognised module type: {module_repr.type}")

        modules.append(module)
    cluster = Cluster(modules, cluster_repr.tailoring_enzymes)

    return cluster

def draw_cluster(cluster_repr: ClusterRepresentation, outfile=None) -> None:
    cluster = build_cluster(cluster_repr)
    cluster.compute_structures(compute_cyclic_products=False)
    cluster.do_tailoring()
    tailored_product = cluster.draw_spaghettis()[-1]
    print (tailored_product)
    if outfile:
        return cluster.draw_cluster(as_string=False, out_file=outfile)
    else:
        cluster.draw_cluster()

def draw_ripp_structure(ripp_cluster: RiPP_Cluster) -> None:
    ripp_cluster.make_peptide()
    ripp_cluster.do_tailoring()
    ripp_cluster.draw_product()
    ripp_cluster.do_macrocyclization()
    ripp_cluster.draw_product()
    ripp_cluster.do_proteolytic_claevage()
    ripp_cluster.draw_product()
    
def get_spaghettis(cluster_repr: ClusterRepresentation) -> List[str]:

    cluster = build_cluster(cluster_repr)
    # for module in cluster.modules:
    #     print(module.domains)
    cluster.compute_structures(compute_cyclic_products=False)
    cluster.do_tailoring()
    cluster.draw_cluster()
    spaghettis = cluster.draw_spaghettis()

    return spaghettis

if __name__ == "__main__":
    ripp_cluster = RiPP_Cluster("best_gene", "mkaekslkayawyiwy", cleavage_sites=[CleavageSiteRepresentation("Y", 10, "leader")], macrocyclisations= [MacrocyclizationRepresentation("O_11","N_161")])
    cluster_repr = ClusterRepresentation([ModuleRepresentation("PKS", "PKS_CIS", "ACETYL_COA",
                                                               [DomainRepresentation("Gene 1", 'AT', None, None, True,
                                                                                     True),
                                                                DomainRepresentation("Gene 1", 'ACP', None, None, True,
                                                                                     True)
                                                                ]),
                                          ModuleRepresentation("PKS", "PKS_CIS", "METHYLMALONYL_COA",
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
                                                                ]),
                                          ModuleRepresentation("NRPS", None, "glycine",
                                                               [DomainRepresentation("Gene 1", 'C', None, None, True,
                                                                                     True),
                                                                DomainRepresentation("Gene 1", 'nMT', None, None, True,
                                                                                     True),
                                                                DomainRepresentation("Gene 1", 'nMT', None, None, True,
                                                                                     False),
                                                                DomainRepresentation("Gene 1", 'A', None, None, True,
                                                                                     True),
                                                                DomainRepresentation("Gene 1", 'PCP', None, None, True,
                                                                                     True),
                                                                DomainRepresentation("Gene 1", 'E', None, None, True,
                                                                                     True)
                                                                ]),
                                          ModuleRepresentation("PKS", "PKS_CIS", "METHYLMALONYL_COA",
                                                               [DomainRepresentation("Gene 1", 'AT', None, None, True,
                                                                                     True),
                                                                DomainRepresentation("Gene 2", 'KR', 'A1', None, True,
                                                                                     True),
                                                                DomainRepresentation("Gene 2", 'DH', None, None, True,
                                                                                     True),
                                                                DomainRepresentation("Gene 2", 'ACP', None, None, True,
                                                                                     True)
                                                                ]),
                                          ModuleRepresentation("PKS", "PKS_CIS", "METHYLMALONYL_COA",
                                                               [DomainRepresentation("Gene 2", 'AT', None, None, True,
                                                                                     True),
                                                                DomainRepresentation("Gene 3", 'KR', 'A1', None, True,
                                                                                     True),
                                                                DomainRepresentation("Gene 3", 'DH', None, None, True,
                                                                                     True),
                                                                DomainRepresentation("Gene 3", 'ACP', None, None, True,
                                                                                     True)
                                                                ]),
                                          ModuleRepresentation("NRPS", None, "tryptophan",
                                                               [DomainRepresentation("Gene 4", 'C', None, None, True,
                                                                                     True),
                                                                DomainRepresentation("Gene 4", 'C', None, None, True,
                                                                                     False),
                                                                DomainRepresentation("Gene 4", 'A', None, None, True,
                                                                                     True),
                                                                DomainRepresentation("Gene 4", 'A', None, None, True,
                                                                                     False),
                                                                DomainRepresentation("Gene 4", 'PCP', None, None, True,
                                                                                     True),
                                                                DomainRepresentation("Gene 4", 'PCP', None, None, True,
                                                                                     False)
                                                                ]),

                                          ], [TailoringRepresentation("gene_7", "P450_OXIDATIVE_BOND_FORMATION", [["C_34","C_36"]])]
                                          )
    draw_cluster(cluster_repr, "nrps_cluster.svg")
   # draw_ripp_structure(ripp_cluster)
