from typing import List, Tuple, Union
from raichu.cluster import Cluster
from raichu.domain.domain import TailoringDomain, CarrierDomain, SynthesisDomain, RecognitionDomain, \
    TerminationDomain, UnknownDomain
from raichu.module import PKSModuleSubtype, NRPSModule, LinearPKSModule, IterativePKSModule, TransATPKSModule,\
    ModuleType
from raichu.domain.domain_types import TailoringDomainType, TerminationDomainType, CarrierDomainType, \
    SynthesisDomainType, RecognitionDomainType, DomainSuperClass, KRDomainSubtype, KSDomainSubtype


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


def build_cluster(cluster_repr: List[Tuple[str, str, Union[str, None], str,
                                           List[Tuple[str, Union[str, None], bool, bool]]]],
                  strict: bool = False) -> Cluster:
    for i, module_repr in enumerate(cluster_repr):
        gene_name, module_type, module_subtype, substrate, domain_reprs = module_repr
        domains = []
        for domain_type, domain_subtype, domain_active, domain_used in domain_reprs:
            domain = make_domain(domain_type, domain_subtype, domain_active, domain_used, substrate, strict=strict)
            domains.append(domain)

        module_type = ModuleType.from_string(module_type)
        if module_type.name == 'PKS':
            module_subtype = PKSModuleSubtype.from_string()



def make_domain(domain_type, domain_subtype, domain_active, domain_used, substrate, strict=False):
    domain_class = DOMAIN_TO_SUPERTYPE.get(domain_type)
    if domain_class:
        if domain_class == RecognitionDomain:
            domain = domain_class(domain_type, substrate, domain_subtype=domain_subtype, active=domain_active,
                                  used=domain_used)
        else:
            domain = domain_class(domain_type, domain_subtype=domain_subtype, active=domain_active,
                                  used=domain_used)
    elif strict:
        raise ValueError(f"Unrecognised domain type: {domain_type}")
    else:
        domain = UnknownDomain(domain_type, domain_subtype=domain_subtype, active=domain_active, used=False)




def get_spaghettis(cluster_repr: List[Tuple[str, str, str, List[str]]]) -> List[str]:
    spaghettis = []
    cluster = build_cluster(cluster_repr)


    return spaghettis



# Each module is a tuple:
# 0: gene name (str)
# 1: module type (str, 'NRPS' or 'PKS')
# 2: substrate (str, SMILES)
# 3: tailoring domains (list of str)