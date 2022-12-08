import random
from raichu.run_raichu import DomainRepresentation, ModuleRepresentation, ClusterRepresentation, draw_cluster
from raichu.domain.domain_types import KRDomainSubtype, KSDomainSubtype, TailoringDomainType
from raichu.substrate import PksStarterSubstrate, PksElongationSubstrate
from raichu.module import PKSDomainType, NRPSDomainType

TAILORING_DOMAINS_PKS = list(PKSDomainType.__members__)
TAILORING_DOMAINS_NRPS = list(NRPSDomainType.__members__)
KR_SUBTYPES = list(KRDomainSubtype.__members__)
KS_SUBTYPES = list(KSDomainSubtype.__members__)
PKS_STARTERS = list(PksStarterSubstrate.__members__)
PKS_ELONGATION = list(PksElongationSubstrate.__members__)
AMINO_ACIDS = ['alanine', 'arginine', 'asparagine',
           'cysteine', 'glutamine', 'glycine',
           'histidine', 'isoleucine', 'leucine', 'lysine',
           'methionine', 'phenylalanine', 'proline', 'serine',
           'threonine', 'tryptophan', 'tyrosine', 'valine', 'glycine']

def make_PKS_module_starter(number_of_domains):
    module_type = random.choice(["PKS_CIS", "PKS_TRANS"])
    domains = []
    if module_type == "PKS_CIS":
        domains += [DomainRepresentation("gene 1", 'AT', None, None, True, True)]
    list_domains = random.sample(TAILORING_DOMAINS_PKS, number_of_domains-2)
    for domain_type in list_domains:
        used = random.choice([True, False])
        active = random.choice([True, False])
        if domain_type in ["ACP", "AT", "TE", "KS", "UNKNOWN", "TD"]:
            continue
        if domain_type == "KR":
            subtype = random.choice(KR_SUBTYPES)
        else:
            subtype = None
        if "DUMMY" in domain_type:
            continue
        domain = DomainRepresentation("gene 1", domain_type, subtype, None, active, used)
        domains += [domain]
    domains += [DomainRepresentation("gene 1", 'ACP', None, None, True, True)]
    substrate = random.choice(PKS_STARTERS)
    module = ModuleRepresentation("PKS", module_type, substrate, domains)
    return module

def make_PKS_module(number_of_domains, TE):
    module_type = random.choice(["PKS_CIS", "PKS_TRANS"])
    domains = []
    if module_type == "PKS_CIS":
        domains += [DomainRepresentation("gene 1", 'AT', None, None, True, True)]
        domains += [DomainRepresentation("gene 1", 'KS', None, None, True, True)]
    if module_type == "PKS_TRANS":
        subtype = random.choice(KS_SUBTYPES)
        if subtype == "CIS":
            subtype = "DB"
        domains += [DomainRepresentation("gene 1", 'KS', subtype, None, True, True)]
    if TE:
        domains += [DomainRepresentation("gene 1", 'TD', None, None, True, True)]
    list_domains = random.sample(TAILORING_DOMAINS_PKS, number_of_domains-2)
    for domain_type in list_domains:
        if domain_type in ["ACP", "AT", "TE", "KS", "UNKNOWN", "TD"]:
            continue
        used = random.choice([True, False])
        active = random.choice([True, False])
        if domain_type == "KR":
            subtype = random.choice(KR_SUBTYPES)
        else:
            subtype = None
        if "DUMMY" in domain_type:
            continue
        domain = DomainRepresentation("gene 1", domain_type, subtype, None, active, used)
        domains += [domain]
    domains += [DomainRepresentation("gene 1", 'ACP', None, None, True, True)]
    substrate = random.choice(PKS_ELONGATION)
    module = ModuleRepresentation("PKS", module_type, substrate, domains)
    return module

def make_trans_AT_PKS_module(number_of_domains, TE):
    module_type = "PKS_TRANS"
    domains = []
    if module_type == "PKS_CIS":
        domains += [DomainRepresentation("gene 1", 'AT', None, None, True, True)]
        domains += [DomainRepresentation("gene 1", 'KS', None, None, True, True)]
    if module_type == "PKS_TRANS":
        subtype = random.choice(KS_SUBTYPES)
        print(subtype)
        if subtype == "CIS":
            subtype = "DB"
        domains += [DomainRepresentation("gene 1", 'KS', subtype, None, True, True)]
    if TE:
        domains += [DomainRepresentation("gene 1", 'TD', None, None, True, True)]
    list_domains = random.sample(TAILORING_DOMAINS_PKS, number_of_domains-2)
    for domain_type in list_domains:
        if domain_type in ["ACP", "AT", "TE", "KS", "UNKNOWN", "TD"]:
            continue
        used = random.choice([True, False])
        active = random.choice([True, False])
        if domain_type == "KR":
            subtype = random.choice(KR_SUBTYPES)
        else:
            subtype = None
        if "DUMMY" in domain_type:
            continue
        domain = DomainRepresentation("gene 1", domain_type, subtype, None, active, used)
        domains += [domain]
    domains += [DomainRepresentation("gene 1", 'ACP', None, None, True, True)]
    substrate = random.choice(PKS_ELONGATION)
    module = ModuleRepresentation("PKS", module_type, substrate, domains)
    return module

def make_NRPS_module_starter(number_of_domains):
    module_type = None
    domains = []
    domains += [DomainRepresentation("gene 1", 'A', None, None, True, True)]
    list_domains = random.sample(TAILORING_DOMAINS_NRPS, number_of_domains-2)
    for domain_type in list_domains:
        if domain_type in ["C", "A", "TE", "PCP", "UNKNOWN", "TD"]:
            continue
        used = random.choice([True, False])
        active = random.choice([True, False])
        subtype = None
        domain = DomainRepresentation("gene 1", domain_type, subtype, None, active, used)
        domains += [domain]
    domains += [DomainRepresentation("gene 1", 'PCP', None, None, True, True)]
    substrate = random.choice(AMINO_ACIDS)
    module = ModuleRepresentation("NRPS", module_type, substrate, domains)
    return module

def make_NRPS_module(number_of_domains, TE):
    module_type = None
    domains = []
    domains += [DomainRepresentation("gene 1", 'C', None, None, True, True)]
    domains += [DomainRepresentation("gene 1", 'A', None, None, True, True)]
    if TE:
        domains += [DomainRepresentation("gene 1", 'TE', None, None, True, True)]
    list_domains = random.sample(TAILORING_DOMAINS_NRPS, number_of_domains-2)
    for domain_type in list_domains:
        if domain_type in ["C", "A", "TE", "PCP", "UNKNOWN", "TD"]:
            continue
        used = random.choice([True, False])
        active = random.choice([True, False])
        subtype = None
        domain = DomainRepresentation("gene 1", domain_type, subtype, None, active, used)
        domains += [domain]
    domains += [DomainRepresentation("gene 1", 'PCP', None, None, True, True)]
    substrate = random.choice(AMINO_ACIDS)
    module = ModuleRepresentation("NRPS", module_type, substrate, domains)
    return module


def create_random_cluster(number_of_modules):
    modules=[]
    starter_type = random.choice (["nrps", "pks"])
    if starter_type == "nrps":
        modules += [make_NRPS_module_starter(random.randint(2,7))]
    else:
        modules += [make_PKS_module_starter(random.randint(2,7))]
    for i in range(0,number_of_modules-2):
        type = random.choice (["nrps", "pks"])
        if type == "nrps":
            modules += [make_NRPS_module(random.randint(2,7), False)]
        else:
            modules += [make_PKS_module(random.randint(2,7), False)]
    type_termination = random.choice (["nrps", "pks"])
    if type_termination == "nrps":
        modules += [make_NRPS_module(random.randint(2,7), True)]
    else:
        modules += [make_PKS_module(random.randint(2,7), True)]
    cluster_representation = ClusterRepresentation(modules)
    return cluster_representation

def create_random_cluster_trans_at(number_of_modules):
    modules=[]
    starter_type = random.choice (["nrps", "pks"])
    if starter_type == "nrps":
        modules += [make_NRPS_module_starter(random.randint(2,7))]
    else:
        modules += [make_PKS_module_starter(random.randint(2,7))]
    for i in range(0,number_of_modules-2):
        type = random.choices (["nrps", "pks"], weights=[20,80])
        if type == "nrps":
            modules += [make_NRPS_module(random.randint(2,7), False)]
        else:
            modules += [make_trans_AT_PKS_module(random.randint(4,7), False)]
    type_termination = random.choice (["nrps", "pks"])
    if type_termination == "nrps":
        modules += [make_NRPS_module(random.randint(2,7), True)]
    else:
        modules += [make_trans_AT_PKS_module(random.randint(4,7), True)]
    cluster_representation = ClusterRepresentation(modules)
    return cluster_representation

if __name__ == "__main__":
    for i in range(1, 10):
        cluster_repr = create_random_cluster_trans_at(random.randint(5,7))
        draw_cluster(cluster_repr,f'demo_cluster_{i}.svg')
