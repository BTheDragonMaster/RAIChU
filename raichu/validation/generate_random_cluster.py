import random
from sys import argv

from raichu.run_raichu import DomainRepresentation, ModuleRepresentation, ClusterRepresentation, draw_cluster
from raichu.domain.domain_types import KRDomainSubtype, KSDomainSubtype
from raichu.substrate import PksStarterSubstrate, PksElongationSubstrate
from raichu.module import PKSDomainType, NRPSDomainType

TAILORING_DOMAINS_PKS = list(PKSDomainType.__members__)
TAILORING_DOMAINS_NRPS = list(NRPSDomainType.__members__)
KR_SUBTYPES = list(KRDomainSubtype.__members__)
KS_SUBTYPES = list(KSDomainSubtype.__members__)
PKS_STARTERS = list(PksStarterSubstrate.__members__)
PKS_ELONGATION = list(PksElongationSubstrate.__members__)
AMINO_ACIDS = ['alanine', 'arginine', 'asparagine', 'cysteine', 'glutamine', 'glycine', 'histidine', 'isoleucine',
               'leucine', 'lysine', 'methionine', 'phenylalanine', 'proline', 'serine', 'threonine', 'tryptophan',
               'tyrosine', 'valine', 'glycine']


def make_pks_starter_module(number_of_domains):
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
    domains.append(DomainRepresentation("gene 1", 'ACP', None, None, True, True))
    substrate = random.choice(PKS_STARTERS)
    module = ModuleRepresentation("PKS", module_type, substrate, domains)
    return module


def make_pks_module(number_of_domains, te):
    module_type = random.choice(["PKS_CIS", "PKS_TRANS"])
    domains = []
    if module_type == "PKS_CIS":
        domains.append(DomainRepresentation("gene 1", 'AT', None, None, True, True))
        domains.append(DomainRepresentation("gene 1", 'KS', None, None, True, True))
    if module_type == "PKS_TRANS":
        subtype = random.choice(KS_SUBTYPES)
        if subtype == "CIS":
            subtype = "DB"
        domains.append(DomainRepresentation("gene 1", 'KS', subtype, None, True, True))

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
        domains.append(domain)
    domains.append(DomainRepresentation("gene 1", 'ACP', None, None, True, True))
    if te:
        domains.append(DomainRepresentation("gene 1", 'TE', None, None, True, True))
    substrate = random.choice(PKS_ELONGATION)
    module = ModuleRepresentation("PKS", module_type, substrate, domains)
    return module


def make_trans_at_pks_module(number_of_domains, te):
    module_type = random.choices(["PKS_CIS", "PKS_TRANS"], weights=[20, 80])[0]

    domains = []
    if module_type == "PKS_CIS":
        domains.append(DomainRepresentation("gene 1", 'AT', None, None, True, True))
        domains.append(DomainRepresentation("gene 1", 'KS', None, None, True, True))
    if module_type == "PKS_TRANS":
        subtype = random.choice(KS_SUBTYPES)
        if subtype == "CIS":
            subtype = "DB"
        domains.append(DomainRepresentation("gene 1", 'KS', subtype, None, True, True))
    list_domains = random.sample(TAILORING_DOMAINS_PKS, number_of_domains-2)
    for domain_type in list_domains:
        if domain_type in ["ACP", "AT", "TE", "KS", "UNKNOWN", "TD", "CYC"]:
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
        domains.append(domain)
    domains.append(DomainRepresentation("gene 1", 'ACP', None, None, True, True))
    if te:
        domains.append(DomainRepresentation("gene 1", 'TD', None, None, True, True))
    substrate = random.choice(PKS_ELONGATION)
    module = ModuleRepresentation("PKS", module_type, substrate, domains)
    return module


def make_nrps_starter_module(number_of_domains):
    module_type = None
    domains = [DomainRepresentation("gene 1", 'A', None, None, True, True)]
    list_domains = random.sample(TAILORING_DOMAINS_NRPS, number_of_domains - 2)
    for domain_type in list_domains:
        if domain_type in ["C", "A", "TE", "PCP", "UNKNOWN", "TD"]:
            continue
        used = random.choice([True, False])
        active = random.choice([True, False])
        subtype = None
        domain = DomainRepresentation("gene 1", domain_type, subtype, None, active, used)
        domains.append(domain)
    domains.append(DomainRepresentation("gene 1", 'PCP', None, None, True, True))
    substrate = random.choice(AMINO_ACIDS)
    module = ModuleRepresentation("NRPS", module_type, substrate, domains)
    return module


def make_nrps_module(number_of_domains, te):
    module_type = None
    domains = [DomainRepresentation("gene 1", 'C', None, None, True, True),
               DomainRepresentation("gene 1", 'A', None, None, True, True)]

    list_domains = random.sample(TAILORING_DOMAINS_NRPS, number_of_domains-2)
    for domain_type in list_domains:
        if domain_type in ["C", "A", "TE", "PCP", "UNKNOWN", "TD"]:
            continue
        used = random.choice([True, False])
        active = random.choice([True, False])
        subtype = None
        domain = DomainRepresentation("gene 1", domain_type, subtype, None, active, used)
        domains.append(domain)
    domains.append(DomainRepresentation("gene 1", 'PCP', None, None, True, True))

    if te:
        domains.append(DomainRepresentation("gene 1", 'TE', None, None, True, True))
    substrate = random.choice(AMINO_ACIDS)
    module = ModuleRepresentation("NRPS", module_type, substrate, domains)
    return module


def create_random_cluster(number_of_modules):
    modules = []
    starter_type = random.choice(["nrps", "pks"])
    if starter_type == "nrps":
        modules.append(make_nrps_starter_module(random.randint(2, 7)))
    else:
        modules.append(make_pks_starter_module(random.randint(2, 7)))

    for _ in range(0, number_of_modules - 2):
        module_type = random.choice(["nrps", "pks"])
        if module_type == "nrps":
            modules.append(make_nrps_module(random.randint(2, 7), False))
        else:
            modules.append(make_pks_module(random.randint(2, 7), False))
    type_termination = random.choice(["nrps", "pks"])
    if type_termination == "nrps":
        modules.append(make_nrps_module(random.randint(2, 7), True))
    else:
        modules.append(make_pks_module(random.randint(2, 7), True))
    cluster_representation = ClusterRepresentation(modules)
    return cluster_representation


def create_random_cluster_trans_at(number_of_modules):
    modules = []
    starter_type = random.choice(["nrps", "pks"])
    if starter_type == "nrps":
        modules.append(make_nrps_starter_module(random.randint(2, 7)))
    else:
        modules.append(make_pks_starter_module(random.randint(2, 7)))
    for _ in range(0, number_of_modules - 2):
        module_type = random.choices(["nrps", "pks"], weights=[20, 80])[0]
        if module_type == "nrps":
            modules.append(make_nrps_module(random.randint(2, 7), False))
        else:
            modules.append(make_trans_at_pks_module(random.randint(9, 12), False))
    type_termination = random.choice(["nrps", "pks"])
    if type_termination == "nrps":
        modules.append(make_nrps_module(random.randint(2, 7), True))
    else:
        modules.append(make_trans_at_pks_module(random.randint(5, 9), True))
    cluster_representation = ClusterRepresentation(modules)
    return cluster_representation


if __name__ == "__main__":
    for i in range(1, int(argv[1]) + 1):
        cluster_repr = create_random_cluster_trans_at(random.randint(5, 7))
        draw_cluster(cluster_repr, f'demo_cluster_transat_{i}.svg')

        cluster_repr = create_random_cluster(random.randint(5, 7))
        draw_cluster(cluster_repr, f'demo_cluster_{i}.svg')
