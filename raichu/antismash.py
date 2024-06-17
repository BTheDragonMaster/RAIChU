
import os
from typing import Optional
from itertools import permutations

from sys import argv
from Bio import SeqIO
from raichu.representations import (
    ClusterRepresentation,
    ModuleRepresentation,
    DomainRepresentation,
)
from raichu.run_raichu import draw_cluster
from raichu.substrate import PksStarterSubstrate
from dataclasses import dataclass, field


NRPS_PKS_RULES = {
    "biosynthetic-additional (rule-based-clusters) AMP-binding",
    "biosynthetic-additional (rule-based-clusters) PP-binding",
    "biosynthetic (rule-based-clusters) NRPS: Condensation",
    "biosynthetic (rule-based-clusters) T1PKS: PKS_AT",
    "biosynthetic (rule-based-clusters) T1PKS: mod_KS",
}

AS_TO_NRPS = {
    "ala": "alanine",
    "arg": "arginine",
    "asn": "asparagine",
    "asp": "aspartic acid",
    "cys": "cysteine",
    "gln": "glutamine",
    "glu": "glutamic acid",
    "gly": "glycine",
    "his": "histidine",
    "ile": "isoleucine",
    "leu": "leucine",
    "lys": "lysine",
    "met": "methionine",
    "phe": "phenylalanine",
    "pro": "proline",
    "ser": "serine",
    "thr": "threonine",
    "trp": "tryptophan",
    "tyr": "tyrosine",
    "val": "valine",
    "3-me-glu": "4-methylglutamic acid",
    "4ppro": "**Unknown**",
    "aad": "2-aminoadipic acid",
    "abu": "2-aminobutyric acid",
    "aeo": "2-amino-9,10-epoxy-8-oxodecanoid acid",
    "ala-b": "beta-alanine",
    "ala-d": "d-alanine",
    "allo-thr": "allo-threonine",
    "b-ala": "beta-alanine",
    "beta-ala": "beta-alanine",
    "bmt": "4-butenyl-4-methyl threonine",
    "cap": "capreomycidine",
    "bht": "**Unknown**",
    "dab": "2,4-diaminobutyric acid",
    "dhb": "2,3-dihydroxybenzoic acid",
    "dhpg": "3,5-dihydroxyphenylglycine",
    "dht": "dehydrobutyrine",
    "dpg": "3,5-dihydroxyphenylglycine",
    "hiv": "2-hydroxyisovalerate",
    "hiv-d": "d-2-hydroxyisovalerate",
    "hmp-d": "**Unknown**",
    "horn": "**Unknown**",
    "hpg": "4-hydroxyphenylglycine",
    "hyv": "4-hydroxyvaline",
    "hyv-d": "**Unknown**",
    "iva": "isovaline",
    "lys-b": "beta-lysine",
    "orn": "ornithine",
    "phg": "phenylglycine",
    "pip": "pipecolic acid",
    "sal": "salicylic acid",
    "tcl": "**Unknown**",
    "vol": "valinol",
    "ldap": "**Unknown**",
    "meval": "tert-leu",
    "alle": "allo-isoleucine",
    "alaninol": "alaninol",
    "n-(1,1-dimethyl-1-allyl)trp": "**Unknown**",
    "d-lyserg": "d-lysergic acid",
    "ser-thr": "**Unknown**",
    "mephe": "**Unknown**",
    "haorn": "**Unknown**",
    "hasn": "**Unknown**",
    "hforn": "**Unknown**",
    "s-nmethoxy-trp": "**Unknown**",
    "alpha-hydroxy-isocaproic-acid": "**Unknown**",
    "mehoval": "**Unknown**",
    "2-oxo-isovaleric-acid": "alpha-ketoisovaleric acid",
    "aoda": "**Unknown**",
    "x": "**Unknown**",
}

antiSMASH_DOMAIN_TO_RAICHU_DOMAIN = {
    "PKS_KS": "KS",
    "PKS_AT": "AT",
    "PKS_PP": "CP",
    "PCP": "CP",
    "ACP": "CP",
    "PKS_KR": "KR",
    "PKS_DH": "DH",
    "PKS_DH2": "DH",
    "PKS_DHt": "DH",
    "PKS_ER": "ER",
    "Thioesterase": "TE",
    "AMP-binding": "A",
    "Epimerization": "E",
    "Condensation": "C",
    "nMT": "nMT",
    "PP-binding": "CP",
    "Heterocyclization": "CYC",
    "CAL_domain": "CAL",
    "TD": "TD"
}

domains_to_be_refined = ["KR", "KS", "A", "AT"]

AS_TO_PKS = {
    "mal": "MALONYL_COA",
    "mmal": "METHYLMALONYL_COA",
    "mxmal": "METHOXYMALONYL_COA",
    "emal": "ETHYLMALONYL_COA",
}

DOMAIN_TO_MODULE_TYPE = {"KS": "PKS",
                         "AT": "PKS",
                         "ER": "PKS",
                         "DH": "PKS",
                         "KR": "PKS",
                         "C": "NRPS",
                         "A": "NRPS",
                         "E": "NRPS",
                         "nMT": "NRPS",
                         "CYC": "NRPS",
                         "CP": "EITHER",
                         "TE": "EITHER",
                         "TD": "EITHER",
                         "CAL": "EITHER",
                         "UNKNOWN": "EITHER"}


@dataclass
class AntiSmashDomain:
    id: str
    type: str
    start: int
    end: int
    strand: str
    gene: str

    subtype: Optional[str] = None
    active: bool = True
    substrate: Optional[str] = None
    name: Optional[str] = None
    raichu_type: Optional[str] = None
    raichu_subtype: Optional[str] = None

    def set_raichu_type(self):
        if self.type in antiSMASH_DOMAIN_TO_RAICHU_DOMAIN:
            self.raichu_type = antiSMASH_DOMAIN_TO_RAICHU_DOMAIN[self.type]
        elif self.type == 'MT' and self.subtype and self.subtype == 'nMT':
            self.raichu_type = 'nMT'
        else:
            self.raichu_type = "UNKNOWN"
            self.name = self.type

    def add_domain_information(self, features):
        if "domain_subtypes" in features:
            self.subtype = features["domain_subtypes"][0]
        if "specificity" in features:
            for spec in features["specificity"]:
                if "consensus:" in spec:
                    specificity = spec.split("consensus:")[1].strip().lower()

                    if specificity in AS_TO_NRPS:
                        substrate_name = AS_TO_NRPS[specificity]
                        self.substrate = substrate_name
                    elif specificity in AS_TO_PKS:
                        substrate_name = AS_TO_PKS[specificity]
                        self.substrate = substrate_name
                elif "KR activity:" in spec:
                    self.active = True if spec.split("KR activity:")[1].strip() == "active" else False

                elif "KR stereochemistry:" in spec:
                    domain_subtype = spec.split("KR stereochemistry:")[1].strip()
                    if domain_subtype == '(unknown)':
                        self.raichu_subtype = "UNKNOWN"
                    else:
                        self.raichu_subtype = domain_subtype

                elif "transATor:" in spec:
                    self.raichu_subtype = spec.split("transATor:")[1].strip().replace("-", "_").replace("/", "_").upper().replace("(UNKNOWN)", "MISCELLANEOUS")
        if not self.active and self.raichu_subtype == "C2":
            self.active = True
            
    def __repr__(self):
        return f"{self.raichu_type}_{self.start}_{self.end}"


@dataclass
class AntiSmashModule:
    start: int
    end: int
    strand: int

    complete: bool = False
    starter_complete: bool = False
    type: Optional[str] = None
    subtype: Optional[str] = None
    domains: list[AntiSmashDomain] = field(default_factory=list)

    def __post_init__(self):
        self.set_module_type()
        self.determine_completeness()

    def __eq__(self, other):
        if self.start == other.start and self.end == other.end and self.strand == other.strand:
            return True
        return False

    def __hash__(self):
        return self.start, self.end

    def __repr__(self):
        return f"{self.start}_{self.end}_{'-'.join([d.raichu_type for d in self.domains])}"

    def set_module_type(self):

        for domain in self.domains:
            if domain.raichu_type == 'C':
                if self.type == 'NRPS' or self.type is None:
                    self.type = 'NRPS'
                else:
                    raise ValueError(f"Can't set module type for module {self.__repr__()}. Cannot have C domain in {self.type} module.")

            if domain.raichu_type == 'CYC':
                if self.type == 'NRPS' or self.type is None:
                    self.type = 'NRPS'
                else:
                    raise ValueError(f"Can't set module type for module {self.__repr__()}. Cannot have heterocyclistion domain in {self.type} module.")

            if domain.raichu_type == 'A':
                if self.type == 'NRPS' or self.type is None:
                    self.type = 'NRPS'
                else:
                    raise ValueError(
                        f"Can't set module type for module {self.__repr__()}. Cannot have A domain in {self.type} module.")
            if domain.raichu_type == 'KR':
                if self.type == 'PKS' or self.type is None:
                    self.type = 'PKS'
                else:
                    raise ValueError(
                        f"Can't set module type for module {self.__repr__()}. Cannot have KR domain in {self.type} module.")

            if domain.raichu_type == 'KS':
                if self.type == 'PKS' or self.type is None:
                    self.type = 'PKS'
                else:
                    raise ValueError(
                        f"Can't set module type for module {self.__repr__()}. Cannot have KS domain in {self.type} module.")

            if domain.raichu_type == 'ER':
                if self.type == 'PKS' or self.type is None:
                    self.type = 'PKS'
                else:
                    raise ValueError(
                        f"Can't set module type for module {self.__repr__()}. Cannot have ER domain in {self.type} module.")

            if domain.raichu_type == 'DH':
                if self.type == 'PKS' or self.type is None:
                    self.type = 'PKS'
                else:
                    raise ValueError(
                        f"Can't set module type for module {self.__repr__()}. Cannot have DH domain in {self.type} module.")

            if domain.raichu_type == 'AT':

                if self.type == 'PKS' or self.type is None:
                    self.type = 'PKS'
                else:
                    raise ValueError(
                        f"Can't set module type for module {self.__repr__()}. Cannot have AT domain in {self.type} module.")

        if self.type == 'PKS':
            has_recognition_domain = False
            for domain in self.domains:
                if domain.raichu_type == 'AT' or domain.raichu_type == 'CAL':
                    has_recognition_domain = True

            if has_recognition_domain:
                self.subtype = 'PKS_CIS'
            else:
                self.subtype = 'PKS_TRANS'

        if self.type is None:
            self.type = "NRPS"

    def determine_completeness(self):
        self.complete = False
        self.starter_complete = False

        if self.type == 'PKS':

            has_ks = False
            has_acp = False
            has_at = False
            has_cal = False

            for domain in self.domains:
                if domain.raichu_type == 'KS':
                    has_ks = True
                if domain.raichu_type == 'CP':
                    has_acp = True
                if domain.raichu_type == 'AT':
                    has_at = True
                if domain.raichu_type == 'CAL':
                    has_cal = True

            if has_acp:
                if self.subtype == 'PKS_TRANS':
                    self.starter_complete = True
                elif self.subtype == 'PKS_CIS' and (has_at or has_cal):
                    self.starter_complete = True

                if has_ks:
                    if self.subtype == 'PKS_TRANS':
                        self.complete = True
                    elif self.subtype == 'PKS_CIS' and has_at:
                        self.complete = True

        elif self.type == 'NRPS':
            has_c = False
            has_cyc = False
            has_a = False
            has_pcp = False
            has_cal = False

            for domain in self.domains:
                if domain.raichu_type == 'C':
                    has_c = True
                if domain.raichu_type == 'CP':
                    has_pcp = True
                if domain.raichu_type == 'A':
                    has_a = True
                if domain.raichu_type == 'CAL':
                    has_cal = True
                if domain.raichu_type == 'CYC':
                    has_cyc = True

            if has_pcp and (has_a or has_cal):
                self.starter_complete = True

                if has_c or has_cyc:
                    self.complete = True

    def is_broken(self, module_nr):
        if module_nr == 0:
            if self.complete or self.starter_complete:
                return False
        elif len(self.domains) == 1 and self.domains[0].raichu_type == 'UNKNOWN':
            return False

        elif self.complete:
            return False

        return True

    def has_termination_domain(self):
        for domain in self.domains:
            if domain.raichu_type in ["TE", "TD"]:
                return True
        return False

    def has_c_starter(self):
        if self.domains and self.domains[0].subtype == "Condensation_Starter":
            return True
        return False

    def is_potential_starter(self):
        if self.starter_complete and not self.complete:
            return True

        return False


@dataclass
class ModuleBlock:
    start: int
    end: int
    strand: int
    modules: list[AntiSmashModule] = field(default_factory=list)
    is_starter_candidate: bool = False
    has_c_starter: bool = False
    has_termination_domain: bool = False

    def __post_init__(self):
        self.label_module_block()

    def __repr__(self):
        return '\n'.join([m.__repr__() for m in self.modules])

    def __eq__(self, other):
        if self.start == other.start and self.end == other.end and self.strand == other.strand:
            return True

        return False

    def __hash__(self):
        return self.start, self.end, self.strand

    def label_module_block(self):
        if self.modules:
            if self.modules[0].is_potential_starter():
                self.is_starter_candidate = True
            if self.modules[0].has_c_starter():
                self.has_c_starter = True
            if self.modules[-1].has_termination_domain():
                self.has_termination_domain = True


class ModuleOrder:

    def __init__(self, module_blocks: list[ModuleBlock]):
        self.module_blocks = module_blocks
        self.block_order = module_blocks
        self.modules_broken = None
        self.modules = []
        self.standalone_modules = []

    def sort_collinear_blocks(self):

        starter_candidates = []
        terminator_candidates = []

        for block in self.module_blocks:
            if block.is_starter_candidate:
                starter_candidates.append(block)

        for block in self.module_blocks:
            if block.has_c_starter:
                starter_candidates.append(block)
            if block.has_termination_domain:
                terminator_candidates.append(block)

        found_optimal_solution = False
        best_block_order = self.block_order[:]
        best_solution = self.check_broken_modules()

        if starter_candidates and terminator_candidates:
            for starter_candidate in starter_candidates:
                for terminator_candidate in terminator_candidates:
                    central_blocks = []
                    for module_block in self.module_blocks:
                        if module_block != starter_candidate and module_block != terminator_candidate:
                            central_blocks.append(module_block)

                    for module_order in permutations(central_blocks):
                        self.block_order = [starter_candidate] + list(module_order) + [terminator_candidate]
                        nr_broken_modules = self.check_broken_modules()
                        if nr_broken_modules < best_solution:
                            best_block_order = self.block_order[:]
                            best_solution = nr_broken_modules
                        if nr_broken_modules == 0:
                            found_optimal_solution = True
                            break

                    if found_optimal_solution:
                        break
                if found_optimal_solution:
                    break
        elif starter_candidates:
            for starter_candidate in starter_candidates:
                central_blocks = []
                for module_block in self.module_blocks:
                    if module_block != starter_candidate:
                        central_blocks.append(module_block)

                for module_order in permutations(central_blocks):
                    self.block_order = [starter_candidate] + list(module_order)
                    nr_broken_modules = self.check_broken_modules()
                    if nr_broken_modules < best_solution:
                        best_block_order = self.block_order[:]
                        best_solution = nr_broken_modules
                    if nr_broken_modules == 0:
                        found_optimal_solution = True
                        break
                if found_optimal_solution:
                    break

        elif terminator_candidates:
            for terminator_candidate in terminator_candidates:
                central_blocks = []
                for module_block in self.module_blocks:
                    if module_block != terminator_candidate:
                        central_blocks.append(module_block)

                for module_order in permutations(central_blocks):
                    self.block_order = list(module_order) + [terminator_candidate]
                    nr_broken_modules = self.check_broken_modules()
                    if nr_broken_modules < best_solution:
                        best_block_order = self.block_order[:]
                        best_solution = nr_broken_modules
                    if nr_broken_modules == 0:
                        found_optimal_solution = True
                        break
                if found_optimal_solution:
                    break

        self.block_order = best_block_order
        self.modules_broken = self.check_broken_modules()

        self.move_broken_modules()

    def remove_standalone_domains(self):

        for module_block in self.module_blocks[:]:
            if len(module_block.modules) == 1 and len(module_block.modules[0].domains) == 1 and \
                    module_block.modules[0].domains[0].raichu_type not in ["TE", "TD"]:
                self.module_blocks.remove(module_block)
                self.standalone_modules.append(module_block.modules[0])

    def move_broken_modules(self):
        module_block = self.recreate_modules()
        modules = module_block.modules[:]
        broken_modules = []
        complete_modules = []
        end_cluster = False

        for i, module in enumerate(module_block.modules):
            if module.is_broken(i) or end_cluster:

                broken_modules.append(module)
            else:
                complete_modules.append(module)

            if module.has_termination_domain() and not module.is_broken(i):
                end_cluster = True

        genes_with_complete_modules = set()

        for module in complete_modules:
            for domain in module.domains:
                genes_with_complete_modules.add(domain.gene)

        movable_broken_modules = []

        for broken_module in broken_modules:
            genes = set()
            for domain in broken_module.domains:
                genes.add(domain.gene)

            if not genes.intersection(genes_with_complete_modules):
                movable_broken_modules.append(broken_module)

        gene_to_broken = {}

        for module in movable_broken_modules:
            for domain in module.domains:
                if domain.gene not in gene_to_broken:
                    gene_to_broken[domain.gene] = []
                if module not in gene_to_broken[domain.gene]:
                    gene_to_broken[domain.gene].append(module)

        for gene, gene_modules in gene_to_broken.items():
            if gene_modules:
                modules_start = min([m.start for m in gene_modules] + [m.end for m in gene_modules])

                for gene_module in gene_modules:
                    modules.remove(gene_module)

                inserted = False
                for i, module in enumerate(modules[:]):
                    if module.start > modules_start:
                        modules = modules[:i] + gene_modules + modules[i:]
                        inserted = True
                        break

                if not inserted:
                    modules += gene_modules

        for standalone_module in self.standalone_modules:
            inserted = False
            for i, module in enumerate(modules[:]):
                if module.start > standalone_module.start:
                    modules = modules[:i] + [standalone_module] + modules[i:]
                    inserted = True
                    break

            if not inserted:
                modules.append(standalone_module)

        self.modules = modules

    def make_raichu_cluster(self):
        module_representations = []
        for i, module in enumerate(self.modules):
            domain_representations = []

            substrates = [domain.substrate for domain in module.domains if domain.substrate is not None]
            if len(substrates) > 0:
                substrate = substrates[0]

                if module.type == "PKS" and i == 0:
                    if substrate not in [v.name for v in PksStarterSubstrate]:
                        substrate = "WILDCARD"
            else:
                if module.type == "NRPS":
                    substrate = "**Unknown**"
                elif module.type == "PKS":
                    substrate = "WILDCARD"
                else:
                    raise ValueError("Module type can only be NRPS or PKS")

            for domain in module.domains:
                if domain.raichu_type == 'CP':
                    if module.type == 'NRPS':
                        domain.raichu_type = "PCP"
                    else:
                        domain.raichu_type = "ACP"

                if module.type == "PKS":
                    if domain.raichu_type == "CAL" and substrate == "WILDCARD":
                        substrate = "**Unknown**"

                domain_representation = DomainRepresentation(domain.gene, domain.raichu_type, domain.raichu_subtype,
                                                             domain.active, True)
                if domain_representation.type in [d.type for d in domain_representations]:
                    domain_representation.used = False
                if domain_representation.type == 'UNKNOWN':
                    domain_representation.used = False
                if domain_representation.type in ["AT", "A", "CAL"]:
                    if ["AT"] in [d.type for d in domain_representations] or \
                            ["A"] in [d.type for d in domain_representations] or \
                            ["CAL"] in [d.type for d in domain_representations]:
                        domain_representation.used = False

                domain_representations.append(domain_representation)
            module_representation = ModuleRepresentation(module.type, module.subtype, substrate, domain_representations)
            module_representations.append(module_representation)
        cluster_representation = ClusterRepresentation(module_representations)
        return cluster_representation

    def recreate_modules(self):
        domains = []
        for block in self.block_order:
            for module in block.modules:
                domains += module.domains

        module_block = make_modules(domains)
        return module_block

    def check_broken_modules(self):
        module_block = self.recreate_modules()
        nr_broken = 0
        end_cluster = False

        for i, module in enumerate(module_block.modules):
            if end_cluster:
                nr_broken += 1
                continue
            if module.is_broken(i):
                nr_broken += 1
            if module.has_termination_domain() and not module.is_broken(i):
                end_cluster = True

        return nr_broken


@dataclass
class AntiSmashGene:
    name: str
    start: int
    end: int
    strand: int
    synonyms: list[str] = field(default_factory=list)
    domains: list[AntiSmashDomain] = field(default_factory=list)
    modules: list[AntiSmashModule] = field(default_factory=list)

    def sort_domains(self):
        self.domains.sort(key=lambda d: d.start)
        if self.strand == -1:
            self.domains.reverse()

    def sort_modules(self):
        self.modules.sort(key=lambda m: m.start)
        if self.strand == -1:
            self.modules.reverse()

        for i, module in enumerate(self.modules):
            module.nr = i

    def __eq__(self, other):
        if self.start == other.start and self.end == other.end and self.strand == other.strand:
            return True
        return False

    def __hash__(self):
        return self.start, self.end

    def __repr__(self):
        return f"{self.name}"


def get_nrps_pks_modules(antismash_gbk):
    as_domains = parse_antismash_domains_gbk(antismash_gbk)
    genes = []
    for record in SeqIO.parse(antismash_gbk, "genbank"):
        for feature in record.features:
            if feature.type == 'CDS':

                qualifiers = feature.qualifiers
                synonyms = []

                if 'gene' in qualifiers:
                    synonyms.append(qualifiers['gene'][0])
                if 'protein' in qualifiers:
                    synonyms.append(qualifiers['protein'][0])
                if 'locus_tag' in qualifiers:
                    synonyms.append(qualifiers['locus_tag'][0])
                if 'protein_id' in qualifiers:
                    synonyms.append(qualifiers['protein_id'][0])

                if synonyms:
                    gene = AntiSmashGene(synonyms[0],
                                         int(feature.location.start),
                                         int(feature.location.end),
                                         int(feature.location.strand),
                                         synonyms)
                    if gene not in genes:
                        genes.append(gene)
                else:
                    print(f"Warning: No identifier found for a CDS. Domains on this gene will not be processed.")
                    # for domain in as_domains:
                    #     if domain.gene in synonyms:

    for domain in as_domains:
        gene_found = False
        for gene in genes:
            if domain.gene in gene.synonyms:
                gene_found = True
                gene.domains.append(domain)
                break
        if not gene_found:
            print(f"Warning: No CDS found for {domain.gene}. Domains on this gene will not be processed.")

    for gene in genes:
        gene.sort_domains()

    gene_groups = sort_genes(genes)

    filtered_groups = []
    for gene_group in gene_groups:
        filtered_groups += filter_genes(gene_group)

    modules_per_collinear_block = ModuleOrder(make_modules_from_gene_groups(filtered_groups))
    if len(filtered_groups) >= 6:
        modules_per_collinear_block.remove_standalone_domains()
    modules_per_collinear_block.sort_collinear_blocks()

    module_blocks = modules_per_collinear_block

    # If there are broken modules, also play around with the gene order
    #
    # if modules_per_collinear_block.modules_broken:
    #     filtered_genes = []
    #
    #     for gene_group in filtered_groups:
    #         for gene in gene_group:
    #             filtered_genes.append([gene])
    #
    #     if len(filtered_genes) < 5:
    #
    #         modules_per_gene = ModuleOrder(make_modules_from_gene_groups(filtered_genes))
    #         modules_per_gene.sort_collinear_blocks()
    #
    #         if modules_per_collinear_block.modules_broken > modules_per_gene.modules_broken:
    #             module_blocks = modules_per_gene

    return module_blocks


def make_modules_from_gene_groups(gene_groups):
    module_groups = []
    for gene_group in gene_groups:
        domains = []

        for gene in gene_group:
            domains += gene.domains

        module_block = make_modules(domains)

        module_groups.append(module_block)

    if not module_groups:
        raise Exception("Cluster is empty. This can happen with Type III PKS clusters. Please check the input file.")

    return module_groups


def make_modules(domains):

    modules = []
    module = []
    end_cluster = False

    type_i_te_domains = []
    type_ii_te_domains = []

    for i, domain in enumerate(domains):
        if domain.raichu_type in ["TD", "TE"]:
            if i == 0:
                type_ii_te_domains.append(domain)
            else:
                if domains[i - 1].gene == domain.gene:
                    type_i_te_domains.append(domain)
                else:
                    type_ii_te_domains.append(domain)

    for i, domain in enumerate(domains):
        if not module and not modules:
            module.append(domain)
            if domain.raichu_type in ["TD", "TE"]:
                modules.append(module)
                module = []
                if not type_i_te_domains:
                    end_cluster = True

        # When encountering a synthesis domain, start a new module

        elif (DOMAIN_TO_MODULE_TYPE[domain.raichu_type] == 'PKS' and "NRPS" in [DOMAIN_TO_MODULE_TYPE[d.raichu_type] for d in module]) or \
                (DOMAIN_TO_MODULE_TYPE[domain.raichu_type] == 'NRPS' and "PKS" in [DOMAIN_TO_MODULE_TYPE[d.raichu_type] for d in module]):
            if module:
                modules.append(module)
                module = [domain]
            else:
                raise Exception("Module can't be empty if it contains a domain. Fix code.")

        elif domain.raichu_type in ["C", "KS", "CYC"]:
            if module:
                # Don't start a new module when the previous domain could be a docking domain
                if module[-1].gene == domain.gene and module[-1].raichu_type == "UNKNOWN" and len(module) == 1:
                    module.append(domain)
                else:
                    modules.append(module)
                    module = [domain]
            else:
                module = [domain]

        elif domain.raichu_type in ["TD", "TE"]:
            if not end_cluster:
                if 'CP' in [d.raichu_type for d in module]:
                    if domain in type_i_te_domains or not type_i_te_domains:
                        module.append(domain)
                        modules.append(module)
                        module = []
                    else:
                        if module:
                            modules.append(module)

                        module = [domain]
                        modules.append(module)
                        module = []

                # Dealing with a standalone TE domain
                else:
                    if module:
                        modules.append(module)

                    module = [domain]
                    modules.append(module)
                    module = []

            # If a previous termination domain has been encountered, always consider the domain as a standalone domain

            else:
                if module:
                    modules.append(module)
                module = [domain]
                modules.append(module)
                module = []

            if domain in type_i_te_domains:
                end_cluster = True
            elif not type_i_te_domains:
                end_cluster = True

        # If a domain of known function is encountered already in the active module
        elif (domain.raichu_type in [d.raichu_type for d in module] and domain.raichu_type != "UNKNOWN") or \
                (domain.raichu_type in ["AT", "A", "CAL"] and ("AT" in [d.raichu_type for d in module] or
                                                               "A" in [d.raichu_type for d in module] or
                                                               "CAL" in [d.raichu_type for d in module])):
            if i != 0:
                previous_domain = domains[i - 1]
                # If the domain is on a new gene, always start a new module
                if previous_domain.gene != domain.gene:
                    modules.append(module)
                    module = [domain]

                # If the domain is on the same gene, start a new module if a carrier domain has already been encountered
                else:

                    # Double CP domains in a module CAN occur, so don't start a new module when the duplicated domain
                    # Is a CP domain on the same gene

                    if 'CP' in [d.raichu_type for d in module] and domain.raichu_type != 'CP':
                        if module:
                            modules.append(module)
                        module = [domain]
                    else:
                        module.append(domain)

            else:
                raise Exception("Module can't be empty if it contains a domain. Fix code.")

        elif domain.raichu_type == "UNKNOWN":
            if i != 0:
                previous_domain = domains[i - 1]
                if previous_domain.gene != domain.gene:
                    if module:
                        modules.append(module)

                    module = [domain]

                    if i != len(domains) - 1:
                        next_domain = domains[i + 1]
                        if next_domain.gene != domain.gene:
                            modules.append(module)
                            module = []

                else:
                    module.append(domain)

            else:
                raise Exception("Module can't be empty if it contains a domain. Fix code.")

        else:
            module.append(domain)

    if module:
        modules.append(module)

    antismash_modules = []

    for module in modules:
        start = min([domain.start for domain in module] + [domain.end for domain in module])
        end = max([domain.start for domain in module] + [domain.end for domain in module])
        strand = module[0].strand
        antismash_module = AntiSmashModule(start, end, strand, domains=module)
        antismash_modules.append(antismash_module)

    if antismash_modules:

        block_start = min([module.start for module in antismash_modules] + [module.end for module in antismash_modules])
        block_end = max([module.start for module in antismash_modules] + [module.end for module in antismash_modules])
        block_strand = antismash_modules[0].strand
        module_block = ModuleBlock(block_start, block_end, block_strand, antismash_modules)
    else:
        raise Exception("Cluster is empty.This can happen with Type III PKS clusters. Please check the input file.")

    return module_block


def filter_genes(genes):
    gene_groups = []
    gene_group = []

    for gene in genes:

        if gene.domains:
            gene_group.append(gene)
        else:
            if gene_group:
                gene_groups.append(gene_group)
                gene_group = []

    if gene_group:
        gene_groups.append(gene_group)

    return gene_groups


def sort_genes(genes):
    genes.sort(key=lambda gene: gene.start)
    gene_groups = []
    gene_group = []
    previous_direction = None

    for i, gene in enumerate(genes):
        current_direction = gene.strand
        if current_direction != previous_direction:
            if gene_group:
                if previous_direction == -1:
                    gene_group.reverse()
                gene_groups.append(gene_group)

            gene_group = [gene]
            previous_direction = current_direction

        else:
            gene_group.append(gene)

        if i == len(genes) - 1:
            if current_direction == -1:
                gene_group.reverse()
            gene_groups.append(gene_group)

    return gene_groups


def parse_antismash_domains_gbk(antismash_gbk, version="7.1.0"):
    domains = []

    for record in SeqIO.parse(antismash_gbk, "genbank"):
        for feature in record.features:
            if feature.type == "aSDomain":
                domain = AntiSmashDomain(feature.qualifiers["domain_id"][0],
                                         feature.qualifiers["aSDomain"][0],
                                         int(feature.location.start),
                                         int(feature.location.end),
                                         feature.location.strand,
                                         feature.qualifiers["locus_tag"][0])
                domain.add_domain_information(feature.qualifiers)
                domain.set_raichu_type()
                domains.append(domain)

    return domains


def load_antismash_gbk(gbk_file, version=7.0):
    module_blocks = get_nrps_pks_modules(gbk_file)
    cluster_representation = module_blocks.make_raichu_cluster()
    return cluster_representation


def parse_antismash_to_cluster_file(gbk_file, out_directory=None, version=7.0):
    if not out_directory:
        out_directory = os.path.splitext(gbk_file)[0]
    if not os.path.exists(out_directory):
        os.mkdir(out_directory)
    cluster = load_antismash_gbk(gbk_file)
    cluster.write_cluster(out_directory)


if __name__ == "__main__":
    draw_cluster(load_antismash_gbk(argv[1]), "test.svg")