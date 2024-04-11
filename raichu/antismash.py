import json
import os
import re
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
    "bmt": "4-butenyl-4-methylthreonine",
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
                    self.subtype = spec.split("transATor:")[1].strip().replace("-", "_").replace("/", "_").upper().replace("(UNKNOWN)", "MISCELLANEOUS")

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

    def determine_completeness(self):
        self.complete = False
        self.starter_complete = False

        if self.type == 'PKS':

            has_ks = False
            has_acp = False
            has_at = False

            for domain in self.domains:
                if domain.raichu_type == 'KS':
                    has_ks = True
                if domain.raichu_type == 'CP':
                    has_acp = True
                if domain.raichu_type == 'AT':
                    has_at = True

            if has_acp:
                if self.subtype == 'PKS_TRANS':
                    self.starter_complete = True
                elif self.subtype == 'PKS_CIS' and has_at:
                    self.starter_complete = True

                if has_ks:
                    if self.subtype == 'PKS_TRANS':
                        self.complete = True
                    elif self.subtype == 'PKS_CIS' and has_at:
                        self.complete = True

        elif self.type == 'NRPS':
            has_c = False
            has_a = False
            has_pcp = False

            for domain in self.domains:
                if domain.raichu_type == 'C':
                    has_c = True
                if domain.raichu_type == 'CP':
                    has_pcp = True
                if domain.raichu_type == 'A':
                    has_a = True

            if has_a and has_pcp:
                self.starter_complete = True

                if has_c:
                    self.complete = True

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


class ModuleBlocks:

    def __init__(self, module_blocks: list[ModuleBlock]):
        self.module_blocks = module_blocks
        self.block_order = module_blocks
        self.optimal_solution = False
        self.modules = []

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
                        if self.check_broken_modules() < best_solution:
                            best_block_order = self.block_order[:]
                        if self.check_broken_modules() == 0:
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
                    if self.check_broken_modules() < best_solution:
                        best_block_order = self.block_order[:]
                        best_solution = self.check_broken_modules()
                    if self.check_broken_modules() == 0:
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
                    if self.check_broken_modules() < best_solution:
                        best_block_order = self.block_order[:]
                        best_solution = self.check_broken_modules()
                    if self.check_broken_modules() == 0:
                        found_optimal_solution = True
                        break
                if found_optimal_solution:
                    break

        self.optimal_solution = found_optimal_solution
        self.block_order = best_block_order
        self.move_broken_modules()

    def move_broken_modules(self):
        module_block = self.recreate_modules()
        print(module_block.modules)
        standalone_gene = False
        for i, module_1 in enumerate(module_block.modules):

            for j, module_2 in enumerate(module_block.modules):
                continue

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
            if not module.complete:
                if i != 0:
                    nr_broken += 1
                elif not module.starter_complete:
                    nr_broken += 1
            if module.has_termination_domain():
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


def get_nrps_pks_genes(antismash_gbk):
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

    modules_per_collinear_block = ModuleBlocks(make_modules_from_gene_groups(filtered_groups))
    modules_per_collinear_block.sort_collinear_blocks()
    print(modules_per_collinear_block.optimal_solution)
    print(modules_per_collinear_block.block_order)

    if not modules_per_collinear_block.optimal_solution:
        filtered_genes = []

        for gene_group in filtered_groups:
            for gene in gene_group:
                filtered_genes.append([gene])

        modules_per_gene = ModuleBlocks(make_modules_from_gene_groups(filtered_genes))
        modules_per_gene.sort_collinear_blocks()
        print(modules_per_gene.optimal_solution)
        print(modules_per_gene.block_order)

    return filtered_groups


def make_modules_from_gene_groups(gene_groups):
    module_groups = []
    for gene_group in gene_groups:
        domains = []

        for gene in gene_group:
            domains += gene.domains

        module_block = make_modules(domains)

        module_groups.append(module_block)

    return module_groups


def make_modules(domains):

    modules = []
    module = []
    end_cluster = False

    for i, domain in enumerate(domains):
        if not module and not modules:
            module.append(domain)

        elif domain.raichu_type in ["C", "KS"]:
            if module:
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
                    module.append(domain)
                    modules.append(module)
                    module = []

                # Dealing with a standalone TE domain
                else:
                    if module:
                        modules.append(module)

                    module = [domain]
                    modules.append(module)
                    module = []

            else:
                if module:
                    modules.append(module)
                module = [domain]
                modules.append(module)
                module = []

            end_cluster = True

        elif domain.raichu_type in [d.raichu_type for d in module] and domain.raichu_type != "UNKNOWN":
            if i != 0:
                previous_domain = domains[i - 1]
                if previous_domain.gene != domain.gene:
                    modules.append(module)
                    module = [domain]
                else:

                    if 'CP' in [d.raichu_type for d in module]:
                        if module:
                            modules.append(module)
                        module = [domain]
                    else:
                        module.append(domain)

            else:
                raise Exception("Well done biology, you have managed to access an impossible location in our code.")

        elif domain.raichu_type == "UNKNOWN":
            if i != 0:
                previous_domain = domains[i - 1]
                if previous_domain.gene != domain.gene:
                    if module:
                        modules.append(module)
                    module = [domain]
                else:
                    module.append(domain)
            else:
                raise Exception("Well done biology, you have managed to access an impossible location in our code.")

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

    block_start = min([module.start for module in antismash_modules] + [module.end for module in antismash_modules])
    block_end = max([module.start for module in antismash_modules] + [module.end for module in antismash_modules])
    block_strand = antismash_modules[0].strand
    module_block = ModuleBlock(block_start, block_end, block_strand, antismash_modules)

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


def map_domains_to_modules_gbk(antismash_gbk, domains):
    modules = []
    strands = []
    gene_to_coords = {}
    gene_to_strand = {}
    gene_to_modules = {}
    for record in SeqIO.parse(antismash_gbk, "genbank"):

        nr_all_modules = len([feature for feature in record.features if feature.type == "aSModule"])
        nr_unknown_modules = len([feature for feature in record.features if
                                  feature.type == "aSModule" and "unknown" in feature.qualifiers["type"]])

        nr_nrps_pks_modules = nr_all_modules - nr_unknown_modules

        for feature in record.features:
            if feature.type == "aSModule":
                module_type = feature.qualifiers["type"][0].upper()
                module_subtype = None

                if module_type == "NRPS" or module_type == "PKS":

                    start = feature.location.start
                    end = feature.location.end
                    gene = feature.qualifiers["locus_tags"][0]
                    if gene not in gene_to_coords:
                        gene_to_coords[gene] = start
                    else:
                        if gene_to_coords[gene] > start:
                            gene_to_coords[gene] = start

                    domains_in_module = [
                        domain
                        for domain in domains
                        if domain["start"] >= start and domain["end"] <= end
                    ]
                    domain_representations = [
                        domain["representation"] for domain in domains_in_module
                    ]

                    strand = domains_in_module[0]["strand"]
                    for domain in domains_in_module:
                        if domain["gene"] not in gene_to_strand:
                            gene_to_strand[domain["gene"]] = domain["strand"]

                    strands.append(strand)
                    substrates = [domain["substrate"] for domain in domains_in_module if domain["substrate"] is not None]
                    if len(substrates) > 0:
                        substrate = substrates[0]

                        # Also account for reversing cluster if on different strand
                        if module_type == "PKS" and ((len(modules) == 0 and strand == 1) or (len(modules) == nr_nrps_pks_modules - 1 and strand == -1)):
                            if substrate not in [v.name for v in PksStarterSubstrate]:

                                substrate = "WILDCARD"
                    else:
                        if module_type == "NRPS":
                            substrate = "**Unknown**"
                        elif module_type == "PKS":
                            substrate = "WILDCARD"
                        else:
                            raise ValueError("Module type can only be NRPS or PKS")

                    if module_type == "PKS":
                        module_subtype = "PKS_TRANS"

                    for domain_representation in domain_representations:
                        if domain_representation.type == "CP":
                            domain_representation.type = (
                                "ACP" if module_type == "PKS" else "PCP"
                            )
                        if module_type == "PKS":
                            if domain_representation.type == "AT" or domain_representation.type == 'CAL':
                                module_subtype = "PKS_CIS"
                            if domain_representation.type == "CAL" and substrate == "WILDCARD":
                                substrate = "**Unknown**"
                    if (
                        module_type == "PKS"
                        and module_subtype == "PKS_TRANS"
                        and not any(
                            [
                                (
                                    domain_representation.subtype == "PKS_TRANS"
                                    or domain_representation.type == "Trans-AT_docking"
                                )
                                for domain_representation in [
                                    domain["representation"] for domain in domains
                                ]
                            ]
                        )
                    ):
                        module_subtype = "PKS_CIS"
                    if strand == -1:
                        domain_representations.reverse()
                    # Check if CP in every module

                    if not any(
                        [
                            domain_representation.type in ["CP", "ACP", "PCP"]
                            for domain_representation in domain_representations
                        ]
                    ):
                        # Check if ACP/PCP exists elsewhere
                        if strand == -1:
                            next_domain_index = domains.index(domains_in_module[0]) - 1
                        else:
                            next_domain_index = domains.index(domains_in_module[-1]) + 1
                        if next_domain_index < len(domains):
                            next_domain = domains[next_domain_index]
                            if next_domain["representation"].type in [
                                "ACP",
                                "PCP",
                                "CP",
                            ]:
                                domain_representations.append(
                                    DomainRepresentation(
                                        feature.qualifiers["locus_tags"][0],
                                        "PCP" if module_type == "NRPS" else "ACP",
                                        subtype=None,
                                        name=None,
                                        active=True,
                                        used=True,
                                    )
                                )
                    module_representation = ModuleRepresentation(
                        module_type, module_subtype, substrate, domain_representations
                    )
                    modules.append(module_representation)

    if all(
        strand == -1 for strand in strands
    ):  # reverse whole cluster, if completely on -1 strand
        modules.reverse()
    gene_groups = sort_genes(gene_to_coords, gene_to_strand)

    cluster_representation = ClusterRepresentation(modules)

    return cluster_representation


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
    domains = parse_antismash_domains_gbk(gbk_file, version)
    cluster_representation = map_domains_to_modules_gbk(gbk_file, domains)
    return cluster_representation


def parse_antismash_to_cluster_file(gbk_file, out_directory=None, version=7.0):
    if not out_directory:
        out_directory = os.path.splitext(gbk_file)[0]
    if not os.path.exists(out_directory):
        os.mkdir(out_directory)
    cluster = load_antismash_gbk(gbk_file)
    cluster.write_cluster(out_directory)


if __name__ == "__main__":
    get_nrps_pks_genes(argv[1])
    # draw_cluster(
    #     load_antismash_gbk(argv[1]),
    #     out_file=argv[2],
    # )
