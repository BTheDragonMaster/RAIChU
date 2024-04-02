import json
import os
import re
from sys import argv
from Bio import SeqIO
from raichu.representations import (
    ClusterRepresentation,
    ModuleRepresentation,
    DomainRepresentation,
)
from raichu.run_raichu import draw_cluster, build_cluster
from raichu.substrate import PksStarterSubstrate


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
}

domains_to_be_refined = ["KR", "KS", "A", "AT"]

AS_TO_PKS = {
    "mal": "MALONYL_COA",
    "mmal": "METHYLMALONYL_COA",
    "mxmal": "METHOXYMALONYL_COA",
    "emal": "ETHYLMALONYL_COA",
}


def map_domains_to_modules_gbk(antismash_gbk, domains):
    modules = []
    strands = []
    for record in SeqIO.parse(antismash_gbk, "genbank"):
        for feature in record.features:
            if feature.type == "aSModule":
                module_type = feature.qualifiers["type"][0].upper()
                module_subtype = None
                if module_type == "NRPS" or module_type == "PKS":
                    start = feature.location.start
                    end = feature.location.end
                    domains_in_module = [
                        domain
                        for domain in domains
                        if domain["start"] >= start and domain["end"] <= end
                    ]
                    domain_representations = [
                        domain["representation"] for domain in domains_in_module
                    ]
                    strand = domains_in_module[0]["strand"]
                    strands.append(strand)
                    substrate = [
                        domain["substrate"]
                        for domain in domains_in_module
                        if domain["substrate"] is not None
                    ]
                    if len(substrate) > 0:
                        substrate = substrate[0]
                        # Also account for reversing cluster if on different strand
                        if module_type == "PKS" and (
                            len(modules) == 0
                            and strand == 1
                            or len(modules)
                            == len(
                                [
                                    feature
                                    for feature in record.features
                                    if feature.type == "aSModule"
                                ]
                            )
                            - 1
                            and strand == -1
                        ):
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
                            if domain_representation.type == "AT":
                                module_subtype = "PKS_CIS"
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
                        # Check if ACP/PCP is actually existing
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
    cluster_representation = ClusterRepresentation(modules)
    return cluster_representation


def parse_antismash_domains_gbk(antismash_gbk, version=7.0):
    domains = []

    for record in SeqIO.parse(antismash_gbk, "genbank"):
        for feature in record.features:
            if feature.type == "aSDomain":
                domain = {}
                domain["active"] = True
                domain["id"] = feature.qualifiers["domain_id"][0]
                domain["type"] = (
                    antiSMASH_DOMAIN_TO_RAICHU_DOMAIN[feature.qualifiers["aSDomain"][0]]
                    if feature.qualifiers["aSDomain"][0]
                    in antiSMASH_DOMAIN_TO_RAICHU_DOMAIN
                    else feature.qualifiers["aSDomain"][0]
                )
                domain["start"] = feature.location.start
                domain["end"] = feature.location.end
                domain["subtype"] = None
                domain["substrate"] = None
                domain["strand"] = feature.location.strand
                if "specificity" in feature.qualifiers:
                    for spec in feature.qualifiers["specificity"]:
                        if "consensus:" in spec:
                            specificity = spec.split("consensus:")[1].strip().lower()
                            if specificity in AS_TO_NRPS:
                                substrate_name = AS_TO_NRPS[specificity]
                                domain["substrate"] = substrate_name
                            elif specificity in AS_TO_PKS:
                                substrate_name = AS_TO_PKS[specificity]
                                domain["substrate"] = substrate_name
                        elif "KR activity:" in spec:
                            domain["active"] = (
                                True
                                if spec.split("KR activity:")[1].strip() == "active"
                                else False
                            )

                        elif "transATor:" in spec:
                            domain["subtype"] = (
                                spec.split("transATor:")[1]
                                .strip()
                                .replace("-", "_")
                                .replace("/", "_")
                                .upper()
                                .replace("(UNKNOWN)", "MISCELLANEOUS")
                            )

                domain["representation"] = DomainRepresentation(
                    feature.qualifiers["locus_tag"][0],
                    domain["type"],
                    domain["subtype"],
                    None,
                    domain["active"],
                    domain["active"],
                )
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


def refine_domain_js(start, gene, details_data_region):
    details_data_region_orf = [
        orf for orf in details_data_region["orfs"] if orf["id"] == gene
    ][0]
    details_data_region_orf_domain = [
        domain
        for domain in details_data_region_orf["domains"]
        if domain["start"] == start
    ][0]
    substrate = None
    KS_subtype = None
    KR_subtype = None
    if details_data_region_orf_domain["abbreviation"] == "AT":
        substrate = details_data_region_orf_domain["predictions"][1][1]
        return [None, substrate.replace("-", "_").upper()]

    if details_data_region_orf_domain["abbreviation"] == "A":
        substrate = AS_TO_NRPS[
            details_data_region_orf_domain["predictions"][0][1].lower()
        ]
        substrate = substrate.lower() if substrate != "**Unknown**" else substrate
        return [None, substrate]

    if details_data_region_orf_domain["abbreviation"] == "KS":
        if len(details_data_region_orf_domain["predictions"]) > 0:
            KS_subtype = details_data_region_orf_domain["predictions"][0][1]
            return [KS_subtype, None]
        else:
            return [None, None]

    if details_data_region_orf_domain["abbreviation"] == "KR":
        if len(details_data_region_orf_domain["predictions"]) > 0:
            KR_subtype = details_data_region_orf_domain["predictions"][1][1]
            if KR_subtype == "(unknown)":
                KR_subtype = None
            return [KR_subtype, None]


def get_cluster_representation_js(
    region_visualizer, details_data_region, version=7.0
) -> ClusterRepresentation:
    module_array = []
    # only get first entry
    modules = region_visualizer[next(iter(region_visualizer))]["modules"]

    for module in modules:
        domain_array = []
        module_subtype = None
        domain_substrate = None
        substrate = "X"
        if module["complete"] == True:
            for domain in module["domains"]:
                name = None
                active = True if domain["inactive"] == False else False
                used = active
                domain_type = domain["name"]
                if domain_type == "":
                    domain_type = domain["description"]
                description = domain["description"]
                gene = domain["cds"]
                start = domain["start"]
                subtype = None
                if domain_type in domains_to_be_refined:
                    subtype, domain_substrate = refine_domain_js(
                        start, gene, details_data_region
                    )
                    if domain_substrate:
                        substrate = domain_substrate
                domain_array.append(
                    DomainRepresentation(gene, domain_type, subtype, name, active, used)
                )
            if "A" in [domain.type for domain in domain_array]:
                module_type = "NRPS"
            else:
                module_type = "PKS"
            if module_type == "PKS":
                if "AT" in [domain.type for domain in domain_array]:
                    module_subtype = "PKS_CIS"
                else:
                    module_subtype = "PKS_TRANS"
            for domain in domain_array:
                if domain.type == "CP":
                    domain.type = "ACP" if module_type == "PKS" else "PCP"
            module_array.append(
                ModuleRepresentation(
                    module_type, module_subtype, substrate, domain_array
                )
            )
    return ClusterRepresentation(module_array)


def load_antismash_js(js_file, region, version=7.0):
    with open(js_file) as file:
        complete_data = file.read()
    pattern_results_data = re.compile(r"var\s+resultsData\s*=\s*([\s\S]*?)(?=$)")
    results_data = pattern_results_data.search(complete_data).group(1)
    pattern_details_data = re.compile(r"var\s+details_data\s*=\s*([\s\S]*?)(?=var|$)")
    details_data = pattern_details_data.search(complete_data).group(1)
    if not details_data and results_data:
        raise ValueError("Not correct antiSMASH input.")
    details_data = json.loads(f'{{"results":{details_data[:-2]}}}')["results"]
    results_data = json.loads(results_data[:-1])
    if "nrpspks" in details_data:
        print(details_data["nrpspks"].keys())
        if region in details_data["nrpspks"]:
            details_data_region = details_data["nrpspks"][region]
        else:
            raise ValueError("Cluster is not a modular cluster.")
    else:
        raise ValueError("Cluster is not a modular cluster.")

    region_visualizer = results_data[region][
        "antismash.outputs.html.visualisers.bubble_view"
    ]
    raw_cluster_representation = get_cluster_representation_js(
        region_visualizer, details_data_region, version
    )

    return raw_cluster_representation


if __name__ == "__main__":
    draw_cluster(
        load_antismash_gbk(argv[1]),
        out_file=argv[2],
    )
