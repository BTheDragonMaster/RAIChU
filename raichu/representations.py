from dataclasses import dataclass
from typing import List, Union
import os


@dataclass
class MacrocyclizationRepresentation:
    atom_1: str
    atom_2: str
    type: str = "oxidative"


@dataclass
class CleavageSiteRepresentation:
    amino_acid: str
    amino_acid_index: int
    peptide_fragment_to_keep: str


@dataclass
class TailoringRepresentation:
    gene_name: str
    type: str
    # Some tailoring reactions involve more than one atom
    modification_sites: List[List[str]]
    substrate: Union[str, None] = None


@dataclass
class IsomerizationRepresentation:
    modification_sites: List[List[str]]


@dataclass
class MethylShiftRepresentation:
    modification_sites: List[List[str]]


@dataclass
class DomainRepresentation:
    gene_name: Union[str, None]
    type: str
    subtype: Union[str, None] = None
    name: Union[str, None] = None
    active: bool = True
    used: bool = True


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

    def write_cluster(self, out_dir):
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)
        cluster_out = os.path.join(out_dir, "cluster.txt")
        with open(cluster_out, "w") as out:
            out.write(
                "gene_name\tmodule_nr\tmodule_type\tmodule_subtype\tmodule_substrate\tnr_iterations\tdomain_type\tdomain_subtype\tdomain_name\tdomain_used\tdomain_active\n"
            )
            for i, module in enumerate(self.modules):
                for domain in module.domains:
                    out.write(
                        f"{domain.gene_name}\t{i}\t{module.type}\t{module.subtype}\t{module.substrate}\t{module.iterations}\t{domain.type}\t{domain.subtype}\t{domain.name}\t{domain.used}\t{domain.active}\n"
                    )
        if self.tailoring_enzymes:
            tailoring_out = os.path.join(out_dir, "tailoring.txt")
            with open(tailoring_out, "w") as tailoring:
                tailoring.write("gene_name\ttype\tsubstrate\tmodification_sites\n")

                for enzyme in self.tailoring_enzymes:
                    if enzyme.modification_sites:
                        site_reprs = []
                        for modification_site in enzyme.modification_sites:
                            site_repr = "|".join(list(map(str, modification_site)))
                            site_reprs.append(site_repr)

                        site_str = ":".join(site_reprs)
                    else:
                        site_str = str(None)

                    tailoring.write(
                        f"{enzyme.gene_name}\t{enzyme.type}\t{enzyme.substrate}\t{site_str}\n"
                    )

    @classmethod
    def from_file(cls, in_dir):
        in_cluster = os.path.join(in_dir, "cluster.txt")
        in_tailoring = os.path.join(in_dir, "tailoring.txt")
        module_representations = []
        modules = {}
        with open(in_cluster, "r") as cluster:
            cluster.readline()
            for line in cluster:
                line = line.strip()
                if line:
                    line_info = line.split("\t")
                    line_info_cleaned = []
                    for i, entry in enumerate(line_info):
                        if line_info[i] == "False":
                            line_info_cleaned.append(False)
                        elif line_info[i] == "True":
                            line_info_cleaned.append(True)
                        elif line_info[i] == "None":
                            line_info_cleaned.append(None)
                        else:
                            line_info_cleaned.append(line_info[i])
                    (
                        gene_name,
                        module_nr,
                        module_type,
                        module_subtype,
                        substrate,
                        iterations,
                        domain_type,
                        domain_subtype,
                        domain_name,
                        domain_used,
                        domain_active,
                    ) = line_info_cleaned
                    module_nr = int(module_nr)
                    domain = DomainRepresentation(
                        gene_name,
                        domain_type,
                        domain_subtype,
                        domain_name,
                        domain_active,
                        domain_used,
                    )
                    if module_nr not in modules:
                        modules[module_nr] = {
                            "type": module_type,
                            "subtype": module_subtype,
                            "substrate": substrate,
                            "iterations": iterations,
                            "domains": [domain],
                        }
                    else:
                        modules[module_nr]["domains"].append(domain)

            for module_nr, module_data in modules.items():
                module = ModuleRepresentation(
                    module_data["type"],
                    module_data["subtype"],
                    module_data["substrate"],
                    module_data["domains"],
                    module_data["iterations"],
                )
                module_representations.append(module)

        if os.path.exists(in_tailoring):
            tailoring_enzymes = []
            with open(in_tailoring, "r") as tailoring:
                tailoring.readline()
                for line in tailoring:
                    line_info = line.split("\t")
                    line_info_cleaned = []
                    for i, entry in enumerate(line_info):
                        if line_info[i] == "None":
                            line_info_cleaned.append(None)
                        else:
                            line_info_cleaned.append(line_info[i])

                    gene_name, enzyme_type, substrate, site_str = line_info_cleaned

                    sites = []

                    if site_str:
                        site = site_str.split(":")
                        sites.append(site.split("|"))
                    tailoring_representation = TailoringRepresentation(
                        gene_name, enzyme_type, sites, substrate
                    )
                    tailoring_enzymes.append(tailoring_representation)

        else:
            tailoring_enzymes = None

        return cls(module_representations, tailoring_enzymes)
