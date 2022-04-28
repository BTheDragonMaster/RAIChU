import os
from raichu.visualize_cluster import *
from raichu.modules_to_structure import cluster_to_structure
from raichu.thioesterase_reactions import thioesterase_all_products
matplotlib.use("Agg")
import matplotlib.backends.backend_agg as agg


def get_raichu_cluster(genes):
    cluster = []

    module_nr = 0

    if genes:

        for i, gene in enumerate(genes):
            last_gene = False

            if gene.name:
                file_name = f"{gene}.png"
            else:
                file_name = f"gene_{i}.png"

            for module in gene.modules:
                substrate_specificity = None
                module_nr += 1

                module_name = f"module {module_nr}"

                if module.type == 'NRPS' and module_nr == 1:
                    module_type = 'starter_module_nrps'

                elif module.type == 'NRPS':
                    module_type = 'elongation_module_nrps'
                elif module.type == 'PKS' and module_nr == 1:
                    module_type = 'starter_module_pks'
                    substrate_specificity = 'SC(=O)CC'
                    for domain in module.domains:
                        if domain.type == 'AT':
                            if domain.substrate:
                                substrate_specificity = domain.substrate.smiles

                else:
                    module_type = 'elongation_module_pks'

                    substrate_specificity = 'pk'
                    for domain in module.domains:
                        if domain.type == 'AT':
                            if domain.substrate:
                                substrate_specificity = domain.substrate.name

                    if substrate_specificity == 'methoxymalonylcoa':
                        substrate_specificity = 'methoxymalonylacp'
                    if substrate_specificity == 'wildcard':
                        substrate_specificity = 'pk'

                if module.type == 'NRPS':
                    substrate_specificity = 'nrp'
                    for domain in module.domains:
                        if domain.type == 'A':
                            if domain.substrate:
                                substrate_specificity = domain.substrate.name

                assert substrate_specificity

                last_domain = False

                tailoring_domains = []
                for domain in module.domains:
                    if domain.type not in ['A', 'C', 'AT', 'KR', 'KS', 'ACP', 'PCP', 'TE']:
                        tailoring_domains.append(domain.type)
                    elif domain.type == 'KR':
                        if domain.subtype:
                            tailoring_domains = [f'KR_{domain.subtype}'] + tailoring_domains
                        else:
                            tailoring_domains = ["KR"] + tailoring_domains
                    elif domain.type == 'TE':
                        last_domain = True
                        if module.type == 'NRPS':
                            module_type = 'terminator_module_nrps'
                        elif module.type == 'PKS':

                            module_type = 'terminator_module_pks'

                if module_type.startswith('starter'):
                    cluster.append([module_name, module_type, substrate_specificity])
                else:
                    cluster.append([module_name, module_type, substrate_specificity, tailoring_domains])

                if last_domain:
                    last_gene = True
                    break

            if last_gene:
                break

    return cluster


def render_products(genes, dir_name="products_test"):

    cluster = get_raichu_cluster(genes)
    if cluster:
        if not os.path.exists(dir_name):
            os.mkdir(dir_name)
        attached_product = cluster_to_structure(cluster, attach_to_acp=True)
        thioesterase_all_products(attached_product, out_folder=dir_name)


def render_cluster(genes, file_name="cluster_test.png"):
    cluster = get_raichu_cluster(genes)
    if cluster:
        draw_cluster(cluster, save_fig=file_name)

