import random
from raichu.visualize_cluster import *

starter_units = ['SC(=O)CC', 'SC(CC(O)=O)=O', 'SC(CC(O)=O)=O',
                 'SC(C(C(O)=O)CC)=O', 'SC(C(C(O)=O)OC)=O', 'SC(C*)=O',
                 'SC(C(C)CC)=O', 'SC(C1C(CCC1)C(=O)O)=O', 'SC(C)=O',
                 'SC(C1=CC=CC=C1)=O', 'SC(CC(C)C)=O', 'SC(C(C(=O)O)CC[Cl])=O']
tailoring_domain_combinations = [[], ['KR'], ['KR', 'DH'], ['KR', 'DH', 'ER']]
elongation_units = ['malonylcoa', 'methylmalonylcoa', 'pk',
                    'methoxymalonylacp', 'ethylmalonylcoa']
nrps_tailoring_domain_combinations = [[], ['E'], ['nMT'], ['E', 'nMT'], ['nMT', 'E']]

def generate_random_pks_cluster():
    """Generates and returns a pure PKS cluster in the standard RAIChU format,
    as [[module_1], [module_2], ...], comprising 2 to 8 modules.
    """
    starter_module = ['module 1', 'starter_module_pks',
                      random.choice(starter_units)]
    cluster = [starter_module]
    last_elongation_domain_nr = random.randint(2, 8)
    for i in range(2, last_elongation_domain_nr):
        elongation_module = [f'module {i}', 'elongation_module_pks',
                             random.choice(elongation_units),
                             random.choice(tailoring_domain_combinations)]
        cluster.append(elongation_module)
    terminator_module =  [f'module {last_elongation_domain_nr}',
                          'terminator_module_pks', random.choice(elongation_units),
                          random.choice(tailoring_domain_combinations)]
    cluster.append(terminator_module)

    return cluster

def generate_random_nrps_cluster():
    """Generates and returns a pure NRPS cluster in the standard RAIChU format,
    as [[module_1], [module_2], ...], comprising 2 to 8 modules.
    """
    list_aa = ['alanine', 'arginine', 'asparagine', 'asparticacid',
               'cysteine', 'glutamine', 'glutamicacid', 'glycine',
               'histidine', 'isoleucine', 'leucine', 'lysine',
               'methionine', 'phenylalanine', 'proline', 'serine',
               'threonine', 'tryptophan', 'tyrosine', 'valine', 'glycine']
    starter_module = ['module 1', 'starter_module_nrps',
                      random.choice(list_aa)]
    cluster = [starter_module]
    last_elongation_domain_nr = random.randint(2, 8)
    for i in range(2, last_elongation_domain_nr):
        elongation_module = [f'module {i}', 'elongation_module_nrps',
                             random.choice(list_aa),
                             random.choice(nrps_tailoring_domain_combinations)]
        cluster.append(elongation_module)
    terminator_module = [f'module {last_elongation_domain_nr}',
                         'terminator_module_nrps', random.choice(list_aa),
                         random.choice(nrps_tailoring_domain_combinations)]
    cluster.append(terminator_module)

    return cluster

def generate_random_hybrid_cluster():
    """Generates and returns a hybrid PKS/NRPS cluster in the standard RAIChU
    format, as [[module_1], [module_2], ...], comprising 2 to 8 modules.
    """
    list_aa = ['alanine', 'arginine', 'asparagine', 'asparticacid', 'cysteine',
               'glutamine', 'glutamicacid', 'glycine', 'histidine',
               'isoleucine', 'leucine', 'lysine', 'methionine',
               'phenylalanine', 'proline', 'serine', 'threonine', 'tryptophan',
               'tyrosine', 'valine', 'glycine']
    # Choose between NPRS or PKS starter module
    if random.randint(0,1) == 1:
        starter_module = ['module 1', 'starter_module_pks',
                          random.choice(starter_units)]
    else:
        starter_module = ['module 1', 'starter_module_nrps',
                          random.choice(list_aa)]
    cluster = [starter_module]
    last_elongation_domain_nr = random.randint(2, 8)
    # For each elongation module, choose between PKS or NRPS module
    for i in range(2, last_elongation_domain_nr):
        if random.randint(0,1) == 1:
            elongation_module = [f'module {i}', 'elongation_module_pks',
                                 random.choice(elongation_units),
                                 random.choice(tailoring_domain_combinations)]
        else:
            elongation_module = [f'module {i}', 'elongation_module_nrps',
                                 random.choice(list_aa),
                                 random.choice(nrps_tailoring_domain_combinations)]
        cluster.append(elongation_module)
    # Choose between PKS or NRPS terminator module
    if random.randint(0,1) == 1:
        terminator_module =  [f'module {last_elongation_domain_nr}',
                              'terminator_module_pks',
                              random.choice(elongation_units),
                              random.choice(tailoring_domain_combinations)]
    else:
        terminator_module = [f'module {last_elongation_domain_nr}',
                             'terminator_module_nrps', random.choice(list_aa),
                             random.choice(nrps_tailoring_domain_combinations)]
    cluster.append(terminator_module)

    return cluster

if __name__ == "__main__":
    # Draw and save 500 pure NRPS clusters
    # for i in range(1, 501):
    #     cluster = generate_random_nrps_cluster()
    #     print(cluster)
    #     draw_cluster(cluster, save_fig=f'NRPS_cluster_{i}.png')
    #     plt.close('all')
    #
    # # Draw and save 500 pure PKS clusters
    # for i in range(1, 501):
    #     cluster = generate_random_pks_cluster()
    #     print(cluster)
    #     draw_cluster(cluster, save_fig=f'PKS_cluster_{i}.png')
    #     plt.close('all')

    # Draw and save 500 hybrid PKS/NRPS clusters
    for i in range(1, 501):
        cluster = generate_random_hybrid_cluster()
        print(cluster)
        draw_cluster(cluster, save_fig=f'hybrid_PKS_NRPS_cluster_{i}.png')
        plt.close('all')


