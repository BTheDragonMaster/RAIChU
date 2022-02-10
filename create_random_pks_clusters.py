import random
from visualize_pks_cluster import *

starter_units = ['SC(=O)CC', 'SC(CC(O)=O)=O', 'SC(CC(O)=O)=O',
                 'SC(C(C(O)=O)CC)=O', 'SC(C(C(O)=O)OC)=O', 'SC(C*)=O',
                 'SC(C(C)CC)=O', 'SC(C1C(CCC1)C(=O)O)=O', 'SC(C)=O',
                 'SC(C1=CC=CC=C1)=O', 'SC(CC(C)C)=O', 'SC(C(C(=O)O)CC[Cl])=O']
tailoring_domain_combinations = [[], ['KR'], ['KR', 'DH'], ['KR', 'DH', 'ER']]
elongation_units = ['malonylcoa', 'methylmalonylcoa', 'pk', 'methoxymalonylacp', 'ethylmalonylcoa']

def generate_random_cluster():
    """
    
    """
    starter_module = ['module 1', 'starter_module', random.choice(starter_units)]
    cluster = [starter_module]
    last_elongation_domain_nr = random.randint(2, 8)
    for i in range(2, last_elongation_domain_nr):
        elongation_module = [f'module {i}', 'elongation_module', random.choice(elongation_units), random.choice(tailoring_domain_combinations)]
        cluster.append(elongation_module)
    terminator_module =  [f'module {random.randint(2, 8)+1}', 'terminator_module', random.choice(elongation_units), random.choice(tailoring_domain_combinations)]
    cluster.append(terminator_module)

    return cluster

if __name__ == "__main__":
    #Pick a number (how many random clusters you want to generate and draw)
    times = 100
    for i in range(times):
        generated_cluster = generate_random_cluster()
        print(generated_cluster)
        pks_cluster_to_structure(generated_cluster, attach_to_acp=True)
