import random
from visualize_pks_cluster import *
from pikachu.smiles.smiles import *
from pikachu.chem.structure import *

starter_units = ['SC(=O)CC', 'SC(CC(O)=O)=O', 'SC(CC(O)=O)=O',
                 'SC(C(C(O)=O)CC)=O', 'SC(C(C(O)=O)OC)=O', 'SC(C*)=O',
                 'SC(C(C)CC)=O', 'SC(C1C(CCC1)C(=O)O)=O', 'SC(C)=O',
                 'SC(C1=CC=CC=C1)=O', 'SC(CC(C)C)=O', 'SC(C(C(=O)O)CC[Cl])=O']
tailoring_domain_combinations = [[], ['KR'], ['KR', 'DH'], ['KR', 'DH', 'ER']]
elongation_units = ['malonylcoa', 'methylmalonylcoa', 'pk', 'methoxymalonylacp', 'ethylmalonylcoa']

def generate_random_pks_cluster():
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

def generate_random_nrps_cluster():
    """

    """
    #list_aa = ['2-AMINOADIPICACID', '2-AMINOBUTYRICACID', 'ANTICAPSIN', '2-AMINOISOBUTYRICACID', 'ALLO-ISOLEUCINE', 'ALANINE', 'ALPHA-METHYLCYSTEINE', 'ALPHA-AMINO-PHENYL-VALERICACID', 'ARGININE', 'ASPARAGINE', 'ASPARTICACID', 'ALLO-THREONINE', 'CAPREOMYCIDINE', 'CITRULLINE', 'CORONAMICACID', 'CYSTEINE', 'CYSTEICACID', 'D-2-AMINOADIPICACID', 'D-2-AMINOBUTYRICACID', 'D-ALLO-ISOLEUCINE', 'D-ALANINE', 'D-ALPHA-METHYLCYSTEINE', 'D-ALPHA-AMINO-PHENYL-VALERICACID', 'D-ARGININE', 'D-ASPARAGINE', 'D-ASPARTICACID', 'D-ALLO-THREONINE', 'D-CITRULLINE', 'D-CYSTEINE', 'D-CYSTEICACID', 'D-2,4-DIAMINOBUTYRICACID', '(2R,3R)-2,3-DIAMINOBUTYRICACID', 'D-3,5-DIHYDROXYPHENYLGLYCINE', 'D-2,3-DIAMINOPROPIONICACID', 'D-ENDURACIDIDINE', 'D-ALPHA-ETHYLNORVALINE', 'D-GLUTAMINE', 'D-GLUTAMICACID', 'D-HYDROXYPHENYLGLYCINE', 'D-HOMOPROLINE', 'D-HOMOSERINE', 'D-HOMOSERINELACTONE', 'D-HOMOTYROSINE', 'D-ISOLEUCINE', 'D-ISOVALINE', 'D-KYNURENINE', 'D-LEUCINE', 'D-LYSINE', 'D-METHYL-2-AMINOOCTANOICACID', 'D-ORNITHINE', 'D-PHENYLGLYCINE', 'D-PHENYLALANINE', 'D-PROLINE', 'D-PHOSPHINOTHRICIN', 'D-SERINE', 'D-TERT-LEU', 'D-THREONINE', 'D-TRYPTOPHAN', 'D-TYROSINE', 'D-VALINE', '2,4-DIAMINOBUTYRICACID', '(2S)-2,3-DIAMINOBUTYRICACID', 'DEHYDROALANINE', '2,3-DEHYDRO-TRYPTOPHAN', 'DEHYDRO-CYSTEINE', '3,5-DIHYDROXYPHENYLGLYCINE', 'D-DIHYDROXYPHENYLTHIAZOLGROUP', 'DIHYDROXYPHENYLTHIAZOLGROUP', '2,3-DIAMINOPROPIONICACID', 'ENDURACIDIDINE', 'ALPHA-ETHYLNORVALINE', 'GLUTAMINE', 'GLUTAMICACID', 'GLYCINE', 'HOMOARGININE', 'HOMOISOLEUCINE', 'HISTIDINE', 'HYDROXYPHENYLGLYCINE', 'HOMOPHENYLALANINE', 'HOMOPROLINE', 'HOMOSERINE', 'HOMOSERINELACTONE', 'HOMOTYROSINE', 'ISOLEUCINE', 'ISOVALINE', 'KYNURENINE', '(2S,3S)-2,3-DIAMINOBUTYRICACID', 'LEUCINE', 'LYSINE', 'METHYL-2-AMINOOCTANOICACID', 'METHIONINE', 'NORCORONAMICACID', 'NORVALINE', 'ORNITHINE', 'PHENYLGLYCINE', 'PHENYLALANINE', 'PROLINE', 'PHOSPHINOTHRICIN', 'R-CORONAMICACID', 'R-NORCORONAMICACID', 'SERINE', 'TERT-LEU', 'THREONINE', 'TRYPTOPHAN', 'TYROSINE', 'VALINE', 'DEHYDROBUTYRINE', 'BETA-LYSINE', '(2S,3R)-2-AMINO-3-HYDROXY-4-(4-NITROPHENYL)BUTANOATE', '2,3-DIAMINO-3-METHYL-PROPANOICACID', '4,5-DIHYDROXYORNITHINE', '4,5-DIHYDROXYHOMOTYROSINE', '3-HYDROXYGLUTAMINE', '4-BUTENYL-4-METHYLTHREONINE', '4-HYDROXYVALINE', '2-METHYLSERINE', '4,5-DEHYDROARGININE', '3-(3-PYRIDYL)-ALANINE', '3-METHYLTYROSINE', 'D-LYSERGICACID', 'N-ACETYLPHENYLALANINE', '3-HYDROXY-4-O-METHYL-5-METHYL-TYROSINE', '4-METHYLPROLINE', 'ISOLEUCINE', '3-[(1R,2R)-2-NITROCYCLOPROPYL]-L-ALANINE', 'DEHYDROBUTYRINE', '(R)BETA-TYROSINE', 'PIPECOLICACID', '2-AMINO-9,10-EPOXY-8-OXODECANOICACID', 'L-PIPERAZICACID', 'DEHYDROTYROSINE', 'BETA-HYDROXYCYCLOHEX‐2‐ENYLALANINE', '2-METHYL-BETA-ALANINE', 'GLYCOCYAMINE']
    list_aa = ['alanine', 'arginine', 'asparagine', 'asparticacid', 'cysteine', 'glutamine', 'glutamicacid', 'glycine', 'histidine', 'isoleucine', 'leucine', 'lysine', 'methionine', 'phenylalanine', 'proline', 'serine', 'threonine', 'tryptophan', 'tyrosine', 'valine', 'glycine']
    starter_module = ['module 1', 'starter_module_nrps', random.choice(list_aa)]
    cluster = [starter_module]
    last_elongation_domain_nr = random.randint(2, 8)
    for i in range(2, last_elongation_domain_nr):
        elongation_module = [f'module {i}', 'elongation_module_nrps', random.choice(list_aa)]
        cluster.append(elongation_module)
    terminator_module =  [f'module {random.randint(2, 8)+1}', 'terminator_module_nrps', random.choice(list_aa)]
    cluster.append(terminator_module)

    return cluster

def generate_random_hybrid_cluster():
    """
    
    """
    #list_aa = ['2-AMINOADIPICACID', '2-AMINOBUTYRICACID', 'ANTICAPSIN', '2-AMINOISOBUTYRICACID', 'ALLO-ISOLEUCINE', 'ALANINE', 'ALPHA-METHYLCYSTEINE', 'ALPHA-AMINO-PHENYL-VALERICACID', 'ARGININE', 'ASPARAGINE', 'ASPARTICACID', 'ALLO-THREONINE', 'CAPREOMYCIDINE', 'CITRULLINE', 'CORONAMICACID', 'CYSTEINE', 'CYSTEICACID', 'D-2-AMINOADIPICACID', 'D-2-AMINOBUTYRICACID', 'D-ALLO-ISOLEUCINE', 'D-ALANINE', 'D-ALPHA-METHYLCYSTEINE', 'D-ALPHA-AMINO-PHENYL-VALERICACID', 'D-ARGININE', 'D-ASPARAGINE', 'D-ASPARTICACID', 'D-ALLO-THREONINE', 'D-CITRULLINE', 'D-CYSTEINE', 'D-CYSTEICACID', 'D-2,4-DIAMINOBUTYRICACID', '(2R,3R)-2,3-DIAMINOBUTYRICACID', 'D-3,5-DIHYDROXYPHENYLGLYCINE', 'D-2,3-DIAMINOPROPIONICACID', 'D-ENDURACIDIDINE', 'D-ALPHA-ETHYLNORVALINE', 'D-GLUTAMINE', 'D-GLUTAMICACID', 'D-HYDROXYPHENYLGLYCINE', 'D-HOMOPROLINE', 'D-HOMOSERINE', 'D-HOMOSERINELACTONE', 'D-HOMOTYROSINE', 'D-ISOLEUCINE', 'D-ISOVALINE', 'D-KYNURENINE', 'D-LEUCINE', 'D-LYSINE', 'D-METHYL-2-AMINOOCTANOICACID', 'D-ORNITHINE', 'D-PHENYLGLYCINE', 'D-PHENYLALANINE', 'D-PROLINE', 'D-PHOSPHINOTHRICIN', 'D-SERINE', 'D-TERT-LEU', 'D-THREONINE', 'D-TRYPTOPHAN', 'D-TYROSINE', 'D-VALINE', '2,4-DIAMINOBUTYRICACID', '(2S)-2,3-DIAMINOBUTYRICACID', 'DEHYDROALANINE', '2,3-DEHYDRO-TRYPTOPHAN', 'DEHYDRO-CYSTEINE', '3,5-DIHYDROXYPHENYLGLYCINE', 'D-DIHYDROXYPHENYLTHIAZOLGROUP', 'DIHYDROXYPHENYLTHIAZOLGROUP', '2,3-DIAMINOPROPIONICACID', 'ENDURACIDIDINE', 'ALPHA-ETHYLNORVALINE', 'GLUTAMINE', 'GLUTAMICACID', 'GLYCINE', 'HOMOARGININE', 'HOMOISOLEUCINE', 'HISTIDINE', 'HYDROXYPHENYLGLYCINE', 'HOMOPHENYLALANINE', 'HOMOPROLINE', 'HOMOSERINE', 'HOMOSERINELACTONE', 'HOMOTYROSINE', 'ISOLEUCINE', 'ISOVALINE', 'KYNURENINE', '(2S,3S)-2,3-DIAMINOBUTYRICACID', 'LEUCINE', 'LYSINE', 'METHYL-2-AMINOOCTANOICACID', 'METHIONINE', 'NORCORONAMICACID', 'NORVALINE', 'ORNITHINE', 'PHENYLGLYCINE', 'PHENYLALANINE', 'PROLINE', 'PHOSPHINOTHRICIN', 'R-CORONAMICACID', 'R-NORCORONAMICACID', 'SERINE', 'TERT-LEU', 'THREONINE', 'TRYPTOPHAN', 'TYROSINE', 'VALINE', 'DEHYDROBUTYRINE', 'BETA-ALANINE', 'BETA-LYSINE', '(2S,3R)-2-AMINO-3-HYDROXY-4-(4-NITROPHENYL)BUTANOATE', '2,3-DIAMINO-3-METHYL-PROPANOICACID', '4,5-DIHYDROXYORNITHINE', '4,5-DIHYDROXYHOMOTYROSINE', '3-HYDROXYGLUTAMINE', '4-BUTENYL-4-METHYLTHREONINE', '4-HYDROXYVALINE', '2-METHYLSERINE', '4,5-DEHYDROARGININE', '3-(3-PYRIDYL)-ALANINE', '3-METHYLTYROSINE', 'D-LYSERGICACID', 'N-ACETYLPHENYLALANINE', '3-HYDROXY-4-O-METHYL-5-METHYL-TYROSINE', '4-METHYLPROLINE', 'ISOLEUCINE', '3-[(1R,2R)-2-NITROCYCLOPROPYL]-L-ALANINE', 'DEHYDROBUTYRINE', '(R)BETA-TYROSINE', 'PIPECOLICACID', '2-AMINO-9,10-EPOXY-8-OXODECANOICACID', 'L-PIPERAZICACID', 'DEHYDROTYROSINE', 'BETA-HYDROXYCYCLOHEX‐2‐ENYLALANINE', '2-METHYL-BETA-ALANINE', 'GLYCOCYAMINE']
    list_aa = ['alanine', 'arginine', 'asparagine', 'asparticacid', 'cysteine',
               'glutamine', 'glutamicacid', 'glycine', 'histidine',
               'isoleucine', 'leucine', 'lysine', 'methionine',
               'phenylalanine', 'proline', 'serine', 'threonine', 'tryptophan',
               'tyrosine', 'valine', 'glycine']
    if random.randint(0,1) == 1:
        starter_module = ['module 1', 'starter_module', random.choice(starter_units)]
    else:
        starter_module = ['module 1', 'starter_module_nrps', random.choice(list_aa)]
    cluster = [starter_module]
    last_elongation_domain_nr = random.randint(2, 8)
    for i in range(2, last_elongation_domain_nr):
        if random.randint(0,1) == 1:
            elongation_module = [f'module {i}', 'elongation_module', random.choice(elongation_units), random.choice(tailoring_domain_combinations)]
        else:
            elongation_module = [f'module {i}', 'elongation_module_nrps', random.choice(list_aa)]
        cluster.append(elongation_module)
    if random.randint(0,1) == 1:
        terminator_module =  [f'module {random.randint(2, 8)+1}', 'terminator_module', random.choice(elongation_units), random.choice(tailoring_domain_combinations)]
    else:
        terminator_module = [f'module {random.randint(2, 8) + 1}', 'terminator_module_nrps', random.choice(list_aa)]
    cluster.append(terminator_module)

    return cluster

if __name__ == "__main__":
    # for i in range(472, 501):
    #     cluster = generate_random_nrps_cluster()
    #     print(cluster)
    #     #draw_pks_cluster(cluster, save_fig=f'NRPS_cluster_{i}.png')
    #     draw_pks_cluster(cluster)
    #     plt.close('all')


    # for i in range(1, 501):
    #     generated_cluster = generate_random_pks_cluster()
    #     print(generated_cluster)
    #     draw_pks_cluster(generated_cluster, save_fig=f'cluster_{i}.png')
    #     plt.close('all')

    for i in range(1, 501):
        generated_cluster = generate_random_hybrid_cluster()
        print(generated_cluster)
        #draw_pks_cluster(generated_cluster)
        draw_pks_cluster(generated_cluster, save_fig=f'hybrid_PKS_NRPS_cluster_{i}.png')
        plt.close('all')

