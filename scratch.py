from pks_modules_to_structure import *
from pks_thioesterase_reactions import *
from pikachu.general import structure_to_smiles
from visualize_pks_cluster import draw_pks_cluster
from pikachu.general import read_smiles, structure_to_smiles
from NRPS_condensation import sulphur_to_hydroxyl

from pikachu.drawing.drawing import Drawer as Drawer_pik
erythromycin_cluster = [['module_1', 'starter_module', 'SC(=O)CC'],
                       ['module_2', 'elongation_module', 'methylmalonylcoa', ['KR_B2']],
                       ['module_3', 'elongation_module', 'methylmalonylcoa', ['KR_A1']],
                       ['module_4', 'elongation_module', 'methylmalonylcoa', ['KR_C2']],
                       ['module_5', 'elongation_module', 'methylmalonylcoa', ['KR', 'DH', 'ER']],
                       ['module_6', 'elongation_module', 'methylmalonylcoa', ['KR_A1']],
                       ['module_7', 'terminator_module', 'methylmalonylcoa', ['KR_A1']]]
daptomycin_cluster = [['module_1', 'starter_module_nrps', 'tryptophan'],
                      ['module_2', 'elongation_module_nrps', 'd-asparagine'],
                      ['module_3', 'elongation_module_nrps', 'asparticacid'],
                      ['module_4', 'elongation_module_nrps', 'threonine'],
                      ['module_5', 'elongation_module_nrps', 'glycine'],
                      ['module_6', 'elongation_module_nrps', 'ornithine'],
                      ['module_7', 'elongation_module_nrps', 'asparticacid'],
                      ['module_8', 'elongation_module_nrps', 'd-alanine'],
                      ['module_9', 'elongation_module_nrps', 'asparticacid'],
                      ['module_10', 'elongation_module_nrps', 'glycine'],
                      ['module_11', 'elongation_module_nrps', 'd-serine'],
                      ['module_12', 'elongation_module_nrps', '4-methylglutamicacid'],
                      ['module_13', 'terminator_module_nrps', 'kynurenine']]
baf_cluster = [['pks module 1 [5148:6510]', 'starter_module', 'SC(C(C(O)=O)C)=O'],
               ['pks module 2 [6591:11160]', 'elongation_module', 'methylmalonylcoa', ['KR_B1']],
               ['pks module 3 [11277:15915]', 'elongation_module', 'malonylcoa', ['KR_A1']],
               ['pks module 4 [16029:19218]', 'elongation_module', 'methylmalonylcoa', []],
               ['pks module 5 [19606:24247]', 'elongation_module', 'methylmalonylcoa', ['KR_A2']],
               ['pks module 6 [24325:29008]', 'elongation_module', 'ethylmalonylcoa', ['KR_B1']],
               ['pks module 7 [29086:34684]', 'elongation_module', 'malonylcoa', ['KR', 'DH']],
               ['pks module 8 [35081:40406]', 'elongation_module', 'pk', ['KR_B1', 'DH']],
               ['pks module 9 [40502:46613]', 'elongation_module', 'methylmalonylcoa', ['KR_B1', 'DH', 'ER']],
               ['pks module 10 [47112:51861]', 'elongation_module', 'methylmalonylcoa', ['KR_A2']],
               ['pks module 11 [51984:57198]', 'elongation_module', 'methylmalonylcoa', ['KR_B1', 'DH']],
               ['pks module 12 [57607:63796]', 'terminator_module', 'methylmalonylcoa', ['KR_B1', 'DH']]]
bor = [['pks module 1 [17152:18484]', 'starter_module', 'SC(C1C(CCC1)C(=O)O)=O'],
       ['pks module 2 [18979:23350]', 'elongation_module', 'malonylcoa', ['KR_B1']],
       ['pks module 3 [23790:28788]', 'elongation_module', 'malonylcoa', ['KR_B1', 'DH']],
       ['pks module 4 [28860:33924]', 'elongation_module', 'methylmalonylcoa', ['KR_B1', 'DH']],
       ['pks module 5 [34289:38780]', 'elongation_module', 'methylmalonylcoa', ['KR_A1']],
       ['pks module 6 [39226:45262]', 'elongation_module', 'methylmalonylcoa', ['KR_B1', 'DH', 'ER']],
       ['pks module 7 [45627:50703]', 'terminator_module', 'malonylcoa', ['KR_A1']]]
pik = [['pks module 1 [3181:6154]', 'starter_module', 'SC(CC)=O'],
       ['pks module 2 [6241:10615]', 'elongation_module', 'methylmalonylcoa', ['KR_B2']],
       ['pks module 3 [10726:16357]', 'elongation_module', 'malonylcoa', ['KR', 'DH']],
       ['pks module 4 [16968:21414]', 'elongation_module', 'methylmalonylcoa', ['KR_C2']],
       ['pks module 5 [21480:27792]', 'elongation_module', 'methylmalonylcoa', ['KR', 'DH', 'ER']],
       ['pks module 6 [28268:32591]', 'elongation_module', 'methylmalonylcoa', ['KR_A1']],
       ['pks module 7 [33062:36956]', 'terminator_module', 'methylmalonylcoa', []]]
tyr = [['module_1', 'starter_module_nrps', 'd-phenylalanine'],
                      ['module_2', 'elongation_module_nrps', 'proline'],
                      ['module_3', 'elongation_module_nrps', 'phenylalanine'],
                      ['module_4', 'elongation_module_nrps', 'd-phenylalanine'],
                      ['module_5', 'elongation_module_nrps', 'asparagine'],
                      ['module_6', 'elongation_module_nrps', 'glutamine'],
                      ['module_7', 'elongation_module_nrps', 'tyrosine'],
                      ['module_8', 'elongation_module_nrps', 'valine'],
                      ['module_9', 'elongation_module_nrps', 'ornithine'],
                      ['module_10', 'elongation_module_nrps', 'leucine']]
tylactone = [['pks module 1 [947:3932]', 'starter_module', 'SC(CC)=O'],
             ['pks module 2 [3998:8474]', 'elongation_module', 'methylmalonylcoa', ['KR_B1']],
             ['pks module 3 [8546:13976]', 'elongation_module', 'methylmalonylcoa', ['KR_B1', 'DH']],
             ['pks module 4 [31342:35743]', 'elongation_module', 'methylmalonylcoa', ['KR_A1']],
             ['pks module 5 [20141:24608]', 'elongation_module', 'methylmalonylcoa', ['KR_C2']],
             ['pks module 6 [24680:30905]', 'elongation_module', 'ethylmalonylcoa', ['KR_B1', 'DH', 'ER']],
             ['pks module 7 [14416:19672]', 'elongation_module', 'malonylcoa', ['KR_B1', 'DH']],
             ['pks module 8 [36365:41705]', 'terminator_module', 'malonylcoa', ['KR_A1']]]
