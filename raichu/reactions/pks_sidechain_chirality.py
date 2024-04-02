from pikachu.reactions.functional_groups import find_atoms
from pikachu.chem.chirality import same_chirality
from raichu.data.molecular_moieties import METHYLMALONYL_SIDECHAIN_CARBON, \
    METHYLMALONYL_SULPHUR_CARBON, METHYLMALONYL_MAINCHAIN_CARBON, METHYLMALONYL_CHIRAL_CARBON


def set_sidechain_chirality(structure):
    chiral_carbons = find_atoms(METHYLMALONYL_CHIRAL_CARBON, structure)
    if len(chiral_carbons) == 1:
        chiral_carbon = chiral_carbons[0]
        sulphur_carbon = find_atoms(METHYLMALONYL_SULPHUR_CARBON, structure)[0]
        sidechain_carbon = find_atoms(METHYLMALONYL_SIDECHAIN_CARBON, structure)[0]
        mainchain_carbon = find_atoms(METHYLMALONYL_MAINCHAIN_CARBON, structure)[0]
        if chiral_carbon.has_neighbour('H'):
            hydrogen = chiral_carbon.get_neighbour('H')
            counterclockwise_order = [hydrogen, sidechain_carbon, sulphur_carbon, mainchain_carbon]
            if same_chirality(counterclockwise_order, chiral_carbon.neighbours):
                chiral_carbon.chiral = 'counterclockwise'
            else:
                chiral_carbon.chiral = 'clockwise'