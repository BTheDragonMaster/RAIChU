lines = open('aa_smiles.txt', 'r').readlines()
from pks_elongation_reactions import *
from pikachu.reactions.functional_groups import *
from pikachu.general import structure_to_smiles
HYDROXYL_CARBOXYLIC_ACID = GroupDefiner('Hydroxyl group carboxylic acid amino acid', 'NCC(O)=O', 3)



outputfile = open('smiles_aa_attached.txt', 'w')
for line in lines:
    line = line.strip()
    name, threeletter, oneletter, formula, smiles, polarity, charge, hydropathy_index, essentiality, pi, pk1, pk2 = line.split()
    structure_amino_acid = Smiles(smiles).smiles_to_structure()
    o_to_remove = find_atoms(HYDROXYL_CARBOXYLIC_ACID, structure_amino_acid)[0]
    for neighbour in o_to_remove.neighbours:
        if neighbour.type == 'C':
            structure_amino_acid.break_bond_between_atoms(o_to_remove, neighbour)
            one, two = structure_amino_acid.split_disconnected_structures()
            if len(one.graph) == 2:
                hydroxyl = one
                structure_amino_acid = two
            elif len(two.graph) == 2:
                hydroxyl = two
                structure_amino_acid = one
            structure_amino_acid.add_atom('S', [neighbour])
            structure_amino_acid.refresh_structure()
            structure_amino_acid.find_cycles()
            print(name)
            new_smiles = structure_to_smiles(structure_amino_acid)
            outputfile.write(name)
            outputfile.write('\t')
            outputfile.write(new_smiles)
            outputfile.write('\n')
outputfile.close()
