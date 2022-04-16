from sys import argv

from pikachu.general import read_smiles, svg_from_structure
from pikachu.drawing.colours import RANDOM_PALETTE_1


def draw_paras_smiles(smiles_file):
    with open(smiles_file, 'r') as paras_smiles:
        for i, line in enumerate(paras_smiles):
            name, smiles = line.strip().split()
            structure = read_smiles(smiles)
            for atom in structure.graph:
                atom.draw.colour = RANDOM_PALETTE_1[i % len(RANDOM_PALETTE_1)]
            svg_from_structure(structure, f'{name}.svg')


if __name__ == "__main__":
    smiles_file = argv[1]
    draw_paras_smiles(smiles_file)


