from sys import argv

from pikachu.general import svg_from_smiles


def parse_smiles(smiles_file):
    name_to_smiles = {}
    with open(smiles_file, 'r') as smiles_in:
        for line in smiles_in:
            line = line.strip()
            name, smiles = line.split('\t')
            name_to_smiles[name] = smiles

    return name_to_smiles

name_to_smiles = parse_smiles(argv[1])

for name, smiles in name_to_smiles.items():
    svg_from_smiles(smiles, f"{name}.svg")
