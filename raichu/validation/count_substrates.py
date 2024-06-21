from raichu.smiles_handling import _METADATA, _SMILES

amino_acids = 0
beta_amino_acids = 0
carboxylic_acids = 0
c_term = 0
total = 0
d_substrates = 0
l_substrates = 0

for substrate_name in _SMILES:
    if substrate_name not in _METADATA:
        print(substrate_name, _SMILES[substrate_name])
    else:
        if _METADATA[substrate_name].type in ["amino_acid", "branched_amino_acid"]:
            amino_acids += 1
        elif _METADATA[substrate_name].type == "beta_amino_acid":
            beta_amino_acids += 1
        elif _METADATA[substrate_name].type in ["double_acid", "acid"]:
            carboxylic_acids += 1
        elif _METADATA[substrate_name].type == "ter_amino":
            c_term += 1
        else:
            print(_METADATA[substrate_name].type)

        if _METADATA[substrate_name].chirality == 'D':
            d_substrates += 1
        elif _METADATA[substrate_name].chirality == 'L':
            l_substrates += 1

        total += 1

print(f"Amino acid substrates: {amino_acids}")
print(f"Beta-amino acid substrates: {beta_amino_acids}")
print(f"Carboxylic acid substrates: {carboxylic_acids}")
print(f"Amine substrates: {c_term}")
print(f"D-amino acid substrates: {d_substrates}")
print(f"L-amino acid substrates: {l_substrates}")
print(f"Total substrates: {total}")
