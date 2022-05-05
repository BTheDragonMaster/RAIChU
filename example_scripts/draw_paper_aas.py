from pikachu.general import read_smiles, Drawer

threonine = read_smiles("N[C@@H]([C@@H](C)O)C(=O)SI")

drawing = Drawer(threonine)

for atom in drawing.structure.graph:
    if atom.type == 'I':
        atom.type = 'PCP'

drawing.structure.refresh_structure()
drawing.draw()
drawing.save_svg("threonine.svg")