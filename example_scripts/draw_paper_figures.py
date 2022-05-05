from pikachu.general import read_smiles, Drawer

malonylcoa = read_smiles("OC(=O)CC(=O)SI")
methylmalonylcoa = read_smiles("OC(=O)C(C)C(=O)SI")
ethylmalonylcoa = read_smiles("OC(=O)C(CC)C(=O)SI")
methoxymalonylcoa = read_smiles("OC(=O)C(OC)C(=O)SI")
wildcard = read_smiles("OC(=O)C(*)C(=O)SI")

malonylcoa_drawing = Drawer(malonylcoa)
methylmalonylcoa_drawing = Drawer(methylmalonylcoa)
ethylmalonylcoa_drawing = Drawer(ethylmalonylcoa)
methoxymalonylcoa_drawing = Drawer(methoxymalonylcoa)
wildcard_drawing = Drawer(wildcard)

drawings = [malonylcoa_drawing,
            methylmalonylcoa_drawing,
            ethylmalonylcoa_drawing,
            methoxymalonylcoa_drawing,
            wildcard_drawing]

for i, drawing in enumerate(drawings):

    for atom in drawing.structure.graph:
        if atom.type == 'I':
            print("here")
            atom.type = 'ACP'

    drawing.structure.refresh_structure()
    drawing.draw()
    drawing.save_svg("malonylcoa.svg")
    break
