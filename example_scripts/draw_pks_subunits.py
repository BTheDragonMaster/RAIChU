from pikachu.general import read_smiles
from pikachu.drawing.drawing import Drawer
from pikachu.drawing.colours import RANDOM_PALETTE_1

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

names = ["malonylcoa",
         "methylmalonylcoa",
         "ethylmalonylcoa",
         "methoxymalonylcoa",
         "wildcard"]

for i, drawing in enumerate(drawings):
    colour = RANDOM_PALETTE_1[i]
    print(colour)
    for atom in drawing.structure.graph:
        if atom.type == 'I':
            print("here")
            atom.type = 'CoA'
        atom.draw.colour = colour

    drawing.structure.refresh_structure()
    drawing.draw()

    name = names[i]

    drawing.save_svg(f"{name}.svg")




