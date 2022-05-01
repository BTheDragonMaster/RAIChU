from pikachu.general import read_smiles
from pikachu.drawing.drawing import Drawer, Options
from pikachu.drawing.colours import RANDOM_PALETTE_1

# malonylcoa = read_smiles("OC(=O)CC(=O)SI")
# methylmalonylcoa = read_smiles("OC(=O)C(C)C(=O)SI")
# ethylmalonylcoa = read_smiles("OC(=O)C(CC)C(=O)SI")
# methoxymalonylcoa = read_smiles("OC(=O)C(OC)C(=O)SI")
# wildcard = read_smiles("OC(=O)C(*)C(=O)SI")

prop = read_smiles("ISC(=O)CC")
ac = read_smiles("ISC(C)=O")
benz = read_smiles("ISC(C1=CC=CC=C1)=O")
but3 = read_smiles("ISC(=O)CC(C)C")
but2 = read_smiles("ISC(=O)C(C)CC")
but2s = read_smiles("ISC(=O)[C@@H](C)CC")
but2r = read_smiles("ISC(=O)[C@H](C)CC")
tcpda = read_smiles("ISC([C@H]1[C@@H](CCC1)C(=O)O)=O")
chc = read_smiles("ISC(=O)C1CCCCC1")
hmal = read_smiles("ISC(=O)C(C(=O)O)O")
hmals = read_smiles("ISC(=O)[C@@H](C(=O)O)O")
hmalr = read_smiles("ISC(=O)[C@H](C(=O)O)O")
cemal = read_smiles("ISC(C(C(=O)O)CC[Cl])=O")
isobut = read_smiles("ISC(=O)C(C)C")

structures = [isobut, ac, benz, but3, but2, but2s, but2r, tcpda, chc, hmal, hmals, hmalr, cemal]
names = ["isoBut", "Ac", "Benz", "3mBut", "2mBut", "2S-mBut", "2R-mBut", "tCPDA", "CHC", "hMal", "2S-hMal", "2R-hMal", "ceMal"]

# malonylcoa_drawing = Drawer(malonylcoa)
# methylmalonylcoa_drawing = Drawer(methylmalonylcoa)
# ethylmalonylcoa_drawing = Drawer(ethylmalonylcoa)
# methoxymalonylcoa_drawing = Drawer(methoxymalonylcoa)
# wildcard_drawing = Drawer(wildcard)
#
# drawings = [malonylcoa_drawing,
#             methylmalonylcoa_drawing,
#             ethylmalonylcoa_drawing,
#             methoxymalonylcoa_drawing,
#             wildcard_drawing]
#
# names = ["malonylcoa",
#          "methylmalonylcoa",
#          "ethylmalonylcoa",
#          "methoxymalonylcoa",
#          "wildcard"]

for i, structure in enumerate(structures):
    name = names[i]
    print(name)
    colour = RANDOM_PALETTE_1[i]
    options = Options()
    options.finetune = False
    drawing = Drawer(structure, options=options)

    for atom in drawing.structure.graph:
        if atom.type == 'I':
            atom.type = 'CoA'
        atom.draw.colour = colour

    drawing.structure.refresh_structure()
    drawing.draw()



    drawing.save_svg(f"{name}.svg")




