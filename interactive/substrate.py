import os
import pygame

import interactive.flatfiles
import interactive.images.nrps_substrates
import interactive.images.pks_substrates

from interactive.style import SUBSTRATE_BUTTON_SIZE, SUBSTRATE_BUTTON_PADDING, \
    SUBSTRATE_BUTTONS_PER_LINE, WIDTH, HEIGHT

SUBSTRATE_DIR = os.path.dirname(interactive.flatfiles.__file__)
NRPS_SUBSTRATE_IMG = os.path.dirname(interactive.images.nrps_substrates.__file__)
PKS_SUBSTRATE_IMG = os.path.dirname(interactive.images.pks_substrates.__file__)
SUBSTRATE_TO_ABBREVIATION_FILE = os.path.join(SUBSTRATE_DIR, 'substrate_to_abbreviation.txt')

PROTEINOGENIC_SUBSTRATES = {"Ala": ['Alanine',
                                    'D-Alanine',
                                    'dehydroalanine',
                                    'beta-Alanine',
                                    '2-methyl-beta-alanine'],
                            "Arg": ["Arginine",
                                    "D-Arginine",
                                    "4,5-dehydroarginine",
                                    "homoarginine"],
                            "Asn": ["Asparagine",
                                    "D-Asparagine"],
                            "Asp": ["AsparticAcid",
                                    "D-AsparticAcid"],
                            "Cys": ["Cysteine",
                                    "D-Cysteine",
                                    "dehydro-cysteine",
                                    "cysteicacid",
                                    "D-cysteicacid",
                                    "alpha-methylcysteine",
                                    "D-alpha-methylcysteine"],
                            "Gln": ["Glutamine",
                                    "D-Glutamine",
                                    "3-Hydroxyglutamine"],
                            "Glu": ["GlutamicAcid",
                                    "D-GlutamicAcid",
                                    "4-methylglutamicacid"],
                            "Gly": ["Glycine"],
                            "His": ["Histidine",
                                    "D-Histidine"],
                            "Ile": ["Isoleucine",
                                    "D-Isoleucine",
                                    "allo-Isoleucine",
                                    "D-allo-Isoleucine",
                                    "homoisoleucine"],
                            "Leu": ["Leucine",
                                    "D-Leucine",
                                    "tert-Leu",
                                    "D-tert-Leucine"],
                            "Lys": ["Lysine",
                                    "D-Lysine",
                                    "beta-lysine"],
                            "Met": ["Methionine",
                                    "D-Methionine"],
                            "Phe": ["Phenylalanine",
                                    "D-Phenylalanine",
                                    "Homophenylalanine",
                                    "N-acetylphenylalanine",
                                    "alpha-amino-phenyl-valericacid",
                                    "D-alpha-amino-phenyl-valericacid"],
                            "Pro": ["Proline",
                                    "D-Proline",
                                    "4-methylproline"],
                            "Ser": ["Serine",
                                    "D-Serine",
                                    "Homoserine",
                                    "D-Homoserine",
                                    "2-Methylserine"],
                            "Thr": ["Threonine",
                                    "D-Threonine",
                                    "allo-Threonine",
                                    "D-allo-Threonine",
                                    "4-butenyl-4-methylthreonine"],
                            "Trp": ["Tryptophan",
                                    "D-Tryptophan",
                                    "2,3-Dehydro-Tryptophan"],
                            "Tyr": ["Tyrosine",
                                    "D-Tyrosine",
                                    "dehydrotyrosine",
                                    "(R)beta-tyrosine",
                                    "Homotyrosine",
                                    "D-Homotyrosine",
                                    "3-methyltyrosine",
                                    "4,5-Dihydroxyhomotyrosine",
                                    "3-hydroxy-4-O-methyl-5-methyltyrosine"],
                            "Val": ["Valine",
                                    "D-Valine",
                                    "4-Hydroxyvaline",
                                    "isovaline",
                                    "D-isovaline",
                                    "alpha-ethylnorvaline",
                                    "D-alpha-ethylnorvaline",
                                    "Norvaline"]}

NONPROTEINOGENIC_SUBSTRATES = {"Aad": ["2-Aminoadipicacid",
                                       "D-2-Aminoadipicacid"],
                               "Cit": ["Citrulline",
                                       "D-Citrulline"],
                               "Cor": ["coronamicacid",
                                       "R-coronamicacid",
                                       "norcoronamicacid",
                                       "R-norcoronamicacid"],
                               "Dab": ["2,4-diaminobutyricacid",
                                       "D-2,4-diaminobutyricacid",
                                       "(2S)-2,3-diaminobutyricacid",
                                       "(2S,3S)-2,3-diaminobutyricacid",
                                       "(2R,3R)-2,3-diaminobutyricacid",
                                       "2,3-diamino-3-methyl-propanoicacid",
                                       "D-2,3-Diaminopropionicacid",
                                       "2-Aminobutyricacid",
                                       "D-2-Aminobutyricacid",
                                       "2-Aminoisobutyricacid",
                                       "Dehydrobutyrine"],
                               "End": ["enduracididine",
                                       "D-enduracididine"],
                               "Kyn": ["kynurenine",
                                       "D-kynurenine"],
                               "Orn": ["Ornithine",
                                       "D-Ornithine",
                                       "4,5-Dihydroxyornithine"],
                               "PheGly": ["phenylglycine",
                                          "D-PhenylGlycine",
                                          "HydroxyPhenylGlycine",
                                          "D-HydroxyPhenylGlycine",
                                          "3,5-dihydroxyphenylglycine",
                                          "D-3,5-dihydroxyphenylglycine"],
                               "Pip": ["Pipecolicacid",
                                       "D-Pipecolicacid",
                                       "L-piperazicacid"],
                               "Other": ["Anticapsin",
                                         "Beta-hydroxycyclohex-2-enylalanine",
                                         "capreomycidine",
                                         "phosphinothricin",
                                         "D-phosphinothricin",
                                         "(2S,3R)-2-amino-3-hydroxy-4-(4-nitrophenyl)butanoate",
                                         "3-(3-Pyridyl)-Alanine",
                                         "3-[(1R,2R)-2-Nitrocyclopropyl]-L-alanine"]}

NON_AMINO_ACIDS = {'Aromatic': ["2-carboxyquinoxaline",
                                "Quinoline-2-carboxylicacid",
                                "2,3-dihydroxybenzoicacid",
                                "4-Hydroxyindole-3-carboxylicacid",
                                "dihydroxyphenylthiazolgroup",
                                "D-dihydroxyphenylthiazolgroup",
                                "D-lysergicacid",
                                "hydrocinnamicacid",
                                "hydroxypicolinicacid",
                                "phenylaceticacid",
                                "salicylicacid"],
                   'Cyclic': ["cycloOrnithine",
                              "D-cycloOrnithine",
                              "homoserinelactone",
                              "D-homoserinelactone"],
                   'Other': ["2-hydroxyisovalerate",
                             "3-aminopropanimidamide",
                             "Alaninol",
                             "alpha-ketoisocaproicacid",
                             "alpha-ketoisovalericacid",
                             "D-2-hydroxyisovalerate",
                             "Glycocyamine",
                             "Glycolicacid",
                             "valinol"]}

FATTY_ACIDS = {'Amino group': ['2-Amino-9,10-epoxy-8-oxodecanoidacid',
                               'D-methyl-2-aminooctanoicacid',
                               'methyl-2-aminooctanoicacid'],
                               
               'No amino group': []}


def parse_substrate_to_abbr():
    substrate_to_abbreviation = {}
    with open(SUBSTRATE_TO_ABBREVIATION_FILE, 'r') as abbreviation_file:
        for line in abbreviation_file:
            line = line.strip()
            if line:
                substrate, abbreviation = line.split()
                substrate_to_abbreviation[substrate] = abbreviation

    return substrate_to_abbreviation


SUBSTRATE_TO_ABBREVIATION = parse_substrate_to_abbr()


class SubstrateGroup:
    def __init__(self, name, substrates):
        self.name = name
        self.substrates = substrates
        self.coords = []
        self.create_raster()
        self.set_rectangles()

    def create_raster(self):

        x_coord = int(0.23 * WIDTH)
        y_coord = int(0.52 * HEIGHT)

        for i, substrate in enumerate(self.substrates):
            if i % SUBSTRATE_BUTTONS_PER_LINE == 0:
                if i != 0:
                    y_coord += SUBSTRATE_BUTTON_SIZE + SUBSTRATE_BUTTON_PADDING
                x_coord = int(0.23 * WIDTH)
            else:
                x_coord += SUBSTRATE_BUTTON_SIZE + SUBSTRATE_BUTTON_PADDING

            self.coords.append((x_coord, y_coord))

    def set_rectangles(self):

        for i, substrate in enumerate(self.substrates):
            x, y = self.coords[i]
            substrate.x = x
            substrate.y = y
            substrate.rectangle = pygame.Rect(x, y, SUBSTRATE_BUTTON_SIZE, SUBSTRATE_BUTTON_SIZE)


class Substrate:
    def __init__(self, name, smiles, module_type):
        self.name = name
        self.abbr = SUBSTRATE_TO_ABBREVIATION[self.name]
        self.smiles = smiles
        self.image = None
        self.highlight_image = None
        if module_type == 'NRPS':
            self.image = os.path.join(NRPS_SUBSTRATE_IMG, f'{name}.png')
            self.highlight_image = os.path.join(NRPS_SUBSTRATE_IMG, f'{name}_highlight.png')
        if module_type == 'PKS':
            self.image = os.path.join(PKS_SUBSTRATE_IMG, f'{name}.png')
            self.highlight_image = os.path.join(PKS_SUBSTRATE_IMG, f'{name}_highlight.png')
        self.rectangle = None
        self.x = 0
        self.y = 0
        self.width = SUBSTRATE_BUTTON_SIZE
        self.height = SUBSTRATE_BUTTON_SIZE
