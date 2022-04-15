import pygame
import os

from interactive.style import BUTTON_TEXT_COLOUR, FONT, BUTTON_HIGHLIGHT_COLOUR, BLACK, BUTTON_COLOUR, \
    BUTTON_PANEL_COLOUR, HEIGHT, WIDTH, DOMAIN_BUTTON_SIZE, SUBSTRATE_GROUP_BUTTONS_PER_LINE, \
    SUBSTRATE_GROUP_BUTTON_PADDING, SUBSTRATE_GROUP_BUTTON_SIZE, SUBSTRATE_BUTTON_SIZE
import interactive.images.domains
import interactive.flatfiles
from interactive.gene import Gene
from interactive.textbox import TextBox
from interactive.parsers import parse_smiles
from interactive.substrate import Substrate, SubstrateGroup, \
    PROTEINOGENIC_SUBSTRATES, NONPROTEINOGENIC_SUBSTRATES, FATTY_ACIDS, NON_AMINO_ACIDS

DOMAIN_IMAGE_DIR = os.path.dirname(interactive.images.domains.__file__)
FLATFILES = os.path.dirname(interactive.flatfiles.__file__)
PARAS_SMILES = os.path.join(FLATFILES, "PARAS_smiles.txt")


class Button:
    def __init__(self, text, position, dimensions, font_size=None):
        self.text = text

        self.x, self.y = position
        self.width, self.height = dimensions
        self.rectangle = pygame.Rect(position, dimensions)
        self.font_size = self.height - 10

        self.text_x = 0
        self.text_y = 0

        self.set_text_position()
        self.font = pygame.font.SysFont(FONT, self.font_size, bold=False)

        self.rendered_text = self.font.render(text, True, BUTTON_TEXT_COLOUR)
        self.selected = False

    def __hash__(self):
        return hash(self.text)

    def __eq__(self, other):
        return self.text == other.text

    def __repr__(self):
        return self.text

    def highlight(self, screen):
        pygame.draw.rect(screen, BUTTON_HIGHLIGHT_COLOUR, self.rectangle)
        pygame.draw.rect(screen, BLACK, self.rectangle, 2)
        screen.blit(self.rendered_text, (self.text_x, self.text_y))

    def draw(self, screen):
        pygame.draw.rect(screen, BUTTON_COLOUR, self.rectangle)
        pygame.draw.rect(screen, BLACK, self.rectangle, 2)
        screen.blit(self.rendered_text, (self.text_x, self.text_y))

    def erase_button(self, screen):
        pygame.draw.rect(screen, BUTTON_PANEL_COLOUR, self.rectangle)
        pygame.draw.rect(screen, BUTTON_PANEL_COLOUR, self.rectangle, 2)

    def set_text_position(self):
        self.text_x = self.x + 5
        self.text_y = self.y + int((self.height - self.font_size) / 2)

    def set_font(self):
        self.font = pygame.font.SysFont(FONT, self.font_size, bold=True)


class DomainButton(Button):
    def __init__(self, position, text):
        dimensions = (DOMAIN_BUTTON_SIZE, DOMAIN_BUTTON_SIZE)
        super().__init__(text, position, dimensions)
        self.image = None
        self.highlight_image = None
        self.domain_type = text

        self.image = os.path.join(DOMAIN_IMAGE_DIR, f'{text.lower()}.png')
        self.highlight_image = os.path.join(DOMAIN_IMAGE_DIR, f'{text.lower()}_highlight.png')

    def draw(self, screen):
        domain_image = pygame.image.load(self.image)
        domain_image_scaled = pygame.transform.smoothscale(domain_image, (DOMAIN_BUTTON_SIZE, DOMAIN_BUTTON_SIZE))
        screen.blit(domain_image_scaled, self.rectangle)

    def highlight(self, screen):
        domain_image = pygame.image.load(self.highlight_image)
        domain_image_scaled = pygame.transform.smoothscale(domain_image, (DOMAIN_BUTTON_SIZE, DOMAIN_BUTTON_SIZE))
        screen.blit(domain_image_scaled, self.rectangle)

    def do_action(self, module, mouse, screen, active_buttons):
        module.gene.erase()
        module.add_domain(self.domain_type)
        module.gene.draw(mouse)
        reset_buttons(screen, active_buttons)


class AddGeneButton(Button):
    def __init__(self):
        position = (int(0.02 * WIDTH), int(0.82 * HEIGHT))
        dimensions = (int(0.2 * WIDTH), int(HEIGHT / 25))

        super().__init__("Add gene", position, dimensions)

    def do_action(self, screen, genes, text_box, active_buttons, mouse):
        gene_number = len(genes)
        gene = Gene(screen, gene_number)
        gene.name = text_box.text
        text_box.erase(screen)
        gene.draw(mouse)
        genes.append(gene)
        reset_buttons(screen, active_buttons)


class RemoveModuleButton(Button):

    def __init__(self):
        position = (int(0.02 * WIDTH), int(0.82 * HEIGHT))
        dimensions = (int(0.2 * WIDTH), int(HEIGHT / 25))

        super().__init__("Remove module", position, dimensions)

    def do_action(self, screen, current_module, active_buttons, mouse):
        gene = current_module.gene
        gene.remove_module(current_module.id)

        gene.erase()
        gene.draw(mouse)

        reset_buttons(screen, active_buttons)


class RemoveGeneButton(Button):
    def __init__(self):
        position = (int(0.02 * WIDTH), int(0.82 * HEIGHT))
        dimensions = (int(0.2 * WIDTH), int(HEIGHT / 25))

        super().__init__("Remove gene", position, dimensions)

    def do_action(self, screen, genes, current_gene, active_buttons, mouse):
        gene_number = current_gene.gene_number
        index_to_remove = None
        for i, gene in enumerate(genes):
            if gene.gene_number == gene_number:
                index_to_remove = i
            elif gene.gene_number > gene_number:
                gene.gene_number -= 1

        genes.pop(index_to_remove)
        current_gene.erase()

        for gene in genes:
            gene.erase()
            gene.draw(mouse)

        reset_buttons(screen, active_buttons)




class CreateGeneButton(Button):
    def __init__(self):
        position = (int(0.02 * WIDTH), int(0.77 * HEIGHT))
        dimensions = (int(0.2 * WIDTH), int(HEIGHT / 25))

        super().__init__("Create gene", position, dimensions)

    def do_action(self, screen, active_buttons):
        position = (int(0.02 * WIDTH), int(0.77 * HEIGHT))
        dimensions = (int(0.2 * WIDTH), int(HEIGHT / 25))
        text_box = TextBox(position, dimensions)
        text_box.draw(screen)
        hide_buttons(ALL_BUTTONS, screen, active_buttons)
        show_button(ADD_GENE_BUTTON, screen, active_buttons)
        return text_box


class AddModuleButton(Button):
    def __init__(self):
        position = (int(0.02 * WIDTH), int(0.87 * HEIGHT))
        dimensions = (int(0.2 * WIDTH), int(HEIGHT / 25))

        super().__init__("Add module", position, dimensions)


class NRPSModuleButton(Button):
    def __init__(self):
        position = (int(0.02 * WIDTH), int(0.92 * HEIGHT))
        dimensions = (int(0.09 * WIDTH), int(HEIGHT / 25))

        super().__init__('NRPS', position, dimensions)

    def do_action(self, screen, current_gene, active_buttons, insertion_point, mouse):
        module_nr = insertion_point.insertion_point

        current_gene.add_module(module_nr, 'NRPS')
        current_gene.erase()
        current_gene.insertion_point = None
        current_gene.draw(mouse)
        reset_buttons(screen, active_buttons)


class PKSModuleButton(Button):
    def __init__(self):
        position = (int(0.13 * WIDTH), int(0.92 * HEIGHT))
        dimensions = (int(0.09 * WIDTH), int(HEIGHT / 25))

        super().__init__('PKS', position, dimensions)

    def do_action(self, screen, current_gene, active_buttons, insertion_point, mouse):
        module_nr = insertion_point.insertion_point

        current_gene.add_module(module_nr, 'PKS')
        current_gene.erase()
        current_gene.insertion_point = None
        current_gene.draw(mouse)
        reset_buttons(screen, active_buttons)


class AddDomainButton(Button):
    def __init__(self):
        position = (int(0.02 * WIDTH), int(0.87 * HEIGHT))
        dimensions = (int(0.2 * WIDTH), int(HEIGHT / 25))

        super().__init__("Add domain", position, dimensions)


class ACPDomainButton(DomainButton):
    def __init__(self):
        position = (int(0.35 * WIDTH), int(0.77 * HEIGHT))

        super().__init__(position, 'ACP')


class PCPDomainButton(DomainButton):
    def __init__(self):
        position = (int(0.35 * WIDTH), int(0.77 * HEIGHT))
        super().__init__(position, 'PCP')


class ADomainButton(DomainButton):
    def __init__(self):
        position = (int(0.30 * WIDTH), int(0.77 * HEIGHT))
        super().__init__(position, 'A')


class CDomainButton(DomainButton):
    def __init__(self):
        position = (int(0.25 * WIDTH), int(0.77 * HEIGHT))
        super().__init__(position, 'C')


class KSDomainButton(DomainButton):
    def __init__(self):
        position = (int(0.25 * WIDTH), int(0.77 * HEIGHT))
        super().__init__(position, 'KS')


class KRDomainButton(DomainButton):
    def __init__(self):
        position = (int(0.25 * WIDTH), int(0.87 * HEIGHT))
        super().__init__(position, 'KR')


class DHDomainButton(DomainButton):
    def __init__(self):
        position = (int(0.30 * WIDTH), int(0.87 * HEIGHT))
        super().__init__(position, 'DH')


class ERDomainButton(DomainButton):
    def __init__(self):
        position = (int(0.35 * WIDTH), int(0.87 * HEIGHT))
        super().__init__(position, 'ER')


class EDomainButton(DomainButton):
    def __init__(self):
        position = (int(0.25 * WIDTH), int(0.87 * HEIGHT))
        super().__init__(position, 'E')


class NMTDomainButton(DomainButton):
    def __init__(self):
        position = (int(0.30 * WIDTH), int(0.87 * HEIGHT))
        super().__init__(position, 'nMT')


class ATDomainButton(DomainButton):
    def __init__(self):
        position = (int(0.30 * WIDTH), int(0.77 * HEIGHT))
        super().__init__(position, 'AT')


class RemoveDomainButton(Button):
    def __init__(self):
        position = (int(0.02 * WIDTH), int(0.82 * HEIGHT))
        dimensions = (int(0.2 * WIDTH), int(HEIGHT / 25))

        super().__init__("Remove domain", position, dimensions)

    def do_action(self, current_domain, screen, active_buttons, mouse):
        gene = current_domain.module.gene
        current_domain.module.remove_domain(current_domain)
        gene.erase()
        gene.draw(mouse)
        reset_buttons(screen, active_buttons)


class SelectSubstrateButton(Button):
    def __init__(self):
        position = (int(0.75 * WIDTH), int(0.77 * HEIGHT))
        dimensions = (int(0.2 * WIDTH), int(HEIGHT / 25))

        super().__init__("Select substrate", position, dimensions)


class SelectDomainTypeButton(Button):
    def __init__(self):
        position = (int(0.75 * WIDTH), int(0.82 * HEIGHT))
        dimensions = (int(0.2 * WIDTH), int(HEIGHT / 25))

        super().__init__("Select subtype", position, dimensions)


class SetDomainInactiveButton(Button):

    def __init__(self):
        position = (int(0.75 * WIDTH), int(0.87 * HEIGHT))
        dimensions = (int(0.2 * WIDTH), int(HEIGHT / 25))

        super().__init__("Set to inactive", position, dimensions)


class SubstrateButton(Button):
    def __init__(self, substrate, position, dimensions):

        super().__init__(substrate, position, dimensions)
        self.substrate = substrate
        self.image = substrate.image
        self.highlight_image = substrate.highlight_image

    def draw(self, screen):
        substrate_image = pygame.image.load(self.image)
        substrate_image_scaled = pygame.transform.smoothscale(substrate_image, (SUBSTRATE_BUTTON_SIZE,
                                                                                SUBSTRATE_BUTTON_SIZE))
        screen.blit(substrate_image_scaled, self.rectangle)

    def highlight(self, screen):
        substrate_image = pygame.image.load(self.highlight_image)
        substrate_image_scaled = pygame.transform.smoothscale(substrate_image, (SUBSTRATE_BUTTON_SIZE,
                                                                                SUBSTRATE_BUTTON_SIZE))
        screen.blit(substrate_image_scaled, self.rectangle)

    def do_action(self):
        pass


class SubstrateSupergroupButton(Button):
    def __init__(self, text, position, group):
        dimensions = (int(0.2 * HEIGHT), int(HEIGHT / 25))
        self.group = group

        super().__init__(text, position, dimensions)


class ProteinogenicButton(SubstrateSupergroupButton):
    def __init__(self):
        position = (int(0.30 * WIDTH), int(0.77 * HEIGHT))

        super().__init__("Proteinogenic", position,
                         PROTEINOGENIC_SUBSTRATES)


class NonProteinogenicButton(SubstrateSupergroupButton):
    def __init__(self):
        position = (int(0.30 * WIDTH), int(0.82 * HEIGHT))

        super().__init__("Non-proteinogenic", position,
                         NONPROTEINOGENIC_SUBSTRATES)


class NonAminoAcidButton(SubstrateSupergroupButton):
    def __init__(self):
        position = (int(0.30 * WIDTH), int(0.87 * HEIGHT))

        super().__init__("Non-amino acids", position,
                         NON_AMINO_ACIDS)


class FattyAcidButton(SubstrateSupergroupButton):
    def __init__(self):
        position = (int(0.30 * WIDTH), int(0.87 * HEIGHT))

        super().__init__("Fatty acids", position,
                         FATTY_ACIDS)


CREATE_GENE_BUTTON = CreateGeneButton()
ADD_GENE_BUTTON = AddGeneButton()
REMOVE_GENE_BUTTON = RemoveGeneButton()

ADD_MODULE_BUTTON = AddModuleButton()
REMOVE_MODULE_BUTTON = RemoveModuleButton()

NRPS_MODULE_BUTTON = NRPSModuleButton()
PKS_MODULE_BUTTON = PKSModuleButton()

ADD_DOMAIN_BUTTON = AddDomainButton()
REMOVE_DOMAIN_BUTTON = RemoveDomainButton()

SELECT_SUBSTRATE_BUTTON = SelectSubstrateButton()
SELECT_DOMAIN_TYPE_BUTTON = SelectDomainTypeButton()
SET_DOMAIN_INACTIVE_BUTTON = SetDomainInactiveButton()


KS_DOMAIN_BUTTON = KSDomainButton()
AT_DOMAIN_BUTTON = ATDomainButton()
ACP_DOMAIN_BUTTON = ACPDomainButton()
DH_DOMAIN_BUTTON = DHDomainButton()
ER_DOMAIN_BUTTON = ERDomainButton()
KR_DOMAIN_BUTTON = KRDomainButton()

C_DOMAIN_BUTTON = CDomainButton()
A_DOMAIN_BUTTON = ADomainButton()
PCP_DOMAIN_BUTTON = PCPDomainButton()
E_DOMAIN_BUTTON = EDomainButton()
NMT_DOMAIN_BUTTON = NMTDomainButton()

PKS_DOMAIN_TO_BUTTON = {'KS': KS_DOMAIN_BUTTON,
                        'AT': AT_DOMAIN_BUTTON,
                        'DH': DH_DOMAIN_BUTTON,
                        'ER': ER_DOMAIN_BUTTON,
                        'KR': KR_DOMAIN_BUTTON,
                        'ACP': ACP_DOMAIN_BUTTON}

NRPS_DOMAIN_TO_BUTTON = {'C': C_DOMAIN_BUTTON,
                         'A': A_DOMAIN_BUTTON,
                         'nMT': NMT_DOMAIN_BUTTON,
                         'PCP': PCP_DOMAIN_BUTTON,
                         'E': E_DOMAIN_BUTTON}

PROTEINOGENIC_BUTTON = ProteinogenicButton
NON_PROTEINOGENIC_BUTTON = NonProteinogenicButton
NON_AMINO_ACID_BUTTON = NonAminoAcidButton
FATTY_ACID_BUTTON = FattyAcidButton

ALL_BUTTONS = [CREATE_GENE_BUTTON,
               ADD_GENE_BUTTON,
               ADD_MODULE_BUTTON,
               NRPS_MODULE_BUTTON,
               PKS_MODULE_BUTTON,
               KS_DOMAIN_BUTTON,
               AT_DOMAIN_BUTTON,
               ACP_DOMAIN_BUTTON,
               DH_DOMAIN_BUTTON,
               ER_DOMAIN_BUTTON,
               KR_DOMAIN_BUTTON,
               C_DOMAIN_BUTTON,
               A_DOMAIN_BUTTON,
               PCP_DOMAIN_BUTTON,
               E_DOMAIN_BUTTON,
               NMT_DOMAIN_BUTTON,
               ADD_DOMAIN_BUTTON,
               REMOVE_DOMAIN_BUTTON,
               REMOVE_MODULE_BUTTON,
               REMOVE_GENE_BUTTON,
               SELECT_SUBSTRATE_BUTTON,
               SELECT_DOMAIN_TYPE_BUTTON,
               SET_DOMAIN_INACTIVE_BUTTON,
               PROTEINOGENIC_BUTTON,
               NON_PROTEINOGENIC_BUTTON,
               NON_AMINO_ACID_BUTTON,
               FATTY_ACID_BUTTON]


def show_domain_buttons(module, screen, active_buttons):
    reset_buttons(screen, active_buttons)
    domain_types = set()

    for domain in module.domains:
        domain_types.add(domain.type)

    if module.type == 'NRPS':

        for domain, button in NRPS_DOMAIN_TO_BUTTON.items():
            if domain not in domain_types:
                show_button(button, screen, active_buttons)

    elif module.type == 'PKS':

        for domain, button in PKS_DOMAIN_TO_BUTTON.items():
            if domain not in domain_types:
                show_button(button, screen, active_buttons)


def hide_domain_buttons(screen, active_buttons):
    for button in active_buttons:
        if type(button) == DomainButton:
            hide_button(button, screen, active_buttons)


def highlight_buttons(active_buttons, screen, mouse):
    for button in active_buttons:
        if button.rectangle.collidepoint(mouse) or button.selected:
            button.highlight(screen)
        else:
            button.draw(screen)


def get_mouse_button(active_buttons, mouse):
    for button in active_buttons:
        if button.rectangle.collidepoint(mouse):
            return button

    return None


def make_buttons(screen):
    buttons = set()
    buttons.add(CREATE_GENE_BUTTON)
    CREATE_GENE_BUTTON.draw(screen)

    return buttons


def draw_buttons(active_buttons, screen, mouse):
    for button in active_buttons:
        if button.rectangle.collidepoint(mouse) or button.selected:
            button.highlight(screen)
        else:
            button.draw(screen)


def reset_buttons(screen, active_buttons):
    hide_buttons(ALL_BUTTONS, screen, active_buttons)
    show_button(CREATE_GENE_BUTTON, screen, active_buttons)


def hide_button(button, screen, active_buttons):
    try:
        active_buttons.remove(button)
    except KeyError:
        pass

    button.erase_button(screen)


def hide_buttons(buttons, screen, active_buttons):
    for button in buttons:
        hide_button(button, screen, active_buttons)


def show_button(button, screen, active_buttons):
    button.draw(screen)
    active_buttons.add(button)


def show_buttons(buttons, screen, active_buttons):
    for button in buttons:
        show_button(button, screen, active_buttons)


def make_raster(substrate_groups, buttons):
    x_coord = int(0.25 * WIDTH)
    y_coord = int(0.77 * HEIGHT)

    for i, group in enumerate(substrate_groups):
        if i % SUBSTRATE_GROUP_BUTTONS_PER_LINE == 0:
            if i != 0:
                y_coord += int(0.05 * HEIGHT)
            x_coord = int(0.25 * WIDTH)
        else:
            x_coord += SUBSTRATE_GROUP_BUTTON_SIZE + SUBSTRATE_GROUP_BUTTON_PADDING

        button = Button(group, (x_coord, y_coord), SUBSTRATE_GROUP_BUTTON_SIZE)
        buttons.add(button)


def make_nrps_substrate_group_buttons():
    buttons = set()

    make_raster(PROTEINOGENIC_SUBSTRATES, buttons)
    make_raster(NONPROTEINOGENIC_SUBSTRATES, buttons)
    make_raster(FATTY_ACIDS, buttons)
    make_raster(NON_AMINO_ACIDS, buttons)

    return buttons


def make_nrps_substrate_buttons():
    buttons = set()

    name_to_smiles = parse_smiles(PARAS_SMILES)
    for group, substrate_names in PROTEINOGENIC_SUBSTRATES.items():
        substrates = []

        for substrate_name in substrate_names:
            smiles = name_to_smiles[substrate_name]
            substrate = Substrate(substrate_name, smiles)
            substrates.append(substrate)

        substrate_group = SubstrateGroup(group, substrates)
        for substrate in substrate_group.substrates:
            buttons.add(SubstrateButton(substrate.name,
                                        (substrate.x, substrate.y),
                                        (substrate.width, substrate.height)))

    return buttons

NRPS_SUBSTRATE_BUTTONS = make_nrps_substrate_buttons()
NRPS_SUBSTRATE_GROUP_BUTTONS = make_nrps_substrate_group_buttons()

ALL_BUTTONS += NRPS_SUBSTRATE_BUTTONS
ALL_BUTTONS += NRPS_SUBSTRATE_GROUP_BUTTONS
