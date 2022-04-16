import os

import pygame

from interactive.buttons import make_buttons, get_mouse_button, Button, DomainButton, hide_button, \
    AddGeneButton, CreateGeneButton, show_buttons, ADD_MODULE_BUTTON, show_button, \
    AddModuleButton, NRPSModuleButton, NRPS_MODULE_BUTTON, PKSModuleButton, PKS_MODULE_BUTTON, \
    show_domain_buttons, reset_buttons, AddDomainButton, ADD_DOMAIN_BUTTON, \
    REMOVE_MODULE_BUTTON, REMOVE_DOMAIN_BUTTON, RemoveDomainButton, RemoveModuleButton, RemoveGeneButton, \
    REMOVE_GENE_BUTTON, SELECT_SUBSTRATE_BUTTON, SelectSubstrateButton, SELECT_DOMAIN_TYPE_BUTTON, \
    SelectDomainTypeButton, SET_DOMAIN_INACTIVE_BUTTON, SetDomainInactiveButton, \
    PKS_SUBSTRATE_BUTTONS, NRPS_SUPERGROUP_BUTTONS, SubstrateSupergroupButton,\
    PROTEINOGENIC_BUTTONS, NON_PROTEINOGENIC_BUTTONS, FATTY_ACID_BUTTONS, NON_AMINO_ACID_BUTTONS, \
    ProteinogenicButton, NonProteinogenicButton, FattyAcidButton, NonAminoAcidButton, \
    PROTEINOGENIC_GROUP_TO_BUTTONS, NONPROTEINOGENIC_GROUP_TO_BUTTONS, FATTY_ACID_GROUP_TO_BUTTONS, \
    NON_AMINOACID_GROUP_TO_BUTTONS, SubstrateGroupButton, SubstrateButton, RenderClusterButton, \
    KR_BUTTONS, KRSubtypeButton
from interactive.domain import Domain
from interactive.module import Module
from interactive.gene import Gene
from interactive.insertion_point import InsertionPoint
from interactive.render_cluster import render_cluster
from interactive.style import RENDER_WINDOW_SIZE, WHITE


class RaichuManager:
    states = {'add gene',
              'remove_gene',
              'add module',
              'remove module',
              'add domain',
              'remove domain',
              'set domain subtype',

              }

    def __init__(self, screen):
        self.genes = []
        self.state = None
        self.active_buttons = make_buttons(screen)
        self.screen = screen

        self.selected_module = None
        self.selected_gene = None
        self.selected_domain = None
        self.selected_substrate_group = None

        self.text_box = None
        self.insertion_point = None
        self.cluster_image = None

    def get_mouse_domain(self, mouse):
        for gene in self.genes:
            for module in gene.modules:
                for domain in module.domains:
                    if domain.rectangle.collidepoint(mouse):
                        return domain

        return None

    def get_mouse_module(self, mouse):
        for gene in self.genes:
            for module in gene.modules:
                if module.rectangle.collidepoint(mouse):
                    return module

        return None

    def get_mouse_gene(self, mouse):
        for gene in self.genes:
            if gene.rectangle.collidepoint(mouse):
                return gene

        return None

    def detect_click_location(self, mouse):

        button = get_mouse_button(self.active_buttons, mouse)
        domain = self.get_mouse_domain(mouse)
        module = self.get_mouse_module(mouse)
        gene = self.get_mouse_gene(mouse)

        if button:
            return button
        elif domain:
            self.reset_selections()
            domain.selected = True
            self.selected_domain = domain
            return domain
        elif module:
            self.reset_selections()
            module.selected = True
            self.selected_module = module
            return module
        elif gene:
            self.reset_selections()
            gene.selected = True
            self.selected_gene = gene
            return gene
        else:
            return None

    def reset_selections(self):
        if self.selected_gene:
            self.selected_gene.insertion_point = None
            self.selected_gene.selected = False
            self.selected_gene = None
        if self.selected_module:
            self.selected_module.selected = False
            self.selected_module = None
        if self.selected_domain:
            self.selected_domain.selected = False
            self.selected_domain = None
        self.selected_substrate_group = None

    def do_button_action(self, button, mouse):
        if type(button) == DomainButton or issubclass(type(button), DomainButton):
            button.do_action(self.selected_module, mouse, self.screen, self.active_buttons)
            self.selected_module.selected = False
            self.selected_module = None
            hide_button(button, self.screen, self.active_buttons)
        elif issubclass(type(button), SubstrateSupergroupButton):
            reset_buttons(self.screen, self.active_buttons)
            if type(button) == ProteinogenicButton:
                show_buttons(PROTEINOGENIC_BUTTONS, self.screen, self.active_buttons)
                self.selected_substrate_group = PROTEINOGENIC_GROUP_TO_BUTTONS
            elif type(button) == NonProteinogenicButton:
                show_buttons(NON_PROTEINOGENIC_BUTTONS, self.screen, self.active_buttons)
                self.selected_substrate_group = NONPROTEINOGENIC_GROUP_TO_BUTTONS
            elif type(button) == FattyAcidButton:
                show_buttons(FATTY_ACID_BUTTONS, self.screen, self.active_buttons)
                self.selected_substrate_group = FATTY_ACID_GROUP_TO_BUTTONS
            elif type(button) == NonAminoAcidButton:
                show_buttons(NON_AMINO_ACID_BUTTONS, self.screen, self.active_buttons)
                self.selected_substrate_group = NON_AMINOACID_GROUP_TO_BUTTONS
        elif type(button) == SubstrateGroupButton:
            reset_buttons(self.screen, self.active_buttons)
            show_buttons(self.selected_substrate_group[button.text], self.screen, self.active_buttons)
            self.selected_substrate_group = None

        elif type(button) == AddGeneButton:
            button.do_action(self.screen, self.genes, self.text_box, self.active_buttons, mouse)
            self.text_box = None
        elif type(button) == CreateGeneButton:
            self.reset_selections()
            self.text_box = button.do_action(self.screen, self.active_buttons)
        elif type(button) == AddModuleButton:
            show_buttons([PKS_MODULE_BUTTON, NRPS_MODULE_BUTTON], self.screen, self.active_buttons)
        elif type(button) == PKSModuleButton or type(button) == NRPSModuleButton:
            button.do_action(self.screen, self.selected_gene, self.active_buttons, self.insertion_point, mouse)
            self.reset_selections()
            self.insertion_point = None
        elif type(button) == RemoveGeneButton:
            button.do_action(self.screen, self.genes, self.selected_gene, self.active_buttons, mouse)
        elif type(button) == AddDomainButton:
            reset_buttons(self.screen, self.active_buttons)
            show_domain_buttons(self.selected_module, self.screen, self.active_buttons)
        elif type(button) == RemoveModuleButton:
            button.do_action(self.screen, self.selected_module, self.active_buttons, mouse)
            self.reset_selections()
        elif type(button) == RemoveDomainButton:
            button.do_action(self.selected_domain, self.screen, self.active_buttons, mouse)
            self.reset_selections()
        elif type(button) == RenderClusterButton:
            render_cluster(self.genes)
            reset_buttons(self.screen, self.active_buttons)
            self.reset_selections()
            cluster_image = pygame.image.load(os.path.join(os.getcwd(), 'cluster_test.png'))
            image_width = cluster_image.get_width()
            image_height = cluster_image.get_height()

            if image_width <= RENDER_WINDOW_SIZE[0] and image_height > RENDER_WINDOW_SIZE[1]:
                ratio = RENDER_WINDOW_SIZE[1] / image_height
                image_width = RENDER_WINDOW_SIZE[0] * ratio
                image_height = RENDER_WINDOW_SIZE[1]
            elif image_width > RENDER_WINDOW_SIZE[0] and image_height <= RENDER_WINDOW_SIZE[1]:
                ratio = RENDER_WINDOW_SIZE[0] / image_width
                image_height = RENDER_WINDOW_SIZE[1] * ratio
                image_width = RENDER_WINDOW_SIZE[0]
            elif image_width > RENDER_WINDOW_SIZE[0] and image_height > RENDER_WINDOW_SIZE[1]:
                ratio_1 = RENDER_WINDOW_SIZE[0] / image_width
                ratio_2 = RENDER_WINDOW_SIZE[1] / image_height
                if ratio_1 <= ratio_2:
                    ratio = ratio_1
                else:
                    ratio = ratio_2

                image_height *= ratio
                image_width *= ratio

            render_window_size = (image_width, image_height)

            self.cluster_image = pygame.transform.smoothscale(cluster_image, render_window_size)
        elif type(button) == SelectDomainTypeButton:
            reset_buttons(self.screen, self.active_buttons)
            show_buttons(KR_BUTTONS, self.screen, self.active_buttons)
        elif type(button) == SelectSubstrateButton:
            if self.selected_domain.type == 'A':
                reset_buttons(self.screen, self.active_buttons)
                show_buttons(NRPS_SUPERGROUP_BUTTONS, self.screen, self.active_buttons)
            elif self.selected_domain.type == 'AT':
                reset_buttons(self.screen, self.active_buttons)
                show_buttons(PKS_SUBSTRATE_BUTTONS, self.screen, self.active_buttons)
        elif type(button) == KRSubtypeButton:
            button.do_action(self.selected_domain, mouse, self.screen, self.active_buttons)
            self.reset_selections()

        elif type(button) == SubstrateButton:
            self.selected_domain.substrate = button.substrate
            self.selected_domain.module.gene.erase()
            self.selected_domain.module.gene.draw(mouse)
            reset_buttons(self.screen, self.active_buttons)
            self.reset_selections()

    def do_click_action(self, mouse):

        selected_entity = self.detect_click_location(mouse)
        if selected_entity:
            for gene in self.genes:
                gene.draw(mouse)
            if type(selected_entity) == Button or issubclass(type(selected_entity), Button):
                self.do_button_action(selected_entity, mouse)
            elif type(selected_entity) == Gene:
                reset_buttons(self.screen, self.active_buttons)
                self.selected_gene = selected_entity
                self.insertion_point = InsertionPoint(selected_entity, mouse)
                selected_entity.insertion_point = self.insertion_point
                selected_entity.draw(mouse)
                show_buttons([ADD_MODULE_BUTTON, REMOVE_GENE_BUTTON], self.screen, self.active_buttons)
            elif type(selected_entity) == Module:
                self.selected_module = selected_entity
                reset_buttons(self.screen, self.active_buttons)
                show_buttons([ADD_DOMAIN_BUTTON, REMOVE_MODULE_BUTTON], self.screen,
                             self.active_buttons)
            elif type(selected_entity) == Domain:
                self.selected_domain = selected_entity
                reset_buttons(self.screen, self.active_buttons)
                show_buttons([REMOVE_DOMAIN_BUTTON], self.screen, self.active_buttons)

                if selected_entity.type == 'A' or selected_entity.type == 'AT':
                    show_button(SELECT_SUBSTRATE_BUTTON, self.screen, self.active_buttons)
                elif selected_entity.type == 'KR':
                    show_button(SELECT_DOMAIN_TYPE_BUTTON, self.screen, self.active_buttons)

        else:
            self.random_click(mouse)

    def random_click(self, mouse):
        self.screen.fill(WHITE)

        reset_buttons(self.screen, self.active_buttons)
        for gene in self.genes:
            gene.draw(mouse)

        for gene in self.genes:
            gene.selected = False
            gene.insertion_point = None
            for module in gene.modules:
                module.selected = False
                for domain in module.domains:
                    domain.selected = False

        self.selected_gene = None
        self.selected_module = None
        self.selected_domain = None
        self.selected_substrate_group = None
        self.cluster_image = False






