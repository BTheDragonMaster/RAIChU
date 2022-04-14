import os
import interactive.images.domains
from interactive.style import DOMAIN_SIZE, MODULE_SPACING

import pygame

DOMAIN_IMAGE_DIR = os.path.dirname(interactive.images.domains.__file__)


class Domain:
    def __init__(self, screen, domain_type, module, domain_id):
        self.screen = screen
        self.module = module
        self.type = domain_type
        self.image = os.path.join(DOMAIN_IMAGE_DIR, f"{self.type.lower()}.png")
        self.highlight_image = os.path.join(DOMAIN_IMAGE_DIR, f"{self.type.lower()}_highlight.png")
        self.subtype = None
        self.id = domain_id
        self.rectangle = None

        self.x = 0
        self.y = 0

        self.width = DOMAIN_SIZE
        self.height = DOMAIN_SIZE

        self.set_rectangle()

        self.selected = False

    def set_rectangle(self):

        x = self.module.x + MODULE_SPACING / 2 + self.id * DOMAIN_SIZE
        y = self.module.y
        self.rectangle = pygame.Rect(x, y, self.width, self.height)

    def set_domain_subtype(self, subtype):
        self.subtype = subtype

    def draw(self, mouse):
        self.set_rectangle()

        if not self.selected and not self.rectangle.collidepoint(mouse):
            domain_image = pygame.image.load(self.image)
        else:
            domain_image = pygame.image.load(self.highlight_image)

        domain_image_scaled = pygame.transform.smoothscale(domain_image, (DOMAIN_SIZE, DOMAIN_SIZE))
        self.screen.blit(domain_image_scaled, self.rectangle)



