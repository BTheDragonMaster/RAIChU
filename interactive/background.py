import pygame

from interactive.style import BUTTON_PANEL_COLOUR, BUTTON_PANEL, BLACK


class ButtonPanel:
    def __init__(self, screen):
        self.screen = screen

    def draw(self):
        pygame.draw.rect(self.screen, BUTTON_PANEL_COLOUR, BUTTON_PANEL)
        pygame.draw.line(self.screen, BLACK, (BUTTON_PANEL.x, BUTTON_PANEL.y),
                         (BUTTON_PANEL.x + BUTTON_PANEL.width, BUTTON_PANEL.y))