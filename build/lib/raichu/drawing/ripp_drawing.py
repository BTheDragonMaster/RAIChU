from raichu.drawing.colours import AMINO_ACID_FILL_COLOURS, AMINO_ACID_OUTLINE_COLOURS


def make_circle(x_coord, y_coord, size, amino_acid):
    """Easy function to draw circle for the domain visualization. Returns
    matplotlib.patches.circle object
    x_coord: int, x-coordinate of the center of the circle to be drawn
    amino_acid : amino acid to be visualized
    size: radius of circle
    """
    colour = AMINO_ACID_FILL_COLOURS[amino_acid]
    outline_colour = AMINO_ACID_OUTLINE_COLOURS[amino_acid]
    circle = f"""<circle r="{size}" cx="{int(x_coord)}" cy="{int(y_coord)}" stroke="{outline_colour}"
        fill="{colour}" stroke-width="1"/>"""
    return circle