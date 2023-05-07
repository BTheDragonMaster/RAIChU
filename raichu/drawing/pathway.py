def draw_line(start_x, end_x, start_y, end_y):
    line = f"""<line x1="{start_x}" x2="{end_x}" y1="{start_y}" y2="{end_y}" stroke-width="2" stroke="black" />"""
    return line


def draw_arrow(start_x, end_x, y, arrow_nr):
    svg = ''
    line_1 = draw_line(start_x, end_x, y, y)
    line_2 = draw_line(end_x - 5, end_x, y - 5, y)
    line_3 = draw_line(end_x - 5, end_x, y + 5, y)

    svg += f"""<g id="arrow_{arrow_nr}">\n"""
    svg += f"{line_1}\n"
    svg += f"{line_2}\n"
    svg += f"{line_3}\n"

    svg += "</g>\n"

    return svg


def make_reaction_text(start_x, end_x, y, text):

    text_x = start_x + (end_x - start_x)/2
    text = f"""<text x="{text_x}" y="{y}" fill="'black" text-anchor="middle" font-family="verdana" font-size = "{12}">\
<tspan y="{y}" dy="0.35em">{text}</tspan></text>"""

    return text


def draw_arrow_and_text(arrow_start, arrow_end, arrow_height, arrow_nr, text):
    arrow_svg = draw_arrow(arrow_start, arrow_end, arrow_height, arrow_nr)
    reaction_text = make_reaction_text(arrow_start, arrow_end, arrow_height - 11, text)

    svg = ''

    svg += f"""<g id="labelled_arrow_{arrow_nr}">\n"""
    svg += arrow_svg
    svg += f"{reaction_text}\n"
    svg += "</g>\n"

    return svg

