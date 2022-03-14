from matplotlib.patches import Circle
from matplotlib.widgets import Button
from pks_thioesterase_reactions import *
import matplotlib.lines as lines
import matplotlib.image as mpimg
from copy import deepcopy

# Colour dicts to match antiSMASH domain colouring
colour_fill_dict = {'ACP':'#81bef7', 'AT':'#f78181', 'KS':'#81f781', \
                    'KR':'#80f680', 'DH':'#f7be81', 'ER':'#81f7f3', \
                    'TE':'#f5c4f2', 'KR*':'#80f680', 'C':'#8181f7', \
                    'A':'#bc7ff5', 'PCP':'#81bef7'}
colour_outline_dict = {'ACP':'#3d79d6', 'AT':'#df5d5d', 'KS':'#5fc65f', \
                       'KR':'#5fbb87', 'DH':'#ca9862', 'ER':'#61bbad',  \
                       'TE':'#a25ba0', 'KR*':'#5fbb87', 'C':'#5858b6',\
                       'A':'#74399b', 'PCP':'#306dd2'}



# Define some global variables that I had to use because matplotlib is not
# necessarily the greatest package to make a GUI
images = []
filenames_list = []
global_figure = ''
global_nr_elongation_modules = 0
global_final_polyketide_Drawer_object = 0

def draw_pks_cluster(pks_cluster, interactive=False, save_fig = False):
    """
    Displays a visualization of the module- and domain architecture of the
    input PKS cluster

    pks_cluster: [[PKS module]] representation of the PKS cluster as:
    Starter module: ['module name','starter_module','SMILES starter unit']
    Elongation modules: ['module name', 'elongation_module',
    'malonylcoa'/'methylmalonylcoa', ['KR', 'DH', 'ER']]

    If interactive=True, buttons are added to the canvas for each elongation
    module, which when pressed cause a new window to pop up displaying the
    reaction mechanism carried out by the respective elongation module
    """

    #Draw (don't show!) and save png's of quick reaction mechanisms per module
    global global_final_polyketide_Drawer_object
    if interactive:
        global_final_polyketide_Drawer_object = \
            pks_cluster_to_structure(pks_cluster, attach_to_acp=True,
                                     visualization_mechanism=True,
                                     draw_mechanism_per_module=True)

    if save_fig:
        filename = save_fig

    # Save (don't show!) drawings of the chain intermediate per module
    list_drawings_per_module = pks_cluster_to_structure(pks_cluster,
    attach_to_acp=True, draw_structures_per_module=True)

    # Close all matplotlib windows that were still open when generating
    # the chain intermediate Drawer objects
    plt.close('all')

    # Reset colour of all atoms to black and remove 'ACP' from drawing
    for drawing_list in list_drawings_per_module:
        for drawing in drawing_list:
            for atom in drawing.structure.graph:
                atom.draw.colour = 'black'
                if hasattr(atom, 'domain_type'):
                    atom.domain_type = ''

    # Build list of all modules comprised in the cluster, used to draw
    # module/domain architecture later
    list_all_domains = []
    elongation_modules_with_mechanisms = []
    for module in pks_cluster:
        module_name = module[0]
        module_list_domains = [module[0]]
        module_type = module[1]
        if module_type == 'starter_module':
            module_list_domains += ['AT', 'ACP']
        elif module_type == 'starter_module_nrps':
            module_list_domains += ['A', 'PCP']
        elif module_type == 'elongation_module' or \
                module_type == 'terminator_module':
            elongation_modules_with_mechanisms.append([module_name, \
            f'{module_name}_quick_mechanism.png'])
            module_list_domains += ['KS', 'AT', 'ACP']
            for tailoring_domain in module[3]:
                if tailoring_domain.startswith('KR') and \
                        tailoring_domain != 'KR_inactive':
                    tailoring_domain = 'KR'
                elif tailoring_domain == 'KR_inactive':
                    tailoring_domain = 'KR*'
                module_list_domains.insert(-1, tailoring_domain)
            if module_type == 'terminator_module':
                module_list_domains.append('TE')
        elif module_type == 'elongation_module_nrps' or\
                module_type == 'terminator_module_nrps':
            elongation_modules_with_mechanisms.append([module_name, \
            f'{module_name}_quick_mechanism.png'])
            module_list_domains += ['C', 'A', 'PCP']
            if module_type == 'terminator_module_nrps':
                module_list_domains.append('TE')
        list_all_domains += [module_list_domains]

    # Define height of window, based on max y coordinates largest structure
    last_drawing = list_drawings_per_module[-1][0]
    min_y = 100000000
    max_y = -100000000
    for atom in last_drawing.structure.graph:
        if atom.draw.position.y < min_y:
            min_y = atom.draw.position.y
        if atom.draw.position.y > max_y:
            max_y = atom.draw.position.y
    delta_y = max_y - min_y

    # Find line length to define the width of the window
    x = 30
    index = 0
    list_all_domains_copy = deepcopy(list_all_domains)
    for module in list_all_domains_copy:
        del module[0]
        for domain in module:
            if domain == 'ACP':
                index += 1
            if domain == 'PCP':
                index += 1
            if domain == module[0]:
                x_min = x
            if domain == module[-1]:
                x_max = x
                if domain == module[-1]:
                    length_line = x_max + 30
            x += 60
        x += 30

    # Make fig
    fig, ax = plt.subplots(figsize=((length_line / 70), (delta_y / 27)))
    ax.set_aspect('equal', adjustable='box')
    thismanager = plt.get_current_fig_manager()
    thismanager.window.wm_iconbitmap("raichu_r_icon.ico")
    thismanager.set_window_title('RAIChU - Visualization PKS/NRPS cluster')
    global global_figure
    global_figure = fig

    #Draw domains per module
    domain_text = []
    x = 30
    index = 0
    for module in list_all_domains:
        module_name = module[0]
        del module[0]
        for domain in module:
            domain_circle = make_circle(x, domain)
            ax.add_patch(domain_circle)
            domain_txt = domain
            if domain == 'KR*':
                domain_txt = 'KR'
                # If the KR domain is inactive, draw cross through domain
                # in interactive mode
                first_line = lines.Line2D([x-9.899494937, x+9.899494937],
                                          [-9.899494937, 9.899494937],
                                          axes = ax, color = 'mediumseagreen')
                second_line = lines.Line2D([x - 9.899494937, x + 9.899494937],
                                           [9.899494937, -9.899494937],
                                           axes = ax, color = 'mediumseagreen')
                ax.add_line(first_line)
                ax.add_line(second_line)
            domain_text.append([domain_txt, x])
            if domain == 'ACP':
                list_drawings_per_module[index].append(x)
                index += 1
            if domain == 'PCP':
                list_drawings_per_module[index].append(x)
                index += 1
            if domain == module[0]:
                x_min = x
            if domain == module[-1]:
                x_max = x
                if domain == module[-1]:
                    length_line = x_max + 30
            x += 30
        x_module_name = (x_min + x_max) / 2
        module_txt_size = 15
        font_modules = {'family': 'verdana', 'size': module_txt_size}
        plt.text(x_module_name, 30, module_name, ha = 'center', va = 'center',\
        fontdict = font_modules)
        x += 60

    plt.axis('equal')
    plt.axis('off')
    fig.tight_layout()

    # Draw horizontal line
    ax.plot([0, length_line], [0, 0], color='black', zorder=1)

    # Add module names
    domain_txt_size = (13)
    for domain_x in domain_text:
        domain, x = domain_x
        font_domains = {'family': 'verdana', 'size': domain_txt_size}
        plt.text(x, 0, domain, ha='center', va='center', \
                 fontdict=font_domains)

    # Change coordinates of all atoms of all structures to match cluster drawing
    list_drawings_correct_coord = []
    for drawing_coord in list_drawings_per_module:
        drawer_obj, x_coord = drawing_coord
        set_domain_to_origin(drawer_obj)
        push_drawing_to_right(drawer_obj, x_coord)
        list_drawings_correct_coord.append(drawer_obj)

    # Add molecules to pks cluster visualization
    height = 400
    draw_structures(list_drawings_correct_coord, fig, ax, height)

    # Add buttons to view reaction mechanisms per module if interactive
    if interactive:
        buttons = []
        nr_elongation_modules = len(elongation_modules_with_mechanisms)
        global global_nr_elongation_modules
        global_nr_elongation_modules = nr_elongation_modules
        x_bottomleft = 0
        for elongation_module in elongation_modules_with_mechanisms:
            module_name, mechanism_filename = elongation_module
            global current_filename
            current_filename = mechanism_filename
            global filenames_list
            filenames_list.append(mechanism_filename)
            rel_width_buttons = 1 / (nr_elongation_modules + 1)
            ax_button = plt.axes([x_bottomleft, 0, rel_width_buttons,\
            0.075], anchor = 'C')
            module_button = Button(ax_button, label = module_name, \
            color='#b3d8fa', hovercolor='#74abde')
            module_button.label.set_fontsize(15)
            module_button.label.set_family('verdana')
            module_button.on_clicked(button_action)
            buttons.append(module_button)
            x_bottomleft += (1 / (nr_elongation_modules + 1))

        # Button for thioesterase reaction products
        ax_button_thioesterase = plt.axes([x_bottomleft, 0, rel_width_buttons, 0.075], anchor='C')
        thioesterase_button = Button(ax_button_thioesterase, label='View thioesterase\nreaction products', color='#b3d8fa', hovercolor='#74abde')
        thioesterase_button.label.set_fontsize(15)
        thioesterase_button.label.set_family('verdana')
        thioesterase_button.on_clicked(button_action_thioesterase)
        buttons.append(thioesterase_button)

    # Show plot, or directly save drawing if save_fig argument is given
    if not save_fig:
        plt.show()
    else:
        if filename.endswith('.png'):
            filename = filename
        else:
            filename = filename + '.png'
        plt.savefig(filename)

    # Delete image files of quick reaction mechanisms
    for module in pks_cluster:
        module_name = module[0]
        if path.exists(f'{module_name}_quick_mechanism.png'):
            os.remove(f'{module_name}_quick_mechanism.png')


def button_action(event):
    """Accessory function that is called when one of the module buttons is
    pressed. It opens another matplotlib window and displays the quick reaction
    mechanism of that particular module

    This function cannot take any arguments, so I had to go use a lot of global
    variables to circumvent that
    """
    # Get the max x-coordinate of each button, add to list
    global global_figure
    size = global_figure.get_size_inches() * global_figure.dpi  #size in pixels
    list_max_x_button = []
    max_x = 0
    global global_nr_elongation_modules
    for i in range(global_nr_elongation_modules):
        max_x += ((1/(global_nr_elongation_modules + 1)) * (size[0]))
        list_max_x_button.append(max_x)

    # If the mouse clicks the button within the x-coordinates of the button,
    # that means that particular button is being pressed
    for i in range(len(list_max_x_button)):
        if event.x < list_max_x_button[i]:
            module_nr = i
            break

    # Connect the module nr to the right filename of the quick reaction
    # mechanism and open that image
    global filenames_list
    img = mpimg.imread(filenames_list[module_nr])
    images.append(img)

    # Open a new matplotlib window and display the quick reaction mechanism
    fig2, ax2 = plt.subplots(figsize=(20, 10))
    ax2.imshow(img)
    plt.axis('equal')
    plt.axis('off')
    thismanager = plt.get_current_fig_manager()
    thismanager.window.wm_iconbitmap("raichu_r_icon.ico")
    window_title = filenames_list[module_nr][:-4]
    thismanager.set_window_title(window_title)
    fig2.tight_layout()
    plt.show()

def button_action_thioesterase(event):
    """Accessory function that is called when the 'View thioesterase reaction
     products' is pressed. It opens  matplotlib windows displaying each of the
     thioesterase reaction products

    This function cannot take any arguments, so I had to go use a  global
    variable to circumvent that
    """
    # Fetch the structure of the molecule tethered to the last ACP/PCP domain
    global global_final_polyketide_Drawer_object
    final_polyketide = global_final_polyketide_Drawer_object

    # If the last elongation module was an NRPS module, attach the product to
    # the domain
    for atom in final_polyketide.graph:
        atom.hybridise()
    final_polyketide.refresh_structure()
    final_polyketide.set_connectivities()
    final_polyketide.set_atom_neighbours()
    final_polyketide.find_cycles()

    # Pass Structure to the thioesterase function, this function will open
    # a separate window for each TE product
    thioesterase_all_products(final_polyketide)

def make_circle(x_coord, domain_type):
    """Easy function to draw circle for the domain visualization. Returns
    matplotlib.patches.circle object of a circle with radius 4

    x_coord: int, x-coordinate of the center of the circle to be drawn
    domain_type: str, PKS domain type (ACP, KS, AT, KR, DH or ER)
    """
    circle = Circle((x_coord, 0.0), 14,facecolor=colour_fill_dict[domain_type]\
    , edgecolor = colour_outline_dict[domain_type], zorder = 2)
    return circle

def set_domain_to_origin(drawer_object):
    """Accessory function that alters the coordinates of all atoms in the input
    Drawer object as such that the ACP/PCP domain has coordinates [0,0], after
    which it will return the Drawer object

    drawer_object: PIKAChU Drawer object, structure that needs to be moved as
    such that the ACP/PCP domain is located in the origin of the canvas
    """
    for atom in drawer_object.structure.graph:
        if hasattr(atom, 'domain_type'):
            domain = atom
            domain_x = atom.draw.position.x
            domain_y = atom.draw.position.y
            atom.draw.position.x = 0
            atom.draw.position.y = -11
    for atom in drawer_object.structure.graph:
        if atom != domain:
            atom.draw.position.x -= domain_x
            atom.draw.position.y -= domain_y
            atom.draw.position.y -= 11
    new_drawer_object = drawer_object

    return new_drawer_object

def push_drawing_to_right(drawer_object, shift_to_right):
    """Accessory function that pushes the input Drawer object to the right by
    changing the x coordinate of each atom according to the shift_to_right
    argument

    drawer_object: PIKAChU Drawer object, structure that needs to be moved to
    the right
    shift_to_right: float, magnitude of the shift to the right
    """
    for atom in drawer_object.structure.graph:
        atom.draw.position.x += shift_to_right
    return drawer_object



def draw_structures(drawer_objects, fig, ax, height):
    """Copied from the PIKAChU drawing.py script, used to draw the reaction
    intermediates in the canvas

    drawer_objects: list of PIKAChU Drawer objects
    fig, ax: matplotlib canvas to draw structures in
    height: height of the canvas
    """
    min_x = 100000000
    max_x = -100000000
    min_y = 100000000
    max_y = -100000000
    for i in range(len(drawer_objects)):
        drawer_object = drawer_objects[i]
        #  fig, ax = plt.subplots()
        ax.set_aspect('equal', adjustable='box')
        ax.axis('off')

        font_size = 3500 / height
        drawer_object.line_width = 600 / height


        plt.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0,
                            hspace=0)

        #  figure, ax = plt.subplots(figsize=(8, 8))
        # ax.patch.set_face_color(drawer_object.options.background_color)
        #  ax.set_aspect()('equal', adjustable='box')
        #  plt.gca().set_aspect('equal', adjustable='box')

        params = {'mathtext.default': 'regular',
                  'font.size': font_size}
        plt.rcParams.update(params)

        ring_centers_x = []
        ring_centers_y = []

        for ring in drawer_object.rings:
            drawer_object.set_ring_center(ring)

            ring_centers_x.append(ring.center.x)
            ring_centers_y.append(ring.center.y)

        #   ax.scatter(ring_centers_x, ring_centers_y, color='blue')

        for bond_nr, bond in drawer_object.structure.bonds.items():
            if bond.atom_1.draw.positioned and bond.atom_2.draw.positioned:
                line = Line(bond.atom_1.draw.position,
                            bond.atom_2.draw.position, bond.atom_1,
                            bond.atom_2)
                midpoint = line.get_midpoint()
                truncated_line = line.get_truncated_line(
                    drawer_object.options.short_bond_length)
                if bond.type == 'single':
                    if bond in drawer_object.chiral_bonds:
                        orientation, chiral_center = \
                            drawer_object.chiral_bond_to_orientation[bond]
                        drawer_object.plot_chiral_bond(orientation,
                                                       chiral_center, line,
                                                       ax, midpoint)
                    else:
                        if (bond.atom_1.type == 'S' and hasattr(bond.atom_2, 'domain_type') or (bond.atom_2.type == 'S' and hasattr(bond.atom_1, 'domain_type'))):
                            drawer_object.plot_halflines_s_domain(line, ax, midpoint)
                        else:
                            drawer_object.plot_halflines(line, ax, midpoint)

                elif bond.type == 'double':
                    if not drawer_object.is_terminal(
                            bond.atom_1) and not drawer_object.is_terminal(
                        bond.atom_2):
                        drawer_object.plot_halflines(line, ax, midpoint)

                        common_ring_numbers = drawer_object.get_common_rings(
                            bond.atom_1, bond.atom_2)

                        if common_ring_numbers:
                            common_rings = []
                            for ring_nr in common_ring_numbers:
                                common_rings.append(
                                    drawer_object.get_ring(ring_nr))

                            common_rings.sort(key=lambda x: len(x.members))
                            common_ring = common_rings[0]
                            ring_centre = common_ring.center
                            second_line = line.double_line_towards_center(
                                ring_centre,
                                drawer_object.options.bond_spacing,
                                drawer_object.options.double_bond_length)
                            second_line_midpoint = second_line.get_midpoint()
                            drawer_object.plot_halflines_double(second_line,
                                                    ax, second_line_midpoint)

                        else:
                            bond_neighbours = bond.atom_1.drawn_neighbours \
                                              + bond.atom_2.drawn_neighbours
                            if bond_neighbours:
                                vectors = [atom.draw.position for atom in
                                           bond_neighbours]
                                gravitational_point = Vector.get_average(
                                    vectors)
                                second_line = line.double_line_towards_center(
                                    gravitational_point,
                                    drawer_object.options.bond_spacing,
                                    drawer_object.options.double_bond_length)
                                second_line_midpoint = second_line.get_midpoint()
                                drawer_object.plot_halflines_double(
                                    second_line, ax,
                                    second_line_midpoint)
                            else:
                                print("Shouldn't happen!")
                    else:
                        if drawer_object.is_terminal(
                                bond.atom_1) and drawer_object.is_terminal(
                                bond.atom_2):
                            dummy_1 = Vector(bond.atom_1.draw.position.x + 1,
                                             bond.atom_1.draw.position.y + 1)
                            dummy_2 = Vector(bond.atom_1.draw.position.x - 1,
                                             bond.atom_1.draw.position.y - 1)
                            double_bond_line_1 = line.double_line_towards_center(
                                dummy_1,
                                drawer_object.options.bond_spacing / 2.0,
                                drawer_object.options.double_bond_length)
                            double_bond_line_1_midpoint = double_bond_line_1.get_midpoint()
                            double_bond_line_2 = line.double_line_towards_center(
                                dummy_2,
                                drawer_object.options.bond_spacing / 2.0,
                                drawer_object.options.double_bond_length)
                            double_bond_line_2_midpoint = double_bond_line_2.get_midpoint()

                            drawer_object.plot_halflines_double(
                                double_bond_line_1, ax,
                                double_bond_line_1_midpoint)
                            drawer_object.plot_halflines_double(
                                double_bond_line_2, ax,
                                double_bond_line_2_midpoint)

                        else:

                            if drawer_object.is_terminal(bond.atom_1):
                                terminal_atom = bond.atom_1
                                branched_atom = bond.atom_2
                            else:
                                terminal_atom = bond.atom_2
                                branched_atom = bond.atom_1

                            if len(branched_atom.drawn_neighbours) >= 3:
                                closest_two = drawer_object.get_sorted_distances_from_list(
                                    terminal_atom,
                                    branched_atom.drawn_neighbours)
                                closest_atom_1 = closest_two[0][1]
                                closest_atom_2 = closest_two[1][1]

                                line = Line(terminal_atom.draw.position,
                                            branched_atom.draw.position,
                                            terminal_atom, branched_atom)

                                bond_1_line = Line(branched_atom.draw.position,
                                                   closest_atom_1.draw.position,
                                                   branched_atom,
                                                   closest_atom_1)
                                bond_2_line = Line(branched_atom.draw.position,
                                                   closest_atom_2.draw.position,
                                                   branched_atom,
                                                   closest_atom_1)

                                double_bond_line_1 = line.double_line_towards_center(
                                    closest_atom_1.draw.position,
                                    drawer_object.options.bond_spacing / 2.0,
                                    drawer_object.options.double_bond_length)
                                double_bond_line_2 = line.double_line_towards_center(
                                    closest_atom_2.draw.position,
                                    drawer_object.options.bond_spacing / 2.0,
                                    drawer_object.options.double_bond_length)

                                double_bond_line_1_midpoint = double_bond_line_1.get_midpoint()
                                double_bond_line_2_midpoint = double_bond_line_2.get_midpoint()

                                intersection_1 = double_bond_line_1.find_intersection(
                                    bond_1_line)
                                intersection_2 = double_bond_line_2.find_intersection(
                                    bond_2_line)

                                if terminal_atom.draw.position.x > branched_atom.draw.position.x:
                                    double_bond_line_1.point_1 = intersection_1
                                    double_bond_line_2.point_1 = intersection_2
                                else:
                                    double_bond_line_1.point_2 = intersection_1
                                    double_bond_line_2.point_2 = intersection_2

                                drawer_object.plot_halflines(
                                    double_bond_line_1, ax,
                                    double_bond_line_1_midpoint)
                                drawer_object.plot_halflines(
                                    double_bond_line_2, ax,
                                    double_bond_line_2_midpoint)

                            else:
                                pass
        # Changed by Sophie
        for atom in drawer_object.structure.graph:
            if atom.type != 'C' and atom.draw.positioned:
                text = atom.type
                if hasattr(atom, 'domain_type'):
                    text = atom.domain_type
                horizontal_alignment = 'center'
                if atom.type == '*':
                    for neighbour in atom.neighbours:
                        if neighbour.type == 'C':
                            neighbouring_c = neighbour
                    # In order to not let the number of the sidechain overlap
                    # with the bond, move Rx symbol along hor. axis depending on
                    # if group is added to the left or right of the chain
                    delta_x_r = drawer_object.get_delta_x_sidechain(atom,
                                                           neighbouring_c)
                    atom.draw.position.x += delta_x_r
                    text = fr'$R_{atom.unknown_index}$'
                if atom.draw.has_hydrogen:
                    # if len(atom.drawn_neighbours) == 1 and atom.draw.has_hydrogen:
                    hydrogen_count = 0
                    for neighbour in atom.neighbours:
                        if neighbour.type == 'H' and not neighbour.draw.is_drawn:
                            hydrogen_count += 1

                    if hydrogen_count:

                        orientation = drawer_object.get_hydrogen_text_orientation(
                            atom)
                        if hydrogen_count > 1:
                            if orientation == 'H_before_atom':
                                text = r'$H_{hydrogens}${atom_type}'.format(
                                    hydrogens=hydrogen_count,
                                    atom_type=atom.type)
                                horizontal_alignment = 'right'
                                atom.draw.position.x += 3
                            else:
                                text = r'${atom_type}H_{hydrogens}$'.format(
                                    hydrogens=hydrogen_count,
                                    atom_type=atom.type)
                                horizontal_alignment = 'left'
                                atom.draw.position.x -= 3
                        elif hydrogen_count == 1:
                            if orientation == 'H_before_atom':
                                text = f'H{atom.type}'
                                horizontal_alignment = 'right'
                                atom.draw.position.x += 3
                            else:
                                text = f'{atom.type}H'
                                horizontal_alignment = 'left'
                                atom.draw.position.x -= 3

                plt.text(atom.draw.position.x, atom.draw.position.y,
                         text,
                         horizontalalignment=horizontal_alignment,
                         verticalalignment='center',
                         color=atom.draw.colour, zorder = 1)

