from matplotlib.patches import Circle
from matplotlib.widgets import  Button
from pks_modules_to_structure import *
import matplotlib.image as mpimg

#Colour dicts to match antiSMASH domain colouring
colour_fill_dict = {'ACP':'#81bef7', 'AT':'#f78181', 'KS':'#81f781', \
                    'KR':'#80f680', 'DH':'#f7be81', 'ER':'#81f7f3', \
                    'TE':'#f5c4f2'}
colour_outline_dict = {'ACP':'#3d79d6', 'AT':'#df5d5d', 'KS':'#5fc65f', \
                       'KR':'#5fbb87', 'DH':'#ca9862', 'ER':'#61bbad',  \
                       'TE':'#a25ba0'}
font_domains = {'family' : 'verdana', 'size' : 12}
font_modules = {'family' : 'verdana', 'size' : 15}

#Define some global variables that I had to use because matplotlib is not
#necessarily the greatest package to make a GUI
images = []
filenames_list = []
global_figure = ''
global_nr_elongation_modules = 0

def draw_pks_cluster(pks_cluster, interactive=False):
    """

    pks_cluster: [[PKS module]] representation of the PKS cluster as:
    Starter module: ['module name','starter_module','SMILES starter unit']
    Elongation modules: ['module name', 'elongation_module',
    'malonylcoa'/'methylmalonylcoa', ['KR', 'DH', 'ER']]
    """

    #Draw (don't show!) and save png's of quick reaction mechanisms per module
    if interactive:
        pks_cluster_to_structure(pks_cluster, attach_to_acp=True, \
                                visualization_mechanism=True, \
                                draw_mechanism_per_module=True)


    #Start length is length at ends
    length_line = 6
    nr_modules = len(pks_cluster)
    #Account for space between modules (8)
    length_line += ((nr_modules - 1) * 10)
    list_all_domains = []
    elongation_modules_with_mechanisms = []
    for module in pks_cluster:
        module_name = module[0]
        module_list_domains = [module[0]]
        module_type = module[1]
        if module_type == 'starter_module':
            module_list_domains += ['ACP']
            length_line += 8
        elif module_type == 'elongation_module':
            elongation_modules_with_mechanisms.append([module_name, \
            f'{module_name}_quick_mechanism.png'])
            module_list_domains += ['KS', 'AT', 'ACP']
            length_line += 26
            for tailoring_domain in module[3]:
                length_line += 8
                module_list_domains.insert(-1, tailoring_domain)
        list_all_domains += [module_list_domains]


    #Make fig
    fig, ax = plt.subplots(figsize=(20, 10))
    thismanager = plt.get_current_fig_manager()
    thismanager.window.wm_iconbitmap("raichu_r_icon.ico")
    thismanager.set_window_title('RAIChU - Visualization PKS cluster')
    global global_figure
    global_figure = fig


    #Draw horizontal line
    ax.plot([0,length_line],[3,3], color = 'black', zorder = 1)
    #Add text above buttons:
    ax.text((length_line / 2), -40, \
            'Click module name to view reaction mechanism', ha = 'center', \
            fontdict = font_modules)
    #Add title
    plt.title('Visualization PKS cluster')
    #Draw domains per module
    x = 6
    for module in list_all_domains:
        module_name = module[0]
        del module[0]
        for domain in module:
            domain_circle = make_circle(x, domain)
            ax.add_patch(domain_circle)
            plt.text(x, 3, domain, ha = 'center', va = 'center', \
            fontdict = font_domains)
            if domain == 'ACP':
                plt.plot([x, x], [-3, 3], color = 'black', zorder = 1)
            if domain == module[0]:
                x_min = x
            if domain == module[-1]:
                x_max = x
            x += 9
        x_module_name = (x_min + x_max) / 2
        plt.text(x_module_name, 15, module_name, ha = 'center', va = 'center',\
        fontdict = font_modules)
        x += 8
    plt.axis('equal')
    plt.axis('off')
    fig.tight_layout()




    #Add buttons to view reaction mechanisms per module if interactive
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
            rel_width_buttons = 1 / (nr_elongation_modules)
            ax_button = plt.axes([x_bottomleft, 0.1, rel_width_buttons, 0.075], \
            anchor = 'C')
            module_button = Button(ax_button, label = module_name, color='#b3d8fa'\
            , hovercolor='#74abde')
            module_button.label.set_fontsize(15)
            module_button.label.set_family('verdana')
            module_button.on_clicked(button_action)
            buttons.append(module_button)
            x_bottomleft += (1 / (nr_elongation_modules))



    #Show plot
    plt.show()


    #Delete image files of quick reaction mechanisms
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
    #Get the max x-coordinate of each button, add to list
    global global_figure
    size = global_figure.get_size_inches() * global_figure.dpi  #size in pixels
    list_max_x_button = []
    max_x = 0
    global global_nr_elongation_modules
    for i in range(global_nr_elongation_modules):
        max_x += ((1/global_nr_elongation_modules) * (size[0]))
        list_max_x_button.append(max_x)
    #If the mouse clicks the button within the x-coordinates of the button,
    #That means that particular button is being pressed
    for i in range(len(list_max_x_button)):
        if event.x < list_max_x_button[i]:
            module_nr = i
            break
    #Connect the module nr to the right filename of the quick reactio mechanism
    #and open that image
    global filenames_list
    img = mpimg.imread(filenames_list[module_nr])
    images.append(img)
    #Open a new matplotlib window and display the quick reaction mechanism
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


def make_circle(x_coord, domain_type):
    """Easy function to draw circle for the domain visualization. Returns
    matplotlib.patches.circle object of a circle with radius 4

    x_coord: int, x-coordinate of the center of the circle to be drawn
    domain_type: str, PKS domain type (ACP, KS, AT, KR, DH or ER)
    """
    circle = Circle((x_coord, 3.0), 4, facecolor=colour_fill_dict[domain_type]\
    , edgecolor = colour_outline_dict[domain_type], zorder = 2)
    return circle



if __name__ == "__main__":
    ex_module = [['module_1', 'starter_module', 'SC(CC(=O)O)=O'],
                 ['module_2', 'elongation_module', 'methoxymalonylacp', ['KR', 'DH']],
                 ['module_3', 'elongation_module', 'malonylcoa', ['KR']],
                 ['module_4', 'elongation_module', 'methylmalonylcoa', ['KR', 'DH', 'ER']],
                 ['module_5', 'elongation_module', 'methoxymalonylacp', []],
                 ['module_6', 'elongation_module', 'methylmalonylcoa', ['KR', 'DH']]]
    draw_pks_cluster(ex_module, interactive=True)

