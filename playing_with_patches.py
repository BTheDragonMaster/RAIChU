import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.widgets import  Button
from pks_modules_to_structure import *
import matplotlib.image as mpimg
from PIL import Image

colour_fill_dict = {'ACP':'#81bef7', 'AT':'#f78181', 'KS':'#81f781', 'KR':'#80f680', 'DH':'#f7be81', 'ER':'#81f7f3', 'TE':'#f5c4f2'}
colour_outline_dict = {'ACP':'#3d79d6', 'AT':'#df5d5d', 'KS':'#5fc65f', 'KR':'#5fbb87', 'DH':'#ca9862', 'ER':'#61bbad', 'TE':'#a25ba0'}
font_domains = {'family' : 'verdana', 'size' : 12}
font_modules = {'family' : 'verdana', 'size' : 15}

def draw_domains(pks_cluster):
    """

    """
    #Draw (don't show!) and save png's of structures after each module
    #pks_cluster_to_structure(pks_cluster, draw_structures_per_module=True, attach_to_acp=True)


    #Load in png's of polyketide structures after each module
    img = Image.open('module_2.png')


    #Start length is length at ends
    length_line = 6
    nr_modules = len(pks_cluster)
    #Account for space between modules (8)
    length_line += ((nr_modules - 1) * 10)
    list_all_domains = []
    for module in pks_cluster:
        module_list_domains = [module[0]]
        module_type = module[1]
        if module_type == 'starter_module':
            module_list_domains += ['ACP']
            length_line += 8
        elif module_type == 'elongation_module':
            module_list_domains += ['KS', 'AT', 'ACP']
            length_line += 26
            for tailoring_domain in module[3]:
                length_line += 8
                module_list_domains.insert(-1, tailoring_domain)
        list_all_domains += [module_list_domains]
    #Make fig
    fig, ax = plt.subplots(figsize=(20, 10))


    #try to add structure
    list = []
    list.append(img)
    img = img.resize((50,100))
    #plt.imshow(img)


    #Draw horizontal line
    ax.plot([0,length_line],[3,3], color = 'black', zorder = 1)
    #Draw domains per module
    x = 6
    for module in list_all_domains:
        module_name = module[0]
        del module[0]
        for domain in module:
            domain_circle = make_circle(x, domain)
            ax.add_patch(domain_circle)
            plt.text(x, 3, domain, ha = 'center', va = 'center', fontdict = font_domains)
            if domain == 'ACP':
                plt.plot([x, x], [-3, 3], color = 'black', zorder = 1)
            if domain == module[0]:
                x_min = x
            if domain == module[-1]:
                x_max = x
            x += 9
        x_module_name = (x_min + x_max) / 2
        plt.text(x_module_name, 15, module_name, ha = 'center', va = 'center', fontdict = font_modules)
        x += 8
    plt.axis('equal')
    plt.axis('off')
    fig.tight_layout()


    # axButn2 = plt.axes([0.4,0.5,0.05,0.075], anchor = 'C')
    # btn2 = Button(
    #     axButn2, label="Prev", color='pink', hovercolor='tomato')


    #Show plot
    plt.show()


def make_circle(x_coord, domain_type):
    """
    """
    circle = Circle((x_coord, 3.0), 4, facecolor = colour_fill_dict[domain_type], edgecolor = colour_outline_dict[domain_type], zorder = 2)
    return circle



if __name__ == "__main__":
    ex_module = [['module_1', 'starter_module', 'SC(CC(=O)O)=O'],
                 ['module_2', 'elongation_module', 'methoxymalonylacp', ['KR', 'DH']],
                 ['module_3', 'elongation_module', 'malonylcoa', ['KR']],
                 ['module_4', 'elongation_module', 'methylmalonylcoa', ['KR', 'DH', 'ER']],
                 ['module_5', 'elongation_module', 'methoxymalonylacp', []],
                 ['module_6', 'elongation_module', 'methylmalonylcoa', ['KR', 'DH']]]
    draw_domains(ex_module)

