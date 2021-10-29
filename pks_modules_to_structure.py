from pks_tailoring_reactions import *
from pks_elongation_reactions import *
import os
from os import path
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from copy import copy


def pks_cluster_to_structure(modules, visualization_mechanism = False, draw_structures_per_module = False, attach_to_acp = False):
    """If visualization_mechanism == False: returns the final polyketide
    structure resulting from the input PKS cluster
    If visualization_mechanism == True: returns a visualization of the elong-
    ation and/or tailoring reactions performed by each PKS module
    If attach_to_acp == True: the polyketide will be attached to an ACP domain
    object of type Domain. Can be combined with previous optional flags.
    If draw_structures_per_module: saves the polyketide intermediate as
    'module_name'.png, doesn't display any visualzisation

    modules: [[PKS module]] representation of the PKS cluster as:
    Starter module: ['module name','starter_module','SMILES starter unit']
    Elongation modules: ['module name', 'elongation_module',
    'malonylcoa'/'methylmalonylcoa', ['KR', 'DH', 'ER']]
    """
    modules = copy(modules)
    for module in modules:
        #Find starting unit SMILES, convert to structure
        if module == modules[0]:
            print(module)
            assert len(module) == 3
            assert module[1] == 'starter_module'
            starter_unit = Smiles(module[2]).smiles_to_structure()
            if attach_to_acp:
                chain_intermediate = attach_to_domain(starter_unit, 'ACP')
            else:
                chain_intermediate = starter_unit
            del modules[0]
    for module in modules:
        module_name, module_type, elongation_unit, list_domains = module
        assert elongation_unit == 'malonylcoa' or elongation_unit == 'methylmalonylcoa' or elongation_unit == 'methoxymalonylacp'
        #Reset atom colours in the first structure depicted to black
        for atom in chain_intermediate.graph:
            atom.draw.colour = 'black'
        if elongation_unit == 'malonylcoa':
            if visualization_mechanism == True:
                #If a previous process was stopped prematurely, make sure the
                #image files are removed before continuing
                if path.exists('1.png'):
                    os.remove('1.png')
                Drawer(chain_intermediate, save_png='1.png')
            new_chain_intermediate = add_malonylunit(chain_intermediate)
            chain_intermediate = new_chain_intermediate
            if visualization_mechanism == True:
                if path.exists('2.png'):
                    os.remove('2.png')
                Drawer(chain_intermediate, save_png='2.png')
        elif elongation_unit == 'methylmalonylcoa':
            if visualization_mechanism == True:
                if path.exists('1.png'):
                    os.remove('1.png')
                Drawer(chain_intermediate, save_png='1.png')
            new_chain_intermediate = add_methylmalonylunit(chain_intermediate)
            chain_intermediate = new_chain_intermediate
            if visualization_mechanism == True:
                if path.exists('2.png'):
                    os.remove('2.png')
                Drawer(chain_intermediate, save_png='2.png')
        elif elongation_unit == 'methoxymalonylacp':
            if visualization_mechanism == True:
                if path.exists('1.png'):
                    os.remove('1.png')
                Drawer(chain_intermediate, save_png='1.png')
            new_chain_intermediate = add_methoxymalonylunit(chain_intermediate)
            chain_intermediate = new_chain_intermediate
            if visualization_mechanism == True:
                if path.exists('2.png'):
                    os.remove('2.png')
                Drawer(chain_intermediate, save_png='2.png')
        #If the module does not contain any tailoring domains:
        if len(list_domains) == 0:
            if visualization_mechanism == True:
                if elongation_unit == 'methylmalonylcoa':
                    display_reactions(['1.png', '2.png'],list_domains,'methylmalonylcoa', module_name)
                elif elongation_unit == 'malonylcoa':
                    display_reactions(['1.png', '2.png'],list_domains,'malonylcoa', module_name)
                elif elongation_unit == 'methoxymalonylacp':
                    display_reactions(['1.png', '2.png'], list_domains, 'methoxymalonylacp', module_name)
        for domain in list_domains:
            if domain == 'KR':
                new_chain_intermediate = ketoreductase(chain_intermediate)
                chain_intermediate = new_chain_intermediate
                if visualization_mechanism == True:
                    if path.exists('3.png'):
                        os.remove('3.png')
                    Drawer(chain_intermediate, save_png='3.png')
                    #Check if the module also contains a DH domain
                    if 'DH' not in list_domains:
                        display_reactions(['1.png', '2.png', '3.png'], list_domains, elongation_unit, module_name)
            if domain == 'DH':
                assert 'KR' in list_domains
                new_chain_intermediate = dehydratase(chain_intermediate)
                chain_intermediate = new_chain_intermediate
                if visualization_mechanism == True:
                    if path.exists('4.png'):
                        os.remove('4.png')
                    Drawer(chain_intermediate, save_png='4.png')
                    if 'ER' not in list_domains:
                        display_reactions(['1.png', '2.png', '3.png','4.png'], list_domains, elongation_unit, module_name)
            if domain == 'ER':
                assert 'KR' in list_domains
                assert 'DH' in list_domains
                new_chain_intermediate = enoylreductase(chain_intermediate)
                chain_intermediate = new_chain_intermediate
                if visualization_mechanism == True:
                    if path.exists('5.png'):
                        os.remove('5.png')
                    Drawer(chain_intermediate, save_png='5.png')
                    display_reactions(['1.png', '2.png', '3.png','4.png','5.png'], list_domains, elongation_unit, module_name)
        if draw_structures_per_module:
            Drawer(chain_intermediate, save_png=f'{module_name}.png')
    # Reset the atom color in the final structure to black
    for atom in chain_intermediate.graph:
        atom.draw.colour = 'black'
    if not visualization_mechanism and not draw_structures_per_module:
        Drawer(chain_intermediate)
    # Remove all intermediate structure files
    if path.exists('1.png'):
        os.remove('1.png')
    if path.exists('2.png'):
        os.remove('2.png')
    if path.exists('3.png'):
        os.remove('3.png')
    if path.exists('4.png'):
        os.remove('4.png')
    if path.exists('5.png'):
        os.remove('5.png')
    return chain_intermediate



def display_reactions(structures, tailoring_domains, elongation_unit, module_name):
    """Helper function pks_cluster_to structure. Returns the visualization of
    the quick reaction mechanism per module, based on four types of situation:
    1) module catalyzes only elongation 2) module catalyzes elongation and KR
    reaction 3) module catalyzes elongation, KR and DH reaction 4) module
    catalyzes elongation, KR, DH and ER reaction

    structures: ordered list of names of .png files depicting the structures
    tailoring_domains: list of names of tailoring domains as string
    elongation_unit: 'methylmalonylcoa' or 'malonylcoa'
    module_name: module name as string
    """
    #Situation 1 when the module contains no tailoring domains:
    if len(tailoring_domains) == 0:
        before_elongation, after_elongation = structures
        images = []
        #Read in images
        img1 = mpimg.imread(before_elongation)
        img2 = mpimg.imread(after_elongation)
        imgpijl_malonylcoa = mpimg.imread('pijl_malonylcoa.png')
        imgpijl_methylmalonylcoa = mpimg.imread('pijl_methylmalonylcoa.png')
        imgpijl_methoxymalonylacp = mpimg.imread('pijl_methoxymalonylacp.png')

        #Add images to list to prevent collection by the garbage collector
        images.append(img1)
        images.append(img2)
        if elongation_unit == 'malonylcoa':
            images.append(imgpijl_malonylcoa)
        elif elongation_unit == 'methylmalonylcoa':
            images.append(imgpijl_methylmalonylcoa)
        elif elongation_unit == 'methoxymalonylacp':
            images.append(imgpijl_methoxymalonylacp)

        #Create figure
        fig = plt.figure(constrained_layout=True, figsize=[15, 8],frameon=False, num=f"Quick reaction mechanism {module_name}")
        ax = fig.subplots(1, 5)
        ax = fig.add_gridspec(1, 5)

        #Add first structure before elongation reaction to plot
        ax1 = fig.add_subplot(ax[0, 0:2])
        ax1.imshow(img1)

        #Add correct arrow to plot
        ax1 = fig.add_subplot(ax[0, 2])
        if elongation_unit == 'malonylcoa':
            ax1.imshow(imgpijl_malonylcoa)
        elif elongation_unit == 'methylmalonylcoa':
            ax1.imshow(imgpijl_methylmalonylcoa)
        elif elongation_unit == 'methoxymalonylacp':
            ax1.imshow(imgpijl_methoxymalonylacp)
        ax1.set_title('KS, AT, ACP', pad=20, fontdict={'size': 17})
        ax1.text(140,-550,f"{module_name}", fontdict={'size' : 20, 'weight': 'bold', 'color':'darkred'})

        #Add elongation product to plot
        ax1 = fig.add_subplot(ax[0, 3:5])
        ax1.imshow(img2)

        #Remove all axes from the figure
        for ax in fig.get_axes():
            ax.axis('off')
        plt.show()


    #Situation 2 when the module contains only the KR tailoring domain
    elif len(tailoring_domains) == 1:
        assert 'KR' in tailoring_domains
        before_elongation, after_elongation, after_kr = structures
        images = []
        #Read in images
        img1 = mpimg.imread(before_elongation)
        img2 = mpimg.imread(after_elongation)
        img3 = mpimg.imread(after_kr)
        imgpijl = mpimg.imread('pijl.png')
        imgpijl_KR = mpimg.imread('pijl_KR.png')
        imgpijl_malonylcoa = mpimg.imread('pijl_malonylcoa.png')
        imgpijl_methylmalonylcoa = mpimg.imread('pijl_methylmalonylcoa.png')
        imgpijl_methoxymalonylacp = mpimg.imread('pijl_methoxymalonylacp.png')

        # Add images to list to prevent collection by the garbage collector
        images.append(img1)
        images.append(img2)
        images.append(img3)
        images.append(imgpijl_KR)
        if elongation_unit == 'malonylcoa':
            images.append(imgpijl_malonylcoa)
        elif elongation_unit == 'methylmalonylcoa':
            images.append(imgpijl_methylmalonylcoa)
        elif elongation_unit == 'methoxymalonylacp':
            images.append(imgpijl_methoxymalonylacp)

        # Create figure
        fig = plt.figure(constrained_layout=True, figsize=[15, 8], frameon=False, num=f"Quick reaction mechanism {module_name}")
        ax = fig.subplots(1, 8)
        ax = fig.add_gridspec(1, 8)

        # Add first structure before elongation reaction to plot
        ax1 = fig.add_subplot(ax[0, 0:2])
        ax1.imshow(img1)

        # Add correct arrow to plot (depending on elongation unit)
        ax1 = fig.add_subplot(ax[0, 2])
        if elongation_unit == 'malonylcoa':
            ax1.imshow(imgpijl_malonylcoa)
        elif elongation_unit == 'methylmalonylcoa':
            ax1.imshow(imgpijl_methylmalonylcoa)
        elif elongation_unit == 'methoxymalonylacp':
            ax1.imshow(imgpijl_methoxymalonylacp)
        ax1.set_title('KS, AT, ACP', pad=20, fontdict={'size': 17})

        # Add elongation product to plot
        ax1 = fig.add_subplot(ax[0, 3:5])
        ax1.set_title(f"{module_name}", pad=120, fontdict={'size' : 20, 'weight': 'bold', 'color':'darkred'})
        ax1.imshow(img2)

        # Add KR arrow to plot
        ax1 = fig.add_subplot(ax[0, 5])
        ax1.imshow(imgpijl)
        ax1.set_title('KR', pad=20, fontdict={'size': 17, 'color':'red'})

        # Add ketoreductase product to plot
        ax1 = fig.add_subplot(ax[0, 6:8])
        ax1.imshow(img3)

        # Remove all axes from the figure
        for ax in fig.get_axes():
            ax.axis('off')

        plt.show()


    #Situation 3 when the module has both the KR and DH tailoring domains
    elif len(tailoring_domains) == 2:
        assert 'KR' in tailoring_domains
        assert 'DH' in tailoring_domains
        before_elongation, after_elongation, after_kr, after_dh = structures
        images = []
        #Read in images
        img1 = mpimg.imread(before_elongation)
        img2 = mpimg.imread(after_elongation)
        img3 = mpimg.imread(after_kr)
        img4 = mpimg.imread(after_dh)
        imgpijl = mpimg.imread('pijl.png')
        imgpijl_KR = mpimg.imread('pijl_KR.png')
        imgpijl_malonylcoa = mpimg.imread('pijl_malonylcoa.png')
        imgpijl_methylmalonylcoa = mpimg.imread('pijl_methylmalonylcoa.png')
        imgpijl_methoxymalonylacp = mpimg.imread('pijl_methoxymalonylacp.png')

        # Add images to list to prevent collection by the garbage collector
        images.append(img1)
        images.append(img2)
        images.append(img3)
        images.append(img4)
        images.append(imgpijl)
        images.append(imgpijl_KR)
        if elongation_unit == 'malonylcoa':
            images.append(imgpijl_malonylcoa)
        elif elongation_unit == 'methylmalonylcoa':
            images.append(imgpijl_methylmalonylcoa)
        elif elongation_unit == 'methoxymalonylacp':
            images.append(imgpijl_methoxymalonylacp)

        # Create figure
        fig = plt.figure(constrained_layout=True, figsize=[20, 8], frameon=False, num=f"Quick reaction mechanism {module_name}")
        ax = fig.subplots(1, 11)
        ax = fig.add_gridspec(1, 11)

        # Add first structure before elongation reaction to plot
        ax1 = fig.add_subplot(ax[0, 0:2])
        ax1.imshow(img1)

        # Add correct arrow to plot (depending on elongation unit)
        ax1 = fig.add_subplot(ax[0, 2])
        if elongation_unit == 'malonylcoa':
            ax1.imshow(imgpijl_malonylcoa)
        elif elongation_unit == 'methylmalonylcoa':
            ax1.imshow(imgpijl_methylmalonylcoa)
        elif elongation_unit == 'methoxymalonylacp':
            ax1.imshow(imgpijl_methoxymalonylacp)
        ax1.set_title('KS, AT, ACP', pad=20, fontdict={'size': 17})

        # Add elongation product to plot
        ax1 = fig.add_subplot(ax[0, 3:5])
        ax1.imshow(img2)

        # Add KR arrow to plot
        ax1 = fig.add_subplot(ax[0, 5])
        ax1.imshow(imgpijl)
        ax1.set_title('KR', pad=20, fontdict={'size': 17, 'color':'red'})
        ax1.text(0,-1000,f"{module_name}", fontdict={'size' : 20, 'weight': 'bold', 'color':'darkred'})

        # Add ketoreductase product to plot
        ax1 = fig.add_subplot(ax[0, 6:8])
        ax1.imshow(img3)

        # Add DH arrow to plot
        ax1 = fig.add_subplot(ax[0, 8])
        ax1.imshow(imgpijl)
        ax1.set_title('DH', pad=20, fontdict={'size': 17, 'color':'blue'})

        # Add dehydratase product to plot
        ax1 = fig.add_subplot(ax[0, 9:11])
        ax1.imshow(img4)


        # Remove all axes from the figure
        for ax in fig.get_axes():
            ax.axis('off')

        plt.show()


    #Situation 4 when the module has the KR, DH and ER tailoring domains
    elif len(tailoring_domains) == 3:
        assert 'KR' in tailoring_domains
        assert 'DH' in tailoring_domains
        assert 'ER' in tailoring_domains
        before_elongation, after_elongation, after_kr, after_dh, after_er = structures
        images = []
        #Read in images
        img1 = mpimg.imread(before_elongation)
        img2 = mpimg.imread(after_elongation)
        img3 = mpimg.imread(after_kr)
        img4 = mpimg.imread(after_dh)
        img5 = mpimg.imread(after_er)
        imgpijl = mpimg.imread('pijl.png')
        imgpijl_KR = mpimg.imread('pijl_KR.png')
        imgpijl_malonylcoa = mpimg.imread('pijl_malonylcoa.png')
        imgpijl_methylmalonylcoa = mpimg.imread('pijl_methylmalonylcoa.png')
        imgpijl_methoxymalonylacp = mpimg.imread('pijl_methoxymalonylacp.png')

        # Add images to list to prevent collection by the garbage collector
        images.append(img1)
        images.append(img2)
        images.append(img3)
        images.append(img4)
        images.append(img5)
        images.append(imgpijl)
        images.append(imgpijl_KR)
        if elongation_unit == 'malonylcoa':
            images.append(imgpijl_malonylcoa)
        elif elongation_unit == 'methylmalonylcoa':
            images.append(imgpijl_methylmalonylcoa)
        elif elongation_unit == 'methoxymalonylacp':
            images.append(imgpijl_methoxymalonylacp)

        # Create figure
        fig = plt.figure(constrained_layout=True, figsize=[20, 8], frameon=False, num=f"Quick reaction mechanism {module_name}")
        ax = fig.subplots(1, 14)
        ax = fig.add_gridspec(1, 14)

        # Add first structure before elongation reaction to plot
        ax1 = fig.add_subplot(ax[0, 0:2])
        ax1.imshow(img1)

        # Add correct arrow to plot (depending on elongation unit)
        ax1 = fig.add_subplot(ax[0, 2])
        if elongation_unit == 'malonylcoa':
            ax1.imshow(imgpijl_malonylcoa)
        elif elongation_unit == 'methylmalonylcoa':
            ax1.imshow(imgpijl_methylmalonylcoa)
        elif elongation_unit == 'methoxymalonylacp':
            ax1.imshow(imgpijl_methoxymalonylacp)
        ax1.set_title('KS, AT, ACP', pad=20, fontdict={'size': 17})

        # Add elongation product to plot
        ax1 = fig.add_subplot(ax[0, 3:5])
        ax1.imshow(img2)

        # Add KR arrow to plot
        ax1 = fig.add_subplot(ax[0, 5])
        ax1.imshow(imgpijl)
        ax1.set_title('KR', pad=20, fontdict={'size': 17, 'color':'red'})

        # Add ketoreductase product to plot
        ax1 = fig.add_subplot(ax[0, 6:8])
        ax1.imshow(img3)
        ax1.text(100, -1000, f"{module_name}",fontdict={'size': 20, 'weight': 'bold', 'color': 'darkred'})

        # Add DH arrow to plot
        ax1 = fig.add_subplot(ax[0, 8])
        ax1.imshow(imgpijl)
        ax1.set_title('DH', pad=20, fontdict={'size': 17, 'color':'blue'})

        # Add dehydratase product to plot
        ax1 = fig.add_subplot(ax[0, 9:11])
        ax1.imshow(img4)

        # Add ER arrow to plot
        ax1 = fig.add_subplot(ax[0, 11])
        ax1.imshow(imgpijl)
        ax1.set_title('ER', pad=20, fontdict={'size': 17, 'color':'LIME'})

        # Add enoylreductase product to plot
        ax1 = fig.add_subplot(ax[0, 12:14])
        ax1.imshow(img5)

        # Remove all axes from the figure
        for ax in fig.get_axes():
            ax.axis('off')

        plt.show()


if __name__ == "__main__":
    example_cluster = [['module_1', 'starter_module', 'SC(CC(=O)O)=O'],
                    ['module_2', 'elongation_module', 'methylmalonylcoa', ['KR', 'DH']],
                    ['module_3', 'elongation_module', 'malonylcoa', ['KR']],
                    ['module_4', 'elongation_module', 'methylmalonylcoa', ['KR', 'DH', 'ER']],
                    ['module_5', 'elongation_module', 'methoxymalonylacp', []],
                    ['module_6', 'elongation_module', 'methylmalonylcoa', ['KR', 'DH']]]
    pks_cluster_to_structure(example_cluster)
    pks_cluster_to_structure(example_cluster, visualization_mechanism=True)




