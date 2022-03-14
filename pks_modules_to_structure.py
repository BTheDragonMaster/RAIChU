from pks_tailoring_reactions import *
import os
from os import path
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.patches import FancyArrow
from copy import copy
from NRPS_condensation import condensation_nrps
from copy import deepcopy
from central_atoms_pk_starter import find_central_atoms_pk_starter

KR_NO_KETOREDUCTASE_ACTIVITY = ['KR_C1', 'KR_C2', 'KR_inactive']
N_AMINO_ACID = GroupDefiner('Nitrogen atom amino acid', 'NCC(=O)O', 0)
C1_AMINO_ACID = GroupDefiner('C1 atom amino acid', 'NCC(=O)O', 1)
C2_AMINO_ACID = GroupDefiner('C2 atom amino acid', 'NCC(=O)O', 2)

def pks_cluster_to_structure(modules, visualization_mechanism = False, \
    draw_structures_per_module = False, attach_to_acp = False, \
    draw_mechanism_per_module = False):
    """If visualization_mechanism == False: returns the final polyketide
    structure resulting from the input PKS cluster
    If visualization_mechanism == True: returns a visualization of the elong-
    ation and/or tailoring reactions performed by each PKS module
    If attach_to_acp == True: the polyketide will be attached to an ACP domain
    object of type Domain. Can be combined with previous optional flags.
    If draw_structures_per_module == True: saves the polyketide intermediates
    as 'module_name'.png, doesn't display any visualisations

    modules: [Starter module, elongation module, ..., terminator_module] with
    Starter module: ['module name','starter_module','SMILES starter unit'] or
    ['module name','starter_module_nrps','SMILES starter amino acid'],
    elongation modules: ['module name', 'elongation_module',
    'elongation_monomer', ['KR', 'DH', 'ER']] or ['module name',
    'elongation_module_nrps', 'amino acid'],
    terminator module: ['module name', 'terminator_module',
    'elongation_monomer', ['KR', 'DH', 'ER']] or ['module name',
    'terminator_module_nrps', 'amino acid'],
    """
    last_module_nrps = False
    # Check if last module is an NRPS module:
    if modules[-1][1] == 'elongation_module_nrps' or \
            modules[-1][1] == 'starter_module_nrps' or \
            modules[-1][1] == 'terminator_module_nrps':
        last_module_nrps = True

    # Construct dict to find the SMILES of all amino acids recognized by PARAS
    dict_aa_smiles = make_dict_aa_smiles()

    modules = copy(modules)
    list_drawings_per_module = []

    for module in modules:

        # Find starting unit SMILES, convert to structure
        if module == modules[0]:
            assert len(module) == 3
            assert module[1] == 'starter_module' or\
                   module[1] == 'starter_module_nrps'

            # If starter module = PKS module: construct starter unit from
            # supplied SMILES string
            if module[1] == 'starter_module':
                starter_unit = Smiles(module[2]).smiles_to_structure()
                starter_unit = find_central_atoms_pk_starter(starter_unit)
                if attach_to_acp:
                    chain_intermediate = attach_to_domain_pk(starter_unit, 'ACP')
                else:
                    chain_intermediate = starter_unit
                if draw_structures_per_module:
                    drawing = Drawer(chain_intermediate, dont_show=True)
                    list_drawings_per_module.append([drawing])

            # If starter module = NRPS module: find SMILES in PARAS.txt and
            # build starter unit
            elif module[1] == 'starter_module_nrps':
                starter_unit = Smiles(dict_aa_smiles[module[2].upper()])\
                    .smiles_to_structure()
                n_atoms_aa = find_atoms(N_AMINO_ACID, starter_unit)
                c1_atoms_aa = find_atoms(C1_AMINO_ACID, starter_unit)
                c2_atoms_aa = find_atoms(C2_AMINO_ACID, starter_unit)
                assert len(n_atoms_aa) == 1
                assert len(c1_atoms_aa) == 1
                assert len(c2_atoms_aa) == 1
                for atom in starter_unit.graph:
                    if atom == n_atoms_aa[0]:
                        atom.in_central_chain = True
                    elif atom == c1_atoms_aa[0]:
                        atom.in_central_chain = True
                    elif atom == c2_atoms_aa[0]:
                        atom.in_central_chain = True
                    else:
                        atom.in_central_chain = False

                # If attached, attach
                chain_intermediate = starter_unit
                if draw_structures_per_module and attach_to_acp:
                    copy_chain_intermediate = deepcopy(chain_intermediate)
                    copy_chain_intermediate.find_cycles()
                    copy_attached = attach_to_domain_nrp(copy_chain_intermediate, 'PCP')
                    copy_attached.find_cycles()
                    drawing = Drawer(copy_attached, dont_show=True)
                    list_drawings_per_module.append([drawing])

            # Delete starter module and iterate over remaining modules
            del modules[0]
    for module in modules:

        #If module is a PKS module:
        if module[1] == 'elongation_module' or module[1] == 'terminator_module':
            # If last reaction was an NRPS reaction
            locations = chain_intermediate.find_substructures(Smiles('C(=O)(O)CN').smiles_to_structure())
            if len(locations) > 0:
                chain_intermediate = attach_to_domain_nrp(chain_intermediate, 'ACP')

            module_name, module_type, elongation_unit, list_domains = module
            assert elongation_unit == 'malonylcoa' or \
                   elongation_unit == 'methylmalonylcoa' or \
                   elongation_unit == 'methoxymalonylacp' or \
                   elongation_unit == 'ethylmalonylcoa' or \
                   elongation_unit == 'pk'

            # Reset atom colours in the first structure depicted to black
            for atom in chain_intermediate.graph:
                atom.draw.colour = 'black'
            if elongation_unit == 'malonylcoa':
                if visualization_mechanism == True:
                    # If a previous process was stopped prematurely, make sure the
                    # image files are removed before continuing
                    if path.exists('1.png'):
                        os.remove('1.png')
                    Drawer(chain_intermediate, save_png='1.png', dpi_drawer=500)
                new_chain_intermediate = add_malonylunit(chain_intermediate)
                chain_intermediate = new_chain_intermediate
                if visualization_mechanism == True:
                    if path.exists('2.png'):
                        os.remove('2.png')
                    Drawer(chain_intermediate, save_png='2.png', dpi_drawer=500)
            elif elongation_unit == 'methylmalonylcoa':
                if visualization_mechanism == True:
                    if path.exists('1.png'):
                        os.remove('1.png')
                    Drawer(chain_intermediate, save_png='1.png', dpi_drawer=500)
                new_chain_intermediate = add_methylmalonylunit(chain_intermediate)
                chain_intermediate = new_chain_intermediate
                if visualization_mechanism == True:
                    if path.exists('2.png'):
                        os.remove('2.png')
                    Drawer(chain_intermediate, save_png='2.png', dpi_drawer=500)
            elif elongation_unit == 'methoxymalonylacp':
                if visualization_mechanism == True:
                    if path.exists('1.png'):
                        os.remove('1.png')
                    Drawer(chain_intermediate, save_png='1.png', dpi_drawer=500)
                new_chain_intermediate = add_methoxymalonylunit(chain_intermediate)
                chain_intermediate = new_chain_intermediate
                if visualization_mechanism == True:
                    if path.exists('2.png'):
                        os.remove('2.png')
                    Drawer(chain_intermediate, save_png='2.png', dpi_drawer=500)
            elif elongation_unit == 'ethylmalonylcoa':
                if visualization_mechanism == True:
                    if path.exists('1.png'):
                        os.remove('1.png')
                    Drawer(chain_intermediate, save_png='1.png', dpi_drawer=500)
                new_chain_intermediate = add_ethylmalonylunit(chain_intermediate)
                chain_intermediate = new_chain_intermediate
                if visualization_mechanism == True:
                    if path.exists('2.png'):
                        os.remove('2.png')
                    Drawer(chain_intermediate, save_png='2.png', dpi_drawer=500)
            elif elongation_unit == 'pk':
                if visualization_mechanism == True:
                    if path.exists('1.png'):
                        os.remove('1.png')
                    Drawer(chain_intermediate, save_png='1.png', dpi_drawer=500)
                new_chain_intermediate = add_pkunit(chain_intermediate)
                chain_intermediate = new_chain_intermediate
                if visualization_mechanism == True:
                    if path.exists('2.png'):
                        os.remove('2.png')
                    Drawer(chain_intermediate, save_png='2.png', dpi_drawer=500)

            # If the module does not contain any tailoring domains:
            if len(list_domains) == 0:
                if visualization_mechanism == True:
                    if elongation_unit == 'methylmalonylcoa':
                        display_reactions(['1.png', '2.png'],
                                          list_domains,'methylmalonylcoa',
                                          module_name,
                                          draw_mechanism_per_module)
                    elif elongation_unit == 'malonylcoa':
                        display_reactions(['1.png', '2.png'],list_domains,
                                          'malonylcoa', module_name,
                                          draw_mechanism_per_module)
                    elif elongation_unit == 'methoxymalonylacp':
                        display_reactions(['1.png', '2.png'], list_domains,
                                          'methoxymalonylacp', module_name,
                                          draw_mechanism_per_module)
                    elif elongation_unit == 'ethylmalonylcoa':
                        display_reactions(['1.png', '2.png'], list_domains,
                                          'ethylmalonylcoa', module_name,
                                          draw_mechanism_per_module)
                    elif elongation_unit == 'pk':
                        display_reactions(['1.png', '2.png'], list_domains,
                                          'pk', module_name,
                                          draw_mechanism_per_module)

            # If module contains tailoring domains, perform tailoring reactions
            for domain in list_domains:
                inactive_kr_domain = False

                if domain.startswith('KR'):
                    if domain == 'KR':
                        new_chain_intermediate = ketoreductase(chain_intermediate)
                        chain_intermediate = new_chain_intermediate
                    elif domain == 'KR_inactive':
                        inactive_kr_domain = True
                    else:
                        kr_type = domain[-2:]
                        new_chain_intermediate = \
                            ketoreductase(chain_intermediate, kr_type = kr_type)
                        chain_intermediate = new_chain_intermediate
                    if visualization_mechanism == True and inactive_kr_domain == False:
                        if path.exists('3.png'):
                            os.remove('3.png')
                        Drawer(chain_intermediate, save_png='3.png', dpi_drawer=500)
                        # Check if the module also contains a DH domain
                        if 'DH' not in list_domains:
                            display_reactions(['1.png', '2.png', '3.png'],
                                              list_domains, elongation_unit,
                                              module_name,
                                              draw_mechanism_per_module)
                    elif inactive_kr_domain == True and visualization_mechanism == True:
                        list_domains = []
                        if elongation_unit == 'methylmalonylcoa':
                            display_reactions(['1.png', '2.png'], list_domains,
                                              'methylmalonylcoa', module_name,
                                              draw_mechanism_per_module)
                        elif elongation_unit == 'malonylcoa':
                            display_reactions(['1.png', '2.png'], list_domains,
                                              'malonylcoa', module_name,
                                              draw_mechanism_per_module)
                        elif elongation_unit == 'methoxymalonylacp':
                            display_reactions(['1.png', '2.png'], list_domains,
                                              'methoxymalonylacp', module_name,
                                              draw_mechanism_per_module)
                        elif elongation_unit == 'ethylmalonylcoa':
                            display_reactions(['1.png', '2.png'], list_domains,
                                              'ethylmalonylcoa', module_name,
                                              draw_mechanism_per_module)

                if domain == 'DH':
                    assert any(domain.startswith('KR') for domain in list_domains)
                    executable = True
                    for domain_type in list_domains:
                        if domain_type in KR_NO_KETOREDUCTASE_ACTIVITY:
                            executable = False
                    if executable:
                        new_chain_intermediate = dehydratase(chain_intermediate)
                        chain_intermediate = new_chain_intermediate
                        if visualization_mechanism == True:
                            if path.exists('4.png'):
                                os.remove('4.png')
                            Drawer(chain_intermediate, save_png='4.png', dpi_drawer=500)
                            if 'ER' not in list_domains:
                                display_reactions(['1.png', '2.png', '3.png','4.png'],
                                                  list_domains, elongation_unit,
                                                  module_name, draw_mechanism_per_module)

                if domain == 'ER':
                    assert any(domain.startswith('KR') for domain in list_domains)
                    assert 'DH' in list_domains
                    executable = True
                    for domain_type in list_domains:
                        if domain_type in KR_NO_KETOREDUCTASE_ACTIVITY:
                            executable = False
                    if executable:
                        new_chain_intermediate = enoylreductase(chain_intermediate)
                        chain_intermediate = new_chain_intermediate
                        if visualization_mechanism == True:
                            if path.exists('5.png'):
                                os.remove('5.png')
                            Drawer(chain_intermediate, save_png='5.png', dpi_drawer=500)
                            display_reactions(['1.png', '2.png', '3.png','4.png','5.png'],
                                              list_domains, elongation_unit,
                                              module_name,
                                              draw_mechanism_per_module)

            if draw_structures_per_module and attach_to_acp:
                drawing = Drawer(chain_intermediate, dont_show=True, dpi_drawer=500)
                list_drawings_per_module.append([drawing])

        # If the module is an NRPS module:
        elif module[1] == 'elongation_module_nrps' or \
                module[1] == 'terminator_module_nrps':
            module_name, module_type, aa_specifity = module
            if visualization_mechanism == True and attach_to_acp:
                for atom in chain_intermediate.graph:
                    atom.draw.colour = 'black'
                # If a previous process was stopped prematurely, make sure the
                # image files are removed before continuing
                if path.exists('1.png'):
                    os.remove('1.png')
                copy_chain_intermediate = deepcopy(chain_intermediate)
                copy_chain_intermediate.find_cycles()
                if not any(hasattr(atom, 'domain_type') for
                           atom in copy_chain_intermediate.graph):
                    copy_attached = attach_to_domain_nrp(copy_chain_intermediate, 'PCP')
                else:
                    copy_attached = copy_chain_intermediate
                Drawer(copy_attached, save_png='1.png', dpi_drawer=500)

            # Reset atom colours in the first structure depicted to black
            for atom in chain_intermediate.graph:
                atom.draw.colour = 'black'

            # Build amino acid structure from dict containig PARAS SMILES
            aa_specifity = aa_specifity.upper()
            assert aa_specifity in dict_aa_smiles
            aa_structure = Smiles(dict_aa_smiles[aa_specifity]).smiles_to_structure()

            # Check that the structure is an amino acid
            if not (len(aa_structure.find_substructures(
                    Smiles('CN').smiles_to_structure())) > 0 \
                    and len(aa_structure.find_substructures(
                        Smiles('C(O)=O').smiles_to_structure())) > 0):
                raise ValueError(
                    f'The structure: {aa_specifity}, is not an amino acid')

            # Perform condensation reaction
            chain_intermediate = condensation_nrps(aa_structure, chain_intermediate)

            # Save drawings if necessary
            if draw_structures_per_module and attach_to_acp:
                copy_chain_intermediate = deepcopy(chain_intermediate)
                copy_chain_intermediate.find_cycles()
                copy_attached = attach_to_domain_nrp(copy_chain_intermediate, 'PCP')
                copy_attached.refresh_structure()
                copy_attached.set_connectivities()
                copy_attached.find_cycles()
                drawing = Drawer(copy_attached, dont_show=True, dpi_drawer=500)
                list_drawings_per_module.append([drawing])
            elif visualization_mechanism == True and attach_to_acp:
                copy_chain_intermediate = deepcopy(chain_intermediate)
                copy_chain_intermediate.find_cycles()
                copy_attached = attach_to_domain_nrp(copy_chain_intermediate, 'PCP')
                if path.exists('2.png'):
                    os.remove('2.png')
                Drawer(copy_attached, save_png='2.png', dpi_drawer=500)
                list_domains = None
                elongation_unit = aa_specifity.lower()
                display_reactions(['1.png', '2.png'], list_domains, elongation_unit, module_name, draw_mechanism_per_module)
                # Necessary for thioesterase reactions when last module is an
                # NRPS module
                if module == modules[-1]:
                    chain_intermediate = copy_attached

    # Reset the atom color in the final structure to black
    for atom in chain_intermediate.graph:
        atom.draw.colour = 'black'

    if not visualization_mechanism and not draw_structures_per_module:
        if last_module_nrps and attach_to_acp:
            chain_intermediate.find_cycles()
            chain_intermediate = attach_to_domain_nrp(chain_intermediate, 'PCP')
            chain_intermediate.find_cycles()
            Drawer(chain_intermediate)
        else:
            chain_intermediate.find_cycles()
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

    if draw_structures_per_module:
        return(list_drawings_per_module)

    return chain_intermediate



def display_reactions(structures, tailoring_domains, elongation_unit,
                      module_name, draw_mechanism_per_module):
    """Helper function pks_cluster_to structure. Returns the visualization of
    the quick reaction mechanism per module, based on four types of situation:
    1) NRPS elongation module, 2) PKS module that catalyzes only elongation
    3) PKS module that catalyzes elongation and KR reaction
    4) PKS module that catalyzes elongation, KR and DH reaction
    5) PKS module that catalyzes elongation, KR, DH and ER reaction

    structures: ordered list of names of .png files depicting the structures
    tailoring_domains: list of names of tailoring domains as string
    elongation_unit: 'amino_acid_name', 'methylmalonylcoa', 'malonylcoa',
    'methoxymalonylacp', 'ethylmalonylcoa' or 'pk'
    module_name: module name as string
    """
    # Situation 1:  NRPS elongation module
    if tailoring_domains == None:
        before_elongation, after_elongation = structures
        images = []
        # Read in images
        img1 = mpimg.imread(before_elongation)
        img2 = mpimg.imread(after_elongation)

        # Add images to list to prevent collection by the garbage collector
        images.append(img1)
        images.append(img2)

        # Create figure
        fig = plt.figure(constrained_layout=True, figsize=[15, 8],
                         frameon=False,
                         num=f"Quick reaction mechanism {module_name}")
        ax = fig.subplots(1, 5)
        ax = fig.add_gridspec(1, 5)

        # Add first structure before elongation reaction to plot
        ax1 = fig.add_subplot(ax[0, 0:2])
        ax1.imshow(img1)

        # Add correct arrow to plot
        ax1 = fig.add_subplot(ax[0, 2])
        arrow = FancyArrow(0, 0.5, 0.89, 0, overhang=0.3, head_width=0.025,
                           head_length=0.1, color='black')
        images.append(arrow)
        ax1.add_patch(arrow)
        text = elongation_unit
        ax1.text(0.5, 0.55, text, ha='center',
                 fontdict={'size': 12, 'color': 'black'})
        ax1.set_title(f"{module_name}", pad=70,
                      fontdict={'size': 20, 'weight': 'bold',
                                'color': 'darkred'})

        # Add elongation product to plot
        ax1 = fig.add_subplot(ax[0, 3:5])
        ax1.imshow(img2)

        # Remove all axes from the figure
        for ax in fig.get_axes():
            ax.axis('off')

        if draw_mechanism_per_module:
            plt.savefig(f'{module_name}_quick_mechanism.png')
            plt.clf()
            plt.close()
        else:
            plt.show()

    else:
        # Situation 2: PKS module contains no tailoring domains:
        if len(tailoring_domains) == 0:
            before_elongation, after_elongation = structures
            images = []
            # Read in images
            img1 = mpimg.imread(before_elongation)
            img2 = mpimg.imread(after_elongation)

            # Add images to list to prevent collection by the garbage collector
            images.append(img1)
            images.append(img2)

            # Create figure
            fig = plt.figure(constrained_layout=True, figsize=[15, 8],
                             frameon=False,
                             num=f"Quick reaction mechanism {module_name}")
            ax = fig.subplots(1, 5)
            ax = fig.add_gridspec(1, 5)

            # Add first structure before elongation reaction to plot
            ax1 = fig.add_subplot(ax[0, 0:2])
            ax1.imshow(img1)

            # Add correct arrow to plot
            ax1 = fig.add_subplot(ax[0, 2])
            arrow = FancyArrow(0, 0.5,0.89, 0, overhang = 0.3,
                               head_width=0.025, head_length=0.1,
                               color= 'black')
            images.append(arrow)
            ax1.add_patch(arrow)
            if elongation_unit == 'malonylcoa':
                text = 'Malonyl-CoA'
            elif elongation_unit == 'methylmalonylcoa':
                text = 'Methylmalonyl-CoA'
            elif elongation_unit == 'methoxymalonylacp':
                text = 'Methoxymalonyl-ACP'
            elif elongation_unit == 'ethylmalonylcoa':
                text = 'Ethylmalonyl-CoA'
            elif elongation_unit == 'pk':
                text = 'Unknown elongation\nunit'
            ax1.text(0.5, 0.55, text, ha='center',
                     fontdict={'size': 12, 'color': 'black'})
            ax1.set_title(f"{module_name}", pad = 70,
                          fontdict={'size' : 20, 'weight': 'bold', 'color':'darkred'})

            # Add elongation product to plot
            ax1 = fig.add_subplot(ax[0, 3:5])
            ax1.imshow(img2)

            # Remove all axes from the figure
            for ax in fig.get_axes():
                ax.axis('off')

            if draw_mechanism_per_module:
                plt.savefig(f'{module_name}_quick_mechanism.png')
                plt.clf()
                plt.close()
            else:
                plt.show()

        # Situation 3: PKS module contains only the KR tailoring domain
        elif len(tailoring_domains) == 1:
            assert any(domain.startswith('KR') for domain in tailoring_domains)
            kr_domain = tailoring_domains[0]
            before_elongation, after_elongation, after_kr = structures
            images = []
            #Read in images
            img1 = mpimg.imread(before_elongation)
            img2 = mpimg.imread(after_elongation)
            img3 = mpimg.imread(after_kr)

            # Add images to list to prevent collection by the garbage collector
            images.append(img1)
            images.append(img2)
            images.append(img3)

            # Create figure
            fig = plt.figure(constrained_layout=True, figsize=[15, 8],
                             frameon=False, num=f"Quick reaction mechanism {module_name}")
            ax = fig.subplots(1, 8)
            ax = fig.add_gridspec(1, 8)

            # Add first structure before elongation reaction to plot
            ax1 = fig.add_subplot(ax[0, 0:2])
            ax1.imshow(img1)

            # Add correct arrow to plot (depending on elongation unit)
            ax1 = fig.add_subplot(ax[0, 2])
            arrow = FancyArrow(0, 0.5,0.89, 0, overhang = 0.3,
                               head_width=0.025, head_length=0.1, color= 'black')
            images.append(arrow)
            ax1.add_patch(arrow)
            if elongation_unit == 'malonylcoa':
                text = 'Malonyl-CoA'
            elif elongation_unit == 'methylmalonylcoa':
                text = 'Methylmalonyl-CoA'
            elif elongation_unit == 'methoxymalonylacp':
                text = 'Methoxymalonyl-ACP'
            elif elongation_unit == 'ethylmalonylcoa':
                text = 'Ethylmalonyl-CoA'
            elif elongation_unit == 'pk':
                text = 'Unknown elongation unit'
            ax1.text(0.5, 0.55, text, ha='center',
                     fontdict={'size': 12, 'color': 'black'})

            # Add elongation product to plot
            ax1 = fig.add_subplot(ax[0, 3:5])
            ax1.set_title(f"{module_name}", pad=50, fontdict={'size' : 20, 'weight': 'bold', 'color':'darkred'})
            ax1.imshow(img2)

            # Add KR arrow to plot
            ax1 = fig.add_subplot(ax[0, 5])
            arrow = FancyArrow(0, 0.5, 0.89, 0, overhang = 0.3, head_width=0.025, head_length=0.1,
                               color='black')
            images.append(arrow)
            ax1.add_patch(arrow)
            ax1.text(0.5, 0.55, kr_domain, ha='center',
                     fontdict={'size': 17, 'color': 'red'})

            # Add ketoreductase product to plot
            ax1 = fig.add_subplot(ax[0, 6:8])
            ax1.imshow(img3)

            # Remove all axes from the figure
            for ax in fig.get_axes():
                ax.axis('off')

            if draw_mechanism_per_module:
                plt.savefig(f'{module_name}_quick_mechanism.png')
                plt.clf()
                plt.close()
            else:
                plt.show()

        # Situation 4: PKS module has both the KR and DH tailoring domains
        elif len(tailoring_domains) == 2:
            assert any(domain.startswith('KR') for domain in tailoring_domains)
            for domain in tailoring_domains:
                if domain.startswith('KR'):
                    kr_domain = domain
            assert 'DH' in tailoring_domains
            before_elongation, after_elongation, after_kr, after_dh = structures
            images = []
            # Read in images
            img1 = mpimg.imread(before_elongation)
            img2 = mpimg.imread(after_elongation)
            img3 = mpimg.imread(after_kr)
            img4 = mpimg.imread(after_dh)

            # Add images to list to prevent collection by the garbage collector
            images.append(img1)
            images.append(img2)
            images.append(img3)
            images.append(img4)

            # Create figure
            fig = plt.figure(constrained_layout=True, figsize=[20, 8], frameon=False, num=f"Quick reaction mechanism {module_name}")
            ax = fig.subplots(1, 11)
            ax = fig.add_gridspec(1, 11)

            # Add first structure before elongation reaction to plot
            ax1 = fig.add_subplot(ax[0, 0:2])
            ax1.imshow(img1)

            # Add correct arrow to plot (depending on elongation unit)
            ax1 = fig.add_subplot(ax[0, 2])
            arrow = FancyArrow(0, 0.5,0.89, 0, overhang = 0.3, head_width=0.025, head_length=0.1, color= 'black')
            images.append(arrow)
            ax1.add_patch(arrow)
            if elongation_unit == 'malonylcoa':
                text = 'Malonyl-CoA'
            elif elongation_unit == 'methylmalonylcoa':
                text = 'Methylmalonyl-CoA'
            elif elongation_unit == 'methoxymalonylacp':
                text = 'Methoxymalonyl-ACP'
            elif elongation_unit == 'ethylmalonylcoa':
                text = 'Ethylmalonyl-CoA'
            elif elongation_unit == 'pk':
                text = 'Unknown elongation unit'
            ax1.text(0.5, 0.55, text, ha='center',
                     fontdict={'size': 12, 'color': 'black'})

            # Add elongation product to plot
            ax1 = fig.add_subplot(ax[0, 3:5])
            ax1.imshow(img2)

            # Add KR arrow to plot
            ax1 = fig.add_subplot(ax[0, 5])
            arrow = FancyArrow(0, 0.5, 0.89, 0, overhang = 0.3,
                               head_width=0.025, head_length=0.1,
                               color='black')
            images.append(arrow)
            ax1.add_patch(arrow)
            ax1.text(0.5, 0.55, kr_domain, ha='center',
                     fontdict={'size': 17, 'color': 'red'})

            # Add ketoreductase product to plot
            ax1 = fig.add_subplot(ax[0, 6:8])
            ax1.imshow(img3)

            # Add DH arrow to plot
            ax1 = fig.add_subplot(ax[0, 8])
            arrow = FancyArrow(0, 0.5, 0.89, 0, overhang = 0.3,
                               head_width=0.025, head_length=0.1,
                               color='black')
            images.append(arrow)
            ax1.add_patch(arrow)
            ax1.text(0.5, 0.55, 'DH', ha='center',
                     fontdict={'size': 17, 'color': 'blue'})

            # Add dehydratase product to plot
            ax1 = fig.add_subplot(ax[0, 9:11])
            ax1.imshow(img4)

            # Remove all axes from the figure
            for ax in fig.get_axes():
                ax.axis('off')

            if draw_mechanism_per_module:
                plt.savefig(f'{module_name}_quick_mechanism.png')
                plt.clf()
                plt.close()
            else:
                plt.show()

        # Situation 5: PKS module has the KR, DH and ER tailoring domains
        elif len(tailoring_domains) == 3:
            assert any(domain.startswith('KR') for domain in tailoring_domains)
            for domain in tailoring_domains:
                if domain.startswith('KR'):
                    kr_domain = domain
            assert 'DH' in tailoring_domains
            assert 'ER' in tailoring_domains
            before_elongation, after_elongation, after_kr, after_dh, after_er = structures
            images = []
            # Read in images
            img1 = mpimg.imread(before_elongation)
            img2 = mpimg.imread(after_elongation)
            img3 = mpimg.imread(after_kr)
            img4 = mpimg.imread(after_dh)
            img5 = mpimg.imread(after_er)

            # Add images to list to prevent collection by the garbage collector
            images.append(img1)
            images.append(img2)
            images.append(img3)
            images.append(img4)
            images.append(img5)

            # Create figure
            fig = plt.figure(constrained_layout=True, figsize=[20, 8],
                             frameon=False, num=f"Quick reaction mechanism {module_name}")
            ax = fig.subplots(1, 14)
            ax = fig.add_gridspec(1, 14)

            # Add first structure before elongation reaction to plot
            ax1 = fig.add_subplot(ax[0, 0:2])
            ax1.imshow(img1)

            # Add correct arrow to plot (depending on elongation unit)
            ax1 = fig.add_subplot(ax[0, 2])
            arrow = FancyArrow(0, 0.5,0.89, 0, overhang = 0.3,
                               head_width=0.025, head_length=0.1, color= 'black')
            images.append(arrow)
            ax1.add_patch(arrow)
            if elongation_unit == 'malonylcoa':
                text = 'Malonyl-CoA'
            elif elongation_unit == 'methylmalonylcoa':
                text = 'Methylmalonyl-CoA'
            elif elongation_unit == 'methoxymalonylacp':
                text = 'Methoxymalonyl-ACP'
            elif elongation_unit == 'ethylmalonylcoa':
                text = 'Ethylmalonyl-CoA'
            elif elongation_unit == 'pk':
                text = 'Unknown elongation unit'
            ax1.text(0.5, 0.55, text, ha='center',
                     fontdict={'size': 12, 'color': 'black'})

            # Add elongation product to plot
            ax1 = fig.add_subplot(ax[0, 3:5])
            ax1.imshow(img2)

            # Add KR arrow to plot
            ax1 = fig.add_subplot(ax[0, 5])
            arrow = FancyArrow(0, 0.5, 0.89, 0, overhang = 0.3, head_width=0.025, head_length=0.1,
                               color='black')
            images.append(arrow)
            ax1.add_patch(arrow)
            ax1.text(0.5, 0.55, kr_domain, ha='center',
                     fontdict={'size': 17, 'color': 'red'})

            # Add ketoreductase product to plot
            ax1 = fig.add_subplot(ax[0, 6:8])
            ax1.imshow(img3)
            ax1.text(100, -800, f"{module_name}",fontdict={'size': 20, 'weight': 'bold', 'color': 'darkred'})

            # Add DH arrow to plot
            ax1 = fig.add_subplot(ax[0, 8])
            arrow = FancyArrow(0, 0.5, 0.89, 0, overhang = 0.3, head_width=0.025, head_length=0.1,
                               color='black')
            images.append(arrow)
            ax1.add_patch(arrow)
            ax1.text(0.5, 0.55, 'DH', ha='center',
                     fontdict={'size': 17, 'color': 'blue'})

            # Add dehydratase product to plot
            ax1 = fig.add_subplot(ax[0, 9:11])
            ax1.imshow(img4)

            # Add ER arrow to plot
            ax1 = fig.add_subplot(ax[0, 11])
            arrow = FancyArrow(0, 0.5, 0.89, 0, overhang = 0.3,
                               head_width=0.025, head_length=0.1,
                               color='black')
            images.append(arrow)
            ax1.add_patch(arrow)
            ax1.text(0.5, 0.55, 'ER', ha='center',
                     fontdict={'size': 17, 'color': 'lime'})

            # Add enoylreductase product to plot
            ax1 = fig.add_subplot(ax[0, 12:14])
            ax1.imshow(img5)

            # Remove all axes from the figure
            for ax in fig.get_axes():
                ax.axis('off')

            if draw_mechanism_per_module:
                plt.savefig(f'{module_name}_quick_mechanism.png')
                plt.clf()
                plt.close()
            else:
                plt.show()


def make_dict_aa_smiles():
    """
    Returns a dict as {name_aa : SMILES_aa} from the amino acids in the
    PARAS_smiles.txt file
    """
    # Parse list SMILES amino acids attached to PCP to dict name -> Structure
    lines_aa = open('PARAS_smiles.txt', 'r', encoding='utf8').readlines()
    dict_aa_smiles = {}

    for line in lines_aa:
        line = line.strip()
        name, smiles = line.split()
        name = name.upper()
        dict_aa_smiles[name] = smiles

    return dict_aa_smiles







