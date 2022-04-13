from pks_tailoring_reactions import *
from RaichuDrawer import RaichuDrawer
import os
from os import path
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from pikachu.general import read_smiles
from matplotlib.patches import FancyArrow
from attributes import ATTRIBUTES
from NRPS_condensation import condensation_nrps, set_nrps_central_chain
from nrps_tailoring_reactions import nrps_epimerization, nrps_methylation
from central_atoms_pk_starter import find_central_atoms_pk_starter
from attributes import ALL_PKS_ELONGATION_UNITS


KR_NO_KETOREDUCTASE_ACTIVITY = ['KR_C1', 'KR_C2', 'KR_inactive']

ELONGATION_UNIT_TO_TEXT = {'malonylcoa': 'Malonyl-CoA',
                           'methylmalonylcoa': 'Methylmalonyl-CoA',
                           'methoxymalonylacp': 'Methoxymalonyl-ACP',
                           'ethylmalonylcoa': 'Ethylmalonyl-CoA',
                           'pk': 'Unknown elongation unit'}

TAILOR_DOMAIN_SHORT_TO_LONG = {'E': 'Epimerization', 'nMT' : 'N-methylation'}


def cluster_to_structure(modules, visualization_mechanism=False,
                         draw_structures_per_module=False,
                         attach_to_acp=False,
                         draw_mechanism_per_module=False):
    """If visualization_mechanism == False: returns the final polyketide
    structure resulting from the input PKS cluster
    If visualization_mechanism == True: returns a visualization of the elong-
    ation and/or tailoring reactions performed by each PKS module
    If attach_to_acp == True: the polyketide will be attached to an ACP domain
    object of type Domain. Can be combined with previous optional flags.
    If draw_structures_per_module == True: saves the polyketide intermediates
    as 'module_name'.png, doesn't display any visualisations

    modules: [Starter module, elongation module, ..., terminator_module] with
    Starter module: ['module name','starter_module_pks','SMILES starter unit']
    or ['module name','starter_module_nrps','SMILES starter amino acid'],
    elongation modules: ['module name', 'elongation_module_pks',
    'elongation_monomer', ['KR', 'DH', 'ER']] or ['module name',
    'elongation_module_nrps', 'amino acid'],
    terminator module: ['module name', 'terminator_module_pks',
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

    modules = modules[::]
    list_drawings_per_module = []

    # Retrieve and remove starter module

    starter_module = modules.pop(0)

    assert len(starter_module) == 3

    starter_module_type = starter_module[1]
    starter_module_smiles = starter_module[2]

    assert starter_module_type == 'starter_module_pks' or starter_module_type == 'starter_module_nrps'

    # If starter module = PKS module: construct starter unit from
    # supplied SMILES string
    if starter_module_type == 'starter_module_pks':
        starter_unit = read_smiles(starter_module_smiles)
        starter_unit.add_attributes(ATTRIBUTES, boolean=True)
        starter_unit = find_central_atoms_pk_starter(starter_unit)

        if attach_to_acp:
            chain_intermediate = attach_to_domain_pk(starter_unit, 'ACP')
        else:
            chain_intermediate = starter_unit
        if draw_structures_per_module:
            drawing = RaichuDrawer(chain_intermediate, dont_show=True)
            list_drawings_per_module.append([drawing])

    # If starter module = NRPS module: find SMILES in PARAS.txt and
    # build starter unit
    elif starter_module_type == 'starter_module_nrps':

        starter_unit = read_smiles(dict_aa_smiles[starter_module_smiles.upper()])
        starter_unit.add_attributes(ATTRIBUTES, boolean=True)
        set_nrps_central_chain(starter_unit)

        # If attached, attach
        chain_intermediate = starter_unit
        if draw_structures_per_module and attach_to_acp:
            copy_chain_intermediate = chain_intermediate.deepcopy()
            copy_chain_intermediate.find_cycles()
            copy_attached = attach_to_domain_nrp(copy_chain_intermediate, 'PCP')
            copy_attached.refresh_structure(find_cycles=True)
            drawing = RaichuDrawer(copy_attached, dont_show=True)
            list_drawings_per_module.append([drawing])

    # Iterate over remaining modules

    for module in modules:
        print(module)
        list_filenames = []

        # If module is a PKS module:

        if module[1] == 'elongation_module_pks' or module[1] == 'terminator_module_pks':
            # If last reaction was an NRPS reaction
            locations = chain_intermediate.find_substructures(read_smiles('C(=O)(O)CN'))
            if len(locations) > 0:
                chain_intermediate = attach_to_domain_nrp(chain_intermediate, 'ACP')

            module_name, module_type, elongation_unit, list_domains = module
            assert elongation_unit in ALL_PKS_ELONGATION_UNITS

            # Reset atom colours in the first structure depicted to black
            for atom in chain_intermediate.graph:
                atom.draw.colour = 'black'

            elongation_monomer = PKS_SUBUNIT_TO_MONOMER[elongation_unit]

            if visualization_mechanism:
                # If a previous process was stopped prematurely, make sure the
                # image files are removed before continuing
                if path.exists('1.png'):
                    os.remove('1.png')
                RaichuDrawer(chain_intermediate, save_png='1.png', dpi=500)
                list_filenames.append('1.png')

            new_chain_intermediate = pks_elongation(chain_intermediate, elongation_monomer)
            chain_intermediate = new_chain_intermediate

            if visualization_mechanism:
                if path.exists('2.png'):
                    os.remove('2.png')
                RaichuDrawer(chain_intermediate, save_png='2.png', dpi=500)
                list_filenames.append('2.png')

            # If the module does not contain any tailoring domains:
            if len(list_domains) == 0:
                if visualization_mechanism:
                    display_reactions(['1.png', '2.png'],
                                      list_domains, elongation_unit,
                                      module_name,
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
                    if visualization_mechanism and not inactive_kr_domain:
                        if path.exists('3.png'):
                            os.remove('3.png')
                        RaichuDrawer(chain_intermediate, save_png='3.png', dpi=500)
                        list_filenames.append('3.png')
                        # Check if the module also contains a DH domain
                        if 'DH' not in list_domains:
                            display_reactions(['1.png', '2.png', '3.png'],
                                              list_domains, elongation_unit,
                                              module_name,
                                              draw_mechanism_per_module)
                    elif inactive_kr_domain and visualization_mechanism:
                        list_domains = []

                        display_reactions(['1.png', '2.png'], list_domains,
                                          elongation_unit, module_name,
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
                        if visualization_mechanism:
                            if path.exists('4.png'):
                                os.remove('4.png')
                            RaichuDrawer(chain_intermediate, save_png='4.png', dpi=500)
                            list_filenames.append('4.png')
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
                        if visualization_mechanism:
                            if path.exists('5.png'):
                                os.remove('5.png')
                            RaichuDrawer(chain_intermediate, save_png='5.png', dpi=500)
                            list_filenames.append('5.png')
                            display_reactions(['1.png', '2.png', '3.png','4.png','5.png'],
                                              list_domains, elongation_unit,
                                              module_name,
                                              draw_mechanism_per_module)

            if draw_structures_per_module and attach_to_acp:
                drawing = RaichuDrawer(chain_intermediate, dont_show=True)
                list_drawings_per_module.append([drawing])

        # If the module is an NRPS module:
        elif module[1] == 'elongation_module_nrps' or \
                module[1] == 'terminator_module_nrps':

            module_name, module_type, aa_specifity, list_tailoring_domains = module
            if visualization_mechanism and attach_to_acp:
                for atom in chain_intermediate.graph:
                    atom.draw.colour = 'black'
                # If a previous process was stopped prematurely, make sure the
                # image files are removed before continuing
                if path.exists('1.png'):
                    os.remove('1.png')
                copy_chain_intermediate = chain_intermediate.deepcopy()
                copy_chain_intermediate.find_cycles()
                if not any(atom.annotations.domain_type for
                           atom in copy_chain_intermediate.graph):
                    copy_attached = attach_to_domain_nrp(copy_chain_intermediate, 'PCP')
                else:
                    copy_attached = copy_chain_intermediate
                RaichuDrawer(copy_attached, save_png='1.png', dpi=500)
                list_filenames.append('1.png')

            # Reset atom colours in the first structure depicted to black
            for atom in chain_intermediate.graph:
                atom.draw.colour = 'black'

            # Build amino acid structure from dict containig PARAS SMILES
            aa_specifity = aa_specifity.upper()
            assert aa_specifity in dict_aa_smiles
            aa_structure = read_smiles(dict_aa_smiles[aa_specifity])

            # Check that the structure is an amino acid
            if not (len(aa_structure.find_substructures(
                    read_smiles('CN'))) > 0
                    and len(aa_structure.find_substructures(
                        read_smiles('C(O)=O'))) > 0):
                raise ValueError(
                    f'The structure: {aa_specifity}, is not an amino acid')

            # Perform condensation reaction

            chain_intermediate = condensation_nrps(aa_structure, chain_intermediate)
            # if len(list_domains) == 0:
            #     if visualization_mechanism:
            #         display_reactions(['1.png', '2.png'],
            #                           list_domains, elongation_unit,
            #                           module_name,
            #                           draw_mechanism_per_module)
            if len(list_tailoring_domains) == 0:
            # Save drawings if necessary
                if draw_structures_per_module and attach_to_acp:
                    copy_chain_intermediate = chain_intermediate.deepcopy()
                    # copy_chain_intermediate.find_cycles()
                    copy_attached = attach_to_domain_nrp(copy_chain_intermediate, 'PCP')
                    copy_attached.refresh_structure()
                    copy_attached.set_connectivities()
                    copy_attached.find_cycles()
                    drawing = RaichuDrawer(copy_attached, dont_show=True)
                    list_drawings_per_module.append([drawing])

                elif visualization_mechanism and attach_to_acp:
                        copy_chain_intermediate = chain_intermediate.deepcopy()
                        copy_chain_intermediate.find_cycles()
                        copy_attached = attach_to_domain_nrp(copy_chain_intermediate, 'PCP')
                        if path.exists('2.png'):
                            os.remove('2.png')
                        RaichuDrawer(copy_attached, save_png='2.png', dpi=500)
                        list_filenames.append('2.png')
                        elongation_unit = aa_specifity.lower()
                        display_reactions(['1.png', '2.png'], list_tailoring_domains, elongation_unit, module_name, draw_mechanism_per_module)
                        # Necessary for thioesterase reactions when last module is an
                        # NRPS module
                        if module == modules[-1]:
                            chain_intermediate = copy_attached
            else:
                for domain in list_tailoring_domains:
                    if domain == 'E':
                        if draw_structures_per_module and attach_to_acp:
                            chain_intermediate = nrps_epimerization(
                                chain_intermediate)
                            if domain == list_tailoring_domains[-1]:
                                copy_chain_intermediate = chain_intermediate.deepcopy()
                                # copy_chain_intermediate.find_cycles()
                                copy_attached = attach_to_domain_nrp(
                                    copy_chain_intermediate, 'PCP')
                                copy_attached.refresh_structure()
                                copy_attached.set_connectivities()
                                copy_attached.find_cycles()
                                drawing = RaichuDrawer(copy_attached,
                                                       dont_show=True)
                                list_drawings_per_module.append([drawing])
                        elif visualization_mechanism and attach_to_acp:
                            copy_chain_intermediate = chain_intermediate.deepcopy()
                            copy_chain_intermediate.find_cycles()
                            copy_attached = attach_to_domain_nrp(
                                copy_chain_intermediate, 'PCP')
                            if domain == list_tailoring_domains[0]:
                                before_reaction_filename = '2.png'
                                after_reaction_filename = '3.png'
                                if path.exists(before_reaction_filename):
                                    os.remove(before_reaction_filename)
                                RaichuDrawer(copy_attached,
                                             save_png=before_reaction_filename,
                                             dpi=500)
                                list_filenames.append(before_reaction_filename)
                            elif len(list_tailoring_domains) == 2 and domain == list_tailoring_domains[1]:
                                after_reaction_filename = '4.png'

                            chain_intermediate = nrps_epimerization(
                                chain_intermediate)
                            copy_chain_intermediate = chain_intermediate.deepcopy()
                            copy_chain_intermediate.find_cycles()
                            copy_attached = attach_to_domain_nrp(
                                copy_chain_intermediate, 'PCP')

                            if path.exists(after_reaction_filename):
                                os.remove(after_reaction_filename)
                            RaichuDrawer(copy_attached, save_png=after_reaction_filename, dpi=500)
                            list_filenames.append(after_reaction_filename)
                            if domain == list_tailoring_domains[-1]:
                                elongation_unit = aa_specifity.lower()
                                display_reactions(list_filenames,
                                                  list_tailoring_domains,
                                                  elongation_unit, module_name,
                                                  draw_mechanism_per_module)
                                # Necessary for thioesterase reactions when last module is an
                                # NRPS module
                                if module == modules[-1]:
                                    chain_intermediate = copy_attached
                    elif domain == 'nMT':
                        if draw_structures_per_module and attach_to_acp:
                            chain_intermediate = nrps_methylation(
                                chain_intermediate)
                            if domain == list_tailoring_domains[-1]:
                                copy_chain_intermediate = chain_intermediate.deepcopy()
                                # copy_chain_intermediate.find_cycles()
                                copy_attached = attach_to_domain_nrp(
                                    copy_chain_intermediate, 'PCP')
                                copy_attached.refresh_structure()
                                copy_attached.set_connectivities()
                                copy_attached.find_cycles()
                                drawing = RaichuDrawer(copy_attached,
                                                       dont_show=True)
                                list_drawings_per_module.append([drawing])
                        elif visualization_mechanism and attach_to_acp:
                            copy_chain_intermediate = chain_intermediate.deepcopy()
                            copy_chain_intermediate.find_cycles()
                            copy_attached = attach_to_domain_nrp(
                                copy_chain_intermediate, 'PCP')
                            if domain == list_tailoring_domains[0]:
                                before_reaction_filename = '2.png'
                                after_reaction_filename = '3.png'
                                if path.exists(before_reaction_filename):
                                    os.remove(before_reaction_filename)
                                RaichuDrawer(copy_attached,
                                             save_png=before_reaction_filename,
                                             dpi=500)
                                list_filenames.append(before_reaction_filename)
                            elif len(list_tailoring_domains) == 2 and domain == \
                                    list_tailoring_domains[1]:
                                after_reaction_filename = '4.png'

                            chain_intermediate = nrps_methylation(
                                chain_intermediate)
                            copy_chain_intermediate = chain_intermediate.deepcopy()
                            copy_chain_intermediate.find_cycles()
                            copy_attached = attach_to_domain_nrp(
                                copy_chain_intermediate, 'PCP')

                            if path.exists(after_reaction_filename):
                                os.remove(after_reaction_filename)
                            RaichuDrawer(copy_attached,
                                         save_png=after_reaction_filename,
                                         dpi=500)
                            list_filenames.append(after_reaction_filename)
                            if domain == list_tailoring_domains[-1]:
                                elongation_unit = aa_specifity.lower()
                                display_reactions(list_filenames,
                                                  list_tailoring_domains,
                                                  elongation_unit, module_name,
                                                  draw_mechanism_per_module)
                                # Necessary for thioesterase reactions when last module is an
                                # NRPS module
                                if module == modules[-1]:
                                    chain_intermediate = copy_attached

    # Reset the atom color in the final structure to black
    for atom in chain_intermediate.graph:
        atom.draw.colour = 'black'

    if not visualization_mechanism and not draw_structures_per_module:
        if last_module_nrps and attach_to_acp:
            # chain_intermediate.find_cycles()
            chain_intermediate = attach_to_domain_nrp(chain_intermediate, 'PCP')
            chain_intermediate.find_cycles()
            RaichuDrawer(chain_intermediate)
        else:
            # chain_intermediate.find_cycles()
            RaichuDrawer(chain_intermediate)

    # Remove all intermediate structure files
    for i in range(1, 6):
        if path.exists(f"{i}.png"):
            os.remove(f"{i}.png")

    if draw_structures_per_module:
        return(list_drawings_per_module)

    return chain_intermediate


def display_reactions(structures, tailoring_domains, elongation_unit,
                      module_name, draw_mechanism_per_module):
    """Helper function cluster_to structure. Returns the visualization of
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
    if elongation_unit not in ALL_PKS_ELONGATION_UNITS:
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
            if len(tailoring_domains) == 1:
                before_elongation, after_elongation, after_tailoring = structures
                images = []
                # Read in images
                img1 = mpimg.imread(before_elongation)
                img2 = mpimg.imread(after_elongation)
                img3 = mpimg.imread(after_tailoring)

                # Add images to list to prevent collection by the garbage collector
                images.append(img1)
                images.append(img2)
                images.append(img3)

                # Create figure
                fig = plt.figure(constrained_layout=True, figsize=[15, 8],
                                 frameon=False,
                                 num=f"Quick reaction mechanism {module_name}")
                ax = fig.subplots(1, 8)
                ax = fig.add_gridspec(1, 8)

                # Add first structure before elongation reaction to plot
                ax1 = fig.add_subplot(ax[0, 0:2])
                ax1.imshow(img1)

                # Add correct arrow to plot (depending on elongation unit)
                ax1 = fig.add_subplot(ax[0, 2])
                arrow = FancyArrow(0, 0.5, 0.89, 0, overhang=0.3,
                                   head_width=0.025, head_length=0.1,
                                   color='black')
                images.append(arrow)
                ax1.add_patch(arrow)
                text = elongation_unit

                ax1.text(0.5, 0.55, text, ha='center',
                         fontdict={'size': 12, 'color': 'black'})

                # Add elongation product to plot
                ax1 = fig.add_subplot(ax[0, 3:5])
                ax1.set_title(f"{module_name}", pad=50,
                              fontdict={'size': 20, 'weight': 'bold',
                                        'color': 'darkred'})
                ax1.imshow(img2)

                # Add KR arrow to plot
                ax1 = fig.add_subplot(ax[0, 5])
                arrow = FancyArrow(0, 0.5, 0.89, 0, overhang=0.3,
                                   head_width=0.025, head_length=0.1,
                                   color='black')
                images.append(arrow)
                ax1.add_patch(arrow)
                ax1.text(0.5, 0.55, TAILOR_DOMAIN_SHORT_TO_LONG[tailoring_domains[0]], ha='center',
                         fontdict={'size': 12})

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
            elif len(tailoring_domains) == 2:
                print(len(structures), structures)
                before_elongation, after_elongation, after_first_tailoring, after_second_tailoring = structures
                images = []
                # Read in images
                img1 = mpimg.imread(before_elongation)
                img2 = mpimg.imread(after_elongation)
                img3 = mpimg.imread(after_first_tailoring)
                img4 = mpimg.imread(after_second_tailoring)

                # Add images to list to prevent collection by the garbage collector
                images.append(img1)
                images.append(img2)
                images.append(img3)
                images.append(img4)

                # Create figure
                fig = plt.figure(constrained_layout=True, figsize=[20, 8],
                                 frameon=False,
                                 num=f"Quick reaction mechanism {module_name}")
                ax = fig.subplots(1, 11)
                ax = fig.add_gridspec(1, 11)

                # Add first structure before elongation reaction to plot
                ax1 = fig.add_subplot(ax[0, 0:2])
                ax1.imshow(img1)

                # Add correct arrow to plot (depending on elongation unit)
                ax1 = fig.add_subplot(ax[0, 2])
                arrow = FancyArrow(0, 0.5, 0.89, 0, overhang=0.3,
                                   head_width=0.025, head_length=0.1,
                                   color='black')
                images.append(arrow)
                ax1.add_patch(arrow)
                text = elongation_unit

                ax1.text(0.5, 0.55, text, ha='center',
                         fontdict={'size': 12, 'color': 'black'})

                # Add elongation product to plot
                ax1 = fig.add_subplot(ax[0, 3:5])
                ax1.imshow(img2)

                # Add KR arrow to plot
                ax1 = fig.add_subplot(ax[0, 5])
                arrow = FancyArrow(0, 0.5, 0.89, 0, overhang=0.3,
                                   head_width=0.025, head_length=0.1,
                                   color='black')
                images.append(arrow)
                ax1.add_patch(arrow)
                ax1.text(0.5, 0.55, TAILOR_DOMAIN_SHORT_TO_LONG[tailoring_domains[0]], ha='center',
                         fontdict={'size': 12})

                # Add ketoreductase product to plot
                ax1 = fig.add_subplot(ax[0, 6:8])
                ax1.imshow(img3)

                # Add DH arrow to plot
                ax1 = fig.add_subplot(ax[0, 8])
                arrow = FancyArrow(0, 0.5, 0.89, 0, overhang=0.3,
                                   head_width=0.025, head_length=0.1,
                                   color='black')
                images.append(arrow)
                ax1.add_patch(arrow)
                ax1.text(0.5, 0.55, TAILOR_DOMAIN_SHORT_TO_LONG[tailoring_domains[1]], ha='center',
                         fontdict={'size': 12})

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
            text = ELONGATION_UNIT_TO_TEXT[elongation_unit]

            if elongation_unit == 'pk':
                text = 'Unknown elongation\nunit'

            ax1.text(0.5, 0.55, text, ha='center',
                     fontdict={'size': 12, 'color': 'black'})
            ax1.set_title(f"{module_name}", pad=70,
                          fontdict={'size': 20, 'weight': 'bold', 'color': 'darkred'})

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
            # Read in images
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
            text = ELONGATION_UNIT_TO_TEXT[elongation_unit]

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

            kr_domain = None

            for domain in tailoring_domains:
                if domain.startswith('KR'):
                    kr_domain = domain
                    break

            assert kr_domain
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
            text = ELONGATION_UNIT_TO_TEXT[elongation_unit]

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

            kr_domain = None

            for domain in tailoring_domains:
                if domain.startswith('KR'):
                    kr_domain = domain
                    break

            assert kr_domain
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
            text = ELONGATION_UNIT_TO_TEXT[elongation_unit]

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







