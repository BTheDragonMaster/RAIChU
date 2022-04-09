from class_domain import *
from pikachu.reactions.functional_groups import find_atoms, GroupDefiner, combine_structures
from pikachu.reactions.basic_reactions import condensation


POLYKETIDE_S = GroupDefiner('Sulphur atom polyketide', 'SC(C)=O', 0)
NRP_C = GroupDefiner('C atom to attach to PCP domain', 'NCC(O)=O', 2)


def attach_to_domain_pk(polyketide, domain_type):
    """
    Attaches the sulphur atom in the input polyketide to a PKS domain and
    returns the attached structure as a PIKAChU Structure object

    domain_type: Str, domain type
    polyketide: PIKAChU Structure object, to-be attached structure
    """
    #Create domain

    domain = make_scaffold_domain(domain_type)
    sh_bond = domain.bond_lookup[domain.atoms[1]][domain.atoms[2]]
    hydrogen = domain.atoms[2]
    sulphur_1 = domain.atoms[1]

    locations_sulphur = find_atoms(POLYKETIDE_S, polyketide)
    assert len(locations_sulphur) == 1
    sulphur_2 = locations_sulphur[0]
    carbon = sulphur_2.get_neighbour('C')
    sc_bond = sulphur_2.get_bond(carbon)

    structure = combine_structures([domain, polyketide])
    structure.break_bond(sh_bond)
    structure.break_bond(sc_bond)

    structure.make_bond(hydrogen, sulphur_2, structure.find_next_bond_nr())
    structure.make_bond(carbon, sulphur_1, structure.find_next_bond_nr())

    split = structure.split_disconnected_structures()

    tethered_polyketide = None

    for structure in split:
        if carbon in structure.graph:
            tethered_polyketide = structure
            break

    assert tethered_polyketide

    return tethered_polyketide


def attach_to_domain_nrp(nrp, domain_type):
    """
    Attaches the input NRP to a PCP domain and returns the product as a
    PIKAChU Structure object

    domain_type: Str, domain type
    nrp: PIKAChU Structure object, to-be attached NRP
    """
    # Create domain

    domain = make_scaffold_domain(domain_type)

    # Remove OH group from carboxylic acid in NRP, to allow attachment
    # to domain
    locations_c_to_domain = find_atoms(NRP_C, nrp)
    assert len(locations_c_to_domain) == 1
    c_atom_to_domain = locations_c_to_domain[0]
    oxygens = c_atom_to_domain.get_neighbours('O')

    hydroxyl_oxygen = None
    hydroxyl_bond = None

    hydrogen_bond = domain.bond_lookup[domain.atoms[1]][domain.atoms[2]]

    assert hydrogen_bond.has_neighbour('S')
    assert hydrogen_bond.has_neighbour('H')

    for oxygen in oxygens:
        bond = c_atom_to_domain.get_bond(oxygen)
        if bond.type == 'single' and oxygen.has_neighbour('H'):
            hydroxyl_oxygen = oxygen
            hydroxyl_bond = bond
            break

    assert hydroxyl_oxygen and hydroxyl_bond

    structure = condensation(nrp, domain, hydroxyl_bond, hydrogen_bond)[0]

    for atom in structure.graph:
        if not hasattr(atom.annotations, 'in_central_chain'):
            for attribute in ATTRIBUTES:
                atom.annotations.add_annotation(attribute, False)

    return structure

