from pikachu.general import read_smiles, draw_structure
from pikachu.reactions.functional_groups import combine_structures, find_atoms

from raichu.data.attributes import ATTRIBUTES
from raichu.reactions.general_tailoring_reactions import cyclodehydration, single_bond_oxidation
from raichu.data.molecular_moieties import ATTACHED_SERINE_OX, ATTACHED_CYSTEINE_OX, ATTACHED_SERINE_O, \
    ATTACHED_CYSTEINE_S, CYCLIC_CYS_C1_ATTACHED, CYCLIC_CYS_C2_ATTACHED, CYCLIC_SER_C1_ATTACHED, CYCLIC_SER_C2_ATTACHED


def epimerization(chiral_centre):
    new_chirality = None
    if chiral_centre.chiral == 'clockwise':
        new_chirality = 'counterclockwise'
    elif chiral_centre.chiral == 'counterclockwise':
        new_chirality = 'clockwise'

    chiral_centre.chiral = new_chirality


def epimerize(nrp):
    """

    """
    did_reaction = False
    # Define epimerization reaction target
    chiral_c_ep_atoms = []
    for atom in nrp.graph:
        if atom.annotations.chiral_c_ep:
            chiral_c_ep_atoms.append(atom)

    # If the substrate is a (fatty) acid, epimerization is not possible
    # TODO: Make verbose mode that prints this
    # if len(chiral_c_ep_atoms) == 0:
    #    print('Warning: Cannot perform epimerization reaction on non-amino acid substrate!')

    assert len(chiral_c_ep_atoms) < 2

    # Carry out epimerization
    if len(chiral_c_ep_atoms) == 1:
        chiral_c = chiral_c_ep_atoms[0]
        epimerization(chiral_c)
        if chiral_c.chiral:
            did_reaction = True

    return nrp, did_reaction


def methylation(target_atom, structure):

    methyl_group = read_smiles('C')
    methyl_group.add_attributes(ATTRIBUTES, boolean=True)
    carbon = methyl_group.atoms[0]
    hydrogen_1 = methyl_group.atoms[1]
    bond_1 = carbon.get_bond(hydrogen_1)

    hydrogen_2 = target_atom.get_neighbour('H')
    if not hydrogen_2:
        raise Exception("Can't methylate this atom!")

    bond_2 = target_atom.get_bond(hydrogen_2)

    methylated_structure = combine_structures([structure, methyl_group])

    methylated_structure.break_bond(bond_1)
    methylated_structure.break_bond(bond_2)

    methylated_structure.make_bond(carbon, target_atom, methylated_structure.find_next_bond_nr())
    methylated_structure.make_bond(hydrogen_1, hydrogen_2, methylated_structure.find_next_bond_nr())

    structures = methylated_structure.split_disconnected_structures()

    for s in structures:
        if carbon.nr in s.atoms:
            return s


def n_methylate(nrp):
    """

    """
    # Define N-methylation reaction target
    n_meth_locations = []
    for atom in nrp.graph:
        if atom.annotations.n_atom_nmeth:
            n_meth_locations.append(atom)

    # If the substrate is not an amino acid, N-methylation is not possible
    # TODO: Add a verbose mode that prints this
    # if len(n_meth_locations) == 0:
    #     print('Warning: Cannot perform N-methylation on a non-amino acid substrate!')

    assert len(n_meth_locations) < 2

    product = None

    if len(n_meth_locations) == 1:

        n_meth = n_meth_locations[0]

        # Check if the N atom has a hydrogen group necessary for the reaction, and
        # not a cyclic amiono acid such as proline
        executable = True
        if not n_meth.has_neighbour('H'):
            executable = False
            product = nrp

        # Carry out N-methylation
        if executable:
            product = methylation(n_meth, nrp)

    else:
        product = nrp

    assert product

    return product, True


def nrps_oxidation(nrp):
    cyclic_cys_c1s = find_atoms(CYCLIC_CYS_C1_ATTACHED, nrp)
    if len(cyclic_cys_c1s) == 1:
        cyclic_cys_c2s = find_atoms(CYCLIC_CYS_C2_ATTACHED, nrp)
        assert len(cyclic_cys_c2s) == 1
        atom_1 = cyclic_cys_c1s[0]
        atom_2 = cyclic_cys_c2s[0]
    else:
        cyclic_ser_c1s = find_atoms(CYCLIC_SER_C1_ATTACHED, nrp)
        if len(cyclic_ser_c1s) == 1:
            cyclic_ser_c2s = find_atoms(CYCLIC_SER_C2_ATTACHED, nrp)
            assert len(cyclic_ser_c2s) == 1
            atom_1 = cyclic_ser_c1s[0]
            atom_2 = cyclic_ser_c2s[0]
        else:
            return nrp, False

    product = single_bond_oxidation(atom_1, atom_2, nrp)
    assert product
    return product, True


def nrps_cyclodehydration(nrp):
    cysteine_ss = find_atoms(ATTACHED_CYSTEINE_S, nrp)
    if len(cysteine_ss) == 1:
        cysteine_oxs = find_atoms(ATTACHED_CYSTEINE_OX, nrp)
        assert len(cysteine_oxs) == 1
        attacking_atom = cysteine_ss[0]
        keto_group = cysteine_oxs[0]
        carbon = keto_group.get_neighbour("C")
        nitrogen = carbon.get_neighbour("N")
    else:
        serine_os = find_atoms(ATTACHED_SERINE_O, nrp)
        if len(serine_os) == 1:
            serine_oxs = find_atoms(ATTACHED_SERINE_OX, nrp)
            attacking_atom = serine_os[0]
            keto_group = serine_oxs[0]
            carbon = keto_group.get_neighbour("C")
            nitrogen = carbon.get_neighbour("N")

        else:
            return nrp, False

    if nitrogen.has_neighbour('H') and attacking_atom.has_neighbour('H'):

        product = cyclodehydration(nrp, attacking_atom, keto_group)
        assert product

        return product, True
    else:
        return nrp, False


if __name__ == "__main__":
    s = read_smiles(r"NCC(=O)NCC(=O)O")
    target_atom = s.atoms[4]
    structure = methylation(target_atom, s)
    draw_structure(structure)
