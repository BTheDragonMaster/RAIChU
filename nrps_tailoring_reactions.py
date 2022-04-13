from pikachu.general import read_smiles, draw_structure
from pikachu.reactions.functional_groups import combine_structures
from attributes import ATTRIBUTES



def epimerization(chiral_centre):
    new_chirality = None
    if chiral_centre.chiral == 'clockwise':
        new_chirality = 'counterclockwise'
    elif chiral_centre.chiral == 'counterclockwise':
        new_chirality = 'clockwise'

    chiral_centre.chiral = new_chirality

def nrps_epimerization(nrp):
    """

    """
    # Define epimerization reaction target
    chiral_c_ep_atoms = []
    for atom in nrp.graph:
        if atom.annotations.chiral_c_ep:
            chiral_c_ep_atoms.append(atom)
    assert len(chiral_c_ep_atoms) == 1
    chiral_c = chiral_c_ep_atoms[0]

    # Carry out epimerization
    epimerization(chiral_c)

    return nrp


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

def nrps_methylation(nrp):
    """

    """
    # Define N-methylation reaction target
    n_meth_locations = []
    for atom in nrp.graph:
        if atom.annotations.n_atom_nmeth:
            n_meth_locations.append(atom)
    assert len(n_meth_locations) == 1
    n_meth = n_meth_locations[0]

    # Carry out N-methylation
    product = methylation(n_meth, nrp)

    return product

if __name__ == "__main__":
    s = read_smiles(r"NCC(=O)NCC(=O)O")
    target_atom = s.atoms[4]
    structure = methylation(target_atom, s)
    draw_structure(structure)





