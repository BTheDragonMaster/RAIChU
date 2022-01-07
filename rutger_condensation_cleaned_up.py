from pikachu.reactions.functional_groups import find_bonds, find_atoms, BondDefiner
from rutger_groupdefiner import *
from rutger_combine_structure import *
from pikachu.drawing.drawing import *
from copy import deepcopy
from pk_attach_to_domain import attach_to_domain
from pikachu.general import highlight_subsmiles_single, structure_to_smiles

def condensation_nrps(chain_intermediate, amino_acid):
    """
    Adapted from Rutger's domain_condensate() function
    Returns the product of the NRPS condensation reaction as a Structure object

    chain_intermediate: Structure object NRP chain intermediate bound to PCP
    domain (hydroxyl group in carboxylic acid group replaced by sulphur atom)

    amino_acid: Structure of the amino acid that should be attached to the
    chain intermediate, bound to PCP domain (hydroxyl group in carboxylic acid
    group replaced by sulphur atom)
    """

    # Give chain_intermediate and amino_acid unique Bond.nrs and Atom.nrs
    chain_intermediate, amino_acid = correct_structures_index\
        ([chain_intermediate, amino_acid])

    # Find the carbon that is bound to the sulphur atom in the chain interm.
    carbon_to_s = find_group(chain_intermediate, PCP_BOUND_PEPTIDE_CARBON)
    assert len(carbon_to_s) == 1
    carbon_to_s = carbon_to_s[0]

    # Identify and remove S from chain intermediate to allow condensation
    sulphur = find_group(chain_intermediate, PCP_BOUND_PEPTIDE_SULPHUR)
    assert len(sulphur) == 1
    sulphur = sulphur[0]

    chain_intermediate.break_bond_between_atoms(carbon_to_s, sulphur)
    one, two = chain_intermediate.split_disconnected_structures()


    if len(one.graph) == 2:
        single_sulphur = one
        chain_intermediate = two
    elif len(two.graph) == 2:
        single_sulphur = two
        chain_intermediate = one

    # Find N-atom amino group in amino acid, remove one proton from atom
    n_amino_acid = find_group(amino_acid, PEPTIDE_END_N)
    assert len(n_amino_acid) == 1
    n_amino_acid = n_amino_acid[0]

    for neighbour_n in n_amino_acid.neighbours:
        if neighbour_n.type == 'H':
            amino_acid.break_bond_between_atoms(n_amino_acid, neighbour_n)
            h_amino_group = neighbour_n
            amino_acid.remove_atom(h_amino_group)
            break

    # Combine graphs chain intermediate and amino acid
    combined = combine_graphs([chain_intermediate, amino_acid])

    # Form peptide bond
    next_bond_nr = combined.find_next_bond_nr()
    combined.make_bond(carbon_to_s, n_amino_acid, next_bond_nr)
    final_product = combined

    return final_product


def make_nrp(list_amino_acids):
    """
    Returns the NRP product as Structure object

    list_amino_acids: [strings] of amino acid names
    """
    # Parse list SMILES amino acids attached to PCP to dict name -> Structure
    lines_aa = open('smiles_aa_attached.txt', 'r').readlines()
    dict_aa_structure = {}

    for line in lines_aa:
        line = line.strip()
        name, smiles = line.split()
        structure_aa = Smiles(smiles).smiles_to_structure()
        dict_aa_structure[name] = structure_aa

    # Take amino acid Structure object from list and add to growing NRP chain
    nrp_chain_intermediate = deepcopy(dict_aa_structure[list_amino_acids[0]])
    list_amino_acids = list_amino_acids[1:]
    for amino_acid_name in list_amino_acids:
        amino_acid_struct = deepcopy(dict_aa_structure[amino_acid_name])
        nrp_chain_intermediate = condensation_nrps(nrp_chain_intermediate, amino_acid_struct)

    # Refresh chain intermediate
    nrp_chain_intermediate.refresh_structure()
    nrp_chain_intermediate.find_cycles()

    return nrp_chain_intermediate






if __name__ == "__main__":
    peptide = make_nrp(['Tryptophan', 'Arginine', 'Histidine', 'Glutamicacid',\
                        'Proline', 'Lysine', 'Serine', 'Tryptophan', 'Cysteine',\
                        'Valine', 'Methionine'])
    #attach_to_domain(peptide, 'PCP') #Doesn't work
    peptide.refresh_structure()
    peptide.make_bond_lookup()
    Drawer(peptide)
    plt.show()

