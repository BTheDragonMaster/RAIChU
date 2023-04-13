from pikachu.reactions.functional_groups import find_atoms, find_bonds
from pikachu.reactions.basic_reactions import condensation
from raichu.central_chain_detection.label_central_chain import label_nrp_central_chain
from raichu.reactions.general import label_rest_groups, initialise_atom_attributes, reset_nrp_annotations
from raichu.data.molecular_moieties import C_TERMINAL_OH, N_AMINO_ACID, B_N_AMINO_ACID

def ribosomal_elongation(amino_acid, chain_intermediate, amino_acid_number=None):
    """
    Returns the ribosomal condensation product as a PIKAChU Structure object
    
    amino_acid_number: Type and number of amino acid in peptide (e.g. A14)
    amino_acid: PIKAChU Structure object of the amino acid
    chain_intermediate: PIKAChU Structure object of the ripp intermediate
    """

    # Initialize annotations in all atoms in amino acid
    initialise_atom_attributes(amino_acid)
    for atom in amino_acid.atoms.values():
        atom.amino_acid = amino_acid_number
    amino_acid.atoms
    oh_bond = find_bonds(C_TERMINAL_OH, chain_intermediate)
    assert len(oh_bond)==1
    oh_bond = oh_bond[0]

    # Label rest groups
    label_rest_groups(amino_acid, chain_intermediate)


    # Find amino group of amino acid or beta-amino acid
    n_atoms_aa = find_atoms(N_AMINO_ACID, amino_acid)
    if not n_atoms_aa:
        n_atoms_aa = find_atoms(B_N_AMINO_ACID, amino_acid)

    assert len(n_atoms_aa) == 1

    n_atom = n_atoms_aa[0]
    h_bond = None

    for bond in n_atom.bonds:
        for neighbour in bond.neighbours:
            if neighbour.type == 'H':
                h_bond = bond
                break

    assert h_bond

    # Label atoms in amino acid that end up in central peptide chain
    label_nrp_central_chain(amino_acid)

    # Reset chiral_c_ep and n_atom_nmeth AtomAnnotations attribute for all atoms in the intermediate;
    # This way, only the newly incorporated amino acid is labelled for tailoring
    reset_nrp_annotations(chain_intermediate)

    # Carry out condensation reaction using build-in PIKAChU function
    condensation_product = condensation(chain_intermediate, amino_acid, oh_bond, h_bond)[0]

    # Refresh condensation product
    condensation_product.refresh_structure(find_cycles=True)

    # Initialize annotations for atoms that don't have any yet
    initialise_atom_attributes(condensation_product)

    return condensation_product