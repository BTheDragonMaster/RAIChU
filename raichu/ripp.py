from pikachu.reactions.functional_groups import find_bonds
from pikachu.general import read_smiles


from raichu.data.molecular_moieties import PEPTIDE_BOND
from raichu.data.attributes import AMINOACID_ONE_LETTER_TO_NAME, AMINOACID_ONE_LETTER_TO_SMILES
from raichu.reactions.general_tailoring_reactions import proteolytic_cleavage, cyclisation
from raichu.tailoring_enzymes import TailoringEnzyme
from raichu.drawing.drawer import RaichuDrawer

class RiPP_Cluster:
    def __init__(self, gene_name_precursor: str, amino_acid_sequence: str, cleavage_sites: list() = None, macrocyclisations: list() = None, tailoring_enzymes_representation = None) -> None:
        self.gene_name = gene_name_precursor
        self.amino_acid_seqence = amino_acid_sequence.upper()
        self.cleavage_sites = cleavage_sites
        self.cleavage_bonds = []
        self.macrocyclisations = macrocyclisations
        self.tailoring_enzymes_representation = tailoring_enzymes_representation
        self.chain_intermediate = None
        self.linear_product = None
        self.tailored_product = None
        self.cyclised_product = None
        self.final_product = None
        self.initialized_macrocyclization_atoms = []
        
    def make_peptide(self):
        smiles_peptide_chain = ""
        for amino_acid in self.amino_acid_seqence:
            if amino_acid in AMINOACID_ONE_LETTER_TO_SMILES:
                substrate = AMINOACID_ONE_LETTER_TO_SMILES[amino_acid]
            else:
                raise ValueError(f"Unknown amino acid: {amino_acid}")
            smiles_peptide_chain += str(substrate)
        smiles_peptide_chain += "O"
        self.linear_product = read_smiles(smiles_peptide_chain)
        self.chain_intermediate = self.linear_product
        
    
    def initialize_cleavage_sites_on_structure(self) -> list: 
        for cleavage_site in self.cleavage_sites:
            amino_acid_cleavage = cleavage_site.position_amino_acid
            number_cleavage = cleavage_site.position_index
            if self.amino_acid_seqence[number_cleavage-1] == amino_acid_cleavage:
                peptide_bonds = find_bonds(PEPTIDE_BOND, self.linear_product)
                peptide_bonds = sorted(peptide_bonds, key = lambda bond: bond.nr)
                cleavage_bond = peptide_bonds[number_cleavage] 
                self.cleavage_bonds += [[cleavage_bond, cleavage_site.structure_to_keep]]
            else:
                raise ValueError(f"No {AMINOACID_ONE_LETTER_TO_NAME[amino_acid_cleavage]} in position {number_cleavage} for cleavage.")
            
    
    def do_proteolytic_claevage(self):
        self.initialize_cleavage_sites_on_structure()
        for bond, structure_to_keep in self.cleavage_bonds:
            self.chain_intermediate = proteolytic_cleavage(bond, self.chain_intermediate, structure_to_keep=structure_to_keep)
        self.final_product = self.chain_intermediate
     
     
    def initialize_macrocyclization_on_structure(self) -> list(list()):
        if self.macrocyclisations:
            for cyclization in self.macrocyclisations:
                atoms = [atom for atom in self.chain_intermediate.atoms.values() if str(atom) in [cyclization.atom1, cyclization.atom2]]
                self.initialized_macrocyclization_atoms += [atoms]


    def do_macrocyclization(self):
        self.initialize_macrocyclization_on_structure()
        for macrocyclization_atoms in self.initialized_macrocyclization_atoms:
            if [str(atom) for atom in self.initialized_macrocyclization_atoms] != self.macrocyclisations:
                raise ValueError(
                    f'Not all atoms {str(self.initialized_macrocyclization_atoms)} for macrocyclisation exist in the structure.')
            atom1 = self.chain_intermediate.get_atom(macrocyclization_atoms[0])
            atom2 = self.chain_intermediate.get_atom(macrocyclization_atoms[1])
            self.chain_intermediate = cyclisation(self.chain_intermediate, atom1, atom2)
        self.cyclised_product = self.chain_intermediate
     
    def initialize_modification_sites_on_structure(self, modification_sites):
            modification_sites_initialized = []
            for atoms_for_reaction in modification_sites:
                atoms_for_reaction_with_numbers = map(
                    lambda atom: [int(''.join(filter(str.isdigit, atom))), atom], atoms_for_reaction)
                atoms_in_structure = list(map(
                    str, self.chain_intermediate.atoms.values()))
                atoms_for_reaction_initialized = [self.chain_intermediate.atoms[atom[0]]
                                                    for atom in atoms_for_reaction_with_numbers if atom[1] in atoms_in_structure]
                atoms_for_reaction_initialized = list(
                    filter(lambda atom: atom is not None, atoms_for_reaction_initialized))
                modification_sites_initialized += [
                    atoms_for_reaction_initialized]
            return modification_sites_initialized

    def do_tailoring(self):
        if self.tailoring_enzymes_representation:
            for tailoring_enzyme_representation in self.tailoring_enzymes_representation:
                modification_sites = self.initialize_modification_sites_on_structure(
                    tailoring_enzyme_representation.modification_sites)
                if [[str(atom) for atom in atoms_for_reaction]for atoms_for_reaction in modification_sites] != tailoring_enzyme_representation.modification_sites:
                    raise ValueError(
                        f'Not all atoms {tailoring_enzyme_representation.modification_sites} for {tailoring_enzyme_representation.type} exist in the structure.')
                tailoring_enzyme = TailoringEnzyme(
                    tailoring_enzyme_representation.gene_name, tailoring_enzyme_representation.type, modification_sites, tailoring_enzyme_representation.substrate)
                self.tailored_product = tailoring_enzyme.do_tailoring(
                    self.chain_intermediate)
                self.chain_intermediate = self.tailored_product
            

    def draw_product(self, as_string=True, out_file=None):
            assert self.chain_intermediate
            drawing = RaichuDrawer(self.chain_intermediate, dont_show=True, add_url=True, draw_Cs_in_pink=False, draw_straightened=False)
            drawing.draw_structure()
            svg_string = drawing.save_svg_string()
            if as_string:
                return svg_string
            else:
                if out_file is None:
                    raise ValueError("Must provide output svg directory if 'as_string' is set to False.")
                else:
                    with open(out_file, 'w') as svg_out:
                        svg_out.write(svg_string)

