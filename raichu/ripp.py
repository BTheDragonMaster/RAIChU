from pikachu.chem.structure import Structure
from pikachu.reactions.functional_groups import find_bonds
from pikachu.general import read_smiles

from raichu.substrate import RibosomeSubstrate
from raichu.reactions.nrps_elongation_reactions import nrps_elongation
from raichu.data.molecular_moieties import PEPTIDE_BOND
from raichu.data.attributes import AMINOACID_ONE_LETTER_TO_NAME
from raichu.reactions.general_tailoring_reactions import proteolytic_cleavage, cyclisation, ribosomal_elongation
from raichu.tailoring_enzymes import TailoringEnzyme
from raichu.drawing.drawer import RaichuDrawer

class RiPP_Cluster:
    def __init__(self, gene_name_precursor: str, amino_acid_sequence: str, cleavage_sites: list() = None, macrocyclisations: list() = None, tailoring_enzymes_representation = None) -> None:
        self.gene_name = gene_name_precursor
        self.amino_acid_seqence = amino_acid_sequence.upper()
        self.cleavage_sites = cleavage_sites
        self.macrocyclisations = macrocyclisations
        self.tailoring_enzymes_representation = tailoring_enzymes_representation
        self.chain_intermediate = None
        self.linear_product = None
        self.cyclised_product = None
        self.final_product = None
        self.tailoring_enzymes = []
        
    def make_peptide(self):
        for index, amino_acid in enumerate(self.amino_acid_seqence):
            name_amino_acid = f'{amino_acid}{index+1}'
            if amino_acid in AMINOACID_ONE_LETTER_TO_NAME:
                name = AMINOACID_ONE_LETTER_TO_NAME[amino_acid]
            else:
                raise ValueError(f"Unknown amino acid: {amino_acid}")
            substrate = RibosomeSubstrate(name)
            building_block = read_smiles(substrate.smiles)
            if self.chain_intermediate:
                self.chain_intermediate = ribosomal_elongation(building_block, self.chain_intermediate, amino_acid_number=name_amino_acid)
            else:
                self.chain_intermediate = building_block
        self.linear_product = self.chain_intermediate
    
    def initialize_cleavage_sites_on_structure(self) -> list:
        cleavage_bonds = []
        for cleavage_site in self.cleavage_sites:
            amino_acid_cleavage = cleavage_site.position_amino_acid
            number_cleavage = cleavage_site.position_index
            if self.amino_acid_seqence[number_cleavage-1] == amino_acid_cleavage:
                peptide_bonds = find_bonds(PEPTIDE_BOND, self.chain_intermediate)
                peptide_bonds = sorted(peptide_bonds, key = lambda bond: bond.nr, reverse= True)
                cleavage_bond = peptide_bonds[number_cleavage] 
                cleavage_bonds += [[cleavage_bond, cleavage_site.structure_to_keep]]
            else:
                raise ValueError(f"No {AMINOACID_ONE_LETTER_TO_NAME[amino_acid_cleavage]} in position {number_cleavage} for cleavage.")
        return cleavage_bonds
            
    
    def do_proteolytic_claevage(self):
        cleavage_bonds = self.initialize_cleavage_sites_on_structure()
        for bond, structure_to_keep in cleavage_bonds:
            self.chain_intermediate = proteolytic_cleavage(bond, self.chain_intermediate, structure_to_keep=structure_to_keep)
        self.final_product = self.chain_intermediate
     
     
    def initialize_macrocyclization_on_structure(self) -> list(list()):
        initialized_macrocyclization_atoms = []
        if self.macrocyclisations:
            for cyclization in self.macrocyclisations:
                atoms = [atom for atom in self.chain_intermediate.atoms.values() if str(atom) in [cyclization.atom1, cyclization.atom2]]
                if len(atoms)<2:
                    raise ValueError(f"Non-existing atoms for cyclization")
                initialized_macrocyclization_atoms += [atoms]
        return initialized_macrocyclization_atoms


    def do_macrocyclization(self):
        macrocyclizations_atoms = self.initialize_macrocyclization_on_structure()
        for macrocyclization_atoms in macrocyclizations_atoms:
            self.chain_intermediate = cyclisation(self.chain_intermediate, macrocyclization_atoms[0], macrocyclization_atoms[1])
        self.cyclised_product = self.chain_intermediate
     
        
    def initialize_tailoring_enzymes_on_structure(self):
        if self.tailoring_enzymes_representation:
            for tailoring_enzyme_representation in self.tailoring_enzymes_representation:
                atom_array = []
                for atoms_for_reaction in tailoring_enzyme_representation.atoms:
                    atoms_for_reaction_initialized = [atom for atom in self.chain_intermediate.atoms.values() if str(atom) in atoms_for_reaction]
                    if len(atoms_for_reaction_initialized)<len(atoms_for_reaction):
                        raise ValueError(f"Non-existing atoms for tailoring")
                    atom_array += [atoms_for_reaction_initialized]
                self.tailoring_enzymes += [TailoringEnzyme(tailoring_enzyme_representation.gene_name, tailoring_enzyme_representation.type, atom_array)]

    def do_tailoring(self):
        self.initialize_tailoring_enzymes_on_structure()
        for tailoring_enzyme in self.tailoring_enzymes:
            self.chain_intermediate = tailoring_enzyme.do_tailoring(self.chain_intermediate)
            
            
    def draw_product(self):
        assert self.chain_intermediate
        drawing = RaichuDrawer(self.chain_intermediate, dont_show=False)
        drawing.draw_structure()
        svg_string = drawing.save_svg_string()
        
        return svg_string
