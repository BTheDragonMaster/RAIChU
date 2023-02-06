from pikachu.chem.structure import Structure
from pikachu.general import read_smiles
from raichu.substrate import TerpeneCyclaseSubstrate
from raichu.reactions.general_tailoring_reactions import dephosphorylation, oxidative_bond_formation
from raichu.tailoring_enzymes import TailoringEnzyme
from raichu.drawing.drawer import RaichuDrawer

class Terpene_Cluster:
    def __init__(self, gene_name_cyclase: str, precursor: str, terpene_cyclase_type: str = None, macrocyclisations: list() = None, tailoring_enzymes_representation = None) -> None:
        self.gene_name = gene_name_cyclase
        self.terpene_cyclase_type = terpene_cyclase_type
        self.precursor = precursor
        self.macrocyclisations = macrocyclisations
        self.tailoring_enzymes_representation = tailoring_enzymes_representation
        self.chain_intermediate = None
        self.tailored_product = None
        self.cyclised_product = None
        self.final_product = None
        self.tailoring_enzymes = []
        self.initialized_macrocyclization_atoms = []
        
    def create_precursor(self) -> Structure:
        substrate = TerpeneCyclaseSubstrate(self.precursor)
        self.chain_intermediate = read_smiles(substrate.smiles)
        self.initialize_macrocyclization_on_structure()
        self.initialize_tailoring_enzymes_on_structure()

    def initialize_macrocyclization_on_structure(self) -> list(list()):
        if self.macrocyclisations:
            for cyclization in self.macrocyclisations:
                atoms = [atom for atom in self.chain_intermediate.atoms.values() if str(
                    atom) in [cyclization.atom1, cyclization.atom2]]
                if len(atoms)<2:
                    raise ValueError(f"Non-existing atoms for cyclization")
                self.initialized_macrocyclization_atoms += [atoms]


    def do_macrocyclization(self):
        if self.terpene_cyclase_type != "Class_2":
            self.chain_intermediate = dephosphorylation(self.chain_intermediate)
        for macrocyclization_atoms in self.initialized_macrocyclization_atoms:
            atom1 = self.chain_intermediate.get_atom(macrocyclization_atoms[0])
            atom2 = self.chain_intermediate.get_atom(macrocyclization_atoms[1])
            self.chain_intermediate = oxidative_bond_formation(
                atom1, atom2, self.chain_intermediate)
        self.cyclised_product = self.chain_intermediate
     
        
    def initialize_tailoring_enzymes_on_structure(self):
        if self.tailoring_enzymes_representation:
            for tailoring_enzyme_representation in self.tailoring_enzymes_representation:
                atom_array = []
                for atoms_for_reaction in tailoring_enzyme_representation.atoms:
                    atoms_for_reaction_initialized = [
                        atom for atom in self.chain_intermediate.atoms.values() if str(atom) in atoms_for_reaction]
                    atoms_for_reaction_initialized_updated = []
                    for atom in atoms_for_reaction:
                        for atom_initialized in atoms_for_reaction_initialized:
                            if str(atom_initialized) == atom:
                                atoms_for_reaction_initialized_updated.append(
                                    atom_initialized)
                    atoms_for_reaction_initialized = atoms_for_reaction_initialized_updated
                    if len(atoms_for_reaction_initialized)<len(atoms_for_reaction):
                        raise ValueError(f"Non-existing atoms for tailoring")
                    atom_array += [atoms_for_reaction_initialized]
                self.tailoring_enzymes += [TailoringEnzyme(tailoring_enzyme_representation.gene_name, tailoring_enzyme_representation.type, atom_array)]

    def do_tailoring(self):
        for tailoring_enzyme in self.tailoring_enzymes:
            self.tailored_product = tailoring_enzyme.do_tailoring(self.chain_intermediate)
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

