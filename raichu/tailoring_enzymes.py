from enum import Enum, unique
from raichu.reactions.general_tailoring_reactions import dephosphorylation, hydroxylation, addition, oxidative_bond_formation, epoxidation, double_bond_reduction, double_bond_shift
from pikachu.chem.structure import Structure
from pikachu.general import read_smiles, structure_to_smiles
from raichu.substrate import TerpeneCyclaseSubstrate
from raichu.data.attributes import PRENYL_TRANSFERASE_SUBSTRATES_TO_SMILES
@unique
class TailoringEnzymeType(Enum):
    METHYLTRANSFERASE = 1
    C_METHYLTRANSFERASE = 2
    N_METHYLTRANSFERASE = 3
    O_METHYLTRANSFERASE = 4
    P450_HYDROXYLATION = 5
    P450_OXIDATIVE_BOND_FORMATION = 6
    P450_EPOXIDATION = 7
    REDUCTASE_DOUBLE_BOND_REDUCTION = 8
    ISOMERASE_DOUBLE_BOND_SHIFT = 9
    PRENYLTRANSFERASE = 10
    ACETYLTRANSFERASE = 11
    ACYLTRANSFERASE = 12
    @staticmethod
    def from_string(label: str) -> "TailoringEnzymeType":
        for value in TailoringEnzymeType:
            if str(value.name) == label:
                return value
        raise ValueError(f"Unknown tailoring enzyme: {label}")


class TailoringEnzyme:

    def __init__(self, gene_name, enzyme_type, atoms:list = None, substrate:str = None) -> None:
        self.gene_name = gene_name
        self.type = TailoringEnzymeType.from_string(enzyme_type)
        self.atoms = atoms
        self.substrate = substrate

    def do_tailoring(self, structure):
        """
        Performs tailoring reaction
        """
        if self.type.name == "P450_HYDROXYLATION":
            for atom in self.atoms:
                atom = atom[0] #only one atom is hydroxylated at a time
                atom = structure.get_atom(atom)
                structure = addition(atom, "O", structure)
        elif self.type.name in ["METHYLTRANSFERASE", "C_METHYLTRANSFERASE", "N_METHYLTRANSFERASE", "O_METHYLTRANSFERASE"]:
            for atom in self.atoms:
                atom = atom[0] #only one atom is methylated at a time
                atom = structure.get_atom(atom)
                structure = addition(atom, "C", structure)
        elif self.type.name == "PRENYLTRANSFERASE":
            for atom in self.atoms:
                atom = atom[0]  # only one atom is methylated at a time
                atom = structure.get_atom(atom)
                if self.substrate not in PRENYL_TRANSFERASE_SUBSTRATES_TO_SMILES:
                    raise ValueError(
                        f"Not implemented prenyltransferase substrate: {self.substrate}")
                substrate = PRENYL_TRANSFERASE_SUBSTRATES_TO_SMILES[self.substrate]
                structure = addition(
                    atom, substrate, structure)
        elif self.type.name == "ACETYLTRANSFERASE":
            for atom in self.atoms:
                atom = atom[0]  # only one atom is methylated at a time
                atom = structure.get_atom(atom)
                structure = addition(atom, "[H]C(C)=O", structure)
        elif self.type.name == "ACYLTRANSFERASE":
            for atom in self.atoms:
                atom = atom[0]  # only one atom is methylated at a time
                atom = structure.get_atom(atom)
                if self.substrate:
                    structure = addition(atom, self.substrate, structure)
        elif self.type.name == "P450_OXIDATIVE_BOND_FORMATION":
            for atoms in self.atoms:
                atom1 = structure.get_atom(atoms[0])
                atom2 = structure.get_atom(atoms[1])
                structure = oxidative_bond_formation(atom1, atom2, structure)
        elif self.type.name == "P450_EPOXIDATION":
            for atoms in self.atoms:
                atom1 = structure.get_atom(atoms[0])
                atom2 = structure.get_atom(atoms[1])
                structure = epoxidation(atom1, atom2, structure)
        elif self.type.name == "REDUCTASE_DOUBLE_BOND_REDUCTION":
            for atoms in self.atoms:
                atom1 = structure.get_atom(atoms[0])
                atom2 = structure.get_atom(atoms[1])
                structure = double_bond_reduction(atom1, atom2, structure)
        elif self.type.name == "ISOMERASE_DOUBLE_BOND_SHIFT":
            for atoms in self.atoms:
                old_double_bond_atom1 = structure.get_atom(atoms[0])
                old_double_bond_atom2 = structure.get_atom(atoms[1])
                new_double_bond_atom1 = structure.get_atom(atoms[2])
                new_double_bond_atom2 = structure.get_atom(atoms[3])
                structure = double_bond_shift(
                    structure, old_double_bond_atom1, old_double_bond_atom2, new_double_bond_atom1, new_double_bond_atom2)
        
        
        return structure
