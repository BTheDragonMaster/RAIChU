from enum import Enum, unique
from raichu.reactions.general_tailoring_reactions import hydroxylation, methylation, oxidative_bond_formation, epoxidation, double_bond_reduction
from pikachu.chem.structure import Structure

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
    @staticmethod
    def from_string(label: str) -> "TailoringEnzymeType":
        for value in TailoringEnzymeType:
            if str(value.name) == label:
                return value
        raise ValueError(f"Unknown tailoring enzyme: {label}")


class TailoringEnzyme:

    def __init__(self, gene_name, enzyme_type, atoms:list = None) -> None:
        self.gene_name = gene_name
        self.type = TailoringEnzymeType.from_string(enzyme_type)
        self.atoms = atoms

    def do_tailoring(self, structure):
        """
        Performs tailoring reaction
        """
        if self.type.name == "P450_HYDROXYLATION":
            for atom in self.atoms:
                atom = atom[0] #only one atom is hydroxylated at a time
                atom = structure.get_atom(atom)
                structure = hydroxylation(atom, structure)
        elif self.type.name in ["METHYLTRANSFERASE", "C_METHYLTRANSFERASE", "N_METHYLTRANSFERASE", "O_METHYLTRANSFERASE"]:
            for atom in self.atoms:
                atom = atom[0] #only one atom is methylated at a time
                atom = structure.get_atom(atom)
                structure = methylation(atom, structure)
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
        
        
        return structure
