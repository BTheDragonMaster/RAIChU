from enum import Enum, unique
from raichu.reactions.general_tailoring_reactions import hydroxylation, methylation
from pikachu.chem.structure import Structure

@unique
class TailoringEnzymeType(Enum):
    P450 = 1
    METHYLTRANSFERASE = 2
    C_METHYLTRANSFERASE = 3
    N_METHYLTRANSFERASE = 4
    O_METHYLTRANSFERASE = 5
    @staticmethod
    def from_string(label: str) -> "TailoringEnzymeType":
        for value in TailoringEnzymeType:
            if str(value.name) == label:
                return value
        raise ValueError(f"Unknown tailoring enzyme: {label}")


class TailoringEnzyme:

    def __init__(self, gene_name, enzyme_type, atoms = None) -> None:
        self.gene_name = gene_name
        self.type = TailoringEnzymeType.from_string(enzyme_type)
        self.atoms = atoms

    def do_tailoring(self, structure):
        """
        Performs tailoring reaction
        """
        if self.type.name == "P450":
            for atom in self.atoms:
                atom = structure.get_atom(atom)
                structure = hydroxylation(atom, structure)
        elif "METHYLTRANSFERASE" in self.type.name:
            for atom in self.atoms:
                atom = structure.get_atom(atom)
                structure = methylation(atom, structure)
        return structure
