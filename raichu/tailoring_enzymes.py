from enum import Enum, unique
from raichu.reactions.general_tailoring_reactions import hydroxylation, methylation
from pikachu.chem.structure import Structure

@unique
class TailoringEnzymeType(Enum):
    P450_HYDROXYLATION = 1
    METHYLTRANSFERASE = 2

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
        if self.type.name == "P450_HYDROXYLATION":
            for atom in self.atoms:
                atom = structure.get_atom(atom)
                structure = hydroxylation(atom, structure)
        elif self.type.name == "METHYLTRANSFERASE":
            for atom in self.atoms:
                atom = structure.get_atom(atom)
                structure = methylation(atom, structure)
        return structure
