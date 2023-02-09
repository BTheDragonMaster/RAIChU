from enum import Enum, unique
from raichu.reactions.general_tailoring_reactions import remove_atom, remove_group_at_bond, single_bond_oxidation, dephosphorylation, hydroxylation, addition, oxidative_bond_formation, epoxidation, double_bond_reduction, double_bond_shift
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
    AMINOTRANSFERASE = 13
    OXIDASE_DOUBLE_BOND_FORMATION = 14
    REDUCTASE_KETO_REDUCTION = 15
    ALCOHOLE_DEHYDROGENASE = 16
    DEHYDRATASE = 17
    DECARBOXYLASE = 18
    MONOAMINE_OXIDASE = 19
    
    
    @staticmethod
    def from_string(label: str) -> "TailoringEnzymeType":
        for value in TailoringEnzymeType:
            if str(value.name) == label:
                return value
        raise ValueError(f"Unknown tailoring enzyme: {label}")


class TailoringEnzyme:

    def __init__(self, gene_name, enzyme_type, modification_sites:list = None, substrate:str = None) -> None:
        self.gene_name = gene_name
        self.type = TailoringEnzymeType.from_string(enzyme_type)
        self.modification_sites = modification_sites
        self.substrate = substrate

    def do_tailoring(self, structure):
        """
        Performs tailoring reaction
        """
        if len(self.modification_sites)==0:
            return structure
        if self.type.name == "P450_HYDROXYLATION":
            for atoms in self.modification_sites:
                if len(atoms) == 0:
                    continue
                atom = atoms[0] #only one atom is hydroxylated at a time
                atom = structure.get_atom(atom)
                structure = addition(atom, "O", structure)
        elif self.type.name in ["METHYLTRANSFERASE", "C_METHYLTRANSFERASE", "N_METHYLTRANSFERASE", "O_METHYLTRANSFERASE"]:
            for atoms in self.modification_sites:
                if len(atoms) == 0:
                    continue
                atom = atoms[0] #only one atom is methylated at a time
                atom = structure.get_atom(atom)
                structure = addition(atom, "C", structure)
        elif self.type.name == "PRENYLTRANSFERASE":
            for atoms in self.modification_sites:
                if len(atoms) == 0:
                    continue
                atom = atoms[0]  # only one atom is methylated at a time
                atom = structure.get_atom(atom)
                if self.substrate not in PRENYL_TRANSFERASE_SUBSTRATES_TO_SMILES:
                    raise ValueError(
                        f"Not implemented prenyltransferase substrate: {self.substrate}")
                substrate = PRENYL_TRANSFERASE_SUBSTRATES_TO_SMILES[self.substrate]
                structure = addition(
                    atom, substrate, structure)
        elif self.type.name == "ACETYLTRANSFERASE":
            for atoms in self.modification_sites:
                if len(atoms) == 0:
                    continue
                atom = atoms[0]  # only one atom is methylated at a time
                atom = structure.get_atom(atom)
                structure = addition(atom, "[H]C(C)=O", structure)
        elif self.type.name == "ACYLTRANSFERASE":
            for atoms in self.modification_sites:
                if len(atoms) == 0:
                    continue
                atom = atoms[0]  # only one atom is methylated at a time
                atom = structure.get_atom(atom)
                if self.substrate:
                    structure = addition(atom, self.substrate, structure)
        elif self.type.name == "P450_OXIDATIVE_BOND_FORMATION":
            for atoms in self.modification_sites:
                if len(atoms) < 2:
                    continue
                atom1 = structure.get_atom(atoms[0])
                atom2 = structure.get_atom(atoms[1])
                structure = oxidative_bond_formation(atom1, atom2, structure)
        elif self.type.name == "P450_EPOXIDATION":
            for atoms in self.modification_sites:
                if len(atoms) < 2:
                    continue
                atom1 = structure.get_atom(atoms[0])
                atom2 = structure.get_atom(atoms[1])
                structure = epoxidation(atom1, atom2, structure)
        elif self.type.name == "REDUCTASE_DOUBLE_BOND_REDUCTION":
            for atoms in self.modification_sites:
                if len(atoms) < 2:
                    continue
                atom1 = structure.get_atom(atoms[0])
                atom2 = structure.get_atom(atoms[1])
                structure = double_bond_reduction(atom1, atom2, structure)
        elif self.type.name == "OXIDASE_DOUBLE_BOND_FORMATION":
            for atoms in self.modification_sites:
                if len(atoms) < 2:
                    continue
                atom1 = structure.get_atom(atoms[0])
                atom2 = structure.get_atom(atoms[1])
                structure = single_bond_oxidation(atom1, atom2, structure)
        elif self.type.name == "ISOMERASE_DOUBLE_BOND_SHIFT":
            for atoms in self.modification_sites:
                if len(atoms) < 4:
                    continue
                old_double_bond_atom1 = structure.get_atom(atoms[0])
                old_double_bond_atom2 = structure.get_atom(atoms[1])
                new_double_bond_atom1 = structure.get_atom(atoms[2])
                new_double_bond_atom2 = structure.get_atom(atoms[3])
                structure = double_bond_shift(
                    structure, old_double_bond_atom1, old_double_bond_atom2, new_double_bond_atom1, new_double_bond_atom2)
        elif self.type.name == "AMINOTRANSFERASE":
            for atom in self.modification_sites:
                if len(atoms) == 0:
                    continue
                atom1 = atom[0] #only one atom is modified at a time
                atom1 = structure.get_atom(atom1)
                oxygen = atom1.get_neighbour('O')
                structure = double_bond_reduction(atom1, oxygen, structure)
                oxygen = structure.get_atom(oxygen)
                structure = remove_atom(oxygen, structure)
                atom = structure.get_atom(atom)
                structure = addition(atom, "N", structure)
        elif self.type.name == "REDUCTASE_KETO_REDUCTION":
            for atom in self.modification_sites:
                if len(atoms) == 0:
                    continue
                atom1 = atom[0] #only one atom is modified at a time
                atom1 = structure.get_atom(atom1)
                if atom1.type != "O":
                    print(f"Can not perform KETO_REDUCTION on atom {atom1}, since there is no oxygen to be reduced.")
                    continue
                atom2 = atom1.get_neighbour('C')
                structure = double_bond_reduction(atom1, atom2, structure)
        elif self.type.name == "ALCOHOLE_DEHYDROGENASE":
            for atom in self.modification_sites:
                if len(atoms) == 0:
                    continue
                atom1 = atom[0] #only one atom is modified at a time
                atom1 = structure.get_atom(atom1)
                if atom1.type != "O":
                    print(f"Can not perform ALCOHOLE_DEHYDROGENASE on atom {atom1}, since there is no oxygen to be reduced.")
                    continue
                atom2 = atom1.get_neighbour('C')
                structure = single_bond_oxidation(atom1, atom2, structure)
        elif self.type.name == "DECARBOXYLASE":
            for atom in self.modification_sites:
                if len(atoms) == 0:
                    continue
                atom1 = atom[0] #only one atom is modified at a time
                atom1 = structure.get_atom(atom1)
                structure = remove_atom(atom1, structure)     
        elif self.type.name == "DEHYDRATASE":
            for atoms in self.modification_sites:
                if len(atoms) < 2:
                    continue
                atom1 = structure.get_atom(atoms[0])
                atom2 = structure.get_atom(atoms[1])
                oxygen = atom1.get_neighbour('O')
                if not oxygen:
                    oxygen = atom2.get_neighbour('O')
                if not oxygen:
                    print(f"Can not perform DEHYDRATASE on atoms {atom1} and {atom2}, since there is no hydroxygroup to be removed.")
                    continue
                structure = remove_atom(oxygen, structure)
                structure = single_bond_oxidation(atom1, atom2, structure)
        elif self.type.name == "MONOAMINE_OXYDASE":
            for atom in self.modification_sites:
                if len(atoms) == 0:
                    continue
                atom1 = atom[0] #only one atom is modified at a time
                atom1 = structure.get_atom(atom1)
                if atom1.type != "N":
                    print(f"Can not perform MONOAMINE_OXYDASE on atom {atom1}, since there is no nitrogen to be removed.")
                    continue
                structure = remove_atom(atom1, structure) 
        
        return structure
