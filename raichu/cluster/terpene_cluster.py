from typing import List

from pikachu.general import read_smiles

from raichu.substrate import TerpeneCyclaseSubstrate
from raichu.representations import (
    MacrocyclizationRepresentation,
    TailoringRepresentation,
    IsomerizationRepresentation,
    MethylShiftRepresentation,
)
from raichu.reactions.general_tailoring_reactions import (
    dephosphorylation,
    oxidative_bond_formation,
    double_bond_reduction,
)
from raichu.cluster.base_cluster import Cluster
from raichu.reactions.general_tailoring_reactions import (
    double_bond_shift,
    reductive_bond_breakage,
    addition
)


class TerpeneCluster(Cluster):
    def __init__(
        self,
        cyclase_name: str,
        precursor: str,
        cyclase_type: str = None,
        macrocyclisations: List[MacrocyclizationRepresentation] = None,
        double_bond_isomerisations: List[IsomerizationRepresentation] = None,
        methyl_shifts: List[MethylShiftRepresentation] = None,
        tailoring_representations: List[TailoringRepresentation] = None,
    ) -> None:
        super().__init__(tailoring_representations, macrocyclisations)
        self.cyclase_name = cyclase_name
        self.cyclase_type = cyclase_type
        self.precursor = precursor
        self.isomerization_representations = double_bond_isomerisations
        self.methyl_shift_representations = methyl_shifts

        self.chain_intermediate = None
        self.tailored_product = None
        self.final_product = None
        self.terpene = True

    def create_precursor(self) -> None:
        substrate = TerpeneCyclaseSubstrate(self.precursor)
        self.chain_intermediate = read_smiles(substrate.smiles)

    def do_double_bond_isomerization(self):
        if not self.isomerization_representations:
            return
        initialized_isomerization_atoms = self.initialize_modification_sites(
            [
                isomerization_representation.modification_sites
                for isomerization_representation in self.isomerization_representations
            ]
        )
        for atoms in initialized_isomerization_atoms:
            if len(atoms) < 4:
                continue
            old_double_bond_atom1 = self.chain_intermediate.get_atom(atoms[0])
            old_double_bond_atom2 = self.chain_intermediate.get_atom(atoms[1])
            new_double_bond_atom1 = self.chain_intermediate.get_atom(atoms[2])
            new_double_bond_atom2 = self.chain_intermediate.get_atom(atoms[3])
            if len(set(atoms)) == len(atoms):
                raise ValueError(
                    "The bonds need to be adjacent to perform a double bond shift."
                )
            self.chain_intermediate = double_bond_shift(
                self.chain_intermediate,
                old_double_bond_atom1,
                old_double_bond_atom2,
                new_double_bond_atom1,
                new_double_bond_atom2,
            )


    def do_methyl_shift(self):
        if not self.methyl_shift_representations:
            return
        initialized_methyl_shift_atoms = self.initialize_modification_sites(
                [
                    methyl_shift_representation.modification_sites
                    for methyl_shift_representation in self.methyl_shift_representations
                ]
            )
        for atoms in initialized_methyl_shift_atoms:
                if len(atoms) != 2:
                    continue
                transferred_c = atoms[0]
                # Assert its actually a methyl group
                assert [atom.type for atom in transferred_c.neighbours].count("H") == 3
                source = None
                source = [
                    atom for atom in transferred_c.neighbours if atom.type != "H"
                ][0]
                assert source
                destination_c = atoms[1]
                assert destination_c.has_neighbour("H")
                assert transferred_c.type == "C" and destination_c.type == "C"

                self.chain_intermediate = reductive_bond_breakage(source, transferred_c, self.chain_intermediate)
                self.chain_intermediate = addition(destination_c, "C", self.chain_intermediate)
                self.chain_intermediate.refresh_structure(find_cycles=True)



    def do_macrocyclization(self, sequential=True):
        initialized_macrocyclization_atoms = self.initialize_macrocyclization()
        self.cyclic_intermediates.append(self.chain_intermediate.deepcopy())
        self.do_double_bond_isomerization()
        self.do_methyl_shift()
        if self.cyclase_type != "Class_2":
            self.chain_intermediate = dephosphorylation(self.chain_intermediate)

        for (
            macrocyclization_atoms,
            cyclisation_type,
        ) in initialized_macrocyclization_atoms:
            atom1 = self.chain_intermediate.get_atom(macrocyclization_atoms[0])
            atom2 = self.chain_intermediate.get_atom(macrocyclization_atoms[1])
            found_bond = False
            for atom in [atom1, atom2]:
                if found_bond and any(
                    [neighbour.type == "H" for neighbour in atom.neighbours]
                ):
                    break
                for bond in atom.bonds:
                    if bond.type == "double":
                        self.chain_intermediate = double_bond_reduction(
                            *bond.neighbours, self.chain_intermediate
                        )
                        found_bond = True
                        break
            self.chain_intermediate = oxidative_bond_formation(
                atom1, atom2, self.chain_intermediate
            )
            self.cyclic_intermediates.append(self.chain_intermediate.deepcopy())

        self.cyclic_product = self.chain_intermediate
