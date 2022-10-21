from raichu.module import _Module
from typing import List
from raichu.reactions.chain_release import cyclic_release


class Cluster:
    def __init__(self, modules: List[_Module], compute_cyclic_products=True) -> None:
        self.modules = modules
        self.chain_intermediate = None

        self.compute_cyclic_products = compute_cyclic_products

        self.structure_intermediates = []
        self.linear_product = None
        self.cyclised_products = []
        self.module_mechanisms = []

    def compute_structures(self):
        for module in self.modules:
            structure = module.run_module(self.chain_intermediate)
            self.structure_intermediates.append(structure.deepcopy())
            if module.is_termination_module:
                self.linear_product = module.release_chain(structure)
                break
        else:
            raise ValueError("Cluster must contain at least one termination module.")

        if self.compute_cyclic_products:
            self.cyclise_all()

    def cyclise(self, atom):
        cyclic_release(self.linear_product, atom)

    def cyclise_all(self):
        pass

    def draw_cluster(self):
        pass


class Mechanism:
    def __init__(self):
        self.structures = []

    def add_structure(self, structure):
        self.structures.append(structure.deepcopy())

    def draw_mechanism(self):
        pass
