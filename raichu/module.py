from typing import List, Union

from pikachu.general import read_smiles
from pikachu.chem.structure import Structure

from raichu.domain.domain import Domain, TailoringDomain, RecognitionDomain, \
    SynthesisDomain, CarrierDomain, TerminationDomain
from raichu.central_chain_detection.label_central_chain import label_pk_central_chain, label_nrp_central_chain
from raichu.attach_to_domain import attach_to_domain_pk, attach_to_domain_nrp
from enum import Enum, unique


@unique
class ModuleType(Enum):
    NRPS = 1
    PKS = 2

    @staticmethod
    def from_string(label: str) -> "ModuleType":
        for value in ModuleType:
            if str(value.name) == label:
                return value
        raise ValueError(f"Unknown module type: {label}")


@unique
class PKSModuleSubtype(Enum):
    PKS_CIS = 1
    PKS_TRANS = 2
    PKS_ITER = 3

    @staticmethod
    def from_string(label: str) -> "PKSModuleSubtype":
        for value in PKSModuleSubtype:
            if str(value.name) == label:
                return value
        raise ValueError(f"Unknown PKS module subtype: {label}")


@unique
class NRPSDomainType(Enum):
    A = 1
    C = 2
    PCP = 3
    E = 4
    nMT = 5
    TE = 6
    TD = 7
    UNKNOWN = 8

    @staticmethod
    def from_string(label: str) -> "NRPSDomainType":
        for value in NRPSDomainType:
            if str(value.name) == label:
                return value
        raise ValueError(f"Unknown NRPS domain type: {label}")


@unique
class PKSDomainType(Enum):
    AT = 1
    KS = 2
    ACP = 3
    KR = 4
    DH = 5
    ER = 6
    TE = 7
    TD = 8
    UNKNOWN = 9

    @staticmethod
    def from_string(label: str) -> "PKSDomainType":
        for value in PKSDomainType:
            if str(value.name) == label:
                return value
        raise ValueError(f"Unknown PKS domain type: {label}")


class _Module:

    def __init__(self, nr: int, module_type: str, domains: List[Domain], *,
                 module_subtype: Union[str, None] = None,
                 starter: bool = False, terminator: bool = False) -> None:
        self.id = nr
        self.type = ModuleType.from_string(module_type)
        self.subtype = None

        if module_subtype is not None:
            if self.type.name == 'PKS':
                self.subtype = PKSModuleSubtype.from_string(module_subtype)
            else:
                raise ValueError(f"Module subtypes not supported for {self.type.name} modules.")

        self.domains = domains

        self.recognition_domain = None
        self.synthesis_domain = None
        self.carrier_domain = None
        self.termination_domain = None

        self.tailoring_domains = []

        self.is_starter_module = starter
        self.is_termination_module = terminator

        for domain in self.domains:
            if self.type.name == 'NRPS':
                assert domain.type.name in [v.name for v in NRPSDomainType]
            if self.type.name == 'PKS':
                assert domain.type.name in [v.name for v in PKSDomainType]

            if domain.used and domain.active:
                if isinstance(domain, TailoringDomain):
                    if domain.type.name not in [d.type.name for d in self.tailoring_domains]:
                        self.tailoring_domains.append(domain)
                    else:
                        raise ValueError(f"Cannot have two used tailoring domains of type {domain.type.name} in one \
module. Remove domain or set the 'used' or 'active' flag to False")

                if isinstance(domain, RecognitionDomain):
                    if not self.recognition_domain:
                        self.recognition_domain = domain
                    else:
                        raise ValueError(f"Cannot have more than one used recognition domain in one \
module. Remove a domain or set the 'used' or 'active' flag to False")

                if isinstance(domain, SynthesisDomain) and domain.is_elongating:
                    if not self.synthesis_domain:
                        self.synthesis_domain = domain
                    else:
                        raise ValueError("Cannot have more than one used and elongating synthesis domains \
in one module. Remove a domain, set the 'used' or 'active' flag to False, or set the 'is_elongating' flag to False.")

                if isinstance(domain, CarrierDomain):
                    if not self.carrier_domain:
                        self.carrier_domain = domain
                    else:
                        raise ValueError("Cannot have more than one used carrier domain in one \
module. Remove a domain or set the 'used' or 'active' flag to False")

                if isinstance(domain, TerminationDomain):
                    if not self.termination_domain:
                        self.termination_domain = domain
                    else:
                        raise ValueError("Cannot have more than one used termination domain in one \
module. Remove a domain or set the 'used' or 'active' flag to False")

        if not self.is_starter_module and not self.synthesis_domain:
            if self.type.name == 'NRPS':
                self.synthesis_domain = SynthesisDomain('C')
            elif self.type.name == 'PKS':
                self.synthesis_domain = SynthesisDomain('KS')

        if self.is_termination_module and not self.termination_domain:
            self.termination_domain = TerminationDomain("DUMMY_TE")

    def run_module(self, structure: Union[Structure, None] = None):
        raise NotImplementedError

    def get_tailoring_domain(self, domain_name: str) -> Union[TailoringDomain, None]:
        for domain in self.tailoring_domains:
            if domain.type.name == domain_name:
                return domain

        return

    def do_pks_tailoring(self, structure: Structure) -> Structure:
        kr_domain = self.get_tailoring_domain("KR")
        dh_domain = self.get_tailoring_domain("DH")
        er_domain = self.get_tailoring_domain("ER")

        if kr_domain and kr_domain.active:
            assert kr_domain.subtype is not None
            structure = kr_domain.do_tailoring(structure)
            if not kr_domain.subtype.name == 'C1' and not kr_domain.subtype.name == 'C2':
                if dh_domain and dh_domain.active:
                    structure = dh_domain.do_tailoring(structure)
                    if er_domain and er_domain.active:
                        structure = er_domain.do_tailoring(structure)

        return structure

    def do_nrps_tailoring(self, structure: Structure) -> Structure:
        e_domain = self.get_tailoring_domain('E')
        n_mt_domain = self.get_tailoring_domain('nMT')

        if e_domain and e_domain.active:
            structure = e_domain.do_tailoring(structure)
        if n_mt_domain and n_mt_domain.active:
            structure = n_mt_domain.do_tailoring(structure)

        return structure

    def release_chain(self, structure: Structure) -> Structure:
        assert structure
        assert self.termination_domain

        released_structure = self.termination_domain.release_chain(structure)

        return released_structure


class LinearPKSModule(_Module):

    def __init__(self, nr, domains, starter=False, terminator=False):
        super().__init__(nr, "PKS", domains, starter=starter, terminator=terminator, module_subtype="PKS_CIS")

    def run_module(self, structure: Union[Structure, None] = None) -> Structure:
        if structure is None:
            assert self.is_starter_module
            starter_unit = read_smiles(self.recognition_domain.substrate.smiles)
            label_pk_central_chain(starter_unit)
            structure = attach_to_domain_pk(starter_unit)
        else:
            structure = self.synthesis_domain.do_elongation(structure, self.recognition_domain.substrate)

        structure = self.do_pks_tailoring(structure)

        return structure


class IterativePKSModule(_Module):
    def __init__(self, nr, domains, starter=False, terminator=False, iterations: int = 1) -> None:
        super().__init__(nr, "PKS", domains, starter=starter, terminator=terminator, module_subtype="PKS_ITER")
        self.iterations = iterations

    def run_module(self, structure=None):
        self.do_pks_tailoring(structure)
        raise NotImplementedError


class TransATPKSModule(_Module):
    def __init__(self, nr, domains, substrate_name="MALONYL_COA", starter=False, terminator=False) -> None:
        super().__init__(nr, "PKS", domains, starter=starter, terminator=terminator, module_subtype="PKS_TRANS")
        if not self.recognition_domain:
            self.recognition_domain = RecognitionDomain("DUMMY_AT", substrate_name)

    def run_module(self, structure=None):
        self.do_pks_tailoring(structure)
        raise NotImplementedError


class NRPSModule(_Module):

    def __init__(self, nr, domains, starter=False, terminator=False):
        super().__init__(nr, "NRPS", domains, starter=starter, terminator=terminator)

    def run_module(self, structure=None) -> Structure:
        if structure is None:
            assert self.is_starter_module
            starter_unit = read_smiles(self.recognition_domain.substrate.smiles)
            label_nrp_central_chain(starter_unit)
            structure = attach_to_domain_nrp(starter_unit)
        else:
            structure = self.synthesis_domain.do_elongation(structure, self.recognition_domain.substrate)

        structure = self.do_nrps_tailoring(structure)

        return structure
