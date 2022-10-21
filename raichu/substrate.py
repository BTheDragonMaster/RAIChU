from paras.features import get_smiles
from enum import Enum, unique

_PKS_TO_SMILES = {"WILDCARD": "SC(C([*])C(O)=O)=O",
                  "MALONYL_COA": "SC(CC(O)=O)=O",
                  "METHYLMALONYL_COA": "SC(C(C)C(O)=O)=O",
                  "METHOXYMALONYL_COA": "SC(C(C(O)=O)OC)=O",
                  "METHYLBUTYRYL_COA_2S": "SC(=O)[C@@H](C)CC",
                  "METHYLBUTYRYL_COA_2R": "SC(=O)[C@H](C)CC",
                  "ETHYLMALONYL_COA": "SC(C(C(O)=O)CC)=O"}


@unique
class PksElongationSubstrate(Enum):
    WILDCARD = 1
    MALONYL_COA = 2
    METHYLMALONYL_COA = 3
    METHOXYMALONYL_COA = 4
    METHYLBUTYRYL_COA_2S = 5
    METHYLBUTYRYL_COA_2R = 6
    ETHYLMALONYL_COA = 7

    @staticmethod
    def from_string(label: str):
        for value in PksElongationSubstrate:
            if str(value) == label:
                return value
        raise ValueError(f"Unknown PKS elongation substrate: {label}")

    @staticmethod
    def get_smiles(label: str):
        for value in PksElongationSubstrate:
            if str(value) == label:
                return _PKS_TO_SMILES[value]

        raise ValueError(f"Unknown PKS elongation substrate: {label}")


class Substrate:
    def __init__(self, name: str, smiles: str) -> None:
        self.name = name
        self.smiles = smiles


class NRPSSubstrate(Substrate):
    def __init__(self, name: str) -> None:
        smiles = get_smiles(name)
        super().__init__(name, smiles)


class PKSSubstrate(Substrate):
    def __init__(self, name: str) -> None:
        smiles = _PKS_TO_SMILES.get(name, None)
        if smiles is None:
            raise ValueError(f"Cannot fetch SMILES string for PKS substrate {name}.")
        super().__init__(name, smiles)
