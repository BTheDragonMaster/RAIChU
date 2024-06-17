import os
from typing import Dict, Optional
import raichu.data
from dataclasses import dataclass

SMILES_FILE = os.path.join(os.path.dirname(raichu.data.__file__), 'smiles.txt')
METADATA_FILE = os.path.join(os.path.dirname(raichu.data.__file__), 'metadata.txt')


@dataclass
class ParasSubstrate:
    name: str = ''
    smiles: str = ''
    abbreviation: str = ''
    type: str = ''
    chirality: Optional[str] = None
    norine: Optional[str] = None

    @classmethod
    def from_line(cls, line):
        line_data = line.split('\t')
        for i in range(len(line_data)):
            line_data[i] = line_data[i].strip()
            if not line_data[i] or line_data[i] == 'None':
                line_data[i] = None
        name, smiles, abbreviation, substrate_type, chirality, norine = line_data
        assert name and smiles and abbreviation and substrate_type
        if not norine:
            norine = 'X'

        return cls(name, smiles, abbreviation, substrate_type, chirality, norine)


def load_smiles() -> Dict[str, str]:
    name_to_smiles = {}
    with open(SMILES_FILE, 'r') as smiles_file:
        smiles_file.readline()
        for line in smiles_file:
            name, smiles = line.strip().split('\t')
            name_to_smiles[name] = smiles
    return name_to_smiles


def load_metadata() -> Dict[str, ParasSubstrate]:
    name_to_metadata = {}
    with open(METADATA_FILE, 'r') as metadata:
        metadata.readline()
        for line in metadata:
            substrate = ParasSubstrate.from_line(line)
            name_to_metadata[substrate.name] = substrate
    return name_to_metadata


UNKNOWN_SUBSTRATE = "**Unknown**"

_SMILES = load_smiles()
_SMILES[UNKNOWN_SUBSTRATE] = "NC(*)C(=O)O"
_METADATA = load_metadata()
_METADATA[UNKNOWN_SUBSTRATE] = ParasSubstrate(UNKNOWN_SUBSTRATE, _SMILES[UNKNOWN_SUBSTRATE], "X", "amino_acid",
                                              "L", "X")


def get_smiles(substrate_name: str) -> str:
    """
    Return SMILES string from substrate name

    substrate_name: str. Unknown substrate: use UNKNOWN_SUBSTRATE
    """
    if not isinstance(substrate_name, str):
        raise TypeError(f"Expected str, got {type(substrate_name)}.")

    if substrate_name in _SMILES:
        return _SMILES[substrate_name]

    raise ValueError(f"Cannot fetch SMILES string for {substrate_name}.")


def get_metadata(substrate_name: str) -> ParasSubstrate:
    """
    Return PARAS metadata from substrate name

    substrate_name: str. Unknown substrate: use UNKNOWN_SUBSTRATE
    """
    if not isinstance(substrate_name, str):
        raise TypeError(f"Expected str, got {type(substrate_name)}.")

    if substrate_name in _METADATA:
        return _METADATA[substrate_name]

    raise ValueError(f"Cannot fetch metadata for {substrate_name}.")

