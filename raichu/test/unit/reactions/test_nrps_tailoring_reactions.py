from pathlib import Path
import unittest
from raichu.reactions.nrps_tailoring_reactions import epimerize, n_methylate

TEST_FOLDER = Path(__file__).resolve().parent

NRPS_DATA = TEST_FOLDER / "test_data/nrps/marformycin/module 2.structure"


class TestNrpsTailoring(unittest.TestCase):
    def test_epimerize_glycine(self):
        pass

    def test_epimerize_tryptophan(self):
        pass
