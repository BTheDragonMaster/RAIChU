from pikachu.chem.structure import *


class Domain(Atom):
    def __init__(self, atom_type, atom_nr, chiral, charge, aromatic):
        super().__init__(atom_type, atom_nr, chiral, charge, aromatic)
        self.type = atom_type
        self.domain_type = None
        self.domain_nr = atom_nr


    def __repr__(self):
        return f'{self.domain_type}_{self.domain_nr}'

    def __hash__(self):
        return self.nr

    def __eq__(self, domain):
        if type(domain) == Domain:
            return self.nr == domain.nr
        else:
            return False

    def print_domaintype(self):
        print(self.domain_type)

    def set_domain_type(self, domain_type):
        assert domain_type == 'ACP' or domain_type == 'AT' or\
               domain_type == 'KS' or domain_type == 'KR' or\
               domain_type == 'DH' or domain_type == 'ER' or\
               domain_type == 'TE' or domain_type == 'PCP'
        self.domain_type = domain_type


def make_domain(domain_type, next_atom_nr):
    """
    Returns a Domain of the input type containing the appropriate nr in the
    Structure object

    domain_type: Str, domain type
    next_atom_nr: Int, 'atom nr' of the domain in the structure
    """
    domain = Domain('I', next_atom_nr, None, 0, False)
    domain.set_domain_type(domain_type)
    return domain




