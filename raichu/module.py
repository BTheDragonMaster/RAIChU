class Module:
    def __init__(self, nr, cluster_type, tailoring_domains, smiles, starter=False):
        self.id = nr
        self.type = cluster_type
        self.tailoring_domains = tailoring_domains
        self.smiles = smiles
        self.starter = starter
