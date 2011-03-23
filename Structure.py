from ProteinEntity import ProteinEntity

class Structure(ProteinEntity):
    """
    The Structure class contains a collection of Model instances.
    """

    def __init__(self, filename):
        ProteinEntity.__init__(self)
        # get PDBParser to parse out information here
        self.level = 'STRUCTURE'
        self.models

    def __repr__(self):
        #resname=self.get_resname()
        #hetflag, resseq, icode=self.get_id()
        #full_id=(resname, hetflag, resseq, icode)
        #return "<Structure %s het=%s resseq=%s icode=%s>" % full_id


    def chains(self):
        chains = []
        for m in self.models:
            chains += m.chains
        return chains

    def residues(self):
        residues = []
        for c in self.chains():
            residues += c.residues
        return residues

    def atoms(self):
        atoms = []
        for r in self.residues():
            atoms += r.atoms
        return atoms

    def set_childrens_parents(self):
        for m in self.models:
            m.set_parent(self)