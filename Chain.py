from ProteinEntity import ProteinEntity

class Chain(ProteinEntity):
    """
    The object representing a chain of a poly-chain protein.
    """
    def __init__(self, biopython_chain):
        ProteinEntity.__init__(self)
        self.level = 'CHAIN'
        self.residues

    def __repr__(self):
        #resname=self.get_resname()
        #hetflag, resseq, icode=self.get_id()
        #full_id=(resname, hetflag, resseq, icode)
        #return "<Residue %s het=%s resseq=%s icode=%s>" % full_id

    # Identity Methods
    def residues(self, hash_key=None):
        return self.__children(hash_key)

    def atoms(self):
        atoms = []
        for r in self.residues():
            atoms += r.atoms()
        return atoms

    def model(self):
        return self.__parent

    def structure(self):
        return self.model().structure()
