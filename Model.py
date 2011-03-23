from ProteinEntity import ProteinEntity

class Model(ProteinEntity):
    """
    The object representing a model in a structure. In a structure
    derived from an X-ray crystallography experiment, only a single
    model will be present (with some exceptions). NMR structures
    normally contain many different models.

    Each Model contains one or more protein sub-chains that are
    accessed alphabetically, i.e. model.chains['A']
    """

    def __init__(self, biopython_model):
        ProteinEntity.__init__(self)
        self.level = 'MODEL'
        
    def __repr__(self):
        #resname=self.get_resname()
        #hetflag, resseq, icode=self.get_id()
        #full_id=(resname, hetflag, resseq, icode)
        #return "<Residue %s het=%s resseq=%s icode=%s>" % full_id

    # Identity Methods
    def chains(self, hash_key=None):
        return self.__children(hash_key)

    def residues(self):
        residues = []
        for c in self.chains():
            residues += c.residues()
        return residues

    def atoms(self):
        atoms = []
        for r in self.residues():
            atoms += r.atoms()
        return atoms

    def structure(self):
        return self.__parent
