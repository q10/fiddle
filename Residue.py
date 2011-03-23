from ProteinEntity import ProteinEntity

class Residue(ProteinEntity):
    """
    Represents a residue. A Residue contains a list of
    atoms indexed numerically, i.e. residue.atoms[0]
    """
    def __init__(self, ):
        ProteinEntity.__init__(self)
        self.level = 'RESIDUE'
        self.atoms

    def __repr__(self):
        #resname=self.get_resname()
        #hetflag, resseq, icode=self.get_id()
        #full_id=(resname, hetflag, resseq, icode)
        #return "<Residue %s het=%s resseq=%s icode=%s>" % full_id

    def set_childrens_parents(self):
        for atom in self.atoms:
            atom.set_parent(self)

    # Identity Methods
    def chain(self):
        return self.__parent

    def model(self):
        return self.chain().model()

    def structure(self):
        return self.model().structure()
