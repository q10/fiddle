from ProteinEntity import ProteinEntity

class Chain(ProteinEntity):
    """
    The object representing a chain of a poly-chain protein.
    """
    def __init__(self, id):
        ProteinEntity.__init__(self, 'CHAIN', id)

    def __repr__(self):
        return "<Chain id=%s>" % self.id()


## FINISHED METHODS GO BELOW

    # Sub-entity manipulation methods
    def add_residue(self, residue):
        self.add_child(residue, residue.id())

    def has_residue_with_id(self, id):
        return self.has_child_with_id(id)

    def remove_residue_with_id(self, id):
        self.remove_child_with_id(id)


    # Hierarchy Identity Methods
    def residues(self, hash_key=None):
        return self.children(hash_key)

    def atoms(self):
        atoms = []
        for r in self.residues():
            atoms += r.atoms()
        return atoms

    def model(self):
        return self.parent()

    def structure(self):
        try:
            return self.model().structure()
        except:
            return None