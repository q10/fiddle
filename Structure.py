from ProteinEntity import ProteinEntity

class Structure(ProteinEntity):
    """
    The Structure class contains a collection of Model instances.
    """

    def __init__(self, filename):
        # get PDBParser to parse out information here
        ProteinEntity.__init__(self)
        self.level = 'STRUCTURE'








## FINISHED METHODS GO BELOW

    def __repr__(self):
        return "<Structure id=%s>" % self.id()

    # Sub-entity manipulation methods
    def add_model(self, model):
        self.add_child(model)

    def has_model_with_id(self, id):
        return self.has_child_with_id(id)

    def remove_model_with_id(self, id):
        self.remove_child_with_id(id)


    # Hierarchy Identity Methods
    def models(self, hash_key=None):
        return self.children(hash_key)

    def chains(self):
        chains = []
        for m in self.models():
            chains += m.chains()
        return chains

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