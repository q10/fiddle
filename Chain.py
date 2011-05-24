from ProteinEntity import ProteinEntity

class Chain(ProteinEntity):
    """
    The object representing a chain of a poly-chain protein.
    """
    def __init__(self, id):
        ProteinEntity.__init__(self, 'CHAIN', id)

    def __repr__(self):
        return "<Chain id=%s>" % self.id()


    # Sub-entity manipulation methods

    def add_residue(self, residue):
        self.add_child(residue, residue.id())

    def has_residue_with_id(self, id):
        return self.has_child_with_id(id)

    def remove_residue_with_id(self, id):
        self.remove_child_with_id(id)

    def remove_waters(self):
        water_residues = filter(lambda residue: residue.name() == 'HOH', self.residues())
        for water in water_residues:
            self.remove_child_with_id(water.id())


    # Hierarchy Identity Methods

    def residues(self, hash_key=None):
        if hash_key is None:
            return list(sorted(self.children(), key=lambda residue: residue.sequence_number()))
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
            return self.parent().structure()
        except:
            return None