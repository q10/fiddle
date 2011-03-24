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

    def __init__(self, id, serial_num=None):
        ProteinEntity.__init__(self, 'MODEL', id)

        if serial_num is None:
            self.__serial_num = id
        else:
            self.__serial_num = serial_num


    def __repr__(self):
        return "<Model id=%s>" % self.id()


    # Sub-entity manipulation methods
    def add_chain(self, chain):
        self.add_child(chain)

    def has_chain_with_id(self, id):
        return self.has_child_with_id(id)

    def remove_chain_with_id(self, id):
        self.remove_child_with_id(id)


    # Hierarchy Identity Methods
    def chains(self, hash_key=None):
        return self.children(hash_key)

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
        return self.parent()
