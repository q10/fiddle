import os
from Structures.ProteinEntity import ProteinEntity
from Structures.PDBInfo import PDBInfo
from Parsers.PDBCoordsParser import PDBCoordsParser

class Structure(ProteinEntity):
    """
    The Structure class contains a collection of Model instances.
    """

    def __init__(self, filename, id=None):
        # Extract the protein id from file
        assert(isinstance(filename, str))
        if id is None or id is '':
            id = os.path.basename(filename).split('.')[0]
        ProteinEntity.__init__(self, 'STRUCTURE', id)

        # Use new PDBCoordsParser instance to parse out and build protein protein model
        all_models = PDBCoordsParser().get_models(filename)
        for model in all_models:
            self.add_model(model)

        # Remove water "residues" (HOH)
        for chain in self.chains():
            chain.remove_waters()

        # Use PDBInfo to parse out and store protein meta-information
        self.__info = PDBInfo(filename)


    def __repr__(self):
        return "<Structure id=%s>" % self.id()

    def show_info(self):
        print(self.__info)

    def info(self, hash_key=None):
        if hash_key is None:
            return self.__info.full_description()
        else:
            return self.__info[hash_key]

    def info_keys(self):
        return self.__info.keys()


    # Sub-entity manipulation methods

    def add_model(self, model):
        self.add_child(model, model.id())

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