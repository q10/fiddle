from Exceptions import ProteinManipulationException
from ProteinEntity import ProteinEntity, DisorderedProteinEntity

class Residue(ProteinEntity):
    """
    Represents a residue. A Residue contains a list of
    atoms indexed numerically, i.e. residue.atoms[0]
    """
    def __init__(self, id, name, segid):
        ProteinEntity.__init__(self, 'RESIDUE', id)
        self.__disordered = 0
        self.__name = name
        self.__segid = segid


    def __repr__(self):
        hetflag, resseq, icode=self.id()
        full_id=(self.__name, hetflag, resseq, icode)
        return "<Residue %s het=%s resseq=%s icode=%s>" % full_id



    def get_unpacked_list(self):
        """
        Returns the list of all atoms, unpack DisorderedAtoms."
        """
        atom_list=self.get_list()
        undisordered_atom_list=[]
        for atom in atom_list:
            if atom.is_disordered():
                undisordered_atom_list=(undisordered_atom_list+ atom.disordered_get_list())
            else:
                undisordered_atom_list.append(atom)
        return undisordered_atom_list

    def segid(self):
        return self.__segid

    def disorder_identifier(self):
        self.__name


## FINISHED METHODS GO BELOW

    def flag_disordered(self):
        "Set the disordered flag."
        self.__disordered=1

    def is_disordered(self):
        "Return 1 if the residue contains disordered atoms."
        return self.__disordered

    def name(self):
        return self.__name

    # Sub-entity manipulation methods
    def add_atom(self, atom):
        self.add_child(atom, atom.id())

    def has_atom_with_id(self, id):
        return self.has_child_with_id(id)

    def remove_atom_with_id(self, id):
        self.remove_child_with_id(id)


    # Hierarchy Identity Methods
    def atoms(self, hash_key=None):
        return self.children(hash_key)

    def chain(self):
        return self.parent()

    def model(self):
        try:
            return self.chain().model()
        except:
            return None

    def structure(self):
        try:
            return self.model().structure()
        except:
            return None



class DisorderedResidue(DisorderedProteinEntity):
    def __init__(self, id):
        """
        Arguments:
        o id - string, atom name
        """
        self.__last_occupancy=-1
        DisorderedProteinEntity.__init__(self, 'DISORDERED RESIDUE', id)


    def __repr__(self):
        return "<Disordered RESIDUE %s>" % self.id()

    def add_residue(self, residue):
        self.add_child(residue, residue.disorder_identifier())
        self.set_main_disordered_child(residue)

    def add_atom_to_residue(self, atom, residue_id=self.main_disorder_identifier()):
        try:
            assert(residue_id is not None)
            residue = self.children(residue_id)
        except:
            raise ProteinManipulationException( \
                "No such residue with disorder identifier %s " \
                % (residue_id) )
        else:
            residue.add_atom(atom)

    def remove_atom_from_residue(self, id, residue_id=self.main_disorder_identifier()):
        try:
            assert(residue_id is not None)
            residue = self.children(residue_id)
        except:
            raise ProteinManipulationException( \
                "No such residue with disorder identifier %s " \
                % (residue_id) )
        else:
            residue.remove_atom_with_id(id)

    def has_alternate_residue_with_name(self, name):
        return self.has_child_with_id(name)

    def remove_alternate_residue_with_name(self, name):
        self.remove_child_with_id(name)
        # need to check disordered_select on removing child and to check if disordered_children is empty


    # Hierarchy Identity Methods

    def alternate_residues(self, hash_key=None):
        return self.children(hash_key)

    def chain(self):
        return self.parent()

    def model(self):
        try:
            return self.chain().model()
        except:
            return None

    def structure(self):
        try:
            return self.model().structure()
        except:
            return None
    