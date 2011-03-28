from Exceptions import ProteinManipulationException
from ProteinEntity import ProteinEntity, DisorderedProteinEntity

class Residue(ProteinEntity):
    """
    Represents a residue. A Residue contains a list of
    atoms indexed numerically, i.e. residue.atoms[0]
    """
    def __init__(self, id, name, segid):
        ProteinEntity.__init__(self, 'RESIDUE', id)
        self.__disordered = False
        self.__name = name
        self.__segid = segid


    def __repr__(self):
        hetflag, resseq, icode=self.id()
        full_id=(self.__name, hetflag, resseq, icode)
        return "<Residue %s het=%s resseq=%s icode=%s>" % full_id


    def segid(self):
        return self.__segid

    def disorder_identifier(self):
        return self.__name

    def flag_disordered(self):
        "Set the disordered flag."
        self.__disordered = True

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

    def chain(self):
        return self.parent()

    def atoms(self, hash_key=None):
        """
        Returns the list atoms, treating DisorderedAtoms as one."
        """
        return self.children(hash_key)

    def all_atoms(self):
        """
        Returns the list of ALL atoms, unpacking the DisorderedAtoms as well."
        """
        undisordered_atom_list=[]
        for atom in self.atoms():
            if atom.is_disordered():
                undisordered_atom_list = (undisordered_atom_list + atom.disordered_atoms())
            else:
                undisordered_atom_list.append(atom)
        return undisordered_atom_list



class DisorderedResidue(DisorderedProteinEntity):
    def __init__(self, id):
        """
        Arguments:
        o id - string, atom name
        """
        self.__last_occupancy = -1
        DisorderedProteinEntity.__init__(self, 'DISORDERED RESIDUE', id)


    def __repr__(self):
        return "<Disordered RESIDUE %s>" % self.id()


    # Sub-entity manipulation methods
    
    def add_residue(self, residue):
        self.add_child(residue, residue.disorder_identifier())
        self.set_main_disordered_child(residue)

    def add_atom_to_residue(self, atom, residue_id=None):
        residue = self.__get_residue(residue_id)
        residue.add_atom(atom)

    def remove_atom_from_residue(self, atom_id, residue_id=None):
        residue = self.__get_residue(residue_id)
        residue.remove_atom_with_id(atom_id)

    def has_residue_with_name(self, name):
        return self.has_child_with_id(name)

    def remove_residue_with_name(self, name):
        self.remove_child_with_id(name)

    def reset_main_disorder_identifier(self):
        new_id = None
        if len(self.children()) is not 0:
            new_id  = self.children_keys()[-1]
        self.set_main_disorder_identifier(new_id)

    def __get_residue(self, residue_id):
        try:
            if residue_id is None:
                    residue_id = self.main_disorder_identifier()
            residue = self.children(residue_id)
            assert(residue is not None and (not isinstance(residue, list)))
            return residue
        except:
            raise ProteinManipulationException( \
                "No such residue with disorder identifier %s belonging to %s" \
                % (residue_id, self) )


    # Hierarchy Identity Methods

    def disordered_residues(self, hash_key=None):
        return self.children(hash_key)

    def chain(self):
        return self.parent()