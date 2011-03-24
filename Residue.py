from ProteinEntity import ProteinEntity

class Residue(ProteinEntity):
    """
    Represents a residue. A Residue contains a list of
    atoms indexed numerically, i.e. residue.atoms[0]
    """
    def __init__(self, id, resname, segid):
        ProteinEntity.__init__(self, 'RESIDUE', id)
        self.__disordered=0
        self.__resname=resname
        self.__segid=segid


    def __repr__(self):
        #resname=self.get_resname()
        #hetflag, resseq, icode=self.get_id()
        #full_id=(resname, hetflag, resseq, icode)
        #return "<Residue %s het=%s resseq=%s icode=%s>" % full_id


    '''
    def flag_disordered(self):
        "Set the disordered flag."
        self.disordered=1

    def is_disordered(self):
        "Return 1 if the residue contains disordered atoms."
        return self.disordered
    '''
    def get_resname(self):
        return self.resname

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



## FINISHED METHODS GO BELOW

    # Sub-entity manipulation methods
    def add_atom(self, atom):
        self.add_child(atom)

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