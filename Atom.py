from ProteinEntity import ProteinEntity

class Atom(ProteinEntity):
    def __init__(self, name, biopython_atom):
        ProteinEntity.__init__(self)
        self.level = 'ATOM'
        
        self.serial_number = biopython_atom.get_serial_number()
        [self.x, self.y, self.z] = biopython_atom.get_coord()
        self.name = biopython_atom.get_name()

    def __repr__(self):
        #resname=self.get_resname()
        #hetflag, resseq, icode=self.get_id()
        #full_id=(resname, hetflag, resseq, icode)
        #return "<Atom %s het=%s resseq=%s icode=%s>" % full_id

    def coords(self):
        return [self.x, self.y, self.z]
    def set_coords(self, nx, ny, nz):
        self.x, self.y, self.z = nx, ny, nz

    def distance_from(self, other_atom):
        # hold on for numpy use?


    # Identity Methods
    def residue(self):
        return self.__parent

    def chain(self):
        return self.residue().chain()

    def model(self):
        return self.chain().model()

    def structure(self):
        return self.model().structure()
