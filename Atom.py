import numpy, warnings
from Exceptions import PDBConstructionWarning
from ProteinEntity import ProteinEntity, DisorderedProteinEntity

class Atom(ProteinEntity):

    """
    Atom object.

    The Atom object stores atom name (both with and without spaces),
    coordinates, B factor, occupancy, alternative location specifier
    and (optionally) anisotropic B factor and standard deviations of
    B factor and positions.

    @param name: atom name (eg. "CA"). Note that spaces are normally stripped.
    @type name: string

    @param coords: atomic coordinates [x,y,z]
    @type coord: Numeric array (Float0, size 3)

    @param bfactor: isotropic B factor
    @type bfactor: number

    @param occupancy: occupancy (0.0-1.0)
    @type occupancy: number

    @param altloc: alternative location specifier for disordered atoms
    @type altloc: string

    @param fullname: full atom name, including spaces, e.g. " CA ". Normally
    these spaces are stripped from the atom name.
    @type fullname: string

    @param element: atom element, e.g. "C" for Carbon, "HG" for mercury,
    @type fullname: uppercase string (or None if unknown)
    """

    def __init__(self, name, coords, bfactor, occupancy, altloc, fullname, serial_number, element=None):
        ProteinEntity.__init__(self, 'ATOM', name)
        self.__direct_parent = None

        self.__name = name      # eg. CA, spaces are removed from atom name
        self.__coords = coords
        self.__bfactor = bfactor
        self.__occupancy = occupancy
        self.__altloc = altloc
        self.__full_id = None   # (structure id, model id, chain id, residue id, atom id)
        self.__fullname = fullname
        self.__disordered_flag = 0
        self.__anisou_array = None
        self.__siguij_array = None
        self.__sigatm_array = None
        self.__serial_number = serial_number
        self.element = self.__derive_element(name, element)

        # Hash that keeps additional properties
        self.xtra={}


    def __derive_element(self, name, element):
        if not element:
            warnings.warn("Atom object (name=%s) without element" % name, PDBConstructionWarning)
            element = "?"
            print(name, "--> ?")
        elif len(element)>2 or element != element.upper() or element != element.strip():
            raise ValueError(element)
        return element


    def __repr__(self):
        return "<Atom %s>" % self.id()

    # Coordinates methods

    def coords(self):
        return self.__coords

    def x(self):
        return self.__coords[0]

    def y(self):
        return self.__coords[1]

    def z(self):
        return self.__coords[2]
    
    def set_coords(self, new_coords):
        self.__coords = new_coords

    def distance_from(self, other_atom):
        assert(isinstance(other_atom, Atom))
        diff = self.coords() - other_atom.coords()
        return numpy.sqrt(numpy.dot(diff, diff))


   # set methods for other attributes

    def set_serial_number(self, n):
        self.__serial_number=n

    def set_bfactor(self, bfactor):
        self.__bfactor=bfactor

    def set_altloc(self, altloc):
        self.__altloc=altloc

    def set_occupancy(self, occupancy):
        self.__occupancy=occupancy

    def set_sigatm(self, sigatm_array):
        """
        Set standard deviation of atomic parameters.

        The standard deviation of atomic parameters consists
        of 3 positional, 1 B factor and 1 occupancy standard
        deviation.

        @param sigatm_array: standard deviations of atomic parameters.
        @type sigatm_array: Numeric array (length 5)
        """
        self.__sigatm_array=sigatm_array

    def set_siguij(self, siguij_array):
        """
        Set standard deviations of anisotropic temperature factors.

        @param siguij_array: standard deviations of anisotropic temperature factors.
        @type siguij_array: Numeric array (length 6)
        """
        self.__siguij_array=siguij_array

    def set_anisou(self, anisou_array):
        """
        Set anisotropic B factor.

        @param anisou_array: anisotropic B factor.
        @type anisou_array: Numeric array (length 6)
        """
        self.__anisou_array=anisou_array

    def flag_disordered(self):
        """Set the disordered flag to 1.

        The disordered flag indicates whether the atom is disordered or not.
        """
        self.__disordered_flag = 1

    def direct_parent(self):
        return self.__direct_parent

    # Public get methods for attributes

    def set_direct_parent(self, new_parent):
        self.__direct_parent = new_parent

    def detach_direct_parent(self):
        self.set_direct_parent(None)

    def is_disordered(self):
        "Return the disordered flag (1 if disordered, 0 otherwise)."
        return self.__disordered_flag

    def sigatm(self):
        "Return standard deviation of atomic parameters."
        return self.__sigatm_array

    def siguij(self):
        "Return standard deviations of anisotropic temperature factors."
        return self.__siguij_array

    def anisou(self):
        "Return anisotropic B factor."
        return self.__anisou_array

    def serial_number(self):
        return self.__serial_number

    def name(self):
        return self.__name

    def fullname(self):
        return self.__fullname

    def full_id(self):
        """Return the full id of the atom.

        The full id of an atom is the tuple
        (structure id, model id, chain id, residue id, atom name, altloc).
        """
        return self.parent().full_id()+((self.name, self.altloc),)

    def bfactor(self):
        "Return B factor."
        return self.__bfactor

    def occupancy(self):
        "Return occupancy."
        return self.__occupancy

    def fullname(self):
        "Return the atom name, including leading and trailing spaces."
        return self.__fullname

    def altloc(self):
        "Return alternative location specifier."
        return self.__altloc

    def disorder_identifier(self):
        return self.__altloc


    # Hierarchy Identity Methods

    def residue(self):
        try:
            return self.parent()
        except:
            return None

    def chain(self):
        try:
            return self.parent().chain()
        except:
            return None

    def model(self):
        try:
            return self.parent().model()
        except:
            return None

    def structure(self):
        try:
            return self.parent().structure()
        except:
            return None


class DisorderedAtom(DisorderedProteinEntity):
    def __init__(self, id):
        """
        Arguments:
        o id - string, atom name
        """
        self.__last_occupancy=-1.0
        DisorderedProteinEntity.__init__(self, 'DISORDERED ATOM', id)

    def __repr__(self):
        return "<Disordered Atom %s>" % self.id()

    def add_atom(self, atom):
        self.add_child(atom, atom.disorder_identifier())
        occupancy = atom.occupancy()
        if occupancy > self.__last_occupancy:
            self.__last_occupancy = occupancy
            self.set_main_disordered_child(atom)

    def has_atom_with_altloc(self, altloc):
        return self.has_child_with_id(altloc)

    def remove_atom_with_altloc(self, altloc):
        self.remove_child_with_id(altloc)

    # also resets the __last_occupancy
    def reset_main_disorder_identifier(self):
        best_disorder_id = None
        if len(self.children()) is not 0:
            self.__last_occupancy = -1
            for key in self.children_keys():
                occupancy = self.children(key).occupancy()
                if occupancy > self.__last_occupancy:
                    best_disorder_id = key
                    self.__last_occupancy = occupancy
        self.set_main_disorder_identifier(best_disorder_id)


    # Hierarchy Identity Methods

    def main_atom(self):
        return self.children(self.main_disorder_identifier())

    def disordered_atoms(self, hash_key=None):
        return self.children(hash_key)

    def residue(self):
        return self.parent()

    def chain(self):
        try:
            return self.parent().chain()
        except:
            return None

    def model(self):
        try:
            return self.parent().model()
        except:
            return None

    def structure(self):
        try:
            return self.parent().structure()
        except:
            return None