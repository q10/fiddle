# Copyright (C) 2002, Thomas Hamelryck (thamelry@binf.ku.dk)
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Consumer class that builds a Structure object.

This is used by the PDBParser and MMCIFparser classes.
"""

import warnings
from Model import Model
from Chain import Chain
from Residue import Residue, DisorderedResidue
from Atom import Atom, DisorderedAtom
from Exceptions import PDBConstructionException, PDBConstructionWarning


class StructureBuilder:
    """
    Deals with contructing the Structure object. The StructureBuilder class is used
    by the PDBParser classes to translate a file to a Structure object.
    """
    def __init__(self):
        self.line_counter = 0
        self.all_models_generated = []

    # METHOD TO BE CHANGED TO BOOLEAN RETURNS
    def _is_completely_disordered(self, residue):
        "Return 1 if all atoms in the residue have a non blank altloc."
        for atom in residue.all_atoms():
            if atom.altloc() == " ":
                return False
        return True


    # Public methods called by the Parser classes

    def set_line_counter(self, line_counter):
        """
        The line counter keeps track of the line in the PDB file that
        is being parsed.

        Arguments:
        o line_counter - int
        """
        assert(isinstance(line_counter, int))
        self.line_counter=line_counter

    def init_model(self, model_id, serial_num=None):
        """Initiate a new Model object with given id.

        Arguments:
        o id - int
        o serial_num - int
        """
        assert(isinstance(model_id, int) and (isinstance(serial_num, int) or serial_num is None))
        self.all_models_generated.append(Model(model_id,serial_num))
        self.model = self.all_models_generated[-1]

    def init_chain(self, chain_id):
        """Initiate a new Chain object with given id.

        Arguments:
        o chain_id - string
        """
        if self.model.has_chain_with_id(chain_id):
            self.chain = self.model.chains(chain_id)
            warnings.warn("WARNING: Chain %s is discontinuous at line %i."
                          % (chain_id, self.line_counter),
                          PDBConstructionWarning)
        else:
            self.chain = Chain(chain_id)
            self.model.add_chain(self.chain)

    def init_seg(self, segid):
        """Flag a change in segid.

        Arguments:
        o segid - string
        """
        assert(isinstance(segid, str))
        self.segid=segid

    def init_residue(self, resname, field, resseq, icode):
        """
        Initiate a new Residue object.

        Arguments:
        o resname - string, e.g. "ASN"
        o field - hetero flag, "W" for waters, "H" for
            hetero residues, otherwise blank.
        o resseq - int, sequence identifier
        o icode - string, insertion code
        """
        if field == "H":
            # The hetero field consists of H_ + the residue name (e.g. H_FUC)
            field = "H_" + resname
        res_id=(field, resseq, icode)
        if field==" ":
            if self.chain.has_residue_with_id(res_id):
                # There already is a residue with the id (field, resseq, icode).
                # This only makes sense in the case of a point mutation.
                warnings.warn("WARNING: Residue ('%s', %i, '%s') "
                              "redefined at line %i."
                              % (field, resseq, icode, self.line_counter),
                              PDBConstructionWarning)
                duplicate_residue=self.chain.residues(res_id)
                if duplicate_residue.is_disordered() == 2:
                    # The residue in the chain is a DisorderedResidue object.
                    # So just add the last Residue object.
                    if duplicate_residue.has_residue_with_name(resname):
                        # The residue was already made
                        self.residue = duplicate_residue
                        duplicate_residue.set_main_disorder_identifier(resname)
                    else:
                        # Make a new residue and add it to the already
                        # present DisorderedResidue
                        new_residue = Residue(res_id, resname, self.segid)
                        duplicate_residue.add_residue(new_residue)
                        self.residue = duplicate_residue
                        return
                else:
                    # Make a new DisorderedResidue object and put all
                    # the Residue objects with the id (field, resseq, icode) in it.
                    # These residues each should have non-blank altlocs for all their atoms.
                    # If not, the PDB file probably contains an error.
                    if not self._is_completely_disordered(duplicate_residue):
                        # if this exception is ignored, a residue will be missing
                        self.residue = None
                        raise PDBConstructionException(\
                            "Blank altlocs in duplicate residue %s ('%s', %i, '%s')" \
                            % (resname, field, resseq, icode))
                    self.chain.remove_residue_with_id(res_id)
                    new_residue = Residue(res_id, resname, self.segid)
                    disordered_residue = DisorderedResidue(res_id)
                    self.chain.add_residue(disordered_residue)
                    disordered_residue.add_residue(duplicate_residue)
                    disordered_residue.add_residue(new_residue)
                    self.residue = disordered_residue
                    return
        residue = Residue(res_id, resname, self.segid)
        self.chain.add_residue(residue)
        self.residue = residue

    def init_atom(self, name, coord, b_factor, occupancy, altloc, fullname,
                  serial_number=None, element=None):
        """
        Initiate a new Atom object.

        Arguments:
        o name - string, atom name, e.g. CA, spaces should be stripped
        o coord - Numeric array (Float0, size 3), atomic coordinates
        o b_factor - float, B factor
        o occupancy - float
        o altloc - string, alternative location specifier
        o fullname - string, atom name including spaces, e.g. " CA "
        o element - string, upper case, e.g. "HG" for mercury
        """
        residue = self.residue
        # if residue is None, an exception was generated during
        # the construction of the residue
        if residue is None:
            return
        # First check if this atom is already present in the residue.
        # If it is, it might be due to the fact that the two atoms have atom
        # names that differ only in spaces (e.g. "CA.." and ".CA.",
        # where the dots are spaces). If that is so, use all spaces
        # in the atom name of the current atom.
        if residue.has_atom_with_id(name):
                duplicate_atom = residue.atoms(name)
                # atom name with spaces of duplicate atom
                duplicate_fullname = duplicate_atom.fullname()
                if duplicate_fullname != fullname:
                    # name of current atom now includes spaces
                    name = fullname
                    warnings.warn("WARNING: atom names %s and %s differ "
                                  "only in spaces at line %i."
                                  % (duplicate_fullname, fullname,
                                     self.line_counter),
                                  PDBConstructionWarning)
        atom = self.atom = Atom(name, coord, b_factor, occupancy, altloc, fullname, serial_number, element)
        if altloc!=" ":
            # The atom is disordered
            if residue.has_atom_with_id(name):
                # Residue already contains this atom
                duplicate_atom = residue.atoms(name)
                if duplicate_atom.is_disordered() is 2:
                    duplicate_atom.add_atom(atom)
                else:
                    # This is an error in the PDB file:
                    # a disordered atom is found with a blank altloc
                    # Detach the duplicate atom, and put it in a
                    # DisorderedAtom object together with the current
                    # atom.
                    residue.remove_atom_with_id(name)
                    disordered_atom = DisorderedAtom(name)
                    residue.add_atom(disordered_atom)
                    disordered_atom.add_atom(atom)
                    disordered_atom.add_atom(duplicate_atom)
                    residue.flag_disordered()
                    warnings.warn("WARNING: disordered atom found "
                                  "with blank altloc before line %i.\n"
                                  % self.line_counter,
                                  PDBConstructionWarning)
            else:
                # The residue does not contain this disordered atom
                # so we create a new one.
                disordered_atom = DisorderedAtom(name)
                residue.add_atom(disordered_atom)
                # Add the real atom to the disordered atom, and the
                # disordered atom to the residue
                disordered_atom.add_atom(atom)
                residue.flag_disordered()
        else:
            # The atom is not disordered
            residue.add_atom(atom)

    def set_anisou(self, anisou_array):
        "Set anisotropic B factor of current Atom."
        self.atom.set_anisou(anisou_array)

    def set_siguij(self, siguij_array):
        "Set standard deviation of anisotropic B factor of current Atom."
        self.atom.set_siguij(siguij_array)

    def set_sigatm(self, sigatm_array):
        "Set standard deviation of atom position of current Atom."
        self.atom.set_sigatm(sigatm_array)

    def set_symmetry(self, spacegroup, cell):
        pass