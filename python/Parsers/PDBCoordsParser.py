# Copyright (C) 2002, Thomas Hamelryck (thamelry@binf.ku.dk)
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Parser for PDB files."""

import warnings
import numpy
from Utils.Exceptions import PDBConstructionException, PDBConstructionWarning
from Parsers.StructureBuilder import StructureBuilder

# If PDB spec says "COLUMNS 18-20" this means line[17:20]


class PDBCoordsParser:
    """
    Parse a PDB file and return a Structure object.
    """

    def __init__(self, PERMISSIVE=1):
        """
        The PDB parser call a number of standard methods in an aggregated
        StructureBuilder object. Normally this object is instantiated by the
        PDBParser object itself, but if the user provides his own StructureBuilder
        object, the latter is used instead.

        Arguments:
        o PERMISSIVE - int, if this is 0 exceptions in constructing the
        SMCRA data structure are fatal. If 1 (DEFAULT), the exceptions are
        caught, but some residues or atoms will be missing. THESE EXCEPTIONS
        ARE DUE TO PROBLEMS IN THE PDB FILE!.
        o structure_builder - an optional user implemented StructureBuilder class.
        """
        self.structure_builder=StructureBuilder()
        self.trailer=None
        self.line_counter=0
        self.PERMISSIVE=PERMISSIVE

    # Public methods

    def get_models(self, filename):
        """
        Return the models.

        Arguments:
        o file - name of the PDB file
        """
        self.trailer=None

        # Make a StructureBuilder instance (pass id of structure as parameter)
        file = open(filename)
        self._parse(file.readlines())
        return self.structure_builder.all_models_generated

    def get_trailer(self):
        "Return the trailer."
        return self.trailer

    # Private methods

    def _parse(self, header_coords_trailer):
        "Parse the PDB file."
        # Filters the header, parses the atomic data in the coords, then return the rest of the file (PDB file trailer)
        self.trailer = self._parse_coordinates(self._skip_header(header_coords_trailer))

    def _skip_header(self, header_coords_trailer):
        "Skip the header of the PDB file and return the rest."
        for i in range(len(header_coords_trailer)):
            self.structure_builder.set_line_counter(i+1)
            record_type = header_coords_trailer[i][0:6]
            if record_type == 'ATOM  ' or record_type == 'HETATM' or record_type == 'MODEL ':
                break
        # Return the rest of the coords+trailer for further processing
        self.line_counter = i
        return header_coords_trailer[i:]

    def _parse_coordinates(self, coords_trailer):
        "Parse the atomic data in the PDB file."
        structure_builder = self.structure_builder

        # Flag we have an open model
        local_line_counter = current_model_id = model_open = 0
        current_chain_id = current_segid = current_residue_id = current_resname = None

        for i in range(len(coords_trailer)):
            line = coords_trailer[i]
            record_type = line[0:6]
            global_line_counter = self.line_counter + local_line_counter + 1
            structure_builder.set_line_counter(global_line_counter)

            if record_type=='ATOM  ' or record_type=='HETATM':
                # Initialize the Model - there was no explicit MODEL record
                if not model_open:
                    structure_builder.init_model(current_model_id)
                    current_model_id+=1
                    model_open=1

                # Extract all data for one atom
                serial_number = self._get_serial_number(line)
                fullname, name = self._get_name(line)
                altloc, resname = self._get_altloc_and_resname(line)
                chainid = self._get_chainid(line)
                resseq = self._get_residue_sequence_number(line)
                icode = self._get_insertion_code(line)
                coord = self._get_coordinates(line, global_line_counter)
                occupancy, bfactor = self._get_occupancy_and_bfactor(line, global_line_counter)
                segid = self._get_segid(line)
                element = self._get_element(line)
                hetero_flag = self._get_hetero_flag(record_type, resname)
                residue_id=(hetero_flag, resseq, icode)

                if current_segid != segid:
                    current_segid = segid
                    structure_builder.init_seg(current_segid)

                # Build new chain AND residue first
                if current_chain_id != chainid:
                    current_chain_id = chainid
                    structure_builder.init_chain(current_chain_id)
                    current_residue_id=residue_id
                    current_resname = resname
                    try:
                        structure_builder.init_residue(resname, hetero_flag, resseq, icode)
                    except PDBConstructionException as message:
                        self._handle_PDB_exception(message, global_line_counter)

                # Build new residue first (still on same chain)
                elif current_residue_id != residue_id or current_resname != resname:
                    current_residue_id=residue_id
                    current_resname=resname
                    try:
                        structure_builder.init_residue(resname, hetero_flag, resseq, icode)
                    except PDBConstructionException as message:
                        self._handle_PDB_exception(message, global_line_counter)

                # Build atom
                try:
                    structure_builder.init_atom(name, coord, bfactor, occupancy, altloc,
                                                fullname, serial_number, element)
                except PDBConstructionException as message:
                    self._handle_PDB_exception(message, global_line_counter)

            elif record_type == 'MODEL ':
                structure_builder.init_model(current_model_id, self._get_model_serial_num(line))
                current_model_id += 1
                model_open = 1
                current_chain_id = current_residue_id = None

            elif record_type == 'ENDMDL':
                model_open = 0
                current_chain_id = current_residue_id = None

            elif record_type == 'ANISOU':
                structure_builder.set_anisou(self._get_anisou_array(line))
            elif record_type == 'SIGUIJ' :
                structure_builder.set_siguij(self._get_siguij_array(line))
            elif record_type == 'SIGATM':
                structure_builder.set_sigatm(self._get_sigatm_array(line))

            elif record_type == 'END   ' or record_type == 'CONECT':
                # End of atomic data, return the trailer
                self.line_counter += local_line_counter
                return coords_trailer[local_line_counter:]

            local_line_counter += 1

        # EOF (does not end in END or CONECT)
        self.line_counter += local_line_counter
        return []


    def _handle_PDB_exception(self, message, line_counter):
        """
        This method catches an exception that occurs in the StructureBuilder
        object (if PERMISSIVE==1), or raises it again, this time adding the
        PDB line number to the error message.
        """
        message = "%s at line %i." % (message, line_counter)
        if self.PERMISSIVE:
            # just print a warning - some residues/atoms may be missing
            warnings.warn("PDBConstructionException: %s\n"
                          "Exception ignored.\n"
                          "Some atoms or residues may be missing in the data structure."
                          % message, PDBConstructionWarning)
        else:
            # exceptions are fatal - raise again with new message (including line nr)
            raise PDBConstructionException(message)

    def _get_serial_number(self, line):
        try:
            return int(line[6:11])
        except:
            return 0

    def _get_model_serial_num(self, line):
        try:
            return int(line[10:14])
        except:
            self._handle_PDB_exception("Invalid or missing model serial number",
                                       global_line_counter)
            return 0

    def _get_name(self, line):
        fullname = line[12:16]
        # get rid of whitespace in atom names
        split_list = fullname.split()
        if len(split_list)!=1:
            # atom name has internal spaces, e.g. " N B ", so
            # we do not strip spaces
            name = fullname
        else:
            # atom name is like " CA ", so we can strip spaces
            name = split_list[0]
        return fullname, name

    def _get_altloc_and_resname(self, line):
        return line[16:17], line[17:20]

    def _get_chainid(self, line):
        return line[21:22]

    def _get_residue_sequence_number(self, line):
        return int(line[22:26].split()[0])   # sequence identifier
        
    def _get_insertion_code(self, line):
        return line[26:27]           # insertion code

    def _get_coordinates(self, line, global_line_counter):
        # atomic coordinates
        try:
            x, y, z = list(map(float, (line[30:38], line[38:46], line[46:54])))
        except:
            #Should we allow parsing to continue in permissive mode?
            #If so what coordindates should we default to?  Easier to abort!
            raise PDBConstructionException(\
                "Invalid or missing coordinate(s) at line %i." \
                % global_line_counter)
        return numpy.array((x, y, z), 'f')

    def _get_occupancy_and_bfactor(self, line, global_line_counter):
        # occupancy & B factor
        try:
            occupancy = float(line[54:60])
        except:
            self._handle_PDB_exception("Invalid or missing occupancy",
                                       global_line_counter)
            occupancy = 0.0 #Is one or zero a good default?
        try:
            bfactor = float(line[60:66])
        except:
            self._handle_PDB_exception("Invalid or missing B factor",
                                       global_line_counter)
            bfactor = 0.0 #The PDB use a default of zero if the data is missing
        return occupancy, bfactor

    def _get_segid(self, line):
        return line[72:76]

    def _get_element(self, line):
        return line[76:78].strip()
    
    def _get_hetero_flag(self, record_type, resname):
        if record_type == 'HETATM':       # hetero atom flag
            if resname == "HOH" or resname == "WAT":
                return "W"
            else:
                return "H"
        else:
            return " "

    def _get_anisou_array(self, line):
        anisou = list(map(float, (line[28:35], line[35:42], line[43:49], line[49:56], line[56:63], line[63:70])))
        # U's are scaled by 10^4 on the PDB file, so this scales back
        return (numpy.array(anisou, 'f')/10000.0).astype('f')

    def _get_siguij_array(self, line):
        # standard deviation of anisotropic B factor
        siguij = list(map(float, (line[28:35], line[35:42], line[42:49], line[49:56], line[56:63], line[63:70])))
        # U sigma's are scaled by 10^4
        return (numpy.array(siguij, 'f')/10000.0).astype('f')

    def _get_sigatm_array(self, line):
        # standard deviation of atomic positions
        sigatm = list(map(float, (line[30:38], line[38:45], line[46:54], line[54:60], line[60:66])))
        return numpy.array(sigatm, 'f')
