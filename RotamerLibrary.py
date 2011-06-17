import os, shelve
from copy import deepcopy
from Parsers.RotamerLibParser import RotamerLibParser
from grow_side_chain import grow_side_chain
from VDW import VDW

class RotamerLibrary:
    ROTAMER_FOLDER = 'rotamers'

    @staticmethod
    def get_raw_rotamer_library():
        rotamer_library = {}
        for filename in os.listdir(RotamerLibrary.ROTAMER_FOLDER):
            r = RotamerLibParser(RotamerLibrary.ROTAMER_FOLDER + "/" + filename)
            rotamer_library[r.residue_name()] = r.torsions()
        return rotamer_library


    @staticmethod
    def get_and_save_filtered_rotamer_chain_to_file(chain, filename):
        filtered_rotamer_chain = RotamerLibrary.get_filtered_rotamer_chain(chain)
        file = shelve.open(filename)
        file['filtered_rotamer_chain'] = filtered_rotamer_chain
        file.close()


    @staticmethod
    def get_filtered_rotamer_chain(chain):
        RAW_ROTAMER_LIBRARY = RotamerLibrary.get_raw_rotamer_library()
        rotamer_chain = []

        for residue in chain.residues():
            try:
                filtered_rotamer_library = deepcopy(RAW_ROTAMER_LIBRARY[residue.name()])
            except KeyError:
                filtered_rotamer_library = []
            
            for torsion_set in filtered_rotamer_library:
                grow_side_chain(residue, torsion_set)
                if RotamerLibrary.check_rotamer_clash_against_backbone(residue, chain):
                    filtered_rotamer_library.remove(torsion_set)
            rotamer_chain.append(filtered_rotamer_library)

        return rotamer_chain


    @staticmethod
    def check_rotamer_clash_against_backbone(residue, chain):
        for other_residue in chain.residues():
            if other_residue.sequence_number() is not residue.sequence_number():
                for atom in residue.side_chain_atoms():
                    for other_atom in other_residue.backbone_atoms():
                        if atom.distance_from(other_atom) < VDW.ALPHA * (VDW.VDW_RADIUS[atom.element()] + VDW.VDW_RADIUS[other_atom.element()]):
                            return True
        return False
