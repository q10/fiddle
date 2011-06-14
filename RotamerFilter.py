
import os
from copy import deepcopy
from Parsers.RotamerLibParser import RotamerLibParser
from grow_side_chain import grow_side_chain
from VDW import VDW_RADIUS

rotamer_folder = 'rotamers'
pdb_file = 'pdbs/sample.pdb'
ALPHA = 0.8

def get_raw_rotamer_library():
    rotamer_library = {}
    for filename in os.listdir(rotamer_folder):
        r = RotamerLibParser(rotamer_folder + "/" + filename)
        rotamer_library[r.residue_name()] = r.torsions()
    return rotamer_library

def check_rotamer_clash_against_backbone(residue, chain):
    for other_residue in chain:
        if other_residue.id.sequence_number() is not residue.sequence_number():
            for atom in residue.side_chain_atoms():
                for other_atom in other_residue.backbone_atoms():
                    if atom.distance_from(other_atom) < ALPHA * (VDW_RADIUS[atom.element()] + VDW_RADIUS[other_atom.element()]):
                        return True
    return False


RAW_ROTAMER_LIBRARY = get_raw_rotamer_library()

def get_filtered_rotamer_chain(chain):
    rotamer_chain = []
    for residue in chain:
        filtered_rotamer_library = deepcopy(RAW_ROTAMER_LIBRARY[residue.name()])
        for torsion_set in filtered_rotamer_library:
            grow_side_chain(residue, torsion_set, False)
            if check_rotamer_clash_against_backbone(residue, chain):
                filtered_rotamer_library.remove(torsion_set)
        rotamer_chain.append(filtered_rotamer_library)
    return rotamer_chain