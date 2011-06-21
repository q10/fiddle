import os, shelve
from multiprocessing import Pool, cpu_count
from copy import deepcopy
from Parsers.RotamerLibParser import RotamerLibParser
from Methods.grow_side_chain import grow_side_chain
from VDW import VDW

ROTAMER_FOLDER = 'rotamers'

def get_raw_rotamer_library():
    rotamer_library = {}
    for filename in os.listdir(ROTAMER_FOLDER):
        r = RotamerLibParser(ROTAMER_FOLDER + "/" + filename)
        rotamer_library[r.residue_name()] = r.torsions()
    return rotamer_library

RAW_ROTAMER_LIBRARY = get_raw_rotamer_library()


def get_and_save_filtered_rotamer_chain_to_file(chain, filename):
    filtered_rotamer_chain = get_filtered_rotamer_chain(chain)
    file = shelve.open(filename)
    file['filtered_rotamer_chain'] = filtered_rotamer_chain
    file.close()


def get_filtered_rotamer_chain(chain):
    process_pool = Pool(processes = cpu_count()*2)
    filtered_rotamer_libraries = process_pool.map(get_filtered_rotamer_for_residue, chain.residues())
    process_pool.close() # no more tasks
    process_pool.join()  # wrap up current tasks

    rotamer_chain = {}
    for output in filtered_rotamer_libraries:
        rotamer_chain.update(output)

    return rotamer_chain


def get_filtered_rotamer_for_residue(residue):
    try:
        filtered_rotamer_library = deepcopy(RAW_ROTAMER_LIBRARY[residue.name()])
    except KeyError:
        filtered_rotamer_library = []

    for torsion_set in filtered_rotamer_library:
        grow_side_chain(residue, torsion_set)
        if check_rotamer_clash_against_backbone(residue):
            filtered_rotamer_library.remove(torsion_set)

    return {residue.sequence_number() : filtered_rotamer_library}


def check_rotamer_clash_against_backbone(residue):
    for other_residue in residue.chain().residues():
        if other_residue.sequence_number() is not residue.sequence_number():
            for atom in residue.side_chain_atoms():
                for other_atom in other_residue.backbone_atoms():
                    if atom.distance_from(other_atom) < VDW.ALPHA * (VDW.VDW_RADIUS[atom.element()] + VDW.VDW_RADIUS[other_atom.element()]):
                        return True
    return False
