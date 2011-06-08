"""

    ##################################################
    ##                                              ##
    ##  subroutine growscesc  --  Grow Side Chain   ##
    ##                                              ##
    ##################################################

"""

from Structures.Residue import Residue
from Utils.Exceptions import GROWSCESCException

def growscesc(nn,idx,residue):


    # build side chain atoms
    assert(isinstance(residue, Residue))
    resname = residue.name()

    print("Current residue is %s-%s" % (resname, residue.sequence_number()))

    if resname is 'ALA':
        pass

    elif resname is 'ARG':
        catom3(nn,idx, 1, 1.54,109.5,  cf(1), 0) # CG
        catom3(nn,idx, 2, 1.54,109.5,  cf(2), 0) # CD
        catom3(nn,idx, 3, 1.45,109.5,  cf(3), 0) # NE
        catom3(nn,idx, 4, 1.35,120.0,  cf(4), 0) # CZ
        catom3(nn,idx, 5, 1.35,120.0,180.0, 0) # NH1
        catom3(nn,idx, 6, 1.35,120.0,120.0, 1) # NH2
        if not withhydrogen:
            catom3(nn,idx, 7, 1.11,109.4,109.4, 1) # HB2
            catom3(nn,idx, 8, 1.11,109.4,109.4,-1) # HB3
            catom3(nn,idx, 9, 1.11,109.4,109.4, 1) # HG2
            catom3(nn,idx,10, 1.11,109.4,109.4,-1) # HG3
            catom3(nn,idx,11, 1.11,109.4,109.4, 1) # HD2
            catom3(nn,idx,12, 1.11,109.4,109.4,-1) # HD3
            catom3(nn,idx,13, 1.02,120.0,120.0, 1) # HE
            catom3(nn,idx,14, 1.02,120.0,180.0, 0) # 1HH1
            catom3(nn,idx,15, 1.02,120.0,120.0, 1) # 2HH1
            catom3(nn,idx,16, 1.02,120.0,180.0, 0) # 1HH2
            catom3(nn,idx,17, 1.02,120.0,120.0, 1) # 2HH2

    elif resname is 'ASN':
        catom3(nn,idx, 1, 1.51,107.8,  cf(1), 0) # CG
        catom3(nn,idx, 2, 1.22,122.5,  cf(2), 0) # OD1
        catom3(nn,idx, 3, 1.34,112.7,124.0, 1) # ND2
        if not withhydrogen:
            catom3(nn,idx, 4, 1.11,109.4,107.9, 1) # HB2
            catom3(nn,idx, 5, 1.11,109.4,107.9,-1) # HB3
            catom3(nn,idx, 6, 1.02,119.0,  0.0, 0) # 1HD2
            catom3(nn,idx, 7, 1.02,119.0,120.0, 1) # 2HD2

    elif resname is 'ASP':
        catom3(nn,idx, 1, 1.51,107.8,  cf(1), 0) # CG
        catom3(nn,idx, 2, 1.25,117.0,  cf(2), 0) # OD1
        catom3(nn,idx, 3, 1.25,117.0,126.0, 1) # OD2
        if not withhydrogen:
            catom3(nn,idx, 4, 1.11,109.4,107.9, 1) # HB2
            catom3(nn,idx, 5, 1.11,109.4,107.9,-1) # HB3

    elif resname is 'ASH':
        catom3(nn,idx, 1, 1.51,107.8,  cf(1), 0) # CG
        catom3(nn,idx, 2, 1.25,117.0,  cf(2), 0) # OD1
        catom3(nn,idx, 3, 1.25,117.0,126.0, 1) # OD2
        if not withhydrogen:
            catom3(nn,idx, 4, 1.11,109.4,107.9, 1) # HB2
            catom3(nn,idx, 5, 1.11,109.4,107.9,-1) # HB3
            catom3(nn,idx, 6, 0.96,109.5,180.0, 0) # HD2

    elif resname is 'CYS':
        catom3(nn,idx,1, 1.82,109.0,  cf(1), 0) # SG
        if not withhydrogen:
            catom3(nn,idx, 6, 1.11,109.4,112.0, 1) # HB2
            catom3(nn,idx, 7, 1.11,109.4,112.0,-1) # HB3
            catom3(nn,idx, 8, 1.34, 96.0,180.0, 0) # HG

    elif resname is 'CYX':
        catom3(nn,idx,1, 1.82,109.0,  cf(1), 0) # SG
        if not withhydrogen:
            catom3(nn,idx, 6, 1.11,109.4,112.0, 1) # HB2
            catom3(nn,idx, 7, 1.11,109.4,112.0,-1) # HB3

    elif resname is 'GLN':
        catom3(nn,idx, 1, 1.54,109.5,  cf(1), 0) # CG
        catom3(nn,idx, 2, 1.51,107.8,  cf(2), 0) # CD
        catom3(nn,idx, 3, 1.22,122.5,  cf(3), 0) # OE1
        catom3(nn,idx, 4, 1.34,112.7,124.0, 1) # NE2
        if not withhydrogen:
            catom3(nn,idx, 5, 1.11,109.4,109.4, 1) # HB2
            catom3(nn,idx, 6, 1.11,109.4,109.4,-1) # HB3
            catom3(nn,idx, 7, 1.11,109.4,107.9, 1) # HG2
            catom3(nn,idx, 8, 1.11,109.4,107.9,-1) # HG3
            catom3(nn,idx, 9, 1.02,119.0,  0.0, 0) # 1HE2
            catom3(nn,idx,10, 1.02,119.0,120.0, 1) # 2HE2

    elif resname is 'GLU':
        catom3(nn,idx, 1, 1.54,109.5,  cf(1), 0) # CG
        catom3(nn,idx, 2, 1.51,107.8,  cf(2), 0) # CD
        catom3(nn,idx, 3, 1.25,117.0,  cf(3), 0) # OE1
        catom3(nn,idx, 4, 1.25,117.0,126.0, 1) # OE2
        if not withhydrogen:
            catom3(nn,idx, 5, 1.11,109.4,109.4, 1) # HB2
            catom3(nn,idx, 6, 1.11,109.4,109.4,-1) # HB3
            catom3(nn,idx, 7, 1.11,109.4,107.9, 1) # HG2
            catom3(nn,idx, 8, 1.11,109.4,107.9,-1) # HG3

    elif resname is 'GLH':
        catom3(nn,idx, 1, 1.54,109.5,  cf(1), 0) # CG
        catom3(nn,idx, 2, 1.51,107.8,  cf(2), 0) # CD
        catom3(nn,idx, 3, 1.25,117.0,  cf(3), 0) # OE1
        catom3(nn,idx, 4, 1.25,117.0,126.0, 1) # OE2
        if not withhydrogen:
            catom3(nn,idx, 5, 1.11,109.4,109.4, 1) # HB2
            catom3(nn,idx, 6, 1.11,109.4,109.4,-1) # HB3
            catom3(nn,idx, 7, 1.11,109.4,107.9, 1) # HG2
            catom3(nn,idx, 8, 1.11,109.4,107.9,-1) # HG3
            catom3(nn,idx, 9, 0.96,109.5,180.0, 0) # HE2

    elif resname is 'GLY':
        pass
    elif resname is 'HIS':
        catom3(nn,idx, 1, 1.50,109.5,  cf(1), 0) # CG
        catom3(nn,idx, 2, 1.35,126.0,  cf(2), 0) # ND1
        catom3(nn,idx, 3, 1.35,126.0,108.0, 1) # CD2
        catom3(nn,idx, 4, 1.35,108.0,  0.0, 0) # CD1
        catom3(nn,idx, 5, 1.35,108.0,  0.0, 0) # NE2
        if not withhydrogen:
            catom3(nn,idx, 6, 1.11,109.4,109.4, 1) # HB2
            catom3(nn,idx, 7, 1.11,109.4,109.4,-1) # HB3
            catom3(nn,idx, 8, 1.02,126.0,  0.0, 0) # HD1
            catom3(nn,idx, 9, 1.10,126.0,126.0, 1) # HD2
            catom3(nn,idx,10, 1.10,126.0,126.0, 1) # HE1
            catom3(nn,idx,11, 1.02,126.0,126.0, 1) # HE2

    elif resname is 'HID':
        catom3(nn,idx, 1, 1.50,109.5,  cf(1), 0) # CG
        catom3(nn,idx, 2, 1.35,126.0,  cf(2), 0) # ND1
        catom3(nn,idx, 3, 1.35,126.0,108.0, 1) # CD2
        catom3(nn,idx, 4, 1.35,108.0,  0.0, 0) # CD1
        catom3(nn,idx, 5, 1.35,108.0,  0.0, 0) # NE2
        if not withhydrogen:
            catom3(nn,idx, 6, 1.11,109.4,109.4, 1) # HB2
            catom3(nn,idx, 7, 1.11,109.4,109.4,-1) # HB3
            catom3(nn,idx, 8, 1.02,126.0,  0.0, 0) # HD1
            catom3(nn,idx, 9, 1.10,126.0,126.0, 1) # HD2
            catom3(nn,idx,10, 1.10,126.0,126.0, 1) # HE1

    elif resname is 'HIE':
        catom3(nn,idx, 1, 1.50,109.5,  cf(1), 0) # CG
        catom3(nn,idx, 2, 1.35,126.0,  cf(2), 0) # ND1
        catom3(nn,idx, 3, 1.35,126.0,108.0, 1) # CD2
        catom3(nn,idx, 4, 1.35,108.0,  0.0, 0) # CD1
        catom3(nn,idx, 5, 1.35,108.0,  0.0, 0) # NE2
        if not withhydrogen:
            catom3(nn,idx, 6, 1.11,109.4,109.4, 1) # HB2
            catom3(nn,idx, 7, 1.11,109.4,109.4,-1) # HB3
            catom3(nn,idx, 8, 1.10,126.0,126.0, 1) # HD2
            catom3(nn,idx, 9, 1.10,126.0,126.0, 1) # HE1
            catom3(nn,idx,10, 1.02,126.0,126.0, 1) # HE2

    elif resname is 'ILE':
        catom3(nn,idx, 1, 1.54,109.5,  cf(1), 0) # CG1
        catom3(nn,idx, 2, 1.54,109.5,109.5, 1) # CG2
        catom3(nn,idx, 3, 1.54,109.5,  cf(2), 0) # CD
        if not withhydrogen:
            catom3(nn,idx, 4, 1.11,109.4,109.4,-1) # HB
            catom3(nn,idx, 5, 1.11,109.4,109.4, 1) # 2HG1
            catom3(nn,idx, 6, 1.11,109.4,109.4,-1) # 3HG1
            catom3(nn,idx, 7, 1.11,110.0,180.0, 0) # 1HG2
            catom3(nn,idx, 8, 1.11,110.0,109.0, 1) # 2HG2
            catom3(nn,idx, 9, 1.11,110.0,109.0,-1) # 3HG2
            catom3(nn,idx,10, 1.11,110.0,180.0, 0) # 1HD1
            catom3(nn,idx,11, 1.11,110.0,109.0, 1) # 2HD1
            catom3(nn,idx,12, 1.11,110.0,109.0,-1) # 3HD1

    elif resname is 'LEU':
        catom3(nn,idx, 1, 1.54,109.5,  cf(1), 0) # CG
        catom3(nn,idx, 2, 1.54,109.5,  cf(2), 0) # CD1
        catom3(nn,idx, 3, 1.54,109.5,109.4,-1) # CD2
        if not withhydrogen:
            catom3(nn,idx, 4, 1.11,109.4,109.4, 1) # HB2
            catom3(nn,idx, 5, 1.11,109.4,109.4,-1) # HB3
            catom3(nn,idx, 6, 1.11,109.4,109.4, 1) # HG
            catom3(nn,idx, 7, 1.11,109.4,180.0, 0) # 1HD1
            catom3(nn,idx, 8, 1.11,109.4,109.4, 1) # 2HD1
            catom3(nn,idx, 9, 1.11,109.4,109.4,-1) # 3HD1
            catom3(nn,idx,10, 1.11,109.4,180.0, 0) # 1HD2
            catom3(nn,idx,11, 1.11,109.4,109.4, 1) # 2HD2
            catom3(nn,idx,12, 1.11,109.4,109.4,-1) # 3HD2

    elif resname is 'LYS':
        catom3(nn,idx, 1, 1.54,109.5,  cf(1), 0) # CG
        catom3(nn,idx, 2, 1.54,109.5,  cf(2), 0) # CD
        catom3(nn,idx, 3, 1.54,109.5,  cf(3), 0) # CE
        catom3(nn,idx, 4, 1.51,109.5,  cf(4), 0) # NZ
        if not withhydrogen:
            catom3(nn,idx, 5, 1.11,109.4,109.4, 1) # HB2
            catom3(nn,idx, 6, 1.11,109.4,109.4,-1) # HB3
            catom3(nn,idx, 7, 1.11,109.4,109.4, 1) # HG2
            catom3(nn,idx, 8, 1.11,109.4,109.4,-1) # HG3
            catom3(nn,idx, 9, 1.11,109.4,109.4, 1) # HD2
            catom3(nn,idx,10, 1.11,109.4,109.4,-1) # HD3
            catom3(nn,idx,11, 1.11,109.4,108.8, 1) # HE2
            catom3(nn,idx,12, 1.11,109.4,108.8,-1) # HE3
            catom3(nn,idx,13, 1.02,109.5,180.0, 0) # HZ1
            catom3(nn,idx,14, 1.02,109.5,109.5, 1) # HZ2
            catom3(nn,idx,15, 1.02,109.5,109.5,-1) # HZ3

    elif resname is 'LYN':
        catom3(nn,idx, 1, 1.54,109.5,  cf(1), 0) # CG
        catom3(nn,idx, 2, 1.54,109.5,  cf(2), 0) # CD
        catom3(nn,idx, 3, 1.54,109.5,  cf(3), 0) # CE
        catom3(nn,idx, 4, 1.51,109.5,  cf(4), 0) # NZ
        if not withhydrogen:
            catom3(nn,idx, 5, 1.11,109.4,109.4, 1) # HB2
            catom3(nn,idx, 6, 1.11,109.4,109.4,-1) # HB3
            catom3(nn,idx, 7, 1.11,109.4,109.4, 1) # HG2
            catom3(nn,idx, 8, 1.11,109.4,109.4,-1) # HG3
            catom3(nn,idx, 9, 1.11,109.4,109.4, 1) # HD2
            catom3(nn,idx,10, 1.11,109.4,109.4,-1) # HD3
            catom3(nn,idx,11, 1.11,109.4,108.8, 1) # HE2
            catom3(nn,idx,12, 1.11,109.4,108.8,-1) # HE3
            catom3(nn,idx,13, 1.02,109.5, 60.0, 0) # HZ2
            catom3(nn,idx,14, 1.02,109.5,300.0, 0) # HZ3

    elif resname is 'MET':
        catom3(nn,idx, 1, 1.54,109.5,  cf(1), 0) # CG
        catom3(nn,idx, 2, 1.82,109.0,  cf(2), 0) # SD
        catom3(nn,idx, 3, 1.82, 96.3,  cf(3), 0) # CE
        if not withhydrogen:
            catom3(nn,idx, 4, 1.11,109.4,109.4, 1) # HB2
            catom3(nn,idx, 5, 1.11,109.4,109.4,-1) # HB3
            catom3(nn,idx, 6, 1.11,109.4,112.0, 1) # HG2
            catom3(nn,idx, 7, 1.11,109.4,112.0,-1) # HG3
            catom3(nn,idx, 8, 1.11,112.0,180.0, 0) # HE1
            catom3(nn,idx, 9, 1.11,112.0,109.4, 1) # HE2
            catom3(nn,idx,10, 1.11,112.0,109.4,-1) # HE3

    elif resname is 'PHE':
        catom3(nn,idx, 1, 1.50,109.5,  cf(1), 0) # CG
        catom3(nn,idx, 2, 1.39,120.0,  cf(2), 0) # CD1
        catom3(nn,idx, 3, 1.39,120.0,120.0, 1) # CD2
        catom3(nn,idx, 4, 1.39,120.0,180.0, 0) # CE1
        catom3(nn,idx, 5, 1.39,120.0,180.0, 0) # CE2
        catom3(nn,idx, 6, 1.39,120.0,  0.0, 0) # CZ
        if not withhydrogen:
            catom3(nn,idx, 7, 1.11,109.4,109.4, 1) # HB2
            catom3(nn,idx, 8, 1.11,109.4,109.4,-1) # HB3
            catom3(nn,idx, 9, 1.10,120.0,120.0, 1) # HD1
            catom3(nn,idx,10, 1.10,120.0,120.0, 1) # HD2
            catom3(nn,idx,11, 1.10,120.0,120.0, 1) # HE1
            catom3(nn,idx,12, 1.10,120.0,120.0, 1) # HE2
            catom3(nn,idx,13, 1.10,120.0,120.0, 1) # HZ

    elif resname is 'PRO':
        catom3(nn,idx, 1, 1.54,107.0,  cf(1), 0) # CG
        catom3(nn,idx, 2, 1.54,107.0,  cf(2), 0) # CD
        if not withhydrogen:
            catom3(nn,idx, 3, 1.11,109.4,109.4, 1) # HB2
            catom3(nn,idx, 4, 1.11,109.4,109.4,-1) # HB3
            catom3(nn,idx, 5, 1.11,109.4,109.4, 1) # HG2
            catom3(nn,idx, 6, 1.11,109.4,109.4,-1) # HG3
            catom3(nn,idx, 7, 1.11,109.4,109.4, 1) # HD2
            catom3(nn,idx, 8, 1.11,109.4,109.4,-1) # HD3

    elif resname is 'SER':
        catom3(nn,idx, 1, 1.41,107.5,  cf(1), 0) # OG
        if not withhydrogen:
            catom3(nn,idx, 2, 1.11,109.4,106.7, 1) # HB2
            catom3(nn,idx, 3, 1.11,109.4,106.7,-1) # HB3
            catom3(nn,idx, 4, 0.94,106.9,180.0, 0) # HG

    elif resname is 'THR':
        catom3(nn,idx, 1, 1.41,107.5,  cf(1), 0) # OG1
        catom3(nn,idx, 2, 1.54,109.5,107.7, 1) # CG2
        if not withhydrogen:
            catom3(nn,idx, 3, 1.11,109.4,106.7,-1) # HB
            catom3(nn,idx, 4, 0.94,106.9,180.0, 0) # HG1
            catom3(nn,idx, 5, 1.11,110.0,180.0, 0) # 1HG2
            catom3(nn,idx, 6, 1.11,110.0,109.0, 1) # 2HG2
            catom3(nn,idx, 7, 1.11,110.0,109.0,-1) # 3HG2

    elif resname is 'TRP':
        catom3(nn,idx, 1, 1.50,109.5,  cf(1), 0) # CG
        catom3(nn,idx, 2, 1.35,126.0,  cf(2), 0) # CD1
        catom3(nn,idx, 3, 1.35,126.0,108.0, 1) # CD2
        catom3(nn,idx, 4, 1.35,108.0,  0.0, 0) # NE1
        catom3(nn,idx, 5, 1.35,108.0,  0.0, 0) # CE2
        catom3(nn,idx, 6, 1.35,120.0,180.0, 0) # CE3
        catom3(nn,idx, 7, 1.35,120.0,  0.0, 0) # CZ2
        catom3(nn,idx, 8, 1.35,120.0,  0.0, 0) # CZ3
        catom3(nn,idx, 9, 1.35,120.0,  0.0, 0) # CH2
        if not withhydrogen:
            catom3(nn,idx,10, 1.11,109.4,109.4, 1) # HB2
            catom3(nn,idx,11, 1.11,109.4,109.4,-1) # HB3
            catom3(nn,idx,12, 1.10,126.0,126.0, 1) # HD1
            catom3(nn,idx,13, 1.05,126.0,126.0, 1) # HE1
            catom3(nn,idx,14, 1.10,120.0,120.0, 1) # HE3
            catom3(nn,idx,15, 1.10,120.0,120.0, 1) # HZ2
            catom3(nn,idx,16, 1.10,120.0,120.0, 1) # HZ3
            catom3(nn,idx,17, 1.10,120.0,120.0, 1) # HH2

    elif resname is 'TYR':
        catom3(nn,idx, 1, 1.50,109.5,  cf(1), 0) # CG
        catom3(nn,idx, 2, 1.39,120.0,  cf(2), 0) # CD1
        catom3(nn,idx, 3, 1.39,120.0,120.0, 1) # CD2
        catom3(nn,idx, 4, 1.39,120.0,180.0, 0) # CE1
        catom3(nn,idx, 5, 1.39,120.0,180.0, 0) # CE2
        catom3(nn,idx, 6, 1.39,120.0,  0.0, 0) # CZ
        catom3(nn,idx, 7, 1.36,120.0,120.0, 1) # OH
        if not withhydrogen:
            catom3(nn,idx, 8, 1.11,109.4,109.4, 1) # HB2
            catom3(nn,idx, 9, 1.11,109.4,109.4,-1) # HB3
            catom3(nn,idx,10, 1.10,120.0,120.0, 1) # HD1
            catom3(nn,idx,11, 1.10,120.0,120.0, 1) # HD2
            catom3(nn,idx,12, 1.10,120.0,120.0, 1) # HE1
            catom3(nn,idx,13, 1.10,120.0,120.0, 1) # HE2
            catom3(nn,idx,14, 0.97,108.0,  0.0, 0) # HH

    elif resname is 'VAL':
        catom3(nn,idx, 1, 1.54,109.5,  cf(1), 0) # CG1
        catom3(nn,idx, 2, 1.54,109.5,109.5,-1) # CG2
        if not withhydrogen:
            catom3(nn,idx, 3, 1.11,109.4,109.4, 1) # HB
            catom3(nn,idx, 4, 1.11,109.4,180.0, 0) # 1HG1
            catom3(nn,idx, 5, 1.11,109.4,109.4, 1) # 2HG1
            catom3(nn,idx, 6, 1.11,109.4,109.4,-1) # 3HG1
            catom3(nn,idx, 7, 1.11,109.4,180.0, 0) # 1HG2
            catom3(nn,idx, 8, 1.11,109.4,109.4, 1) # 2HG2
            catom3(nn,idx, 9, 1.11,109.4,109.4,-1) # 3HG2

    else:
        raise GROWSCESCException("Residue %s is undefined." % resname)
    
    return
