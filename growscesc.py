"""

    ##################################################
    ##                                              ##
    ##  subroutine growscesc  --  Grow Side Chain   ##
    ##                                              ##
    ##################################################

"""

from Structures.Residue import Residue
from Utils.Exceptions import GROWSCESCException
from pdbatm2 import *

def growscesc(nn,idx,residue, withhydrogens):


    # build side chain atoms
    assert(isinstance(residue, Residue))
    resname = residue.name()

    print("Current residue is %s-%s" % (resname, residue.sequence_number()))

    atoms = residue.atoms_hash()
    [N, CA, CB] = atoms['N'], atoms['CA'], atoms['CB']

    if resname is 'ALA':
        pass

    elif resname is 'ARG':
        [CG, CD, NE, CZ, NH1, NH2] = atoms['CG'], atoms['CD'], atoms['NE'], atoms['CZ'], atoms['NH1'], atoms['NH2']
        pdbatm2(CG, CB, CA, N, 1.54, 109.5, cf(1), 0) # CG
        pdbatm2(CD, CG, CB, CA, 1.54, 109.5, cf(2), 0) # CD
        pdbatm2(NE, CD, CG, CB, 1.45, 109.5, cf(3), 0) # NE
        pdbatm2(CZ, NE, CD, CG, 1.35, 120.0, cf(4), 0) # CZ
        pdbatm2(NH1, CZ, NE, CD, 1.35, 120.0, 180.0, 0) # NH1
        pdbatm2(NH2, CZ, NE, NH1, 1.35, 120.0, 120.0, 1) # NH2


        '''
        pdbatm2(nn,idx, 1, 1.54, 109.5, cf(1), 0) # CG
        pdbatm2(nn,idx, 2, 1.54, 109.5, cf(2), 0) # CD
        pdbatm2(nn,idx, 3, 1.45, 109.5, cf(3), 0) # NE
        pdbatm2(nn,idx, 4, 1.35, 120.0, cf(4), 0) # CZ
        pdbatm2(nn,idx, 5, 1.35, 120.0, 180.0, 0) # NH1
        pdbatm2(nn,idx, 6, 1.35, 120.0, 120.0, 1) # NH2
        '''
        if withhydrogens:
            [HB2, HB3, HG2, HG3, HD2, HD3] = atoms['HB2'], atoms['HB3'], atoms['HG2'], atoms['HG3'], atoms['HD2'], atoms['HD3']
            [HE, IHH1, IIHH1, IHH2, IIHH2] = atoms['HE'], atoms['1HH1'], atoms['2HH1'], atoms['1HH2'], atoms['2HH2']
            pdbatm2(HB2, CB, CA, CG, 1.11, 109.4, 109.4, 1) # HB2
            pdbatm2(HB3, CB, CA, CG, 1.11, 109.4, 109.4, -1) # HB3
            pdbatm2(HG2, CG, CB, CD, 1.11, 109.4, 109.4, 1) # HG2
            pdbatm2(HG3, CG, CB, CD, 1.11, 109.4, 109.4, -1) # HG3
            pdbatm2(HD2, CD, CG, NE, 1.11, 109.4, 109.4, 1) # HD2
            pdbatm2(HD3, CD, CG, NE, 1.11, 109.4, 109.4, -1) # HD3
            pdbatm2(HE, NE, CD, CZ, 1.02, 120.0, 120.0, 1) # HE
            pdbatm2(IHH1, NH1, CZ, NE, 1.02, 120.0, 180.0, 0) # 1HH1
            pdbatm2(IIHH1, NH1, CZ, NE, 1.02, 120.0, 120.0, 1) # 2HH1
            pdbatm2(IHH2, NH2, CZ, NE, 1.02, 120.0, 180.0, 0) # 1HH2
            pdbatm2(IIHH2, NH2, CZ, NE, 1.02, 120.0, 120.0, 1) # 2HH2

            '''
            pdbatm2(nn,idx, 7, 1.11, 109.4, 109.4, 1) # HB2
            pdbatm2(nn,idx, 8, 1.11, 109.4, 109.4, -1) # HB3
            pdbatm2(nn,idx, 9, 1.11, 109.4, 109.4, 1) # HG2
            pdbatm2(nn,idx,10, 1.11, 109.4, 109.4, -1) # HG3
            pdbatm2(nn,idx,11, 1.11, 109.4, 109.4, 1) # HD2
            pdbatm2(nn,idx,12, 1.11, 109.4, 109.4, -1) # HD3
            pdbatm2(nn,idx,13, 1.02, 120.0, 120.0, 1) # HE
            pdbatm2(nn,idx,14, 1.02, 120.0, 180.0, 0) # 1HH1
            pdbatm2(nn,idx,15, 1.02, 120.0, 120.0, 1) # 2HH1
            pdbatm2(nn,idx,16, 1.02, 120.0, 180.0, 0) # 1HH2
            pdbatm2(nn,idx,17, 1.02, 120.0, 120.0, 1) # 2HH2
            '''
    elif resname is 'ASN':
        [CG, OD1, ND2] = atoms['CG'], atoms['OD1'], atoms['ND2']
        pdbatm2(CG, CB, CA, N, 1.51, 107.8, cf(1), 0) # CG
        pdbatm2(OD1, CG, CB, CA, 1.22, 122.5, cf(2), 0) # OD1
        pdbatm2(ND2, CG, CB, OD1, 1.34, 112.7, 124.0, 1) # ND2
        if withhydrogens:
            [HB2, HB3, IHD2, IIHD2] = atoms['HB2'], atoms['HB3'], atoms['1HD2'], atoms['2HD2']
            pdbatm2(HB2, CB, CA, CG, 1.11, 109.4, 107.9, 1) # HB2
            pdbatm2(HB3, CB, CA, CG, 1.11, 109.4, 107.9, -1) # HB3
            pdbatm2(IHD2, ND2, CG, CB, 1.02, 119.0,  0.0, 0) # 1HD2
            pdbatm2(IIHD2, ND2, CG, IHD2, 1.02, 119.0, 120.0, 1) # 2HD2

    elif resname is 'ASP' or resname is 'ASH':  # ASH is protonated ASN
        [CG, OD1, OD2] = atoms['CG'], atoms['OD1'], atoms['OD2']
        pdbatm2(CG, CB, CA, N, 1.51, 107.8, cf(1), 0) # CG
        pdbatm2(OD1, CG, CB, CA, 1.25, 117.0, cf(2), 0) # OD1
        pdbatm2(OD2, CG, CB, OD1, 1.25, 117.0, 126.0, 1) # OD2
        if withhydrogens:
            [HB2, HB3] = atoms['HB2'], atoms['HB3']
            pdbatm2(HB2, CB, CA, CG, 1.11, 109.4, 107.9, 1) # HB2
            pdbatm2(HB3, CB, CA, CG, 1.11, 109.4, 107.9, -1) # HB3
            if resname is 'ASH':
                pdbatm2(atoms['HD2'], OD2, CG, CB, 0.96, 109.5, 180.0, 0) # HD2

    elif resname is 'CYS' or resname is 'CYX':  # CYX is deprotonated CYS
        SG = atoms['SG']
        pdbatm2(SG, CB, CA, N, 1.82, 109.0, cf(1), 0) # SG
        if withhydrogens:
            [HB2, HB3] = atoms['HB2'], atoms['HB3']
            pdbatm2(HB2, CB, CA, SG, 1.11, 109.4, 112.0, 1) # HB2
            pdbatm2(HB3, CB, CA, SG, 1.11, 109.4, 112.0, -1) # HB3
            if resname is 'CYS':
                pdbatm2(atoms['HG'], SG, CB, CA, 1.34, 96.0, 180.0, 0) # HG

    elif resname is 'GLN':
        pdbatm2(nn,idx, 1, 1.54, 109.5, cf(1), 0) # CG
        pdbatm2(nn,idx, 2, 1.51, 107.8, cf(2), 0) # CD
        pdbatm2(nn,idx, 3, 1.22, 122.5, cf(3), 0) # OE1
        pdbatm2(nn,idx, 4, 1.34, 112.7, 124.0, 1) # NE2
        if withhydrogens:
            pdbatm2(nn,idx, 5, 1.11, 109.4, 109.4, 1) # HB2
            pdbatm2(nn,idx, 6, 1.11, 109.4, 109.4, -1) # HB3
            pdbatm2(nn,idx, 7, 1.11, 109.4, 107.9, 1) # HG2
            pdbatm2(nn,idx, 8, 1.11, 109.4, 107.9, -1) # HG3
            pdbatm2(nn,idx, 9, 1.02, 119.0,  0.0, 0) # 1HE2
            pdbatm2(nn,idx,10, 1.02, 119.0, 120.0, 1) # 2HE2

    elif resname is 'GLU' or resname is 'GLH':  # GLH is protonated GLU
        pdbatm2(nn,idx, 1, 1.54, 109.5, cf(1), 0) # CG
        pdbatm2(nn,idx, 2, 1.51, 107.8, cf(2), 0) # CD
        pdbatm2(nn,idx, 3, 1.25, 117.0, cf(3), 0) # OE1
        pdbatm2(nn,idx, 4, 1.25, 117.0, 126.0, 1) # OE2
        if withhydrogens:
            pdbatm2(nn,idx, 5, 1.11, 109.4, 109.4, 1) # HB2
            pdbatm2(nn,idx, 6, 1.11, 109.4, 109.4, -1) # HB3
            pdbatm2(nn,idx, 7, 1.11, 109.4, 107.9, 1) # HG2
            pdbatm2(nn,idx, 8, 1.11, 109.4, 107.9, -1) # HG3
            if resname is 'GLH':
                pdbatm2(nn,idx, 9, 0.96, 109.5, 180.0, 0) # HE2

    elif resname is 'GLY':
        pass

    elif resname is 'HIS' or resname is 'HID' or resname is 'HIE':
        pdbatm2(nn,idx, 1, 1.50, 109.5, cf(1), 0) # CG
        pdbatm2(nn,idx, 2, 1.35, 126.0, cf(2), 0) # ND1
        pdbatm2(nn,idx, 3, 1.35, 126.0, 108.0, 1) # CD2
        pdbatm2(nn,idx, 4, 1.35, 108.0,  0.0, 0) # CD1
        pdbatm2(nn,idx, 5, 1.35, 108.0,  0.0, 0) # NE2
        if withhydrogens:
            pdbatm2(nn,idx, 6, 1.11, 109.4, 109.4, 1) # HB2
            pdbatm2(nn,idx, 7, 1.11, 109.4, 109.4, -1) # HB3
            if resname is 'HIS' or resname is 'HID':
                pdbatm2(nn,idx, 8, 1.02, 126.0,  0.0, 0) # HD1
            pdbatm2(nn,idx, 9, 1.10, 126.0, 126.0, 1) # HD2
            pdbatm2(nn,idx,10, 1.10, 126.0, 126.0, 1) # HE1
            if resname is 'HIS' or resname is 'HIE':
                pdbatm2(nn,idx,11, 1.02, 126.0, 126.0, 1) # HE2

    elif resname is 'ILE':
        pdbatm2(nn,idx, 1, 1.54, 109.5, cf(1), 0) # CG1
        pdbatm2(nn,idx, 2, 1.54, 109.5, 109.5, 1) # CG2
        pdbatm2(nn,idx, 3, 1.54, 109.5, cf(2), 0) # CD
        if withhydrogens:
            pdbatm2(nn,idx, 4, 1.11, 109.4, 109.4, -1) # HB
            pdbatm2(nn,idx, 5, 1.11, 109.4, 109.4, 1) # 2HG1
            pdbatm2(nn,idx, 6, 1.11, 109.4, 109.4, -1) # 3HG1
            pdbatm2(nn,idx, 7, 1.11, 110.0, 180.0, 0) # 1HG2
            pdbatm2(nn,idx, 8, 1.11, 110.0, 109.0, 1) # 2HG2
            pdbatm2(nn,idx, 9, 1.11, 110.0, 109.0, -1) # 3HG2
            pdbatm2(nn,idx,10, 1.11, 110.0, 180.0, 0) # 1HD1
            pdbatm2(nn,idx,11, 1.11, 110.0, 109.0, 1) # 2HD1
            pdbatm2(nn,idx,12, 1.11, 110.0, 109.0, -1) # 3HD1

    elif resname is 'LEU':
        pdbatm2(nn,idx, 1, 1.54, 109.5, cf(1), 0) # CG
        pdbatm2(nn,idx, 2, 1.54, 109.5, cf(2), 0) # CD1
        pdbatm2(nn,idx, 3, 1.54, 109.5, 109.4, -1) # CD2
        if withhydrogens:
            pdbatm2(nn,idx, 4, 1.11, 109.4, 109.4, 1) # HB2
            pdbatm2(nn,idx, 5, 1.11, 109.4, 109.4, -1) # HB3
            pdbatm2(nn,idx, 6, 1.11, 109.4, 109.4, 1) # HG
            pdbatm2(nn,idx, 7, 1.11, 109.4, 180.0, 0) # 1HD1
            pdbatm2(nn,idx, 8, 1.11, 109.4, 109.4, 1) # 2HD1
            pdbatm2(nn,idx, 9, 1.11, 109.4, 109.4, -1) # 3HD1
            pdbatm2(nn,idx,10, 1.11, 109.4, 180.0, 0) # 1HD2
            pdbatm2(nn,idx,11, 1.11, 109.4, 109.4, 1) # 2HD2
            pdbatm2(nn,idx,12, 1.11, 109.4, 109.4, -1) # 3HD2

    elif resname is 'LYS':
        pdbatm2(nn,idx, 1, 1.54, 109.5, cf(1), 0) # CG
        pdbatm2(nn,idx, 2, 1.54, 109.5, cf(2), 0) # CD
        pdbatm2(nn,idx, 3, 1.54, 109.5, cf(3), 0) # CE
        pdbatm2(nn,idx, 4, 1.51, 109.5, cf(4), 0) # NZ
        if withhydrogens:
            pdbatm2(nn,idx, 5, 1.11, 109.4, 109.4, 1) # HB2
            pdbatm2(nn,idx, 6, 1.11, 109.4, 109.4, -1) # HB3
            pdbatm2(nn,idx, 7, 1.11, 109.4, 109.4, 1) # HG2
            pdbatm2(nn,idx, 8, 1.11, 109.4, 109.4, -1) # HG3
            pdbatm2(nn,idx, 9, 1.11, 109.4, 109.4, 1) # HD2
            pdbatm2(nn,idx,10, 1.11, 109.4, 109.4, -1) # HD3
            pdbatm2(nn,idx,11, 1.11, 109.4, 108.8, 1) # HE2
            pdbatm2(nn,idx,12, 1.11, 109.4, 108.8, -1) # HE3
            pdbatm2(nn,idx,13, 1.02, 109.5, 180.0, 0) # HZ1
            pdbatm2(nn,idx,14, 1.02, 109.5, 109.5, 1) # HZ2
            pdbatm2(nn,idx,15, 1.02, 109.5, 109.5, -1) # HZ3

    elif resname is 'LYN':
        pdbatm2(nn,idx, 1, 1.54, 109.5, cf(1), 0) # CG
        pdbatm2(nn,idx, 2, 1.54, 109.5, cf(2), 0) # CD
        pdbatm2(nn,idx, 3, 1.54, 109.5, cf(3), 0) # CE
        pdbatm2(nn,idx, 4, 1.51, 109.5, cf(4), 0) # NZ
        if withhydrogens:
            pdbatm2(nn,idx, 5, 1.11, 109.4, 109.4, 1) # HB2
            pdbatm2(nn,idx, 6, 1.11, 109.4, 109.4, -1) # HB3
            pdbatm2(nn,idx, 7, 1.11, 109.4, 109.4, 1) # HG2
            pdbatm2(nn,idx, 8, 1.11, 109.4, 109.4, -1) # HG3
            pdbatm2(nn,idx, 9, 1.11, 109.4, 109.4, 1) # HD2
            pdbatm2(nn,idx,10, 1.11, 109.4, 109.4, -1) # HD3
            pdbatm2(nn,idx,11, 1.11, 109.4, 108.8, 1) # HE2
            pdbatm2(nn,idx,12, 1.11, 109.4, 108.8, -1) # HE3
            pdbatm2(nn,idx,13, 1.02, 109.5, 60.0, 0) # HZ2
            pdbatm2(nn,idx,14, 1.02, 109.5, 300.0, 0) # HZ3

    elif resname is 'MET':
        pdbatm2(nn,idx, 1, 1.54, 109.5, cf(1), 0) # CG
        pdbatm2(nn,idx, 2, 1.82, 109.0, cf(2), 0) # SD
        pdbatm2(nn,idx, 3, 1.82, 96.3, cf(3), 0) # CE
        if withhydrogens:
            pdbatm2(nn,idx, 4, 1.11, 109.4, 109.4, 1) # HB2
            pdbatm2(nn,idx, 5, 1.11, 109.4, 109.4, -1) # HB3
            pdbatm2(nn,idx, 6, 1.11, 109.4, 112.0, 1) # HG2
            pdbatm2(nn,idx, 7, 1.11, 109.4, 112.0, -1) # HG3
            pdbatm2(nn,idx, 8, 1.11, 112.0, 180.0, 0) # HE1
            pdbatm2(nn,idx, 9, 1.11, 112.0, 109.4, 1) # HE2
            pdbatm2(nn,idx,10, 1.11, 112.0, 109.4, -1) # HE3

    elif resname is 'PHE':
        pdbatm2(nn,idx, 1, 1.50, 109.5, cf(1), 0) # CG
        pdbatm2(nn,idx, 2, 1.39, 120.0, cf(2), 0) # CD1
        pdbatm2(nn,idx, 3, 1.39, 120.0, 120.0, 1) # CD2
        pdbatm2(nn,idx, 4, 1.39, 120.0, 180.0, 0) # CE1
        pdbatm2(nn,idx, 5, 1.39, 120.0, 180.0, 0) # CE2
        pdbatm2(nn,idx, 6, 1.39, 120.0,  0.0, 0) # CZ
        if withhydrogens:
            pdbatm2(nn,idx, 7, 1.11, 109.4, 109.4, 1) # HB2
            pdbatm2(nn,idx, 8, 1.11, 109.4, 109.4, -1) # HB3
            pdbatm2(nn,idx, 9, 1.10, 120.0, 120.0, 1) # HD1
            pdbatm2(nn,idx,10, 1.10, 120.0, 120.0, 1) # HD2
            pdbatm2(nn,idx,11, 1.10, 120.0, 120.0, 1) # HE1
            pdbatm2(nn,idx,12, 1.10, 120.0, 120.0, 1) # HE2
            pdbatm2(nn,idx,13, 1.10, 120.0, 120.0, 1) # HZ

    elif resname is 'PRO':
        pdbatm2(nn,idx, 1, 1.54, 107.0, cf(1), 0) # CG
        pdbatm2(nn,idx, 2, 1.54, 107.0, cf(2), 0) # CD
        if withhydrogens:
            pdbatm2(nn,idx, 3, 1.11, 109.4, 109.4, 1) # HB2
            pdbatm2(nn,idx, 4, 1.11, 109.4, 109.4, -1) # HB3
            pdbatm2(nn,idx, 5, 1.11, 109.4, 109.4, 1) # HG2
            pdbatm2(nn,idx, 6, 1.11, 109.4, 109.4, -1) # HG3
            pdbatm2(nn,idx, 7, 1.11, 109.4, 109.4, 1) # HD2
            pdbatm2(nn,idx, 8, 1.11, 109.4, 109.4, -1) # HD3

    elif resname is 'SER':
        pdbatm2(nn,idx, 1, 1.41, 107.5, cf(1), 0) # OG
        if withhydrogens:
            pdbatm2(nn,idx, 2, 1.11, 109.4, 106.7, 1) # HB2
            pdbatm2(nn,idx, 3, 1.11, 109.4, 106.7, -1) # HB3
            pdbatm2(nn,idx, 4, 0.94, 106.9, 180.0, 0) # HG

    elif resname is 'THR':
        pdbatm2(nn,idx, 1, 1.41, 107.5, cf(1), 0) # OG1
        pdbatm2(nn,idx, 2, 1.54, 109.5, 107.7, 1) # CG2
        if withhydrogens:
            pdbatm2(nn,idx, 3, 1.11, 109.4, 106.7, -1) # HB
            pdbatm2(nn,idx, 4, 0.94, 106.9, 180.0, 0) # HG1
            pdbatm2(nn,idx, 5, 1.11, 110.0, 180.0, 0) # 1HG2
            pdbatm2(nn,idx, 6, 1.11, 110.0, 109.0, 1) # 2HG2
            pdbatm2(nn,idx, 7, 1.11, 110.0, 109.0, -1) # 3HG2

    elif resname is 'TRP':
        pdbatm2(nn,idx, 1, 1.50, 109.5, cf(1), 0) # CG
        pdbatm2(nn,idx, 2, 1.35, 126.0, cf(2), 0) # CD1
        pdbatm2(nn,idx, 3, 1.35, 126.0, 108.0, 1) # CD2
        pdbatm2(nn,idx, 4, 1.35, 108.0,  0.0, 0) # NE1
        pdbatm2(nn,idx, 5, 1.35, 108.0,  0.0, 0) # CE2
        pdbatm2(nn,idx, 6, 1.35, 120.0, 180.0, 0) # CE3
        pdbatm2(nn,idx, 7, 1.35, 120.0,  0.0, 0) # CZ2
        pdbatm2(nn,idx, 8, 1.35, 120.0,  0.0, 0) # CZ3
        pdbatm2(nn,idx, 9, 1.35, 120.0,  0.0, 0) # CH2
        if withhydrogens:
            pdbatm2(nn,idx,10, 1.11, 109.4, 109.4, 1) # HB2
            pdbatm2(nn,idx,11, 1.11, 109.4, 109.4, -1) # HB3
            pdbatm2(nn,idx,12, 1.10, 126.0, 126.0, 1) # HD1
            pdbatm2(nn,idx,13, 1.05, 126.0, 126.0, 1) # HE1
            pdbatm2(nn,idx,14, 1.10, 120.0, 120.0, 1) # HE3
            pdbatm2(nn,idx,15, 1.10, 120.0, 120.0, 1) # HZ2
            pdbatm2(nn,idx,16, 1.10, 120.0, 120.0, 1) # HZ3
            pdbatm2(nn,idx,17, 1.10, 120.0, 120.0, 1) # HH2

    elif resname is 'TYR':
        pdbatm2(nn,idx, 1, 1.50, 109.5, cf(1), 0) # CG
        pdbatm2(nn,idx, 2, 1.39, 120.0, cf(2), 0) # CD1
        pdbatm2(nn,idx, 3, 1.39, 120.0, 120.0, 1) # CD2
        pdbatm2(nn,idx, 4, 1.39, 120.0, 180.0, 0) # CE1
        pdbatm2(nn,idx, 5, 1.39, 120.0, 180.0, 0) # CE2
        pdbatm2(nn,idx, 6, 1.39, 120.0,  0.0, 0) # CZ
        pdbatm2(nn,idx, 7, 1.36, 120.0, 120.0, 1) # OH
        if withhydrogens:
            pdbatm2(nn,idx, 8, 1.11, 109.4, 109.4, 1) # HB2
            pdbatm2(nn,idx, 9, 1.11, 109.4, 109.4, -1) # HB3
            pdbatm2(nn,idx,10, 1.10, 120.0, 120.0, 1) # HD1
            pdbatm2(nn,idx,11, 1.10, 120.0, 120.0, 1) # HD2
            pdbatm2(nn,idx,12, 1.10, 120.0, 120.0, 1) # HE1
            pdbatm2(nn,idx,13, 1.10, 120.0, 120.0, 1) # HE2
            pdbatm2(nn,idx,14, 0.97, 108.0,  0.0, 0) # HH

    elif resname is 'VAL':
        pdbatm2(nn,idx, 1, 1.54, 109.5, cf(1), 0) # CG1
        pdbatm2(nn,idx, 2, 1.54, 109.5, 109.5, -1) # CG2
        if withhydrogens:
            pdbatm2(nn,idx, 3, 1.11, 109.4, 109.4, 1) # HB
            pdbatm2(nn,idx, 4, 1.11, 109.4, 180.0, 0) # 1HG1
            pdbatm2(nn,idx, 5, 1.11, 109.4, 109.4, 1) # 2HG1
            pdbatm2(nn,idx, 6, 1.11, 109.4, 109.4, -1) # 3HG1
            pdbatm2(nn,idx, 7, 1.11, 109.4, 180.0, 0) # 1HG2
            pdbatm2(nn,idx, 8, 1.11, 109.4, 109.4, 1) # 2HG2
            pdbatm2(nn,idx, 9, 1.11, 109.4, 109.4, -1) # 3HG2

    else:
        raise GROWSCESCException("Residue %s is undefined." % resname)
    
    return
