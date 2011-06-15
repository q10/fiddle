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

def grow_side_chain(residue, torsions):

    # build side chain atoms
    assert(isinstance(residue, Residue))
    resname = residue.name()

    print("Current residue is %s-%s" % (resname, residue.sequence_number()))

    atoms = residue.atoms_hash()
    [N, CA, CB] = atoms['N'], atoms['CA'], atoms['CB']

    if resname == 'ALA':
        pass

    elif resname == 'ARG':
        [CG, CD, NE, CZ, NH1, NH2] = atoms['CG'], atoms['CD'], atoms['NE'], atoms['CZ'], atoms['NH1'], atoms['NH2']
        pdbatm2(CG, CB, CA, N, 1.54, 109.5, torsions[0], 0) # CG
        pdbatm2(CD, CG, CB, CA, 1.54, 109.5, torsions[1], 0) # CD
        pdbatm2(NE, CD, CG, CB, 1.45, 109.5, torsions[2], 0) # NE
        pdbatm2(CZ, NE, CD, CG, 1.35, 120.0, torsions[3], 0) # CZ
        pdbatm2(NH1, CZ, NE, CD, 1.35, 120.0, 180.0, 0) # NH1
        pdbatm2(NH2, CZ, NE, NH1, 1.35, 120.0, 120.0, 1) # NH2
        if residue.has_hydrogens():
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
            pdbatm2(IIHH1, NH1, CZ, IHH1, 1.02, 120.0, 120.0, 1) # 2HH1
            pdbatm2(IHH2, NH2, CZ, NE, 1.02, 120.0, 180.0, 0) # 1HH2
            pdbatm2(IIHH2, NH2, CZ, IHH2, 1.02, 120.0, 120.0, 1) # 2HH2

    elif resname == 'ASN':
        [CG, OD1, ND2] = atoms['CG'], atoms['OD1'], atoms['ND2']
        pdbatm2(CG, CB, CA, N, 1.51, 107.8, torsions[0], 0) # CG
        pdbatm2(OD1, CG, CB, CA, 1.22, 122.5, torsions[1], 0) # OD1
        pdbatm2(ND2, CG, CB, OD1, 1.34, 112.7, 124.0, 1) # ND2
        if residue.has_hydrogens():
            [HB2, HB3, IHD2, IIHD2] = atoms['HB2'], atoms['HB3'], atoms['1HD2'], atoms['2HD2']
            pdbatm2(HB2, CB, CA, CG, 1.11, 109.4, 107.9, 1) # HB2
            pdbatm2(HB3, CB, CA, CG, 1.11, 109.4, 107.9, -1) # HB3
            pdbatm2(IHD2, ND2, CG, CB, 1.02, 119.0,  0.0, 0) # 1HD2
            pdbatm2(IIHD2, ND2, CG, IHD2, 1.02, 119.0, 120.0, 1) # 2HD2

    elif resname == 'ASP' or resname == 'ASH':  # ASH is protonated ASN
        [CG, OD1, OD2] = atoms['CG'], atoms['OD1'], atoms['OD2']
        pdbatm2(CG, CB, CA, N, 1.51, 107.8, torsions[0], 0) # CG
        pdbatm2(OD1, CG, CB, CA, 1.25, 117.0, torsions[1], 0) # OD1
        pdbatm2(OD2, CG, CB, OD1, 1.25, 117.0, 126.0, 1) # OD2
        if residue.has_hydrogens():
            [HB2, HB3] = atoms['HB2'], atoms['HB3']
            pdbatm2(HB2, CB, CA, CG, 1.11, 109.4, 107.9, 1) # HB2
            pdbatm2(HB3, CB, CA, CG, 1.11, 109.4, 107.9, -1) # HB3
            if resname == 'ASH':
                pdbatm2(atoms['HD2'], OD2, CG, CB, 0.96, 109.5, 180.0, 0) # HD2

    elif resname == 'CYS' or resname == 'CYX':  # CYX is deprotonated CYS
        SG = atoms['SG']
        pdbatm2(SG, CB, CA, N, 1.82, 109.0, torsions[0], 0) # SG
        if residue.has_hydrogens():
            [HB2, HB3] = atoms['HB2'], atoms['HB3']
            pdbatm2(HB2, CB, CA, SG, 1.11, 109.4, 112.0, 1) # HB2
            pdbatm2(HB3, CB, CA, SG, 1.11, 109.4, 112.0, -1) # HB3
            if resname == 'CYS':
                pdbatm2(atoms['HG'], SG, CB, CA, 1.34, 96.0, 180.0, 0) # HG

    elif resname == 'GLN':
        [CG, CD, OE1, NE2] = atoms['CG'], atoms['CD'], atoms['OE1'], atoms['NE2']
        pdbatm2(CG, CB, CA, N, 1.54, 109.5, torsions[0], 0) # CG
        pdbatm2(CD, CG, CB, CA, 1.51, 107.8, torsions[1], 0) # CD
        pdbatm2(OE1, CD, CG, CB, 1.22, 122.5, torsions[2], 0) # OE1
        pdbatm2(NE2, CD, CG, OE1, 4, 1.34, 112.7, 124.0, 1) # NE2
        if residue.has_hydrogens():
            [HB2, HB3, HG2, HG3, IHE2, IIHE2] = atoms['HB2'], atoms['HB3'], atoms['HG2'], atoms['HG3'], atoms['1HE2'], atoms['2HE2']
            pdbatm2(HB2, CB, CA, CG, 1.11, 109.4, 109.4, 1) # HB2
            pdbatm2(HB3, CB, CA, CG, 1.11, 109.4, 109.4, -1) # HB3
            pdbatm2(HG2, CG, CB, CD, 1.11, 109.4, 107.9, 1) # HG2
            pdbatm2(HG3, CG, CB, CD, 1.11, 109.4, 107.9, -1) # HG3
            pdbatm2(IHE2, NE2, CD, CG, 1.02, 119.0,  0.0, 0) # 1HE2
            pdbatm2(IIHE2, NE2, CD, IHE2, 1.02, 119.0, 120.0, 1) # 2HE2

    elif resname == 'GLU' or resname == 'GLH':  # GLH is protonated GLU
        [CG, CD, OE1, OE2] = atoms['CG'], atoms['CD'], atoms['OE1'], atoms['OE2']
        pdbatm2(CG, CB, CA, N, 1.54, 109.5, torsions[0], 0) # CG
        pdbatm2(CD, CG, CB, CA, 1.51, 107.8, torsions[1], 0) # CD
        pdbatm2(OE1, CD, CG, CB, 1.25, 117.0, torsions[2], 0) # OE1
        pdbatm2(OE2, CD, CG, OE1, 1.25, 117.0, 126.0, 1) # OE2
        if residue.has_hydrogens():
            [HB2, HB3, HG2, HG3] = atoms['HB2'], atoms['HB3'], atoms['HG2'], atoms['HG3']
            pdbatm2(HB2, CB, CA, CG, 1.11, 109.4, 109.4, 1) # HB2
            pdbatm2(HB3, CB, CA, CG, 1.11, 109.4, 109.4, -1) # HB3
            pdbatm2(HG2, CG, CB, CD, 1.11, 109.4, 107.9, 1) # HG2
            pdbatm2(HG3, CG, CB, CD, 1.11, 109.4, 107.9, -1) # HG3
            if resname == 'GLH':
                pdbatm2(atoms['HE2'], OE2, CD, CG, 0.96, 109.5, 180.0, 0) # HE2

    elif resname == 'GLY':
        pass

    elif resname == 'HIS' or resname == 'HID' or resname == 'HIE':
        [CG, ND1, CD2, CE1, NE2] = atoms['CG'], atoms['ND1'], atoms['CD2'], atoms['CE1'], atoms['NE2']
        pdbatm2(CG, CB, CA, N, 1.50, 109.5, torsions[0], 0) # CG
        pdbatm2(ND1, CG, CB, CA, 1.35, 126.0, torsions[1], 0) # ND1
        pdbatm2(CD2, CG, CB, ND1, 1.35, 126.0, 108.0, 1) # CD2
        pdbatm2(CE1, ND1, CG, CD2, 1.35, 108.0,  0.0, 0) # CE1
        pdbatm2(NE2, CD2, CG, ND1, 1.35, 108.0,  0.0, 0) # NE2
        if residue.has_hydrogens():
            pdbatm2(atoms['HB2'], CB, CA, CG, 1.11, 109.4, 109.4, 1) # HB2
            pdbatm2(atoms['HB3'], CB, CA, CG, 1.11, 109.4, 109.4, -1) # HB3
            if resname == 'HID' or resname == 'HIE':
                pdbatm2(atoms['HD2'], CD2, CG, NE2, 1.10, 126.0, 126.0, 1) # HD2
                pdbatm2(atoms['HE1'], CE1, ND1, NE2, 1.10, 126.0, 126.0, 1) # HE1
            if resname == 'HIS' or resname == 'HID':
                pdbatm2(atoms['HD1'], ND1, CG, CD2, 1.02, 126.0,  0.0, 0) # HD1
            if resname == 'HIS' or resname == 'HIE':
                pdbatm2(atoms['HE2'], NE2, CD2, CE1, 1.02, 126.0, 126.0, 1) # HE2

    elif resname == 'ILE':
        [CG1, CG2, CD] = atoms['CG1'], atoms['CG2'], atoms['CD']
        pdbatm2(CG1, CB, CA, N, 1.54, 109.5, torsions[0], 0) # CG1
        pdbatm2(CD, CG1, CB, CA, 1.54, 109.5, torsions[1], 0) # CD
        pdbatm2(CG2, CB, CA, CG1, 1.54, 109.5, 109.5, 1) # CG2
        if residue.has_hydrogens():
            [HB, IIHG1, IIIHG1, IHG2, IIHG2] = atoms['HB'], atoms['2HG1'], atoms['3HG1'], atoms['1HG2'], atoms['2HG2']
            [IIIHG2, IHD1, IIHD1, IIIHD1] = atoms['3HG2'], atoms['1HD1'], atoms['2HD1'], atoms['3HD1']
            pdbatm2(HB, CB, CA, CG1, 1.11, 109.4, 109.4, -1) # HB
            pdbatm2(IIHG1, CG1, CB, CD, 1.11, 109.4, 109.4, 1) # 2HG1
            pdbatm2(IIIHG1, CG1, CB, CD, 1.11, 109.4, 109.4, -1) # 3HG1
            pdbatm2(IHG2, CG2, CB, CA, 1.11, 110.0, 180.0, 0) # 1HG2
            pdbatm2(IIHG2, CG2, CB, IHG2, 1.11, 110.0, 109.0, 1) # 2HG2
            pdbatm2(IIIHG2, CG2, CB, IHG2, 1.11, 110.0, 109.0, -1) # 3HG2
            pdbatm2(IHD1, CD, CG1, CB, 1.11, 110.0, 180.0, 0) # 1HD1
            pdbatm2(IIHD1, CD, CG1, IHD1, 1.11, 110.0, 109.0, 1) # 2HD1
            pdbatm2(IIIHD1, CD, CG1, IHD1, 1.11, 110.0, 109.0, -1) # 3HD1

    elif resname == 'LEU':
        [CG, CD1, CD2] = atoms['CG'], atoms['CD1'], atoms['CD2']
        pdbatm2(CG, CB, CA, N, 1.54, 109.5, torsions[0], 0) # CG
        pdbatm2(CD1, CG, CB, CA, 1.54, 109.5, torsions[1], 0) # CD1
        pdbatm2(CD2, CG, CB, CD1, 1.54, 109.5, 109.4, -1) # CD2
        if residue.has_hydrogens():
            [HB2, HB3, HG, IHD1, IIHD1] = atoms['HB2'], atoms['HB3'], atoms['HG'], atoms['1HD1'], atoms['2HD1']
            [IIIHD1, IHD2, IIHD2, IIIHD2] = atoms['3HD1'], atoms['1HD2'], atoms['2HD2'], atoms['3HD2']
            pdbatm2(HB2, CB, CA, CG, 1.11, 109.4, 109.4, 1) # HB2
            pdbatm2(HB3, CB, CA, CG, 1.11, 109.4, 109.4, -1) # HB3
            pdbatm2(HG, CG, CB, CD1, 1.11, 109.4, 109.4, 1) # HG
            pdbatm2(IHD1, CD1, CG, CB, 1.11, 109.4, 180.0, 0) # 1HD1
            pdbatm2(IIHD1, CD1, CG, IHD1, 1.11, 109.4, 109.4, 1) # 2HD1
            pdbatm2(IIIHD1, CD1, CG, IHD1, 1.11, 109.4, 109.4, -1) # 3HD1
            pdbatm2(IHD2, CD2, CG, CB, 1.11, 109.4, 180.0, 0) # 1HD2
            pdbatm2(IIHD2, CD2, CG, IHD2, 1.11, 109.4, 109.4, 1) # 2HD2
            pdbatm2(IIIHD2, CD2, CG, IHD2, 1.11, 109.4, 109.4, -1) # 3HD2

    elif resname == 'LYS' or resname == 'LYN':  # LYN is deprotonated LYS
        [CG, CD, CE, NZ] = atoms['CG'], atoms['CD'], atoms['CE'], atoms['NZ']
        pdbatm2(CG, CB, CA, N, 1.54, 109.5, torsions[0], 0) # CG
        pdbatm2(CD, CG, CB, CA, 1.54, 109.5, torsions[1], 0) # CD
        pdbatm2(CE, CD, CG, CB, 1.54, 109.5, torsions[2], 0) # CE
        pdbatm2(NZ, CE, CD, CG, 1.51, 109.5, torsions[3], 0) # NZ
        if residue.has_hydrogens():
            [HB2, HB3, HG2, HG3, HD2] = atoms['HB2'], atoms['HB3'], atoms['HG2'], atoms['HG3'], atoms['HD2']
            [HD3, HE2, HE3, HZ2, HZ3] = atoms['HD3'], atoms['HE2'], atoms['HE3'], atoms['HZ2'], atoms['HZ3']
            pdbatm2(HB2, CB, CA, CG, 1.11, 109.4, 109.4, 1) # HB2
            pdbatm2(HB3, CB, CA, CG, 1.11, 109.4, 109.4, -1) # HB3
            pdbatm2(HG2, CG, CB, CD, 1.11, 109.4, 109.4, 1) # HG2
            pdbatm2(HG3, CG, CB, CD, 1.11, 109.4, 109.4, -1) # HG3
            pdbatm2(HD2, CD, CG, CE, 1.11, 109.4, 109.4, 1) # HD2
            pdbatm2(HD3, CD, CG, CE, 1.11, 109.4, 109.4, -1) # HD3
            pdbatm2(HE2, CE, CD, NZ, 1.11, 109.4, 108.8, 1) # HE2
            pdbatm2(HE3, CE, CD, NZ, 1.11, 109.4, 108.8, -1) # HE3
            if resname == 'LYS':
                HZ1 = atoms['HZ1']
                pdbatm2(HZ1, NZ, CE, CD, 1.02, 109.5, 180.0, 0) # HZ1
                pdbatm2(HZ2, NZ, CE, HZ1, 1.02, 109.5, 109.5, 1) # HZ2
                pdbatm2(HZ3, NZ, CE, HZ1, 1.02, 109.5, 109.5, -1) # HZ3
            else: # LYN case
                pdbatm2(HZ2, NZ, CE, CD, 1.02, 109.5, 60.0, 0) # HZ2
                pdbatm2(HZ3, NZ, CE, CD, 1.02, 109.5, 300.0, 0) # HZ3

    elif resname == 'MET':
        [CG, SD, CE] = atoms['CG'], atoms['SD'], atoms['CE']
        pdbatm2(CG, CB, CA, N, 1.54, 109.5, torsions[0], 0) # CG
        pdbatm2(SD, CG, CB, CA, 1.82, 109.0, torsions[1], 0) # SD
        pdbatm2(CE, SD, CG, CB, 1.82, 96.3, torsions[2], 0) # CE
        if residue.has_hydrogens():
            [HB2, HB3, HG2, HG3, HE1, HE2, HE3] = atoms['HB2'], atoms['HB3'], atoms['HG2'], atoms['HG3'], atoms['HE1'], atoms['HE2'], atoms['HE3']
            pdbatm2(HB2, CB, CA, CG, 1.11, 109.4, 109.4, 1) # HB2
            pdbatm2(HB3, CB, CA, CG, 1.11, 109.4, 109.4, -1) # HB3
            pdbatm2(HG2, CG, CB, SD, 1.11, 109.4, 112.0, 1) # HG2
            pdbatm2(HG3, CG, CB, SD, 1.11, 109.4, 112.0, -1) # HG3
            pdbatm2(HE1, CE, SD, CG, 1.11, 112.0, 180.0, 0) # HE1
            pdbatm2(HE2, CE, SD, HE1, 1.11, 112.0, 109.4, 1) # HE2
            pdbatm2(HE3, CE, SD, HE1, 1.11, 112.0, 109.4, -1) # HE3

    elif resname == 'PHE':
        [CG, CD1, CD2, CE1, CE2, CZ] = atoms['CG'], atoms['CD1'], atoms['CD2'], atoms['CE1'], atoms['CE2'], atoms['CZ']
        pdbatm2(CG, CB, CA, N, 1.50, 109.5, torsions[0], 0) # CG
        pdbatm2(CD1, CG, CB, CA, 1.39, 120.0, torsions[1], 0) # CD1
        pdbatm2(CD2, CG, CB, CD1, 1.39, 120.0, 120.0, 1) # CD2
        pdbatm2(CE1, CD1, CG, CB, 1.39, 120.0, 180.0, 0) # CE1
        pdbatm2(CE2, CD2, CG, CB, 1.39, 120.0, 180.0, 0) # CE2
        pdbatm2(CZ, CE1, CD1, CG, 1.39, 120.0,  0.0, 0) # CZ
        if residue.has_hydrogens():
            [HB2, HB3, HD1, HD2, HE1, HE2, HZ] = atoms['HB2'], atoms['HB3'], atoms['HD1'], atoms['HD2'], atoms['HE1'], atoms['HE2'], atoms['HZ']
            pdbatm2(HB2, CB, CA, CG, 1.11, 109.4, 109.4, 1) # HB2
            pdbatm2(HB3, CB, CA, CG, 1.11, 109.4, 109.4, -1) # HB3
            pdbatm2(HD1, CD1, CG, CE1, 1.10, 120.0, 120.0, 1) # HD1
            pdbatm2(HD2, CD2, CG, CE2, 1.10, 120.0, 120.0, 1) # HD2
            pdbatm2(HE1, CE1, CD1, CZ, 1.10, 120.0, 120.0, 1) # HE1
            pdbatm2(HE2, CE2, CD2, CZ, 1.10, 120.0, 120.0, 1) # HE2
            pdbatm2(HZ, CZ, CE1, CE2, 1.10, 120.0, 120.0, 1) # HZ

    elif resname == 'PRO':
        [CG, CD] = atoms['CG'], atoms['CD']
        pdbatm2(CG, CB, CA, N, .54, 107.0, torsions[0], 0) # CG
        pdbatm2(CD, CG, CB, CA, 1.54, 107.0, torsions[1], 0) # CD
        if residue.has_hydrogens():
            [HB2, HB3, HG2, HG3, HD2, HD3] = atoms['HB2'], atoms['HB3'], atoms['HG2'], atoms['HG3'], atoms['HD2'], atoms['HD3']
            pdbatm2(HB2, CB, CA, CG, 1.11, 109.4, 109.4, 1) # HB2
            pdbatm2(HB3, CB, CA, CG, 1.11, 109.4, 109.4, -1) # HB3
            pdbatm2(HG2, CG, CB, CD, 1.11, 109.4, 109.4, 1) # HG2
            pdbatm2(HG3, CG, CB, CD, 1.11, 109.4, 109.4, -1) # HG3
            pdbatm2(HD2, CD, CG, N, 1.11, 109.4, 109.4, 1) # HD2
            pdbatm2(HD3, CD, CG, N, 1.11, 109.4, 109.4, -1) # HD3

    elif resname == 'SER':
        OG = atoms['OG']
        pdbatm2(OG, CB, CA, N, 1.41, 107.5, torsions[0], 0) # OG
        if residue.has_hydrogens():
            [HB2, HB3, HG] = atoms['HB2'], atoms['HB3'], atoms['HG']
            pdbatm2(HB2, CB, CA, OG, 1.11, 109.4, 106.7, 1) # HB2
            pdbatm2(HB3, CB, CA, OG, 1.11, 109.4, 106.7, -1) # HB3
            pdbatm2(HG, OG, CB, CA, 0.94, 106.9, 180.0, 0) # HG

    elif resname == 'THR':
        [OG1, CG2] = atoms['OG1'], atoms['CG2']
        pdbatm2(OG1, CB, CA, N, 1.41, 107.5, torsions[0], 0) # OG1
        pdbatm2(CG2, CB, CA, OG1, 1.54, 109.5, 107.7, 1) # CG2
        if residue.has_hydrogens():
            [HB, HG1, IHG2, IIHG2, IIIHG2] = atoms['HB'], atoms['HG1'], atoms['1HG2'], atoms['2HG2'], atoms['3HG2']
            pdbatm2(HB, CB, CA, OG1, 1.11, 109.4, 106.7, -1) # HB
            pdbatm2(HG1, OG1, CB, CA, 0.94, 106.9, 180.0, 0) # HG1
            pdbatm2(IHG2, CG2, CB, CA, 1.11, 110.0, 180.0, 0) # 1HG2
            pdbatm2(IIHG2, CG2, CB, IHG2, 1.11, 110.0, 109.0, 1) # 2HG2
            pdbatm2(IIIHG2, CG2, CB, IHG2, 1.11, 110.0, 109.0, -1) # 3HG2

    elif resname == 'TRP':
        [CG, CD1, CD2, NE1, CE2] = atoms['CG'], atoms['CD1'], atoms['CD2'], atoms['NE1'], atoms['CE2']
        [CE3, CZ2, CZ3, CH2] = atoms['CE3'], atoms['CZ2'], atoms['CZ3'], atoms['CH2']
        pdbatm2(CG, CB, CA, N, 1.50, 109.5, torsions[0], 0) # CG
        pdbatm2(CD1, CG, CB, CA, 1.35, 126.0, torsions[1], 0) # CD1
        pdbatm2(CD2, CG, CB, CD1, 1.35, 126.0, 108.0, 1) # CD2
        pdbatm2(NE1, CD1, CG, CD2, 1.35, 108.0,  0.0, 0) # NE1
        pdbatm2(CE2, NE1, CD1, CG, 1.35, 108.0,  0.0, 0) # CE2
        pdbatm2(CE3, CD2, CE2, NE1, 1.35, 120.0, 180.0, 0) # CE3
        pdbatm2(CZ2, CE2, CD2, CE3, 1.35, 120.0,  0.0, 0) # CZ2
        pdbatm2(CZ3, CE3, CD2, CE2, 1.35, 120.0,  0.0, 0) # CZ3
        pdbatm2(CH2, CZ2, CE2, CD2, 1.35, 120.0,  0.0, 0) # CH2
        if residue.has_hydrogens():
            [HB2, HB3, HD1, HE1, HE3] = atoms['HB2'], atoms['HB3'], atoms['HD1'], atoms['HE1'], atoms['HE3']
            [HZ2, HZ3, HH2] = atoms['HZ2'], atoms['HZ3'], atoms['HH2']
            pdbatm2(HB2, CB, CA, CG, 1.11, 109.4, 109.4, 1) # HB2
            pdbatm2(HB3, CB, CA, CG, 1.11, 109.4, 109.4, -1) # HB3
            pdbatm2(HD1, CD1, CG, NE1, 1.10, 126.0, 126.0, 1) # HD1
            pdbatm2(HE1, NE1, CD1, CE2, 1.05, 126.0, 126.0, 1) # HE1
            pdbatm2(HE3, CE3, CD2, CZ3, 1.10, 120.0, 120.0, 1) # HE3
            pdbatm2(HZ2, CZ2, CE2, CH2, 1.10, 120.0, 120.0, 1) # HZ2
            pdbatm2(HZ3, CZ3, CE3, CH2, 1.10, 120.0, 120.0, 1) # HZ3
            pdbatm2(HH2, CH2, CZ2, CZ3, 1.10, 120.0, 120.0, 1) # HH2

    elif resname == 'TYR':
        [CG, CD1, CD2, CE1, CE2, CZ, OH] = atoms['CG'], atoms['CD1'], atoms['CD2'], atoms['CE1'], atoms['CE2'], atoms['CZ'], atoms['OH']
        pdbatm2(CG, CB, CA, N, 1.50, 109.5, torsions[0], 0) # CG
        pdbatm2(CD1, CG, CB, CA, 1.39, 120.0, torsions[1], 0) # CD1
        pdbatm2(CD2, CG, CB, CD1, 1.39, 120.0, 120.0, 1) # CD2
        pdbatm2(CE1, CD1, CG, CB, 1.39, 120.0, 180.0, 0) # CE1
        pdbatm2(CE2, CD2, CG, CB, 1.39, 120.0, 180.0, 0) # CE2
        pdbatm2(CZ, CE1, CD1, CG, 1.39, 120.0,  0.0, 0) # CZ
        pdbatm2(OH, CZ, CE1, CE2, 1.36, 120.0, 120.0, 1) # OH
        if residue.has_hydrogens():
            [HB2, HB3, HD1, HD2, HE1, HE2, HH] = atoms['HB2'], atoms['HB3'], atoms['HD1'], atoms['HD2'], atoms['HE1'], atoms['HE2'], atoms['HH']
            pdbatm2(HB2, CB, CA, CG, 1.11, 109.4, 109.4, 1) # HB2
            pdbatm2(HB3, CB, CA, CG, 1.11, 109.4, 109.4, -1) # HB3
            pdbatm2(HD1, CD1, CG, CE1, 1.10, 120.0, 120.0, 1) # HD1
            pdbatm2(HD2, CD2, CG, CE2, 1.10, 120.0, 120.0, 1) # HD2
            pdbatm2(HE1, CE1, CD1, CZ, 1.10, 120.0, 120.0, 1) # HE1
            pdbatm2(HE2, CE2, CD2, CZ, 1.10, 120.0, 120.0, 1) # HE2
            pdbatm2(HH, OH, CZ, CE1, 0.97, 108.0,  0.0, 0) # HH

    elif resname == 'VAL':
        [CG1, CG2] = atoms['CG1'], atoms['CG2']
        pdbatm2(CG1, CB, CA, N, 1.54, 109.5, torsions[0], 0) # CG1
        pdbatm2(CG2, CB, CA, CG1, 1.54, 109.5, 109.5, -1) # CG2
        if residue.has_hydrogens():
            [HB, IHG1, IIHG1, IIIHG1, IHG2, IIHG2, IIIHG2] = atoms['HB'], atoms['1HG1'], atoms['2HG1'], atoms['3HG1'], atoms['1HG2'], atoms['2HG2'], atoms['3HG2']
            pdbatm2(HB, CB, CA, CG1, 1.11, 109.4, 109.4, 1) # HB
            pdbatm2(IHG1, CG1, CB, CA, 1.11, 109.4, 180.0, 0) # 1HG1
            pdbatm2(IIHG1, CG1, CB, IHG1, 1.11, 109.4, 109.4, 1) # 2HG1
            pdbatm2(IIIHG1, CG1, CB, IHG1, 1.11, 109.4, 109.4, -1) # 3HG1
            pdbatm2(IHG2, CG2, CB, CA, 1.11, 109.4, 180.0, 0) # 1HG2
            pdbatm2(IIHG2, CG2, CB, IHG2, 1.11, 109.4, 109.4, 1) # 2HG2
            pdbatm2(IIIHG2, CG2, CB, IHG2, 1.11, 109.4, 109.4, -1) # 3HG2

    else:
        raise GROWSCESCException("Residue %s is undefined." % resname)
    
    return
