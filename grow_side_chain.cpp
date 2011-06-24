#include "common.h"

void grow_side_chain(Residue * residue, int torsions_index) {
    // build side chain atoms
    cout << "Current residue is " << residue << endl;

    string resname = residue->name;
    map <string, Atom *> * atoms = (residue->atoms);
    vector<double *> * rotamer_lib = (*ROTAMER_LIBRARY)[resname];
    double * torsions = (*rotamer_lib)[torsions_index];

    Atom *N = (*atoms)["N"], *CA = (*atoms)["CA"], *CB;
    if (resname.compare("GLY") != 0)
        CB = (*atoms)["CB"];

    if (resname.compare("ALA") == 0)
        return;

    else if (resname.compare("ARG") == 0) {
        Atom *CG = (*atoms)["CG"], *CD = (*atoms)["CD"], *NE = (*atoms)["NE"], *CZ = (*atoms)["CZ"], *NH1 = (*atoms)["NH1"], *NH2 = (*atoms)["NH2"];
        pdbatm2(CG, CB, CA, N, 1.54, 109.5, torsions[0], 0); // CG
        pdbatm2(CD, CG, CB, CA, 1.54, 109.5, torsions[1], 0); // CD
        pdbatm2(NE, CD, CG, CB, 1.45, 109.5, torsions[2], 0); // NE
        pdbatm2(CZ, NE, CD, CG, 1.35, 120.0, torsions[3], 0); // CZ
        pdbatm2(NH1, CZ, NE, CD, 1.35, 120.0, 180.0, 0); // NH1
        pdbatm2(NH2, CZ, NE, NH1, 1.35, 120.0, 120.0, 1); // NH2
        if (residue->has_hydrogens()) {
            Atom *HB2 = (*atoms)["HB2"], *HB3 = (*atoms)["HB3"], *HG2 = (*atoms)["HG2"], *HG3 = (*atoms)["HG3"], *HD2 = (*atoms)["HD2"], *HD3 = (*atoms)["HD3"];
            Atom *HE = (*atoms)["HE"], *IHH1 = (*atoms)["1HH1"], *IIHH1 = (*atoms)["2HH1"], *IHH2 = (*atoms)["1HH2"], *IIHH2 = (*atoms)["2HH2"];
            pdbatm2(HB2, CB, CA, CG, 1.11, 109.4, 109.4, 1); // HB2
            pdbatm2(HB3, CB, CA, CG, 1.11, 109.4, 109.4, -1); // HB3
            pdbatm2(HG2, CG, CB, CD, 1.11, 109.4, 109.4, 1); // HG2
            pdbatm2(HG3, CG, CB, CD, 1.11, 109.4, 109.4, -1); // HG3
            pdbatm2(HD2, CD, CG, NE, 1.11, 109.4, 109.4, 1); // HD2
            pdbatm2(HD3, CD, CG, NE, 1.11, 109.4, 109.4, -1); // HD3
            pdbatm2(HE, NE, CD, CZ, 1.02, 120.0, 120.0, 1); // HE
            pdbatm2(IHH1, NH1, CZ, NE, 1.02, 120.0, 180.0, 0); // 1HH1
            pdbatm2(IIHH1, NH1, CZ, IHH1, 1.02, 120.0, 120.0, 1); // 2HH1
            pdbatm2(IHH2, NH2, CZ, NE, 1.02, 120.0, 180.0, 0); // 1HH2
            pdbatm2(IIHH2, NH2, CZ, IHH2, 1.02, 120.0, 120.0, 1); // 2HH2
        }
    } else if (resname.compare("ASN") == 0) {
        Atom *CG = (*atoms)["CG"], *OD1 = (*atoms)["OD1"], *ND2 = (*atoms)["ND2"];
        pdbatm2(CG, CB, CA, N, 1.51, 107.8, torsions[0], 0); // CG
        pdbatm2(OD1, CG, CB, CA, 1.22, 122.5, torsions[1], 0); // OD1
        pdbatm2(ND2, CG, CB, OD1, 1.34, 112.7, 124.0, 1); // ND2
        if (residue->has_hydrogens()) {
            Atom *HB2 = (*atoms)["HB2"], *HB3 = (*atoms)["HB3"], *IHD2 = (*atoms)["1HD2"], *IIHD2 = (*atoms)["2HD2"];
            pdbatm2(HB2, CB, CA, CG, 1.11, 109.4, 107.9, 1); // HB2
            pdbatm2(HB3, CB, CA, CG, 1.11, 109.4, 107.9, -1); // HB3
            pdbatm2(IHD2, ND2, CG, CB, 1.02, 119.0, 0.0, 0); // 1HD2
            pdbatm2(IIHD2, ND2, CG, IHD2, 1.02, 119.0, 120.0, 1); // 2HD2
        }
    } else if (resname.compare("ASP") == 0 or resname.compare("ASH") == 0) { // ASH is protonated ASN
        Atom *CG = (*atoms)["CG"], *OD1 = (*atoms)["OD1"], *OD2 = (*atoms)["OD2"];
        pdbatm2(CG, CB, CA, N, 1.51, 107.8, torsions[0], 0); // CG
        pdbatm2(OD1, CG, CB, CA, 1.25, 117.0, torsions[1], 0); // OD1
        pdbatm2(OD2, CG, CB, OD1, 1.25, 117.0, 126.0, 1); // OD2
        if (residue->has_hydrogens()) {
            Atom *HB2 = (*atoms)["HB2"], *HB3 = (*atoms)["HB3"];
            pdbatm2(HB2, CB, CA, CG, 1.11, 109.4, 107.9, 1); // HB2
            pdbatm2(HB3, CB, CA, CG, 1.11, 109.4, 107.9, -1); // HB3
            if (resname.compare("ASH") == 0)
                pdbatm2((*atoms)["HD2"], OD2, CG, CB, 0.96, 109.5, 180.0, 0); // HD2
        }
    } else if (resname.compare("CYS") == 0 or resname.compare("CYX") == 0) {// CYX is deprotonated CYS
        Atom *SG = (*atoms)["SG"];
        pdbatm2(SG, CB, CA, N, 1.82, 109.0, torsions[0], 0); // SG
        if (residue->has_hydrogens()) {
            Atom *HB2 = (*atoms)["HB2"], *HB3 = (*atoms)["HB3"];
            pdbatm2(HB2, CB, CA, SG, 1.11, 109.4, 112.0, 1); // HB2
            pdbatm2(HB3, CB, CA, SG, 1.11, 109.4, 112.0, -1); // HB3
            if (resname.compare("CYS") == 0)
                pdbatm2((*atoms)["HG"], SG, CB, CA, 1.34, 96.0, 180.0, 0); // HG
        }
    } else if (resname.compare("GLN") == 0) {
        Atom *CG = (*atoms)["CG"], *CD = (*atoms)["CD"], *OE1 = (*atoms)["OE1"], *NE2 = (*atoms)["NE2"];
        pdbatm2(CG, CB, CA, N, 1.54, 109.5, torsions[0], 0); // CG
        pdbatm2(CD, CG, CB, CA, 1.51, 107.8, torsions[1], 0); // CD
        pdbatm2(OE1, CD, CG, CB, 1.22, 122.5, torsions[2], 0); // OE1
        pdbatm2(NE2, CD, CG, OE1, 1.34, 112.7, 124.0, 1); // NE2
        if (residue->has_hydrogens()) {
            Atom *HB2 = (*atoms)["HB2"], *HB3 = (*atoms)["HB3"], *HG2 = (*atoms)["HG2"], *HG3 = (*atoms)["HG3"], *IHE2 = (*atoms)["1HE2"], *IIHE2 = (*atoms)["2HE2"];
            pdbatm2(HB2, CB, CA, CG, 1.11, 109.4, 109.4, 1); // HB2
            pdbatm2(HB3, CB, CA, CG, 1.11, 109.4, 109.4, -1); // HB3
            pdbatm2(HG2, CG, CB, CD, 1.11, 109.4, 107.9, 1); // HG2
            pdbatm2(HG3, CG, CB, CD, 1.11, 109.4, 107.9, -1); // HG3
            pdbatm2(IHE2, NE2, CD, CG, 1.02, 119.0, 0.0, 0); // 1HE2
            pdbatm2(IIHE2, NE2, CD, IHE2, 1.02, 119.0, 120.0, 1); // 2HE2
        }
    } else if (resname.compare("GLU") == 0 or resname.compare("GLH") == 0) { // GLH is protonated GLU
        Atom *CG = (*atoms)["CG"], *CD = (*atoms)["CD"], *OE1 = (*atoms)["OE1"], *OE2 = (*atoms)["OE2"];
        pdbatm2(CG, CB, CA, N, 1.54, 109.5, torsions[0], 0); // CG
        pdbatm2(CD, CG, CB, CA, 1.51, 107.8, torsions[1], 0); // CD
        pdbatm2(OE1, CD, CG, CB, 1.25, 117.0, torsions[2], 0); // OE1
        pdbatm2(OE2, CD, CG, OE1, 1.25, 117.0, 126.0, 1); // OE2
        if (residue->has_hydrogens()) {
            Atom *HB2 = (*atoms)["HB2"], *HB3 = (*atoms)["HB3"], *HG2 = (*atoms)["HG2"], *HG3 = (*atoms)["HG3"];
            pdbatm2(HB2, CB, CA, CG, 1.11, 109.4, 109.4, 1); // HB2
            pdbatm2(HB3, CB, CA, CG, 1.11, 109.4, 109.4, -1); // HB3
            pdbatm2(HG2, CG, CB, CD, 1.11, 109.4, 107.9, 1); // HG2
            pdbatm2(HG3, CG, CB, CD, 1.11, 109.4, 107.9, -1); // HG3
            if (resname.compare("GLH") == 0)
                pdbatm2((*atoms)["HE2"], OE2, CD, CG, 0.96, 109.5, 180.0, 0); // HE2
        }
    } else if (resname.compare("GLY") == 0)
        return;

    else if (resname.compare("HIS") == 0 or resname.compare("HID") == 0 or resname.compare("HIE") == 0) {
        Atom *CG = (*atoms)["CG"], *ND1 = (*atoms)["ND1"], *CD2 = (*atoms)["CD2"], *CE1 = (*atoms)["CE1"], *NE2 = (*atoms)["NE2"];
        pdbatm2(CG, CB, CA, N, 1.50, 109.5, torsions[0], 0); // CG
        pdbatm2(ND1, CG, CB, CA, 1.35, 126.0, torsions[1], 0); // ND1
        pdbatm2(CD2, CG, CB, ND1, 1.35, 126.0, 108.0, 1); // CD2
        pdbatm2(CE1, ND1, CG, CD2, 1.35, 108.0, 0.0, 0); // CE1
        pdbatm2(NE2, CD2, CG, ND1, 1.35, 108.0, 0.0, 0); // NE2
        if (residue->has_hydrogens()) {
            pdbatm2((*atoms)["HB2"], CB, CA, CG, 1.11, 109.4, 109.4, 1); // HB2
            pdbatm2((*atoms)["HB3"], CB, CA, CG, 1.11, 109.4, 109.4, -1); // HB3
            if (resname.compare("HID") == 0 or resname.compare("HIE") == 0) {
                pdbatm2((*atoms)["HD2"], CD2, CG, NE2, 1.10, 126.0, 126.0, 1); // HD2
                pdbatm2((*atoms)["HE1"], CE1, ND1, NE2, 1.10, 126.0, 126.0, 1); // HE1
            }
            if (resname.compare("HIS") == 0 or resname.compare("HID") == 0)
                pdbatm2((*atoms)["HD1"], ND1, CG, CD2, 1.02, 126.0, 0.0, 0); // HD1
            if (resname.compare("HIS") == 0 or resname.compare("HIE") == 0)
                pdbatm2((*atoms)["HE2"], NE2, CD2, CE1, 1.02, 126.0, 126.0, 1); // HE2
        }
    } else if (resname.compare("ILE") == 0) {
        Atom *CG1 = (*atoms)["CG1"], *CG2 = (*atoms)["CG2"], *CD1 = (*atoms)["CD1"];
        pdbatm2(CG1, CB, CA, N, 1.54, 109.5, torsions[0], 0); // CG1
        pdbatm2(CG2, CB, CA, CG1, 1.54, 109.5, 109.5, 1); // CG2
        pdbatm2(CD1, CG1, CB, CA, 1.54, 109.5, torsions[1], 0); // CD1
        if (residue->has_hydrogens()) {
            Atom *HB = (*atoms)["HB"], *IIHG1 = (*atoms)["2HG1"], *IIIHG1 = (*atoms)["3HG1"], *IHG2 = (*atoms)["1HG2"], *IIHG2 = (*atoms)["2HG2"];
            Atom *IIIHG2 = (*atoms)["3HG2"], *IHD1 = (*atoms)["1HD1"], *IIHD1 = (*atoms)["2HD1"], *IIIHD1 = (*atoms)["3HD1"];
            pdbatm2(HB, CB, CA, CG1, 1.11, 109.4, 109.4, -1); // HB
            pdbatm2(IIHG1, CG1, CB, CD1, 1.11, 109.4, 109.4, 1); // 2HG1
            pdbatm2(IIIHG1, CG1, CB, CD1, 1.11, 109.4, 109.4, -1); // 3HG1
            pdbatm2(IHG2, CG2, CB, CA, 1.11, 110.0, 180.0, 0); // 1HG2
            pdbatm2(IIHG2, CG2, CB, IHG2, 1.11, 110.0, 109.0, 1); // 2HG2
            pdbatm2(IIIHG2, CG2, CB, IHG2, 1.11, 110.0, 109.0, -1); // 3HG2
            pdbatm2(IHD1, CD1, CG1, CB, 1.11, 110.0, 180.0, 0); // 1HD1
            pdbatm2(IIHD1, CD1, CG1, IHD1, 1.11, 110.0, 109.0, 1); // 2HD1
            pdbatm2(IIIHD1, CD1, CG1, IHD1, 1.11, 110.0, 109.0, -1); // 3HD1
        }
    } else if (resname.compare("LEU") == 0) {
        Atom *CG = (*atoms)["CG"], *CD1 = (*atoms)["CD1"], *CD2 = (*atoms)["CD2"];
        pdbatm2(CG, CB, CA, N, 1.54, 109.5, torsions[0], 0); // CG
        pdbatm2(CD1, CG, CB, CA, 1.54, 109.5, torsions[1], 0); // CD1
        pdbatm2(CD2, CG, CB, CD1, 1.54, 109.5, 109.4, -1); // CD2
        if (residue->has_hydrogens()) {
            Atom *HB2 = (*atoms)["HB2"], *HB3 = (*atoms)["HB3"], *HG = (*atoms)["HG"], *IHD1 = (*atoms)["1HD1"], *IIHD1 = (*atoms)["2HD1"];
            Atom *IIIHD1 = (*atoms)["3HD1"], *IHD2 = (*atoms)["1HD2"], *IIHD2 = (*atoms)["2HD2"], *IIIHD2 = (*atoms)["3HD2"];
            pdbatm2(HB2, CB, CA, CG, 1.11, 109.4, 109.4, 1); // HB2
            pdbatm2(HB3, CB, CA, CG, 1.11, 109.4, 109.4, -1); // HB3
            pdbatm2(HG, CG, CB, CD1, 1.11, 109.4, 109.4, 1); // HG
            pdbatm2(IHD1, CD1, CG, CB, 1.11, 109.4, 180.0, 0); // 1HD1
            pdbatm2(IIHD1, CD1, CG, IHD1, 1.11, 109.4, 109.4, 1); // 2HD1
            pdbatm2(IIIHD1, CD1, CG, IHD1, 1.11, 109.4, 109.4, -1); // 3HD1
            pdbatm2(IHD2, CD2, CG, CB, 1.11, 109.4, 180.0, 0); // 1HD2
            pdbatm2(IIHD2, CD2, CG, IHD2, 1.11, 109.4, 109.4, 1); // 2HD2
            pdbatm2(IIIHD2, CD2, CG, IHD2, 1.11, 109.4, 109.4, -1); // 3HD2
        }
    } else if (resname.compare("LYS") == 0 or resname.compare("LYN") == 0) { // LYN is deprotonated LYS
        Atom *CG = (*atoms)["CG"], *CD = (*atoms)["CD"], *CE = (*atoms)["CE"], *NZ = (*atoms)["NZ"];
        pdbatm2(CG, CB, CA, N, 1.54, 109.5, torsions[0], 0); // CG
        pdbatm2(CD, CG, CB, CA, 1.54, 109.5, torsions[1], 0); // CD
        pdbatm2(CE, CD, CG, CB, 1.54, 109.5, torsions[2], 0); // CE
        pdbatm2(NZ, CE, CD, CG, 1.51, 109.5, torsions[3], 0); // NZ
        if (residue->has_hydrogens()) {
            Atom *HB2 = (*atoms)["HB2"], *HB3 = (*atoms)["HB3"], *HG2 = (*atoms)["HG2"], *HG3 = (*atoms)["HG3"], *HD2 = (*atoms)["HD2"];
            Atom *HD3 = (*atoms)["HD3"], *HE2 = (*atoms)["HE2"], *HE3 = (*atoms)["HE3"], *HZ2 = (*atoms)["HZ2"], *HZ3 = (*atoms)["HZ3"];
            pdbatm2(HB2, CB, CA, CG, 1.11, 109.4, 109.4, 1); // HB2
            pdbatm2(HB3, CB, CA, CG, 1.11, 109.4, 109.4, -1); // HB3
            pdbatm2(HG2, CG, CB, CD, 1.11, 109.4, 109.4, 1); // HG2
            pdbatm2(HG3, CG, CB, CD, 1.11, 109.4, 109.4, -1); // HG3
            pdbatm2(HD2, CD, CG, CE, 1.11, 109.4, 109.4, 1); // HD2
            pdbatm2(HD3, CD, CG, CE, 1.11, 109.4, 109.4, -1); // HD3
            pdbatm2(HE2, CE, CD, NZ, 1.11, 109.4, 108.8, 1); // HE2
            pdbatm2(HE3, CE, CD, NZ, 1.11, 109.4, 108.8, -1); // HE3
            if (resname.compare("LYS") == 0) {
                Atom *HZ1 = (*atoms)["HZ1"];
                pdbatm2(HZ1, NZ, CE, CD, 1.02, 109.5, 180.0, 0); // HZ1
                pdbatm2(HZ2, NZ, CE, HZ1, 1.02, 109.5, 109.5, 1); // HZ2
                pdbatm2(HZ3, NZ, CE, HZ1, 1.02, 109.5, 109.5, -1); // HZ3
            } else { // LYN case
                pdbatm2(HZ2, NZ, CE, CD, 1.02, 109.5, 60.0, 0); // HZ2
                pdbatm2(HZ3, NZ, CE, CD, 1.02, 109.5, 300.0, 0); // HZ3
            }
        }
    } else if (resname.compare("MET") == 0) {
        Atom *CG = (*atoms)["CG"], *SD = (*atoms)["SD"], *CE = (*atoms)["CE"];
        pdbatm2(CG, CB, CA, N, 1.54, 109.5, torsions[0], 0); // CG
        pdbatm2(SD, CG, CB, CA, 1.82, 109.0, torsions[1], 0); // SD
        pdbatm2(CE, SD, CG, CB, 1.82, 96.3, torsions[2], 0); // CE
        if (residue->has_hydrogens()) {
            Atom *HB2 = (*atoms)["HB2"], *HB3 = (*atoms)["HB3"], *HG2 = (*atoms)["HG2"], *HG3 = (*atoms)["HG3"], *HE1 = (*atoms)["HE1"], *HE2 = (*atoms)["HE2"], *HE3 = (*atoms)["HE3"];
            pdbatm2(HB2, CB, CA, CG, 1.11, 109.4, 109.4, 1); // HB2
            pdbatm2(HB3, CB, CA, CG, 1.11, 109.4, 109.4, -1); // HB3
            pdbatm2(HG2, CG, CB, SD, 1.11, 109.4, 112.0, 1); // HG2
            pdbatm2(HG3, CG, CB, SD, 1.11, 109.4, 112.0, -1); // HG3
            pdbatm2(HE1, CE, SD, CG, 1.11, 112.0, 180.0, 0); // HE1
            pdbatm2(HE2, CE, SD, HE1, 1.11, 112.0, 109.4, 1); // HE2
            pdbatm2(HE3, CE, SD, HE1, 1.11, 112.0, 109.4, -1); // 
        }
    } else if (resname.compare("PHE") == 0) {
        Atom *CG = (*atoms)["CG"], *CD1 = (*atoms)["CD1"], *CD2 = (*atoms)["CD2"], *CE1 = (*atoms)["CE1"], *CE2 = (*atoms)["CE2"], *CZ = (*atoms)["CZ"];
        pdbatm2(CG, CB, CA, N, 1.50, 109.5, torsions[0], 0); // CG
        pdbatm2(CD1, CG, CB, CA, 1.39, 120.0, torsions[1], 0); // CD1
        pdbatm2(CD2, CG, CB, CD1, 1.39, 120.0, 120.0, 1); // CD2
        pdbatm2(CE1, CD1, CG, CB, 1.39, 120.0, 180.0, 0); // CE1
        pdbatm2(CE2, CD2, CG, CB, 1.39, 120.0, 180.0, 0); // CE2
        pdbatm2(CZ, CE1, CD1, CG, 1.39, 120.0, 0.0, 0); // CZ
        if (residue->has_hydrogens()) {
            Atom *HB2 = (*atoms)["HB2"], *HB3 = (*atoms)["HB3"], *HD1 = (*atoms)["HD1"], *HD2 = (*atoms)["HD2"], *HE1 = (*atoms)["HE1"], *HE2 = (*atoms)["HE2"], *HZ = (*atoms)["HZ"];
            pdbatm2(HB2, CB, CA, CG, 1.11, 109.4, 109.4, 1); // HB2
            pdbatm2(HB3, CB, CA, CG, 1.11, 109.4, 109.4, -1); // HB3
            pdbatm2(HD1, CD1, CG, CE1, 1.10, 120.0, 120.0, 1); // HD1
            pdbatm2(HD2, CD2, CG, CE2, 1.10, 120.0, 120.0, 1); // HD2
            pdbatm2(HE1, CE1, CD1, CZ, 1.10, 120.0, 120.0, 1); // HE1
            pdbatm2(HE2, CE2, CD2, CZ, 1.10, 120.0, 120.0, 1); // HE2
            pdbatm2(HZ, CZ, CE1, CE2, 1.10, 120.0, 120.0, 1); // HZ
        }
    } else if (resname.compare("PRO") == 0) {
        Atom *CG = (*atoms)["CG"], *CD = (*atoms)["CD"];
        pdbatm2(CG, CB, CA, N, 1.54, 107.0, torsions[0], 0); // CG
        pdbatm2(CD, CG, CB, CA, 1.54, 107.0, torsions[1], 0); // CD
        if (residue->has_hydrogens()) {
            Atom *HB2 = (*atoms)["HB2"], *HB3 = (*atoms)["HB3"], *HG2 = (*atoms)["HG2"], *HG3 = (*atoms)["HG3"], *HD2 = (*atoms)["HD2"], *HD3 = (*atoms)["HD3"];
            pdbatm2(HB2, CB, CA, CG, 1.11, 109.4, 109.4, 1); // HB2
            pdbatm2(HB3, CB, CA, CG, 1.11, 109.4, 109.4, -1); // HB3
            pdbatm2(HG2, CG, CB, CD, 1.11, 109.4, 109.4, 1); // HG2
            pdbatm2(HG3, CG, CB, CD, 1.11, 109.4, 109.4, -1); // HG3
            pdbatm2(HD2, CD, CG, N, 1.11, 109.4, 109.4, 1); // HD2
            pdbatm2(HD3, CD, CG, N, 1.11, 109.4, 109.4, -1); // HD3
        }
    } else if (resname.compare("SER") == 0) {
        Atom *OG = (*atoms)["OG"];
        pdbatm2(OG, CB, CA, N, 1.41, 107.5, torsions[0], 0); // OG
        if (residue->has_hydrogens()) {
            Atom *HB2 = (*atoms)["HB2"], *HB3 = (*atoms)["HB3"], *HG = (*atoms)["HG"];
            pdbatm2(HB2, CB, CA, OG, 1.11, 109.4, 106.7, 1); // HB2
            pdbatm2(HB3, CB, CA, OG, 1.11, 109.4, 106.7, -1); // HB3
            pdbatm2(HG, OG, CB, CA, 0.94, 106.9, 180.0, 0); // HG
        }
    } else if (resname.compare("THR") == 0) {
        Atom *OG1 = (*atoms)["OG1"], *CG2 = (*atoms)["CG2"];
        pdbatm2(OG1, CB, CA, N, 1.41, 107.5, torsions[0], 0); // OG1
        pdbatm2(CG2, CB, CA, OG1, 1.54, 109.5, 107.7, 1); // CG2
        if (residue->has_hydrogens()) {
            Atom *HB = (*atoms)["HB"], *HG1 = (*atoms)["HG1"], *IHG2 = (*atoms)["1HG2"], *IIHG2 = (*atoms)["2HG2"], *IIIHG2 = (*atoms)["3HG2"];
            pdbatm2(HB, CB, CA, OG1, 1.11, 109.4, 106.7, -1); // HB
            pdbatm2(HG1, OG1, CB, CA, 0.94, 106.9, 180.0, 0); // HG1
            pdbatm2(IHG2, CG2, CB, CA, 1.11, 110.0, 180.0, 0); // 1HG2
            pdbatm2(IIHG2, CG2, CB, IHG2, 1.11, 110.0, 109.0, 1); // 2HG2
            pdbatm2(IIIHG2, CG2, CB, IHG2, 1.11, 110.0, 109.0, -1); // 3HG2
        }
    } else if (resname.compare("TRP") == 0) {
        Atom *CG = (*atoms)["CG"], *CD1 = (*atoms)["CD1"], *CD2 = (*atoms)["CD2"], *NE1 = (*atoms)["NE1"], *CE2 = (*atoms)["CE2"];
        Atom *CE3 = (*atoms)["CE3"], *CZ2 = (*atoms)["CZ2"], *CZ3 = (*atoms)["CZ3"], *CH2 = (*atoms)["CH2"];
        pdbatm2(CG, CB, CA, N, 1.50, 109.5, torsions[0], 0); // CG
        pdbatm2(CD1, CG, CB, CA, 1.35, 126.0, torsions[1], 0); // CD1
        pdbatm2(CD2, CG, CB, CD1, 1.35, 126.0, 108.0, 1); // CD2
        pdbatm2(NE1, CD1, CG, CD2, 1.35, 108.0, 0.0, 0); // NE1
        pdbatm2(CE2, NE1, CD1, CG, 1.35, 108.0, 0.0, 0); // CE2
        pdbatm2(CE3, CD2, CE2, NE1, 1.35, 120.0, 180.0, 0); // CE3
        pdbatm2(CZ2, CE2, CD2, CE3, 1.35, 120.0, 0.0, 0); // CZ2
        pdbatm2(CZ3, CE3, CD2, CE2, 1.35, 120.0, 0.0, 0); // CZ3
        pdbatm2(CH2, CZ2, CE2, CD2, 1.35, 120.0, 0.0, 0); // CH2
        if (residue->has_hydrogens()) {
            Atom *HB2 = (*atoms)["HB2"], *HB3 = (*atoms)["HB3"], *HD1 = (*atoms)["HD1"], *HE1 = (*atoms)["HE1"], *HE3 = (*atoms)["HE3"];
            Atom *HZ2 = (*atoms)["HZ2"], *HZ3 = (*atoms)["HZ3"], *HH2 = (*atoms)["HH2"];
            pdbatm2(HB2, CB, CA, CG, 1.11, 109.4, 109.4, 1); // HB2
            pdbatm2(HB3, CB, CA, CG, 1.11, 109.4, 109.4, -1); // HB3
            pdbatm2(HD1, CD1, CG, NE1, 1.10, 126.0, 126.0, 1); // HD1
            pdbatm2(HE1, NE1, CD1, CE2, 1.05, 126.0, 126.0, 1); // HE1
            pdbatm2(HE3, CE3, CD2, CZ3, 1.10, 120.0, 120.0, 1); // HE3
            pdbatm2(HZ2, CZ2, CE2, CH2, 1.10, 120.0, 120.0, 1); // HZ2
            pdbatm2(HZ3, CZ3, CE3, CH2, 1.10, 120.0, 120.0, 1); // HZ3
            pdbatm2(HH2, CH2, CZ2, CZ3, 1.10, 120.0, 120.0, 1); // HH2
        }
    } else if (resname.compare("TYR") == 0) {
        Atom *CG = (*atoms)["CG"], *CD1 = (*atoms)["CD1"], *CD2 = (*atoms)["CD2"], *CE1 = (*atoms)["CE1"], *CE2 = (*atoms)["CE2"], *CZ = (*atoms)["CZ"], *OH = (*atoms)["OH"];
        pdbatm2(CG, CB, CA, N, 1.50, 109.5, torsions[0], 0); // CG
        pdbatm2(CD1, CG, CB, CA, 1.39, 120.0, torsions[1], 0); // CD1
        pdbatm2(CD2, CG, CB, CD1, 1.39, 120.0, 120.0, 1); // CD2
        pdbatm2(CE1, CD1, CG, CB, 1.39, 120.0, 180.0, 0); // CE1
        pdbatm2(CE2, CD2, CG, CB, 1.39, 120.0, 180.0, 0); // CE2
        pdbatm2(CZ, CE1, CD1, CG, 1.39, 120.0, 0.0, 0); // CZ
        pdbatm2(OH, CZ, CE1, CE2, 1.36, 120.0, 120.0, 1); // OH
        if (residue->has_hydrogens()) {
            Atom *HB2 = (*atoms)["HB2"], *HB3 = (*atoms)["HB3"], *HD1 = (*atoms)["HD1"], *HD2 = (*atoms)["HD2"], *HE1 = (*atoms)["HE1"], *HE2 = (*atoms)["HE2"], *HH = (*atoms)["HH"];
            pdbatm2(HB2, CB, CA, CG, 1.11, 109.4, 109.4, 1); // HB2
            pdbatm2(HB3, CB, CA, CG, 1.11, 109.4, 109.4, -1); // HB3
            pdbatm2(HD1, CD1, CG, CE1, 1.10, 120.0, 120.0, 1); // HD1
            pdbatm2(HD2, CD2, CG, CE2, 1.10, 120.0, 120.0, 1); // HD2
            pdbatm2(HE1, CE1, CD1, CZ, 1.10, 120.0, 120.0, 1); // HE1
            pdbatm2(HE2, CE2, CD2, CZ, 1.10, 120.0, 120.0, 1); // HE2
            pdbatm2(HH, OH, CZ, CE1, 0.97, 108.0, 0.0, 0); // HH
        }
    } else if (resname.compare("VAL") == 0) {
        Atom *CG1 = (*atoms)["CG1"], *CG2 = (*atoms)["CG2"];
        pdbatm2(CG1, CB, CA, N, 1.54, 109.5, torsions[0], 0); // CG1
        pdbatm2(CG2, CB, CA, CG1, 1.54, 109.5, 109.5, -1); // CG2
        if (residue->has_hydrogens()) {
            Atom *HB = (*atoms)["HB"], *IHG1 = (*atoms)["1HG1"], *IIHG1 = (*atoms)["2HG1"], *IIIHG1 = (*atoms)["3HG1"], *IHG2 = (*atoms)["1HG2"], *IIHG2 = (*atoms)["2HG2"], *IIIHG2 = (*atoms)["3HG2"];
            pdbatm2(HB, CB, CA, CG1, 1.11, 109.4, 109.4, 1); // HB
            pdbatm2(IHG1, CG1, CB, CA, 1.11, 109.4, 180.0, 0); // 1HG1
            pdbatm2(IIHG1, CG1, CB, IHG1, 1.11, 109.4, 109.4, 1); // 2HG1
            pdbatm2(IIIHG1, CG1, CB, IHG1, 1.11, 109.4, 109.4, -1); // 3HG1
            pdbatm2(IHG2, CG2, CB, CA, 1.11, 109.4, 180.0, 0); // 1HG2
            pdbatm2(IIHG2, CG2, CB, IHG2, 1.11, 109.4, 109.4, 1); // 2HG2
            pdbatm2(IIIHG2, CG2, CB, IHG2, 1.11, 109.4, 109.4, -1); // 3HG2
        }
    } else
        ASSERT(false, "Residue " + STRING(resname) + " is undefined.\n");

    return;
}
