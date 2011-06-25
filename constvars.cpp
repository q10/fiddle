#include "common.h"

const double VDW_ALPHA = 0.8;
map <string, double> * VDW_RADIUS = new map <string, double>;
map <string, vector<double *> *> * ROTAMER_LIBRARY = new map <string, vector<double *> *>;

void initialize_constants() {

    (*VDW_RADIUS)["H"] = 1.20;
    (*VDW_RADIUS)["C"] = 1.70;
    (*VDW_RADIUS)["N"] = 1.55;
    (*VDW_RADIUS)["O"] = 1.52;
    (*VDW_RADIUS)["S"] = 1.80;

    initialize_rotamer_library();

    return;
}

void initialize_rotamer_library() {
    string filename, residue, almost_all_residues [] = {"arg", "asn", "asp", "cys", "gln", "glu", "his", "ile", "leu", "lys", "met", "phe", "pro", "ser", "thr", "trp", "tyr", "val"};
    (*ROTAMER_LIBRARY)["GLY"] = new vector<double *>;
    (*ROTAMER_LIBRARY)["ALA"] = new vector<double *>;
    for (int i = 0; i < 18; i++) {
        residue = almost_all_residues[i];
        filename = "rotamers/" + residue + ".lib";
        transform(residue.begin(), residue.end(), residue.begin(), ::toupper); // convert to upper case
        (*ROTAMER_LIBRARY)[residue] = get_rotamer_lib_from_file(filename);
    }
    return;
}

void test_rotamer_library_creation() {
    initialize_constants();
    vector<double *> * ile_rotamer_lib = (*ROTAMER_LIBRARY)["ILE"];
    for (int i = 0; i < ile_rotamer_lib->size(); i++)
        cout << (*ile_rotamer_lib)[i][0] << " " << (*ile_rotamer_lib)[i][1] << endl;
    return;
}
