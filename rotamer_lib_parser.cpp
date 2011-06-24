#include "common.h"

vector<double *> * get_rotamer_lib_from_file(string & filename) {
    ifstream filestream(filename.c_str());
    ASSERT(filestream.is_open(), "Could not open input rotamer file " + STRING(filename));

    vector< vector< double > > tmp_rotamer_list;
    vector< double * > * rotamer_list = new vector< double * >;
    double val, *coords;
    string line;

    while (getline(filestream, line)) {
        istringstream iss(line);
        vector< double > rotamer;

        while (iss >> val)
            rotamer.push_back(val);
        tmp_rotamer_list.push_back(rotamer);
    }

    int size = tmp_rotamer_list[0].size();

    for (int j = 0; j < tmp_rotamer_list.size(); j++) {
        ASSERT(tmp_rotamer_list[j].size() == size, "Not the same number of torsion angles per line in .lib file.");
        coords = new double [size];
        for (int k = 0; k < size; k++)
            coords[k] = tmp_rotamer_list[j][k];
        rotamer_list->push_back(coords);
    }

    return rotamer_list;
}

void test_rotamer_lib_parser() {
    string filename = "rotamers/val.lib";
    vector< double * > * val_rotamer_list = get_rotamer_lib_from_file(filename);
    for (int i = 0; i < val_rotamer_list->size(); i++)
        cout << (*val_rotamer_list)[i][0] << " " << (*val_rotamer_list)[i][1] << endl;
    return;
}
