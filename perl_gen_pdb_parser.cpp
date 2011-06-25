#include "common.h"

Chain * get_chain_from_perl_gen_pdb_file(string & filename) {
    ifstream filestream(filename.c_str());
    ASSERT(filestream.is_open(), "Could not open input PDB file " + STRING(filename));

    Chain * chain = NULL;
    Residue * current_residue = NULL;
    Atom * current_atom = NULL;

    string line, line_key, atom_name, resname, chain_id, element;
    int atom_id, residue_id;
    double * coords = new double [3];

    while (getline(filestream, line)) {
        element = line[13];
        istringstream iss(line);
        iss >> line_key; // get ATOM
        iss >> atom_id; // get atom id
        iss >> atom_name; // get atom name (ex CA)
        iss >> resname; // get residue name (ex HIS)
        iss >> chain_id; // get chain id
        iss >> residue_id; // get residue id
        iss >> coords[0];
        iss >> coords[1];
        iss >> coords[2];

        if (chain == NULL)
            chain = new Chain(chain_id);

        if (current_residue == NULL or current_residue->id != residue_id) {
            current_residue = new Residue(residue_id, resname);
            chain->add_residue(current_residue);
        }

        current_atom = new Atom(atom_id, atom_name, coords, element);
        current_residue->add_atom(current_atom);
    }
    delete coords;
    return chain;
}

void test_get_chain_from_perl_gen_pdb_file() {
    string filename = "pdbs/sample.pdb";
    Chain * chain = get_chain_from_perl_gen_pdb_file(filename);
    Residue * residue = (*chain->residues)[1];

    map <string, Atom *> * atoms = residue->atoms;

    Atom * atom = (*residue->atoms)["CG"];
    Atom * btom = (*atoms)["N"];

    cout << chain << endl << residue << endl << atom
            << " {" << atom->coords[0] << " " << atom->coords[1] << " " << atom->coords[2] << "} "
            << "Element: " << atom->element << endl;
    cout << btom << " {" << btom->coords[0] << " " << btom->coords[1] << " " << btom->coords[2] << "} "
            << "Element: " << btom->element << endl;
    return;
}
