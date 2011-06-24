#include "common.h"

ostream & operator<<(ostream & out, Residue * residue) {
    return out << "<Residue " << residue->name << "-" << residue->id << ">";
}

Residue::Residue(int tmp_id, string & tmp_name) {
    id = tmp_id;
    name = tmp_name;
    chain = NULL;
    atoms = new map <string, Atom *>;
    backbone_atoms = new vector< Atom * >;
    side_chain_atoms = new vector< Atom * >;
    non_clashing_rotamer_ids = new vector< int >;
}

Residue::~Residue() {
    delete atoms;
    delete backbone_atoms;
    delete side_chain_atoms;
    delete non_clashing_rotamer_ids;
}

bool Residue::has_hydrogens() {
    for (int i = 0; i < backbone_atoms->size(); i++)
        if ((*backbone_atoms)[i]->element.compare("H") == 0)
            return true;
    return false;
}

inline bool Residue::is_backbone_atom(Atom * atom) {
    return atom->name.compare("N") == 0 || atom->name.compare("CA") == 0 || atom->name.compare("C") == 0 ||
            atom->name.compare("O") == 0 || atom->name.compare("H") == 0 || atom->name.compare("H1") == 0 ||
            atom->name.compare("H2") == 0 || atom->name.compare("H3") == 0 || atom->name.compare("HA") == 0 ||
            atom->name.compare("HA2") == 0 || atom->name.compare("HA3") == 0 || atom->name.compare("HB") == 0 ||
            atom->name.compare("CB") == 0 || atom->name.compare("HB1") == 0 || atom->name.compare("HB2") == 0 ||
            atom->name.compare("HB3") == 0;
}

void Residue::add_atom(Atom * atom) {
    (*atoms)[atom->name] = atom;
    if (is_backbone_atom(atom))
        backbone_atoms->push_back(atom);
    else
        side_chain_atoms->push_back(atom);
    atom->residue = this;
    return;
}

void Residue::remove_atom_with_name(string & name) {
    map <string, Atom *>::iterator it = atoms->find(name);
    Atom * atom = it->second;

    if (is_backbone_atom(atom)) {
        for (int i = 0; i < backbone_atoms->size(); i++)
            if ((*backbone_atoms)[i]->element.compare(name) == 0)
                backbone_atoms->erase(backbone_atoms->begin() + i);
    } else {
        for (int i = 0; i < side_chain_atoms->size(); i++)
            if ((*side_chain_atoms)[i]->element.compare(name) == 0)
                side_chain_atoms->erase(side_chain_atoms->begin() + i);
    }

    atom->residue = NULL;
    atoms->erase(it);
    return;
}

void test_residue() {
    Residue * residue = sample_residue();
    Atom * atom = sample_atom();
    residue->add_atom(atom);
    cout << residue << endl;
    cout << atom->residue << endl;
    map <string, Atom *> * atoms = (residue->atoms);
    cout << (*atoms)["CA"] << endl;
    Atom * btom = atom, * ctom = atom;
    cout << ctom << endl;
    string resname = residue->name;
    cout << resname << endl;
    return;
}

Residue * sample_residue() {
    string name = "TYR";
    Residue * residue = new Residue(42, name);
    return residue;
}
