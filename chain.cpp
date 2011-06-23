#include "common.h"

ostream & operator<<(ostream & out, Chain * chain) {
    return out << "<Chain id=" << chain->id << ">";
}

Chain::Chain(string & tmp_id) {
    id = tmp_id;
    residues = new map <int, Residue *>;
}

Chain::~Chain() {
    delete residues;
}

void Chain::add_residue(Residue * residue) {
    (*residues)[residue->id] = residue;
    residue->chain = this;
    return;
}

void Chain::remove_residue_with_id(int id) {
    map <int, Residue *>::iterator it = residues->find(id);
    Residue * residue = it->second;
    residue->chain = NULL;
    residues->erase(it);
    return;
}

void test_chain() {
    return;
}

Chain * sample_chain() {
    string id = "A";
    Chain * chain = new Chain(id);
    return chain;
}
