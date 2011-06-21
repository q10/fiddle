#include "common.h"

Atom::Atom(int tmp_id, string tmp_name, double * tmp_coords, string tmp_element) {
    id = tmp_id;
    name = tmp_name;
    element = tmp_element;
    coords = new double[3];
    for (int i = 0; i < 3; i++)
        coords[i] = tmp_coords[i];
}

Atom::~Atom() {
    delete [] coords;
}

ostream & Atom::operator<<(ostream & out, const Atom& atom) {
    out << "<Atom " << atom.name << ">";
    return out;
}

void Atom::set_coords(double * new_coords) {
    for (int i = 0; i < 3; i++)
        coords[i] = new_coords[i];
    return;
}

void Atom::set_residue(Residue * tmp_residue) {
    residue = tmp_residue;
    return;
}

double Atom::distance_from(Atom * other_atom) {
    double * other_coords = other_atom.coords;
    double dx = coords[0] - other_coords[0];
    double dy = coords[1] - other_coords[1];
    double dz = coords[2] - other_coords[2];
    return sqrt(dx * dx + dy * dy + dz * dz);
}

Residue * Atom::residue() {
    return residue;
}

Chain * Atom::chain() {
    return residue.chain;
}