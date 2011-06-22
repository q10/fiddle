#include "common.h"

ostream & operator<<(ostream & out, Atom atom) {
    return out << "<Atom " << atom.name << ">";
}

ostream & operator<<(ostream & out, Atom * atom) {
    return out << "<Atom " << atom->name << ">";
}

Atom::Atom(int tmp_id, string & tmp_name, double * tmp_coords, string & tmp_element) {
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
    double * other_coords = other_atom->coords;
    double dx = coords[0] - other_coords[0];
    double dy = coords[1] - other_coords[1];
    double dz = coords[2] - other_coords[2];
    return sqrt(dx * dx + dy * dy + dz * dz);
}

Chain * Atom::chain() {
    return residue->chain;
}

void test_atom() {
    string s = "CA", t = "C";
    double * coords = new double [3];
    coords[0] = 1.3;
    coords[1] = 2.0;
    coords[2] = 3.2;
    Atom atom(1, s, coords, t);         // bad way of making new classes, b/c not in heap, so data can be overwritten (see cout atom.coords[0])
    Atom * btom = new Atom(1, s, coords, t); // correct/safe way of creating new objects
    cout << btom << " " << atom << endl;
    cout << atom.id << " " << btom->id << endl;
    cout << atom.coords[0] << ", " << atom.coords[1] << ", " << atom.coords[2] << endl;
    cout << btom->coords[0] << ", " << btom->coords[1] << ", " << btom->coords[2] << endl;
    return;
}
