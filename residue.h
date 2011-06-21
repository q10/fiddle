/* 
 * File:   residue.h
 * Author: BENSON J MA
 *
 * Created on June 21, 2011, 2:11 AM
 */

#ifndef RESIDUE_H
#define	RESIDUE_H

using namespace std;

class Residue {
public:
    int id;
    string name;

    Residue(int tmp_id, string tmp_name);
    ~Residue();
    ostream & operator<<(ostream & out, const Residue& residue);
    bool has_hydrogens();
    void add_atom(Atom * atom);
    void remove_atom_with_id(int id);
    Chain * chain();
    List Atom atoms(); // wrong at the moment
    List Atom backbone_atoms(); // wrong at the moment
    List Atom side_chain_atoms(); // wrong at the moment
};

#endif	/* RESIDUE_H */
