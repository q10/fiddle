/* 
 * File:   residue.h
 * Author: BENSON J MA
 *
 * Created on June 21, 2011, 2:11 AM
 */

#ifndef RESIDUE_H
#define	RESIDUE_H

#pragma once

using namespace std;

class Chain;
class Atom;

class Residue {
public:
    int id;
    string name;
    Chain * chain;

    Residue(int tmp_id, string tmp_name);
    ~Residue();

    friend ostream & operator<<(ostream & out, Residue residue);
    friend ostream & operator<<(ostream & out, Residue * residue);
    
    bool has_hydrogens();
    void add_atom(Atom * atom);
    void remove_atom_with_id(int id);
    //List Atom atoms(); // wrong at the moment
    //List Atom backbone_atoms(); // wrong at the moment
    //List Atom side_chain_atoms(); // wrong at the moment
};

void test_residue();

#endif	/* RESIDUE_H */
