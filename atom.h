/* 
 * File:   atom.h
 * Author: BENSON J MA
 *
 * Created on June 21, 2011, 2:11 AM
 */

#pragma once
#ifndef ATOM_H
#define	ATOM_H

class Chain;
class Residue;

class Atom {
public:
    int id;
    std::string name, element;
    double * coords;
    Residue * residue;

    Atom(int tmp_id, std::string & tmp_name, double * tmp_coords, std::string & tmp_element);
    ~Atom();
    
    // Must be friend to access private members.
    friend std::ostream & operator<<(std::ostream & out, Atom atom);
    friend std::ostream & operator<<(std::ostream & out, Atom * atom);
    
    void set_coords(double * new_coords);
    void set_residue(Residue * tmp_residue);
    double distance_from(Atom * other_atom);
    Chain * chain();
};

void test_atom();
Atom * sample_atom();

#endif	/* ATOM_H */
