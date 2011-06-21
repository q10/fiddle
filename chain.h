/* 
 * File:   chain.h
 * Author: BENSON J MA
 *
 * Created on June 21, 2011, 2:12 AM
 */

#ifndef CHAIN_H
#define	CHAIN_H

using namespace std;

class Chain {
public:
    string id;

    Chain(string id);
    ~Chain();
    ostream & operator<<(ostream & out, const Chain& chain);
    void add_residue(Atom * atom);
    void remove_residue_with_id(int id);
    List Residue residues(); // wrong at the moment
    List Atom atoms(); // wrong at the moment
};

#endif	/* CHAIN_H */
