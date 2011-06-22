/* 
 * File:   chain.h
 * Author: BENSON J MA
 *
 * Created on June 21, 2011, 2:12 AM
 */

#ifndef CHAIN_H
#define	CHAIN_H

#pragma once

using namespace std;

class Residue;
class Atom;

class Chain {
public:
    string id;

    Chain(string id);
    ~Chain();
    
    friend ostream & operator<<(ostream & out, Chain chain);
    friend ostream & operator<<(ostream & out, Chain * chain);

    void add_residue(Atom * atom);
    void remove_residue_with_id(int id);
    //List Residue residues(); // wrong at the moment
    //List Atom atoms(); // wrong at the moment
};

void test_chain();

#endif	/* CHAIN_H */
