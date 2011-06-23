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
    map <int, Residue *> * residues;

    Chain(string & tmp_id);
    ~Chain();

    friend ostream & operator<<(ostream & out, Chain * chain);

    void add_residue(Residue * residue);
    void remove_residue_with_id(int id);
    //List Atom atoms(); // wrong at the moment
};

void test_chain();
Chain * sample_chain();

#endif	/* CHAIN_H */
