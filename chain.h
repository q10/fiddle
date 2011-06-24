/* 
 * File:   chain.h
 * Author: BENSON J MA
 *
 * Created on June 21, 2011, 2:12 AM
 */

#pragma once
#ifndef CHAIN_H
#define	CHAIN_H

class Residue;

class Chain {
public:
    std::string id;
    std::map <int, Residue *> * residues;

    Chain(std::string & tmp_id);
    ~Chain();

    friend std::ostream & operator<<(std::ostream & out, Chain * chain);

    void add_residue(Residue * residue);
    void remove_residue_with_id(int id);
    //List Atom atoms(); // wrong at the moment
};

void test_chain();
Chain * sample_chain();

#endif	/* CHAIN_H */
