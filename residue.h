/* 
 * File:   residue.h
 * Author: BENSON J MA
 *
 * Created on June 21, 2011, 2:11 AM
 */

#pragma once
#ifndef RESIDUE_H
#define	RESIDUE_H

class Residue {
public:
    int id;
    std::string name;
    Chain * chain;
    std::map <std::string, Atom *> * atoms;
    std::vector< Atom * > * backbone_atoms;
    std::vector< Atom * > * side_chain_atoms;
    std::vector< int > * non_clashing_rotamer_ids;

    Residue(int tmp_id, std::string & tmp_name);
    ~Residue();

    friend std::ostream & operator<<(std::ostream & out, Residue * residue);

    bool has_hydrogens();
    bool is_backbone_atom(Atom * atom);
    void add_atom(Atom * atom);
    void remove_atom_with_name(std::string & name);
    void set_chain(Chain * tmp_chain);
};

void test_residue();
Residue * sample_residue();

#endif	/* RESIDUE_H */
