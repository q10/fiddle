#include "common.h"

void initialize_rotamers_for_each_residue(Chain * chain) {
    for (map <int, Residue *>::iterator it = chain->residues->begin(); it != chain->residues->end(); it++)
        initialize_rotamer_for_residue(it->second);
    return;
}

void initialize_rotamer_for_residue(Residue * residue) {
    for (int torsions_index = 0; torsions_index < (*ROTAMER_LIBRARY)[residue->name]->size(); torsions_index++) {
        grow_side_chain(residue, torsions_index);
        if (!check_rotamer_clash_against_backbone(residue))
            residue->non_clashing_rotamer_ids->push_back(torsions_index);
    }
    return;
}

bool check_rotamer_clash_against_backbone(Residue * residue) {
    map <int, Residue *> * residues = residue->chain->residues;
    vector <Atom *> *side_chain_atoms, *backbone_atoms;
    Residue * other_residue;
    Atom *atom, *other_atom;
    double vdw_distance;

    for (map <int, Residue *>::iterator other_res = residues->begin(); other_res != residues->end(); other_res++) {
        other_residue = other_res->second;

        if (other_residue->id != residue->id) {
            side_chain_atoms = residue->side_chain_atoms;
            backbone_atoms = other_residue->backbone_atoms;
            
            for (int i = 0; i < side_chain_atoms->size(); i++) {
                atom = (*side_chain_atoms)[i];                

                for (int j = 0; j < backbone_atoms->size(); j++) {
                    other_atom = (*backbone_atoms)[j];
                    vdw_distance = VDW_ALPHA * ((*VDW_RADIUS)[atom->element] + (*VDW_RADIUS)[atom->element]);
                    
                    if (atom->distance_from(other_atom) < vdw_distance)
                        return true;
                }
            }
        }
    }
    return false;
}