/* 
 * File:   backbone_clash_checker.h
 * Author: BENSON J MA
 *
 * Created on June 24, 2011, 12:48 AM
 */

#pragma once
#ifndef BACKBONE_CLASH_CHECKER_H
#define	BACKBONE_CLASH_CHECKER_H

void initialize_rotamers_for_each_residue_in_chain(Chain * chain);
void initialize_rotamer_for_residue(Residue * residue);
bool check_rotamer_clash_against_backbone(Residue * residue);
void test_initialize_rotamer_for_residue();
void test_initialize_rotamers_for_chain();

#endif	/* BACKBONE_CLASH_CHECKER_H */
