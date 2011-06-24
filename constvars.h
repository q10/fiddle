/* 
 * File:   constvars.h
 * Author: BENSON J MA
 *
 * Created on June 21, 2011, 2:13 AM
 */

#pragma once
#ifndef CONSTVARS_H
#define	CONSTVARS_H

#define BUF_SIZE 1024

extern const double VDW_ALPHA;
extern std::map <std::string, double> * VDW_RADIUS;
extern std::map <std::string, std::vector<double *> *> * ROTAMER_LIBRARY;

void initialize_constants();
void initialize_rotamer_library();
void test_rotamer_library_creation();

#endif	/* CONSTVARS_H */
