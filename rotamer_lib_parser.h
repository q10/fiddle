/* 
 * File:   rotamer_lib_parser.h
 * Author: BENSON J MA
 *
 * Created on June 21, 2011, 2:12 AM
 */

#pragma once
#ifndef ROTAMER_LIB_PARSER_H
#define	ROTAMER_LIB_PARSER_H

std::vector<double *> * get_rotamer_lib_from_file(std::string & filename);
void test_rotamer_lib_parser();

#endif	/* ROTAMER_LIB_PARSER_H */

