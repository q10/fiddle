/* 
 * File:   random.h
 * Author: BENSON J MA
 *
 * Created on June 22, 2011, 12:39 AM
 */

#pragma once
#ifndef RANDOM_H
#define	RANDOM_H

using namespace std;

double RAN3();
int RANDINT(int low, int high);
double MIN(double a, double b);
double ROUND(double d);
double RANDGAUSS(double mean, double stdev);
double RANDGAUSS();
double * RANDUNITVECTOR();

string TIMESTAMP();
void ASSERT(bool expression, string error_msg);

template <typename T> string STRING(T tval) {
    stringstream out;
    out << tval;
    return out.str();
}

inline double DEG2RAD(double degrees) {
    return degrees * M_PI / 180.0;
}

#endif	/* RANDOM_H */
