/*
 * CLinAlg.h
 *
 *  Created on: 22 Dec 2014
 *      Author: raiden
 */

#ifndef CLINALG_H_
#define CLINALG_H_

#include <math.h>

namespace std
{

class CLinAlg
{
public:
static double* tridiag_solver(double** A, double* b, int m);
static double* periodic_tridiag_solver(double** A, double* b, int m);
static double* block_tridiag_iterative_solver(double** A, double* b, int Mx, int My);
static double* block_tridiag_solver(double** A, double* b, int Mx, int My);
static double* tridiag_dot_vector(double** A, double* b, int M);
static double* block_tridiag_dot_vector(double** A, double* b, int Mx, int My);
static double** tridiag_inverse(double** A, int M);
static double** tridiag_dot_full_matrix(double** T, double** F, int M);
};

} /* namespace std */

#endif /* CLINALG_H_ */
