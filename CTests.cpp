/*
 * CTests.cpp
 *
 *  Created on: 23 Dec 2014
 *      Author: raiden
 */

#include "CTests.h"

#include "CLinAlg.h"
#include <math.h>
#include <iostream>
#include <fstream>

using namespace std;

void write_array(string filename, double* x, int m)
{
	string filename_x = filename + ".txt";

	ofstream myfile(filename_x.c_str());
	if (myfile.is_open())
	{
		for (int i = 0; i < m; ++i)
			myfile << x[i] << endl;

		myfile.close();
	}
	else
	{
		cout << "Unable to open file";
		throw;
	}
}

void CTests::thomas_algorithm()
{
	int M = 20;

	double **A, *x, *debug_b, *b;

	A = new double *[M];
	x = new double[M];
	debug_b = new double[M];
	b = new double[M];

	for (int i = 0; i < M; ++i)
	{
		A[i] = new double[3];

		A[i][0] = -1.0;
		A[i][1] = 2.0;
		A[i][2] = -1.0;

		double _x = i - M / 2;
		b[i] = exp(-_x * _x);
	}

	x = CLinAlg::tridiag_solver(A, b, M);

	for (int i = 0; i < M; ++i)
	{
		debug_b[i] = 0;

		for (int j = 0; j < 3; ++j)
		{
			if ((i == 0) && (j == 0))
				continue;

			if ((i == M - 1) && (j == 2))
				continue;

			debug_b[i] += A[i][j] * x[i + (j - 1)];
		}

		if (fabs(debug_b[i] - b[i]) >= 1E-10)
		{
			cout << "Failure at " << i << endl;
			write_array("__b", b, M);
			write_array("__debug_b", debug_b, M);
			write_array("__x", x, M);
			throw;
		}
	}
}

void CTests::periodic_thomas_algorithm()
{
	int M = 20;

	double **A, *x, *debug_b, *b;

	A = new double *[M];
	x = new double[M];
	debug_b = new double[M];
	b = new double[M];

	for (int i = 0; i < M; ++i)
	{
		A[i] = new double[3];

		A[i][0] = -1.0;
		A[i][1] = 2.0;
		A[i][2] = -1.0;

		double _x = i - M / 2;
		b[i] = exp(-_x * _x);
	}

	x = CLinAlg::periodic_tridiag_solver(A, b, M);

	for (int i = 0; i < M; ++i)
	{
		debug_b[i] = 0;

		for (int j = 0; j < 3; ++j)
		{
			if (i == 0 && j == 0)
			{
				debug_b[i] += A[0][0] * x[M - 1];
				continue;
			}
			if ((i == M - 1) && (j == 2))
			{
				debug_b[i] += A[M - 1][2] * x[0];
				continue;
			}

			debug_b[i] += A[i][j] * x[i + (j - 1)];
		}

		if (fabs(debug_b[i] - b[i]) >= 1E-10)
		{
			cout << "Failure at " << i << endl;
			write_array("__b", b, M);
			write_array("__debug_b", debug_b, M);
			write_array("__x", x, M);
			throw;
		}
	}
}

void CTests::successive_over_relaxations()
{
	int M = 20;

	double** A = new double*[M * M];
	double* b = new double[M * M];

	for (int i = 0; i < M * M; ++i)
	{
		A[i] = new double[5];
		b[i] = (double) i / M;
		A[i][0] = -1.0;
		A[i][1] = -1.0;
		A[i][2] = 4;
		A[i][3] = -1;
		A[i][4] = -1;
	}

	double* x = new double[M * M];

	x = CLinAlg::block_tridiag_iterative_solver(A, b, M, M);

	double *_b = new double[M * M];

	_b = CLinAlg::block_tridiag_dot_vector(A, x, M, M);

	for (int i = 0; i < M * M; ++i)
	{
		if (fabs(b[i] - _b[i]) > 1E-7)
		{
			cout << "Algorithm did not converge" << endl;
			throw;
		}
	}

	cout << "Done." << endl;
}
