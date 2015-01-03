/*
 * CLinAlg.cpp
 *
 *  Created on: 22 Dec 2014
 *      Author: raiden
 */

#include "CLinAlg.h"

using namespace std;

#define MIN(x,y) ((x)>(y)?(y):(x))
#define MAX(x,y) ((x)>(y)?(x):(y))

double norm_2(double* x, int m)
{
	double norm = 0;

	for (int i = 0; i < m; ++i)
		norm += (x[i] * x[i]);

	return sqrt(norm);
}

double vector_deviation(double* a, double* b, int m)
{
	double* c = new double[m];
	for (int i = 0; i < m; ++i)
		c[i] = a[i] - b[i];

	return (norm_2(c, m) / norm_2(a, m));
}

// Thomas Algorithm
double* CLinAlg::tridiag_solver(double** A, double* d, int m)
{
	double* x = new double[m];

	static double *c_i, *d_i;

	if (!c_i)
	{
		c_i = new double[m - 1];
		d_i = new double[m - 1];
	}

	c_i[0] = A[0][2] / A[0][1];
	d_i[0] = d[0] / A[0][1];

	for (int i = 1; i < m - 1; ++i)
	{
		double denominator = A[i][1] - A[i][0] * c_i[i - 1];

		c_i[i] = A[i][2] / denominator;
		d_i[i] = (d[i] - A[i][0] * d_i[i - 1]) / denominator;
	}

	x[m - 1] = (d[m - 1] - A[m - 1][0] * d_i[m - 2])
			/ (A[m - 1][1] - A[m - 1][0] * c_i[m - 2]);

	for (int i = m - 2; i >= 0; --i)
		x[i] = d_i[i] - c_i[i] * x[i + 1];

	return x;
}

// Modified Thomas Algorithm
double* CLinAlg::periodic_tridiag_solver(double** A, double* d, int m)
{
	double* x = new double[m];

	static double *x_1, *x_2, *q;

	if (!x_1)
	{
		x_1 = new double[m - 1];
		x_2 = new double[m - 1];
		q = new double[m - 1];
	}

	for (int i = 0; i < m - 1; ++i)
		q[i] = 0;

	q[0] = -A[0][0];
	q[m - 2] = -A[m - 1][2];

	x_1 = tridiag_solver(A, d, m - 1);
	x_2 = tridiag_solver(A, q, m - 1);

	double num = d[m - 1] - A[m - 1][2] * x_1[0] - A[m - 1][0] * x_1[m - 2];
	double den = A[m - 1][1] + A[m - 1][2] * x_2[0] + A[m - 1][	0] * x_2[m - 2];

	x[m - 1] = num / den;

	for (int i = 0; i < m - 1; ++i)
		x[i] = x_1[i] + x_2[i] * x[m - 1];

	return x;
}

// Successive OverRelaxation with optimal omega as of the 2D Poisson linear system
// --- Suitable only for Dirichlet boundary conditions ---
double* CLinAlg::block_tridiag_iterative_solver(double** A, double* b, int Mx,
		int My)
{
	int m = Mx * My;

	double* x0 = new double[m];
	double* x = new double[m];
	const double tolerance = 1E-10;
	const int iter_max = 10000;

	double spectral_radius = 1.0 / 4 * (cos(M_PI / Mx) + cos(M_PI / My));
	double omega = 2.0 / (1 + sqrt(1 - spectral_radius * spectral_radius));

	double err = 0;
	int counter = 0;

	double lower_sum = 0;
	double upper_sum = 0;

	for (int i = 0; i < m; ++i)
	{
		// initial guess
		x0[i] = 1;

		// initialisation
		x[i] = 0;
	}

	do
	{
		counter++;
		if (counter > iter_max)
			throw;

		for (int i = 0; i < m; ++i)
		{
			lower_sum = 0;
			upper_sum = 0;

			if (i >= 1)
			{
				lower_sum += A[i][1] * x[i - 1];
				if (i >= Mx)
					lower_sum += A[i][0] * x[i - Mx];
			}

			if (i <= m - 2)
			{
				upper_sum += A[i][3] * x0[i + 1];
				if (i + Mx <= Mx * My - 1)
					upper_sum += A[i][4] * x0[i + Mx];
			}

			x[i] = omega * ((-lower_sum - upper_sum + b[i]) / A[i][2])
					+ (1 - omega) * x0[i];
		}

		err = vector_deviation(x0, x, m);

		for (int i = 0; i < m; ++i)
			x0[i] = x[i];

	} while (err > tolerance);

	delete[] x0;

	return x;
}

double* CLinAlg::block_tridiag_dot_vector(double** A, double* b, int Mx, int My)
{
	double* ret = new double[Mx * My];

	for (int i = 0; i < Mx * My; ++i)
	{
		ret[i] = A[i][2] * b[i];

		if (i >= 1)
		{
			ret[i] += A[i][1] * b[i - 1];
			if (i >= Mx)
				ret[i] += A[i][0] * b[i - Mx];
		}

		if (i <= Mx * My - 2)
		{
			ret[i] += A[i][3] * b[i + 1];
			if (i + Mx <= Mx * My - 1)
				ret[i] += A[i][4] * b[i + Mx];
		}
	}

	return ret;
}

// Extended Thomas Algorithm
double* CLinAlg::block_tridiag_solver(double** A, double* d, int Mx, int My)
{
	double* x = new double[Mx * My];

	static double ***a, ***b, ***c, ***c_prime, ***d_prime;

	int ij = 0;

	if (!a)
	{
		a = new double**[My];
		b = new double**[My];
		c = new double**[My];
		c_prime = new double**[My];
		d_prime = new double**[My];

		for (int j = 0; j < My; ++j)
		{
			a[j] = new double*[Mx];
			b[j] = new double*[Mx];
			c[j] = new double*[Mx];
			c_prime[j] = new double*[Mx];
			d_prime[j] = new double*[Mx];

			for (int i = 0; i < Mx; ++i)
			{
				a[j][i] = new double[3];
				b[j][i] = new double[3];
				c[j][i] = new double[3];
				c_prime[j][i] = new double[3];
				d_prime[j][i] = new double[3];
			}
		}
	}

	for (int j = 0; j < My; ++j)
	{
		for (int i = 0; i < Mx; ++i)
		{
			ij = i + Mx * j;

			a[j][i][0] = a[j][i][2] = 0;
			a[j][i][1] = A[ij][0];

			b[j][i][0] = A[ij][1];
			b[j][i][1] = A[ij][2];
			b[j][i][2] = A[ij][3];

			c[j][i][0] = c[j][i][2] = 0;
			c[j][i][1] = A[ij][4];
		}
	}

//	c_prime[0] = A[0][2] / A[0][1];
//	d_prime[0] = d[0] / A[0][1];
//
//	for (int i = 1; i < m - 1; ++i)
//	{
//		double denominator = A[i][1] - A[i][0] * c_i[i - 1];
//
//		c_i[i] = A[i][2] / denominator;
//		d_i[i] = (d[i] - A[i][0] * d_i[i - 1]) / denominator;
//	}
//
//	x[m - 1] = (d[m - 1] - A[m - 1][0] * d_i[m - 2])
//			/ (A[m - 1][1] - A[m - 1][0] * c_i[m - 2]);
//
//	for (int i = m - 2; i >= 0; --i)
//		x[i] = d_i[i] - c_i[i] * x[i + 1];

	return x;
}

