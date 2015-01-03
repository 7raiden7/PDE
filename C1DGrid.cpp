/*
 * C1DGrid.cpp
 *
 *  Created on: 14 Dec 2014
 *      Author: raiden
 */

#include "C1DGrid.h"

namespace std
{

void C1DGrid::compose()
{
	double dx = (range[1] - range[0]) / (Mx - 1);

	for (int i = 0; i < Mx; ++i)
		x[i] = range[0] + dx * i;
}

void C1DGrid::explicit_euler_evolution_operator()
{
	for (int i = 0; i < Mx; ++i)
	{
		A[i][0] = 1.0;
		A[i][1] = -2.0;
		A[i][2] = 1.0;
	}
}

void C1DGrid::implicit_euler_evolution_operator()
{
	for (int i = 0; i < Mx; ++i)
	{
		A[i][0] = -1.0;
		A[i][1] = 2.0;
		A[i][2] = -1.0;
	}
}

void C1DGrid::crank_nicholson_evolution_operator()
{
	for (int i = 0; i < Mx; ++i)
	{
		A[i][0] = -.5;
		A[i][1] = 1.0;
		A[i][2] = -.5;

		B[i][0] = .5;
		B[i][1] = -1.0;
		B[i][2] = .5;
	}
}

C1DGrid::C1DGrid(double* boundaries, int n, int Mx) :
		CGrid(boundaries, n, Mx)
{
	if (n != 1)
		throw;

	x = new double[Mx];
}

C1DGrid::~C1DGrid()
{
	delete[] x;
}

} /* namespace std */
