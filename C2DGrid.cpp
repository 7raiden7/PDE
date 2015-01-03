/*
 * C2DGrid.cpp
 *
 *  Created on: 25 Dec 2014
 *      Author: raiden
 */

#include "C2DGrid.h"

namespace std
{

void C2DGrid::compose()
{
	double dx = (range[1] - range[0]) / (Mx - 1);
	double dy = (range[3] - range[2]) / (My - 1);

	for (int i = 0; i < Mx; ++i)
		x[i] = range[0] + dx * i;

	for (int i = 0; i < My; ++i)
		y[i] = range[2] + dy * i;
}

void C2DGrid::explicit_euler_evolution_operator()
{
	for (int i = 0; i < Mx * My; ++i)
	{
		A[i][0] = 1.0;
		A[i][1] = 1.0;
		A[i][2] = -2.0; // This will be multiplied times CFL_x + CFL_y
		A[i][3] = 1.0;
		A[i][4] = 1.0;
	}
}

void C2DGrid::implicit_euler_evolution_operator()
{
	for (int i = 0; i < Mx * My; ++i)
	{
		A[i][0] = -1.0;
		A[i][1] = -1.0;
		A[i][2] = 2.0; // This will be multiplied times CFL_x + CFL_y
		A[i][3] = -1.0;
		A[i][4] = -1.0;
	}
}

void C2DGrid::crank_nicholson_evolution_operator()
{
	for (int i = 0; i < Mx * My; ++i)
	{
		A[i][0] = -.5;
		A[i][1] = -.5;
		A[i][2] = 1.0; // This will be multiplied times CFL_x + CFL_y
		A[i][3] = -.5;
		A[i][4] = -.5;

		B[i][0] = .5;
		B[i][1] = .5;
		B[i][2] = -1.0; // This will be multiplied times CFL_x + CFL_y
		B[i][3] = .5;
		B[i][4] = .5;
	}
}

C2DGrid::C2DGrid(double* boundaries, int n, int Mx, int My) :
		CGrid(boundaries, n, Mx, My)
{
	if (n != 2)
		throw;

	x = new double[Mx];
	y = new double[My];
}

C2DGrid::~C2DGrid()
{
	delete[] x;
	delete[] y;
}

} /* namespace std */
