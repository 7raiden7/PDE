/*
 * CGrid.cpp
 *
 *  Created on: 14 Dec 2014
 *      Author: raiden
 */

#include "CGrid.h"

namespace std
{

CGrid::CGrid(double* range, int n, int Mx)
{
	if (n != 1)
		throw;

	this->range = new double[2];
	this->Mx = Mx;
	this->My = My;

	int tot_dim = 1;

	this->A = new double *[Mx];

	for (int i = 0; i < Mx; ++i)
		A[i] = new double[3];

	this->B = new double *[Mx];

	for (int i = 0; i < Mx; ++i)
		B[i] = new double[3];

	for (int i = 0; i < 2; ++i)
		this->range[i] = range[i];
}

CGrid::CGrid(double* range, int n, int Mx, int My)
{
	if (n != 2)
		throw;

	this->range = new double[4];
	this->Mx = Mx;
	this->My = My;

	this->A = new double *[Mx * My];

	for (int i = 0; i < Mx * My; ++i)
		A[i] = new double[5];

	this->B = new double *[Mx * My];

	for (int i = 0; i < Mx * My; ++i)
		B[i] = new double[5];

	for (int i = 0; i < 4; ++i)
		this->range[i] = range[i];
}

CGrid::~CGrid()
{
	delete[] range;
	delete[] A;
	delete[] B;
}

void CGrid::compose()
{
	throw "Not Implemented";
}

void CGrid::explicit_euler_evolution_operator()
{
	throw "Not Implemented";
}

void CGrid::implicit_euler_evolution_operator()
{
	throw "Not Implemented";
}

void CGrid::crank_nicholson_evolution_operator()
{
	throw "Not Implemented";
}

} /* namespace std */
