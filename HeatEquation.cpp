/*
 * HeatEquation.cpp
 *
 *  Created on: 14 Dec 2014
 *      Author: raiden
 */

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "C1DHeat.h"
#include "C2DHeat.h"
#include "CGNUPlot.h"
#include "CTests.h"
#include "CInput.h"

using namespace std;

double* initial_condition(double* inputX, int lenX, double* inputY, int lenY)
{
	double *output = new double[lenX * lenY];
	int ij = 0;
	double square = 0;
	for (int i = 0; i < lenX; ++i)
	{
		for (int j = 0; j < lenY; ++j)
		{
			square = inputX[i] * inputX[i] + inputY[j] * inputY[j];
			ij = i + lenX * j;
			//output[ij] = exp(-square);
			output[ij] = 2 * inputX[i] + 3 * inputY[j];
			//output[ij] = sin(2 * M_PI * (inputX[i] + inputY[j]) / 10);
		}
	}

	return output;
}

double* source_term(double* inputX, int lenX, double* inputY, int lenY,
		double* z_xy)
{
	double* f_xy = new double[lenX * lenY];
	int xy = 0;

	for (int y = 0; y < lenY; ++y)
	{
		for (int x = 0; x < lenX; ++x)
		{
			xy = x + lenX * y;
			f_xy[xy] = 5;
		}
	}

	return f_xy;
}

double* left_boundary_condition_1D(int N)
{
	double *bc = new double[N];

	for (int t = 0; t < N; ++t)
		bc[t] = 0;

	return bc;
}

double* right_boundary_condition_1D(int N)
{
	double *bc = new double[N];

	for (int i = 0; i < N; ++i)
		bc[i] = 0;

	return bc;
}

double** left_boundary_condition_2D(double* inputY, int My, double inputX0,
		int N)
{
	double **bc = new double*[N];

	double* aux_x = new double[1];
	aux_x[0] = inputX0;

	double* aux = initial_condition(aux_x, 1, inputY, My);

	for (int t = 0; t < N; ++t)
	{
		bc[t] = new double[My];
		for (int y = 0; y < My; ++y)
			bc[t][y] = aux[y];
	}

	return bc;
}

double** right_boundary_condition_2D(double* inputY, int My, double inputX0,
		int N)
{
	double **bc = new double*[N];

	double* aux_x = new double[1];
	aux_x[0] = inputX0;

	double* aux = initial_condition(aux_x, 1, inputY, My);

	for (int t = 0; t < N; ++t)
	{
		bc[t] = new double[My];
		for (int y = 0; y < My; ++y)
			bc[t][y] = aux[y];
	}

	return bc;
}

double** upper_boundary_condition_2D(double* inputX, int Mx, double inputY0,
		int N)
{
	double **bc = new double*[N];

	for (int t = 0; t < N; ++t)
	{
		bc[t] = new double[Mx];
		for (int x = 0; x < Mx; ++x)
			bc[t][x] = 2 * inputX[x] + 3 * inputY0;
	}

	return bc;
}

double** lower_boundary_condition_2D(double* inputX, int Mx, double inputY0,
		int N)
{
	double **bc = new double*[N];

	for (int t = 0; t < N; ++t)
	{
		bc[t] = new double[Mx];
		for (int x = 0; x < Mx; ++x)
			bc[t][x] = 2 * inputX[x] + 3 * inputY0;
	}

	return bc;
}

void solve_1D()
{
	double delay = 0.002;
	int Mx = 20;
	double T = 2;
	double kappa = 1;
	CHeat::SolverType solver_type = CHeat::ImplicitEuler;
	CHeat::BoundaryConditionType bc_type = CHeat::Periodic;
	double CFL = 0.45;
	double x_range[] =
	{ -2, 2 };
	CInput::InitialCondition ic = &initial_condition;
	CInput::SourceTerm source = &source_term;
	CInput::BoundaryCondition1D left_bc = &left_boundary_condition_1D;
	CInput::BoundaryCondition1D right_bc = &right_boundary_condition_1D;

	CInput* input = new CInput(Mx, T, kappa, solver_type, bc_type, CFL, x_range,
			ic, source, left_bc, right_bc);

	cout << "Solving 1D Heat Equation..." << endl;
	C1DHeat* heat = new C1DHeat(input);
	heat->evolve();

	cout << "Completed." << endl;

	CGNUPlot *plot = new CGNUPlot();

	cout << "Plotting Results..." << endl;

	double* t = new double[heat->N];

	for (int n = 0; n < heat->N; ++n)
		t[n] = n * T / heat->N;

	plot->plot_2D_animation(heat->u_x_t, (heat->grid)->x, t, Mx, heat->N,
			"Heat Equation", "x", "u(x,t)", solver_type + " Solution", delay);
//	plot->plot_3D(heat->u_x_t, (heat->grid)->x, t, Mx, heat->N, "Heat Equation",
//			"x", "Time", "u(x,t)", solver_type + " Solution");
	cout << "Done." << endl;

	delete[] t;
	delete[] heat;
	delete[] plot;
}

void solve_2D()
{
	double delay = 0.002;
	int Mx = 20;
	int My = 30;
	double T = 2;
	double kappa = 1;
	CHeat::SolverType solver_type = CHeat::ExplicitEuler;
	CHeat::BoundaryConditionType bc_type = CHeat::Periodic;
	double CFL = 0.45;
	double x_range[] =
	{ -5, 5, -5, 5 };
	CInput::InitialCondition ic = &initial_condition;
	CInput::SourceTerm source = &source_term;
	CInput::BoundaryCondition2D left_bc = &left_boundary_condition_2D;
	CInput::BoundaryCondition2D right_bc = &right_boundary_condition_2D;
	CInput::BoundaryCondition2D down_bc = &lower_boundary_condition_2D;
	CInput::BoundaryCondition2D up_bc = &upper_boundary_condition_2D;

	CInput* input = new CInput(Mx, My, T, kappa, solver_type, bc_type, CFL,
			x_range, ic, source, left_bc, right_bc, down_bc, up_bc);

	cout << "Solving 2D Heat Equation..." << endl;
	C2DHeat* heat = new C2DHeat(input);
	heat->evolve();
	cout << "Completed." << endl;

	cout << "Plotting Results..." << endl;
	CGNUPlot *plot = new CGNUPlot();

	double* t = new double[heat->N];
	for (int n = 0; n < heat->N; ++n)
		t[n] = n * T / heat->N;

	plot->plot_3D_animation(heat->u_xy_t, (heat->grid)->x, (heat->grid)->y, t,
			Mx, My, heat->N, "Heat Equation", "x", "y", "u",
			solver_type + " Solution", delay);

	cout << "Done." << endl;

	delete[] heat;
	delete[] plot;
}

int main()
{
	solve_2D();

	return 0;
}
