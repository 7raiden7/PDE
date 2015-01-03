/*
 * Heat.cpp
 *
 *  Created on: 15 Dec 2014
 *      Author: raiden
 */

#include "C1DHeat.h"

using namespace std;

void C1DHeat::initialise()
{
	grid->compose();

	switch (solver_type)
	{
	case ExplicitEuler:
		grid->explicit_euler_evolution_operator();
		break;
	case ImplicitEuler:
		grid->implicit_euler_evolution_operator();
		break;
	case CrankNicholson:
		grid->crank_nicholson_evolution_operator();
		break;
	default:
		cout << "Not Implemented" << endl;
		throw;
	}

	find_dt();
}

void C1DHeat::explicit_euler_evolve(int t)
{
	int delta_kron = 0;
	double left_bc = 0;
	double right_bc = 0;

	for (int x = 2; x < Mx - 2; ++x)
	{
		for (int i = 0; i < 3; ++i)
		{
			delta_kron = i == 1 ? 1 : 0;
			u_x_t[t][x] += (delta_kron + CFL * grid->A[x][i])
					* u_x_t[t - 1][x + (i - 1)];
		}
		u_x_t[t][x] += dt * source_x[x];
	}

	u_x_t[t][1] = dt * source_x[1];
	u_x_t[t][Mx - 2] = dt * source_x[Mx - 2];

	switch (bc_type)
	{
	case Dirichlet:

		u_x_t[t][0] = left_boundary_condition_t[t];
		u_x_t[t][Mx - 1] = right_boundary_condition_t[t];

		for (int i = 0; i < 3; ++i)
		{
			delta_kron = i == 1 ? 1 : 0;
			left_bc =
					i == 0 ?
							left_boundary_condition_t[t - 1] :
							u_x_t[t - 1][1 + (i - 1)];
			right_bc =
					i == 2 ?
							right_boundary_condition_t[t - 1] :
							u_x_t[t - 1][Mx - 2 + (i - 1)];

			u_x_t[t][1] += (delta_kron + CFL * grid->A[1][i]) * left_bc;

			u_x_t[t][Mx - 2] += (delta_kron + CFL * grid->A[Mx - 2][i])
					* right_bc;
		}

		break;
	case Neumann:
	{
		double dx_left = grid->x[1] - grid->x[0];
		double dx_right = grid->x[Mx - 1] - grid->x[Mx - 2];

		u_x_t[t][1] += (1 + CFL * (grid->A[1][1] + 1)) * u_x_t[t - 1][1]
				+ (CFL * grid->A[1][2]) * u_x_t[t - 1][2]
				+ (CFL * grid->A[1][0])
						* (-dx_left * left_boundary_condition_t[t - 1]);

		u_x_t[t][Mx - 2] += (1 + CFL * (grid->A[Mx - 2][1] + 1))
				* u_x_t[t - 1][Mx - 2]
				+ (CFL * grid->A[Mx - 2][0]) * u_x_t[t - 1][Mx - 3]
				+ (CFL * grid->A[Mx - 2][2])
						* (dx_right * right_boundary_condition_t[t - 1]);

		u_x_t[t][0] = u_x_t[t][1] - dx_left * left_boundary_condition_t[t - 1];
		u_x_t[t][Mx - 1] = u_x_t[t][Mx - 2]
				+ dx_right * right_boundary_condition_t[t - 1];

		break;
	}
	case Periodic:

		u_x_t[t][1] += (1 + CFL * grid->A[1][1]) * u_x_t[t - 1][1]
				+ (CFL * grid->A[1][2]) * u_x_t[t - 1][2]
				+ (CFL * grid->A[1][0]) * u_x_t[t - 1][Mx - 2];

		u_x_t[t][Mx - 2] += (1 + CFL * grid->A[Mx - 2][1]) * u_x_t[t - 1][Mx - 2]
				+ (CFL * grid->A[Mx - 2][0]) * u_x_t[t - 1][Mx - 3]
				+ (CFL * grid->A[Mx - 2][2]) * u_x_t[t - 1][1];

		u_x_t[t][0] = u_x_t[t][Mx - 2];
		u_x_t[t][Mx - 1] = u_x_t[t][1];

		break;
	default:
		cout << "Not implemented" << endl;
		throw;
	}
}

void C1DHeat::implicit_method_evolve(int t)
{
	// Without Source Term:
	// Implicit Euler is Fully Implicit
	// Crank-Nicholson is Semi-Implicit

	// With Source Term:
	// Implicit Euler is Implicit-Explicit (IMEX) - source term is treated as explicit
	// Crank-Nicholson is Implicit-Explicit (IMEX) - source term is treated as explicit
	bool is_fully_implicit = false;

	switch (solver_type)
	{
	case ImplicitEuler:
		is_fully_implicit = true;
		break;
	case CrankNicholson:
		// do nothing: CN is semi-implicit
		break;
	default:
		cout << "Not Supported" << endl;
		throw;
	}
	int delta_kron = 0;
	double** l_matrix = new double*[Mx];

	double* r_vector = new double[Mx];
	double* system_solution = new double[Mx];

	for (int x = 1; x < Mx - 1; ++x)
	{
		l_matrix[x] = new double[3];
		r_vector[x] = 0;

		for (int i = 0; i < 3; ++i)
		{
			delta_kron = i == 1 ? 1 : 0;
			l_matrix[x][i] = delta_kron + CFL * grid->A[x][i];

			r_vector[x] += (delta_kron + CFL * grid->B[x][i])
					* u_x_t[t - 1][x + (i - 1)];
		}

		if (is_fully_implicit)
			r_vector[x] = u_x_t[t - 1][x];

		r_vector[x] += dt * source_x[x];
	}

	l_matrix[0] = new double[3];
	l_matrix[Mx - 1] = new double[3];

	switch (bc_type)
	{
	case Dirichlet:

		l_matrix[0][0] = 0;
		l_matrix[0][1] = 1;
		l_matrix[0][2] = 0;

		l_matrix[Mx - 1][0] = 0;
		l_matrix[Mx - 1][1] = 1;
		l_matrix[Mx - 1][2] = 0;

		r_vector[0] = left_boundary_condition_t[t];
		r_vector[Mx - 1] = right_boundary_condition_t[t];

		break;

	case Neumann:
	{
		double dx_left = grid->x[1] - grid->x[0];
		double dx_right = grid->x[Mx - 1] - grid->x[Mx - 2];

		l_matrix[0][0] = 0;
		l_matrix[0][1] = 1;
		l_matrix[0][2] = -1;

		l_matrix[Mx - 1][0] = -1;
		l_matrix[Mx - 1][1] = 1;
		l_matrix[Mx - 1][2] = 0;

		r_vector[0] = -dx_left * left_boundary_condition_t[t];
		r_vector[Mx - 1] = dx_right * right_boundary_condition_t[t];

		break;
	}
	case Periodic:

		l_matrix[0][0] = 0;
		l_matrix[0][1] = 1;
		l_matrix[0][2] = 0;

		l_matrix[Mx - 1][0] = 0;
		l_matrix[Mx - 1][1] = 1;
		l_matrix[Mx - 1][2] = 0;

		r_vector[0] = r_vector[Mx - 2];
		r_vector[Mx - 1] = r_vector[1];

		break;
	default:
		cout << "Not implemented" << endl;
		throw;
	}

	system_solution = CLinAlg::tridiag_solver(l_matrix, r_vector, Mx);

	for (int x = 0; x < Mx; ++x)
		u_x_t[t][x] = system_solution[x];

	delete[] l_matrix;
	delete[] r_vector;
	delete[] system_solution;
}

void C1DHeat::evolve()
{
	u_x_t[0] = ic_x;

	for (int t = 1; t < N; ++t)
	{
		u_x_t[t] = new double[Mx];

		for (int x = 0; x < Mx; ++x)
			u_x_t[t][x] = 0;

		double y0[1] =
		{ 0 };
		source_x = source_term(grid->x, Mx, y0, 1, u_x_t[t - 1]);

		switch (solver_type)
		{
		case ExplicitEuler:
			explicit_euler_evolve(t);
			break;
		case ImplicitEuler:
		case CrankNicholson:
			implicit_method_evolve(t);
			break;
		default:
			cout << "Not Implemented" << endl;
			throw;
		}
	}
}

void C1DHeat::find_dt()
{
	double dx = grid->x[1] - grid->x[0];

	dt = CFL * (dx * dx) / (kappa * kappa);

	N = (int) (T / dt);
}

C1DHeat::C1DHeat(CInput* input) :
		CHeat(input->T, input->kappa, input->CFL, input->solver_type,
				input->bc_type)
{
	Mx = input->Mx;
	grid = new C1DGrid(input->x_range, 1, Mx);

	initialise();

	double y0[1] =
	{ 0 };
	ic_x = input->initial_condition(grid->x, Mx, y0, 1);
	this->source_term = input->source_term;
	left_boundary_condition_t = input->left_boundary_condition_1D(N);
	right_boundary_condition_t = input->right_boundary_condition_1D(N);

	u_x_t = new double*[N];
	source_x = new double[Mx];
}

C1DHeat::~C1DHeat()
{
	delete[] grid;
	delete[] u_x_t;
	delete[] ic_x;

	delete[] left_boundary_condition_t;
	delete[] right_boundary_condition_t;
}
