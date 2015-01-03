/*
 * C2DHeat.cpp
 *
 *  Created on: 25 Dec 2014
 *      Author: raiden
 */

#include "C2DHeat.h"

using namespace std;

double C2DHeat::courant_number(int i)
{
	switch (i)
	{
	case 0:
	case 4:
		return CFL_y;
	case 1:
	case 3:
		return CFL_x;
	case 2:
		return CFL;
	default:
		throw;
	}
}

int delta_kron(int i)
{
	return i == 2 ? 1 : 0;
}

int C2DHeat::shift(int i)
{
	int shift_y = i == 0 ? -1 : (i == 4 ? 1 : 0);
	int shift_x = (i != 0 && i != 4) ? (i - 2) : 0;

	return (shift_x + Mx * shift_y);
}

void C2DHeat::initialise()
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

void C2DHeat::explicit_euler_evolve(int t)
{
	double left_bc = 0;
	double right_bc = 0;
	double upper_bc = 0;
	double lower_bc = 0;

	int xy = 0;
	int xy_left = 0;
	int xy_right = 0;
	int xy_down = 0;
	int xy_up = 0;

	for (int y = 2; y < My - 2; ++y)
	{
		for (int x = 2; x < Mx - 2; ++x)
		{
			xy = x + Mx * y;

			for (int i = 0; i < 5; ++i)
				u_xy_t[t][xy] += (delta_kron(i)
						+ courant_number(i) * grid->A[xy][i])
						* u_xy_t[t - 1][xy + shift(i)];

			u_xy_t[t][xy] += dt * source_xy[xy];
		}
	}

	// y = 1 and y = My - 2 are calculated in the following loop
	for (int y = 2; y < My - 2; ++y)
	{
		xy_left = 1 + Mx * y;
		xy_right = (Mx - 2) + Mx * y;

		u_xy_t[t][xy_left] = dt * source_xy[xy_left];
		u_xy_t[t][xy_right] = dt * source_xy[xy_right];
	}

	for (int x = 1; x < Mx - 1; ++x)
	{
		xy_down = x + Mx * 1;
		xy_up = x + Mx * (My - 2);

		u_xy_t[t][xy_up] = dt * source_xy[xy_up];
		u_xy_t[t][xy_down] = dt * source_xy[xy_down];
	}

	switch (bc_type)
	{
	case Dirichlet:

		for (int y = 0; y < My; ++y)
		{
			xy_left = 0 + Mx * y;
			u_xy_t[t][xy_left] = left_boundary_condition_y_t[t][y];

			xy_right = (Mx - 1) + Mx * y;
			u_xy_t[t][xy_right] = right_boundary_condition_y_t[t][y];
		}

		for (int x = 0; x < Mx; ++x)
		{
			xy_down = x + Mx * 0;
			u_xy_t[t][xy_down] = lower_boundary_condition_x_t[t][x];

			xy_up = x + Mx * (My - 1);
			u_xy_t[t][xy_up] = upper_boundary_condition_x_t[t][x];
		}

		// y = 1 and y = My - 2 are calculated in the following loop
		for (int y = 2; y < My - 2; ++y)
		{
			xy_left = 1 + Mx * y;
			xy_right = (Mx - 2) + Mx * y;

			for (int i = 0; i < 5; ++i)
			{
				left_bc =
						i == 1 ?
								left_boundary_condition_y_t[t - 1][y] :
								u_xy_t[t - 1][xy_left + shift(i)];
				right_bc =
						i == 3 ?
								right_boundary_condition_y_t[t - 1][y] :
								u_xy_t[t - 1][xy_right + shift(i)];

				u_xy_t[t][xy_left] += (delta_kron(i)
						+ courant_number(i) * grid->A[xy_left][i]) * left_bc;

				u_xy_t[t][xy_right] += (delta_kron(i)
						+ courant_number(i) * grid->A[xy_right][i]) * right_bc;
			}
		}

		for (int x = 1; x < Mx - 1; ++x)
		{
			xy_down = x + Mx * 1;
			xy_up = x + Mx * (My - 2);

			for (int i = 0; i < 5; ++i)
			{
				lower_bc =
						i == 0 ?
								lower_boundary_condition_x_t[t - 1][x] :
								u_xy_t[t - 1][xy_down + shift(i)];

				upper_bc =
						i == 4 ?
								upper_boundary_condition_x_t[t - 1][x] :
								u_xy_t[t - 1][xy_up + shift(i)];

				u_xy_t[t][xy_up] += (delta_kron(i)
						+ courant_number(i) * grid->A[xy_up][i]) * upper_bc;

				u_xy_t[t][xy_down] += (delta_kron(i)
						+ courant_number(i) * grid->A[xy_down][i]) * lower_bc;
			}
		}

		break;

	case Neumann:
	{
		double dx_left = grid->y[1] - grid->y[0];
		double dx_right = grid->y[Mx - 1] - grid->y[Mx - 2];

		double dy_down = grid->y[1] - grid->y[0];
		double dy_up = grid->y[Mx - 1] - grid->y[Mx - 2];

		// y = 1 and y = My - 2 are calculated in the following loop
		for (int y = 2; y < My - 2; ++y)
		{
			xy_left = 1 + Mx * y;
			xy_right = (Mx - 2) + Mx * y;

			for (int i = 0; i < 5; ++i)
			{
				left_bc =
						i == 1 ?
								(-dx_left
										* left_boundary_condition_y_t[t - 1][y]) :
								u_xy_t[t - 1][xy_left + shift(i)];
				right_bc =
						i == 3 ?
								(dx_right
										* right_boundary_condition_y_t[t - 1][y]) :
								u_xy_t[t - 1][xy_right + shift(i)];

				u_xy_t[t][xy_left] += (delta_kron(i) * (1 + CFL_x)
						+ courant_number(i) * grid->A[xy_left][i]) * left_bc;

				u_xy_t[t][xy_right] += (delta_kron(i) * (1 + CFL_x)
						+ courant_number(i) * grid->A[xy_right][i]) * right_bc;
			}
		}

		for (int x = 1; x < Mx - 1; ++x)
		{
			xy_down = x + Mx * 1;
			xy_up = x + Mx * (My - 2);

			for (int i = 0; i < 5; ++i)
			{
				lower_bc =
						i == 0 ?
								(-dy_down
										* lower_boundary_condition_x_t[t - 1][x]) :
								u_xy_t[t - 1][xy_down + shift(i)];

				upper_bc =
						i == 4 ?
								(dy_up * upper_boundary_condition_x_t[t - 1][x]) :
								u_xy_t[t - 1][xy_up + shift(i)];

				u_xy_t[t][xy_up] += (delta_kron(i) * (1 + CFL_y)
						+ courant_number(i) * grid->A[xy_up][i]) * upper_bc;

				u_xy_t[t][xy_down] += (delta_kron(i) * (1 + CFL_y)
						+ courant_number(i) * grid->A[xy_down][i]) * lower_bc;
			}
		}

		for (int y = 0; y < My; ++y)
		{
			xy_left = 0 + Mx * y;
			u_xy_t[t][xy_left] = u_xy_t[t][xy_left + 1]
					- dx_left * left_boundary_condition_y_t[t][y];

			xy_right = (Mx - 1) + Mx * y;
			u_xy_t[t][xy_right] = u_xy_t[t][xy_right - 1]
					+ dx_right * right_boundary_condition_y_t[t][y];
		}

		for (int x = 0; x < Mx; ++x)
		{
			xy_down = x + Mx * 0;
			u_xy_t[t][xy_down] = u_xy_t[t][xy_down + Mx]
					- dy_down * lower_boundary_condition_x_t[t][x];

			xy_up = x + Mx * (My - 1);
			u_xy_t[t][xy_up] = u_xy_t[t][xy_up - Mx]
					+ dy_up * upper_boundary_condition_x_t[t][x];
		}

		break;
	}
	case Periodic:

		// y = 1 and y = My - 2 are calculated in the following loop
		for (int y = 2; y < My - 2; ++y)
		{
			xy_left = 1 + Mx * y;
			xy_right = (Mx - 2) + Mx * y;

			for (int i = 0; i < 5; ++i)
			{
				left_bc =
						i == 1 ?
								u_xy_t[t - 1][(Mx - 2) + Mx * y] :
								u_xy_t[t - 1][xy_left + shift(i)];
				right_bc =
						i == 3 ?
								u_xy_t[t - 1][1 + Mx * y] :
								u_xy_t[t - 1][xy_right + shift(i)];

				u_xy_t[t][xy_left] += (delta_kron(i)
						+ courant_number(i) * grid->A[xy_left][i]) * left_bc;

				u_xy_t[t][xy_right] += (delta_kron(i)
						+ courant_number(i) * grid->A[xy_right][i]) * right_bc;
			}
		}

		for (int x = 1; x < Mx - 1; ++x)
		{
			xy_down = x + Mx * 1;
			xy_up = x + Mx * (My - 2);

			for (int i = 0; i < 5; ++i)
			{
				lower_bc =
						i == 0 ?
								u_xy_t[t - 1][x + Mx * (My - 2)] :
								u_xy_t[t - 1][xy_down + shift(i)];

				upper_bc =
						i == 4 ?
								u_xy_t[t - 1][x + Mx * 1] :
								u_xy_t[t - 1][xy_up + shift(i)];

				u_xy_t[t][xy_up] += (delta_kron(i)
						+ courant_number(i) * grid->A[xy_up][i]) * upper_bc;

				u_xy_t[t][xy_down] += (delta_kron(i)
						+ courant_number(i) * grid->A[xy_down][i]) * lower_bc;
			}
		}

		for (int y = 0; y < My; ++y)
		{
			xy_left = 0 + Mx * y;
			u_xy_t[t][xy_left] = u_xy_t[t][(Mx - 2) + Mx * y];

			xy_right = (Mx - 1) + Mx * y;
			u_xy_t[t][xy_right] = u_xy_t[t][1 + Mx * y];
		}

		for (int x = 0; x < Mx; ++x)
		{
			xy_down = x + Mx * 0;
			u_xy_t[t][xy_down] = u_xy_t[t][x + Mx * (My - 2)];

			xy_up = x + Mx * (My - 1);
			u_xy_t[t][xy_up] = u_xy_t[t][x + Mx * 1];
		}

		break;
	default:
		cout << "Not implemented" << endl;
		throw;
	}
}

void C2DHeat::implicit_method_evolve(int t)
{
	bool is_fully_implicit = false;

	int xy = 0;
	int xy_left = 0;
	int xy_right = 0;
	int xy_down = 0;
	int xy_up = 0;

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

	double** l_matrix = new double*[Mx * My];

	double* r_vector = new double[Mx * My];
	double* system_solution = new double[Mx * My];

	for (int y = 1; y < My - 1; ++y)
	{
		for (int x = 1; x < Mx - 1; ++x)
		{
			xy = x + Mx * y;
			l_matrix[xy] = new double[5];
			r_vector[xy] = 0;

			for (int i = 0; i < 5; ++i)
			{
				l_matrix[xy][i] = delta_kron(i)
						+ courant_number(i) * grid->A[x][i];

				r_vector[xy] += (delta_kron(i)
						+ courant_number(i) * grid->B[x][i])
						* u_xy_t[t - 1][xy + shift(i)];
			}

			if (is_fully_implicit)
				r_vector[xy] = u_xy_t[t - 1][xy];

			r_vector[xy] += dt * source_xy[xy];
		}
	}

	for (int x = 0; x < Mx; ++x)
	{
		xy_down = x + Mx * 0;
		l_matrix[xy_down] = new double[5];

		xy_up = x + Mx * (My - 1);
		l_matrix[xy_up] = new double[5];
	}
	for (int y = 0; y < My; ++y)
	{
		xy_left = 0 + Mx * y;
		l_matrix[xy_left] = new double[5];

		xy_right = (Mx - 1) + Mx * y;
		l_matrix[xy_right] = new double[5];
	}

	switch (bc_type)
	{
	case Dirichlet:

		for (int x = 0; x < Mx; ++x)
		{
			xy_down = x + Mx * 0;
			l_matrix[xy_down][0] = 0;
			l_matrix[xy_down][1] = 0;
			l_matrix[xy_down][2] = 1;
			l_matrix[xy_down][3] = 0;
			l_matrix[xy_down][4] = 0;

			r_vector[xy_down] = lower_boundary_condition_x_t[t][x];

			xy_up = x + Mx * (My - 1);

			l_matrix[xy_up][0] = 0;
			l_matrix[xy_up][1] = 0;
			l_matrix[xy_up][2] = 1;
			l_matrix[xy_up][3] = 0;
			l_matrix[xy_up][4] = 0;

			r_vector[xy_up] = upper_boundary_condition_x_t[t][x];
		}

		for (int y = 0; y < My; ++y)
		{
			xy_left = 0 + Mx * y;
			l_matrix[xy_left][0] = 0;
			l_matrix[xy_left][1] = 0;
			l_matrix[xy_left][2] = 1;
			l_matrix[xy_left][3] = 0;
			l_matrix[xy_left][4] = 0;

			r_vector[xy_left] = left_boundary_condition_y_t[t][y];

			xy_right = (Mx - 1) + Mx * y;
			l_matrix[xy_right][0] = 0;
			l_matrix[xy_right][1] = 0;
			l_matrix[xy_right][2] = 1;
			l_matrix[xy_right][3] = 0;
			l_matrix[xy_right][4] = 0;

			r_vector[xy_right] = right_boundary_condition_y_t[t][y];
		}

		break;

	case Neumann:
	{
		double dx_left = grid->y[1] - grid->y[0];
		double dx_right = grid->y[Mx - 1] - grid->y[Mx - 2];

		double dy_down = grid->y[1] - grid->y[0];
		double dy_up = grid->y[Mx - 1] - grid->y[Mx - 2];

		for (int x = 0; x < Mx; ++x)
		{
			xy_down = x + Mx * 0;
			l_matrix[xy_down][0] = -1;
			l_matrix[xy_down][1] = 0;
			l_matrix[xy_down][2] = 1;
			l_matrix[xy_down][3] = 0;
			l_matrix[xy_down][4] = 0;

			r_vector[xy_down] = -dy_down * lower_boundary_condition_x_t[t][x];

			xy_up = x + Mx * (My - 1);

			l_matrix[xy_up][0] = 0;
			l_matrix[xy_up][1] = 0;
			l_matrix[xy_up][2] = -1;
			l_matrix[xy_up][3] = 0;
			l_matrix[xy_up][4] = 1;

			r_vector[xy_up] = dy_up * upper_boundary_condition_x_t[t][x];
		}

		for (int y = 0; y < My; ++y)
		{
			xy_left = 0 + Mx * y;
			l_matrix[xy_left][0] = 0;
			l_matrix[xy_left][1] = -1;
			l_matrix[xy_left][2] = 1;
			l_matrix[xy_left][3] = 0;
			l_matrix[xy_left][4] = 0;

			r_vector[xy_left] = -dx_left * left_boundary_condition_y_t[t][y];

			xy_right = (Mx - 1) + Mx * y;
			l_matrix[xy_right][0] = 0;
			l_matrix[xy_right][1] = 0;
			l_matrix[xy_right][2] = -1;
			l_matrix[xy_right][3] = 1;
			l_matrix[xy_right][4] = 0;

			r_vector[xy_right] = dx_right * right_boundary_condition_y_t[t][y];
		}

		break;
	}
	case Periodic:

		for (int x = 0; x < Mx; ++x)
		{
			xy_down = x + Mx * 0;
			l_matrix[xy_down][0] = 0;
			l_matrix[xy_down][1] = 0;
			l_matrix[xy_down][2] = 1;
			l_matrix[xy_down][3] = 0;
			l_matrix[xy_down][4] = 0;

			r_vector[xy_down] = r_vector[x + Mx * (My - 2)];

			xy_up = x + Mx * (My - 1);

			l_matrix[xy_up][0] = 0;
			l_matrix[xy_up][1] = 0;
			l_matrix[xy_up][2] = 1;
			l_matrix[xy_up][3] = 0;
			l_matrix[xy_up][4] = 0;

			r_vector[xy_up] = r_vector[x + Mx * 1];
		}

		for (int y = 0; y < My; ++y)
		{
			xy_left = 0 + Mx * y;
			l_matrix[xy_left][0] = 0;
			l_matrix[xy_left][1] = 0;
			l_matrix[xy_left][2] = 1;
			l_matrix[xy_left][3] = 0;
			l_matrix[xy_left][4] = 0;

			r_vector[xy_left] = r_vector[(Mx - 2) + Mx * y];

			xy_right = (Mx - 1) + Mx * y;
			l_matrix[xy_right][0] = 0;
			l_matrix[xy_right][1] = 0;
			l_matrix[xy_right][2] = 1;
			l_matrix[xy_right][3] = 0;
			l_matrix[xy_right][4] = 0;

			r_vector[xy_right] = r_vector[1 + Mx * y];
		}

		break;
	default:
		cout << "Not implemented" << endl;
		throw;
	}

	switch (bc_type)
	{
	case Dirichlet:
	case Periodic:
		system_solution = CLinAlg::block_tridiag_iterative_solver(l_matrix,
				r_vector, Mx, My);
		break;
	default:
		cout
				<< "Poisson Solver with Neumann BC to be implemented."
				<< endl;
		throw;
	}

	for (int xy = 0; xy < Mx * My; ++xy)
		u_xy_t[t][xy] = system_solution[xy];

	delete[] l_matrix;
	delete[] r_vector;
	delete[] system_solution;
}

void C2DHeat::evolve()
{
	u_xy_t[0] = ic_xy;

	for (int t = 1; t < N; ++t)
	{
		u_xy_t[t] = new double[Mx * My];

		for (int xy = 0; xy < Mx * My; ++xy)
			u_xy_t[t][xy] = 0;

		source_xy = source_term(grid->x, Mx, grid->y, My, u_xy_t[t - 1]);

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

void C2DHeat::find_dt()
{
	double dx = grid->x[1] - grid->x[0];
	double dy = grid->y[1] - grid->y[0];

	double dx_2 = dx * dx;
	double dy_2 = dy * dy;
	double kappa_2 = kappa * kappa;

	dt = CFL * (dx_2 * dy_2) / (kappa_2 * (dx_2 + dy_2));

	CFL_x = kappa_2 * dt / dx_2;
	CFL_y = kappa_2 * dt / dy_2;

	N = (int) (T / dt);
}

C2DHeat::C2DHeat(CInput* input) :
		CHeat(input->T, input->kappa, input->CFL, input->solver_type,
				input->bc_type)
{
	Mx = input->Mx;
	My = input->My;
	grid = new C2DGrid(input->xy_range, 2, Mx, My);

	initialise();

	ic_xy = input->initial_condition(grid->x, Mx, grid->y, My);

	this->source_term = input->source_term;

	left_boundary_condition_y_t = input->left_boundary_condition_2D(grid->y, My,
			grid->x[0], N);
	right_boundary_condition_y_t = input->right_boundary_condition_2D(grid->y,
			My, grid->x[Mx - 1], N);
	lower_boundary_condition_x_t = input->lower_boundary_condition_2D(grid->x,
			Mx, grid->y[0], N);
	upper_boundary_condition_x_t = input->upper_boundary_condition_2D(grid->x,
			Mx, grid->y[My - 1], N);

	u_xy_t = new double*[N];
}

C2DHeat::~C2DHeat()
{
	delete[] grid;
	delete[] u_xy_t;
	delete[] ic_xy;

	delete[] left_boundary_condition_y_t;
	delete[] right_boundary_condition_y_t;
	delete[] lower_boundary_condition_x_t;
	delete[] upper_boundary_condition_x_t;
}

