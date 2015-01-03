/*
 * CInput.h
 *
 *  Created on: 24 Dec 2014
 *      Author: raiden
 */

#ifndef CINPUT_H_
#define CINPUT_H_

#include "CHeat.h"

namespace std
{

class CInput
{
public:
	typedef double* (*InitialCondition)(double* inputX, int lenX,
			double* inputY, int lenY);
	typedef double* (*BoundaryCondition1D)(int N);
	typedef double** (*BoundaryCondition2D)(double* inputXY, int Mxy,
			double inputXY0, int N);
	typedef double* (*SourceTerm)(double* inputX, int lenX, double* inputY,
			int lenY, double* z_xy);

	CInput(int Mx, double T, double kappa, CHeat::SolverType solver_type,
			CHeat::BoundaryConditionType bc_type, double CFL,
			double (&x_range)[2], InitialCondition initial_condition,
			SourceTerm source, BoundaryCondition1D left_boundary_condition,
			BoundaryCondition1D right_boundary_condition);

	CInput(int Mx, int My, double T, double kappa,
			CHeat::SolverType solver_type, CHeat::BoundaryConditionType bc_type,
			double CFL, double (&x_range)[4],
			InitialCondition initial_condition, SourceTerm source,
			BoundaryCondition2D left_boundary_condition,
			BoundaryCondition2D right_boundary_condition,
			BoundaryCondition2D lower_boundary_condition,
			BoundaryCondition2D upper_boundary_condition);

	virtual ~CInput();

	int Mx;
	int My;
	double T;
	double kappa;
	CHeat::SolverType solver_type;
	CHeat::BoundaryConditionType bc_type;

	InitialCondition initial_condition;
	SourceTerm source_term;

	BoundaryCondition1D left_boundary_condition_1D;
	BoundaryCondition1D right_boundary_condition_1D;

	BoundaryCondition2D left_boundary_condition_2D;
	BoundaryCondition2D right_boundary_condition_2D;

	BoundaryCondition2D upper_boundary_condition_2D;
	BoundaryCondition2D lower_boundary_condition_2D;
	double CFL;
	double x_range[2];
	double xy_range[4];
};

} /* namespace std */

#endif /* CINPUT_H_ */
