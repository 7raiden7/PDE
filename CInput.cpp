/*
 * CInput.cpp
 *
 *  Created on: 24 Dec 2014
 *      Author: raiden
 */

#include "CInput.h"

namespace std
{

CInput::CInput(int Mx, double T, double kappa, CHeat::SolverType solver_type,
		CHeat::BoundaryConditionType bc_type, double CFL, double (&x_range)[2],
		InitialCondition initial_condition, SourceTerm source_term,
		BoundaryCondition1D left_boundary_condition,
		BoundaryCondition1D right_boundary_condition)
{
	this->Mx = Mx;
	this->T = T;
	this->kappa = kappa;
	this->solver_type = solver_type;
	this->bc_type = bc_type;
	this->initial_condition = initial_condition;
	this->source_term = source_term;
	this->left_boundary_condition_1D = left_boundary_condition;
	this->right_boundary_condition_1D = right_boundary_condition;
	this->CFL = CFL;

	this->x_range[0] = x_range[0];
	this->x_range[1] = x_range[1];
}

CInput::CInput(int Mx, int My, double T, double kappa,
		CHeat::SolverType solver_type, CHeat::BoundaryConditionType bc_type,
		double CFL, double (&x_range)[4], InitialCondition initial_condition,
		SourceTerm source_term,
		BoundaryCondition2D left_boundary_condition,
		BoundaryCondition2D right_boundary_condition,
		BoundaryCondition2D lower_boundary_condition,
		BoundaryCondition2D upper_boundary_condition)
{
	this->Mx = Mx;
	this->My = My;
	this->T = T;
	this->kappa = kappa;
	this->solver_type = solver_type;
	this->bc_type = bc_type;

	this->initial_condition = initial_condition;
	this->source_term = source_term;
	this->left_boundary_condition_2D = left_boundary_condition;
	this->right_boundary_condition_2D = right_boundary_condition;
	this->lower_boundary_condition_2D = lower_boundary_condition;
	this->upper_boundary_condition_2D = upper_boundary_condition;

	this->CFL = CFL;

	this->xy_range[0] = x_range[0];
	this->xy_range[1] = x_range[1];
	this->xy_range[2] = x_range[2];
	this->xy_range[3] = x_range[3];
}

CInput::~CInput()
{
}

} /* namespace std */
