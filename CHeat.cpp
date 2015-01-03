/*
 * CHeat.cpp
 *
 *  Created on: 19 Dec 2014
 *      Author: raiden
 */

#include "CHeat.h"

namespace std
{

CHeat::CHeat(double T, double kappa, double CFL, SolverType solver_type,
		BoundaryConditionType bc_type)
{
	this->T = T;
	this->kappa = kappa;
	this->CFL = CFL;
	this->solver_type = solver_type;
	this->bc_type = bc_type;
}

CHeat::~CHeat()
{
}

}

/* namespace std */
