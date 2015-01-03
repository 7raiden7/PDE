/*
 * CHeat.h
 *
 *  Created on: 19 Dec 2014
 *      Author: raiden
 */

#include <math.h>
#include <iostream>
#include "CGrid.h"
#include "CLinAlg.h"
#include "CInput.h"

#ifndef CHEAT_H_
#define CHEAT_H_

namespace std
{

class CHeat
{
public:
	enum SolverType{ExplicitEuler, ImplicitEuler, CrankNicholson};
	enum BoundaryConditionType{Dirichlet, Neumann, Periodic};
	CHeat(double kappa, double T, double CFL, SolverType st, BoundaryConditionType bc);
	virtual ~CHeat();
	int N;
	SolverType solver_type;
	BoundaryConditionType bc_type;
	virtual void evolve() = 0;
	virtual void explicit_euler_evolve(int t) = 0;
	virtual void implicit_method_evolve(int t) = 0;
protected:
	double kappa;
	double T;
	double CFL;
	double dt;
};

} /* namespace std */

#endif /* CHEAT_H_ */
