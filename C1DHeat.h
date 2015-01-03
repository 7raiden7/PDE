/*
 * Heat.h
 *
 *  Created on: 15 Dec 2014
 *      Author: raiden
 */

#include "C1DGrid.h"
#include "CHeat.h"

#ifndef HEAT1D_H_
#define HEAT1D_H_

namespace std
{

class C1DHeat : public CHeat
{
public:
	C1DHeat(CInput* input);
	double** u_x_t;
	C1DGrid* grid;
	virtual ~C1DHeat();
	virtual void evolve();
	virtual void explicit_euler_evolve(int t);
	virtual void implicit_method_evolve(int t);
private:
	int Mx;
	double* ic_x;
	double* source_x;
	void find_dt();
	void initialise();
	double* left_boundary_condition_t;
	double* right_boundary_condition_t;
	double* upper_boundary_condition_t;
	double* lower_boundary_condition_t;
	CInput::SourceTerm source_term;
};

} /* namespace std */

#endif /* HEAT_H_ */
