/*
 * C2DHeat.h
 *
 *  Created on: 25 Dec 2014
 *      Author: raiden
 */

#ifndef C2DHEAT_H_
#define C2DHEAT_H_

#include "C2DGrid.h"
#include "CHeat.h"

namespace std
{

class C2DHeat : public CHeat
{
public:
	C2DHeat(CInput* input);
	double** u_xy_t;
	C2DGrid* grid;
	virtual ~C2DHeat();
	virtual void evolve();
	virtual void explicit_euler_evolve(int t);
	virtual void implicit_method_evolve(int t);
private:
	int Mx;
	int My;
	double* ic_xy;
	double* source_xy;
	double CFL_x;
	double CFL_y;
	void find_dt();
	void initialise();
	double courant_number(int i);
	int shift(int i);
	double** left_boundary_condition_y_t;
	double** right_boundary_condition_y_t;
	double** upper_boundary_condition_x_t;
	double** lower_boundary_condition_x_t;
	CInput::SourceTerm source_term;
};

} /* namespace std */

#endif /* C2DHEAT_H_ */
