/*
 * CGrid.h
 *
 *  Created on: 14 Dec 2014
 *      Author: raiden
 */

#ifndef CGRID_H_
#define CGRID_H_

namespace std
{

class CGrid
{
public:
	CGrid(double* boundaries, int n, int Mx, int My);
	CGrid(double* boundaries, int n, int Mx);
	virtual void compose() = 0;
	virtual void explicit_euler_evolution_operator() = 0;
	virtual void implicit_euler_evolution_operator() = 0;
	virtual void crank_nicholson_evolution_operator() = 0;
	double** A;
	double** B;
	virtual ~CGrid();
	int Mx;
	int My;
	double* x;
	double* y;
protected:
	double* range;
};

} /* namespace std */

#endif /* CGRID_H_ */
