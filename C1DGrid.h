/*
 * C1DGrid.h
 *
 *  Created on: 14 Dec 2014
 *      Author: raiden
 */

#ifndef C1DGRID_H_
#define C1DGRID_H_

#include "CGrid.h"

namespace std
{

class C1DGrid : public CGrid
{
public:
	C1DGrid(double* range, int n, int Mx);
	virtual ~C1DGrid();
	virtual void compose();
	virtual void explicit_euler_evolution_operator();
	virtual void implicit_euler_evolution_operator();
	virtual void crank_nicholson_evolution_operator();
};

} /* namespace std */

#endif /* C2DGRID_H_ */
