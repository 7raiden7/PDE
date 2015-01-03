/*
 * C2DGrid.h
 *
 *  Created on: 25 Dec 2014
 *      Author: raiden
 */

#ifndef C2DGRID_H_
#define C2DGRID_H_

#include "CGrid.h"

namespace std
{

class C2DGrid : public CGrid
{
public:
	C2DGrid(double* range, int n, int Mx, int My);
	virtual ~C2DGrid();
	virtual void compose();
	virtual void explicit_euler_evolution_operator();
	virtual void implicit_euler_evolution_operator();
	virtual void crank_nicholson_evolution_operator();
};

} /* namespace std */

#endif /* C2DGRID_H_ */
