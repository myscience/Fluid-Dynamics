#pragma once

#include "grid.h"
#include "levelset.h"

class Velocity
{
public:
	Velocity(int& w, int& h);
	Velocity(const Grid& u, const Grid& v);
	~Velocity();

	/* Assigment operator */
	Velocity& operator= (const Velocity& V);

	/* Boundary condition operator */
	void operator&= (const LevelSet& ls);

	/* Operation on velocity */
	Grid ndiv(const LevelSet& fls, const LevelSet& sls, const double& dx);
	void advect(const double& dt);

	friend class Solver;
public:
	int w, h;

	Grid u;
	Grid v;
};

