#include "velocity.h"

Velocity::Velocity(int& w, int& h) : w(w), h(h), u(w + 1, h, 0., 0.5, 0.), v(w, h + 1, 0.5, 0., 0.)
{
}

Velocity::Velocity(const Grid& u, const Grid& v) : w(v.getw()), h(u.geth()), u(u), v(v)
{
}

Velocity::~Velocity()
{
}

Velocity& Velocity::operator=(const Velocity& V)
{
	u = V.u;
	v = V.v;

	return *this;
}

void Velocity::operator&=(const LevelSet& sls)
{
	for (int j = 0; j < h; j++)
		for (int i = 0; i < w + 1; i++)
			if (sls(i - 1, j) == IN or sls(i, j) == IN) u[j * (w + 1) + i] = 0.;

	for (int j = 0; j < h + 1; j++)
		for (int i = 0; i < w; i++)
			if (sls(i, j - 1) == IN or sls(i, j) == IN) v[j * w + i] = 0.;
}

Grid Velocity::ndiv(const LevelSet& fls, const LevelSet& sls, const double& dx)
{
	/* Here we prepare the negative divergence Grid */
	Grid ndiv(w, h, 0.5, 0.5, 0.);
	
	double invdx = 1. / dx;

	for (int j = 0; j < h; j++) {
		for (int i = 0; i < w; i++) {
			/* We exclude non-fluid region from computation */
			if (fls(i, j) != IN) continue;
			
			Id idx(i, j);
			Id r(i + 1, j);
			Id t(i, j + 1);

			ndiv[idx] = -invdx * ((1. - sls.Fr(idx)) * u[r] - (1. - sls.Fl(idx)) * u[idx] + 
								  (1. - sls.Ft(idx)) * v[t] - (1. - sls.Fb(idx)) * v[idx]);
			//ndiv[idx] = -invdx * (u[e] - u[idx] + v[s] - v[idx]);
		}
	}

	return ndiv;
}

void Velocity::advect(const double& dt)
{
	/* Here we advect velocity field via itself */
	Velocity tmp(u, v);

	u.advect(dt, tmp);
	v.advect(dt, tmp);
}
