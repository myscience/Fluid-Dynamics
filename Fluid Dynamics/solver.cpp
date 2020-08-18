#include "solver.h"

Solver::Solver(const int& w, const int& h, const double& dx, const double& dt) : 
	w(w), h(h), dx(dx), dt(dt), rho(1.), A(w, h), P(w, h, 0.5, 0.5, 0.)
{
}

Solver::~Solver()
{
}

Velocity Solver::project(const LevelSet& fls, const LevelSet& sls, Velocity& V)
{
	if (w != fls.getw() || w != sls.getw() || h != fls.geth() || h != sls.geth())
		throw std::invalid_argument("Solver project: Shape mismatch.");

	/* Create new Velocity for output */
	Velocity Vout(V);

	/* First we solve for pressure */
	setA(fls, sls);

	Grid ndiv = V.ndiv(fls, sls, dx);

	/* Code for adjustment for solid velocity goes here */
	// TODO: Add code for solid velocity

	Grid p = PCG(ndiv, fls, sls, 1e-5);

	/* Then we project velocity to a solenoidal field */
	double scale = dt / (rho * dx);

	/* Update u component */
	for (int j = 0; j < h; j++) {
		for (int i = 0; i < w + 1; i++) {
			int idx = j * (w + 1) + i;

			/* WARNING: Not so sure about the construct fls != MIS. However nieghter fis == IN
						Is correct as it mark as UNKWN solid cell whose normal velocity is instead
						fixed by the border condition. A possible issue that arises is that (due to
						extrapolation) in-solid velocities arises, thus breaking the simulation.
			*/
			if (fls(i - 1, j) == IN or fls(i, j) == IN) {
				if (sls(i - 1, j) == IN or sls(i, j) == IN) Vout.u[idx] = 0.;
				else Vout.u[idx] -= scale * (p[Id (i, j)] - p[Id (i - 1, j)]);
			}
			else Vout.u[idx] = UNKWN;
		}
	}

	for (int j = 0; j < h + 1; j++) {
		for (int i = 0; i < w; i++) {
			int idx = j * w + i;

			/* Update v component */
			if (fls(i, j - 1) == IN or fls(i, j) == IN) {
				if (sls(i, j - 1) == IN or sls(i, j) == IN) Vout.v[idx] = 0.;
				else Vout.v[idx] -= scale * (p[idx] - p[idx - w]);
			}
			else Vout.v[idx] = UNKWN;
		}
	}

	/* Fill unknwon velocity value with extrapolation */
	Vout.u.extrapolate();
	Vout.v.extrapolate();

	return Vout;
}

void Solver::setA(const LevelSet& fls, const LevelSet& sls)
{
	/* Here we populate the interation matrix (stored as a compressed Stencil) */
	double scale = dt / (rho * dx * dx);

	for (int j = 0; j < h; j++) {
		for (int i = 0; i < w; i++) {
			int idx = j * w + i;

			/* Here we clean the Stencil */
			A.diag[idx] = A.plusX[idx] = A.plusY[idx] = 0.;

			/* Only FLUID cell are relevant for interaction */
			if (fls(i, j) != IN) continue;

			Id h(i, j);
			
			/* Handling negative X */
			A.diag[idx] += (1. - sls.Fl(h)) * scale;

			/* Handling positive X */
			A.diag[idx]  += (1. - sls.Fr(h)) * scale;
			A.plusX[idx] = -(1. - sls.Fr(h)) * scale;

			/* Handling negative Y */
			A.diag[idx] += (1. - sls.Fb(h)) * scale;

			/* Handling positive Y */
			A.diag[idx]  += (1. - sls.Ft(h)) * scale;
			A.plusY[idx] = -(1. - sls.Ft(h)) * scale;
		}
	}
}

void Solver::setP(const LevelSet& fls, const LevelSet& sls, const double& tau, const double& sigma)
{
	for (int j = 0; j < h; j++) {
		for (int i = 0; i < w; i++) {
			int idx = j * w + i;

			if (fls(i, j) != IN) continue;

			double e = A.diag[idx];
			double pxx = A.plusX[idx - 1] * P[idx - 1];
			double pyy = A.plusY[idx - w] * P[idx - w];
			double pxy = A.plusX[idx - w] * P[idx - w];
			double pyx = A.plusY[idx - 1] * P[idx - 1];

			e -= (pxx * pxx + pyy * pyy + tau * (pxx * pyx + pyy * pxy));

			if (e < sigma * A.diag[idx]) e = A.diag[idx];

			P[idx] = 1. / std::sqrt(e);
		}
	}
}

Grid Solver::applyP(const Grid& r, const LevelSet& fls, const LevelSet& sls)
{
	Grid z(w, h, 0.5, 0.5, 0.);
	Grid q(w, h, 0.5, 0.5, 0.);

	for (int j = 0; j < h; j++) {
		for (int i = 0; i < w; i++) {
			int idx = j * w + i;

			/* We skip non-fluid cell */
			if (fls(i, j) != IN) continue;

			double t = r[idx] - A.plusX[idx - 1] * P[idx - 1] * q[idx - 1]
							  - A.plusY[idx - w] * P[idx - w] * q[idx - w];

			q[idx] = t * P[idx];
		}
	}

	for (int j = h - 1; j >= 0; j--) {
		for (int i = w - 1; i >= 0; i--) {
			int idx = j * w + i;

			/* We skip non-fluid cell */
			if (fls(i, j) != IN) continue;

			double t = q[idx] - A.plusX[idx] * P[idx] * z[idx + 1]
							  - A.plusY[idx] * P[idx] * z[idx + w];

			z[idx] = t * P[idx];
		}
	}

	return z;
}

Grid Solver::PCG(const Grid& rhs, const LevelSet& fls, const LevelSet& sls, const double& eps)
{
	/* We build the MIC(0) Preconditioner */
	setP(fls, sls, 0.97, 0.25);

	/* We prepare the solution and the residual vectors */
	Grid p(w, h, 0.5, 0.5, 0.);
	Grid r(rhs);

	/* Exit if rhs is trivial */
	if (r == 0) return p;

	/* We set auxiliary and search vectors */
	Grid z = applyP(r, fls, sls);	
	Grid s(z);

	/* Init sigma ad the dot-product */
	double sigma = (z, r);

	for (size_t k = 0; k < maxIt; k++) {
		z = A * s;

		double alpha = sigma / (z, s);

		/* Here we update our estimate for the pressure */
		p += s * alpha;
		r -= z * alpha;

		/* Check if computed pressure is close enough */
		if (r.norm_inf(fls) < eps) return p;

		z = applyP(r, fls, sls);

		double sigma_new = (z, r);
		double beta = sigma_new / sigma;

		s = z + s * beta;
		sigma = sigma_new;
	}

	/* Report exceeded maximum iterations */
	std::cerr << "PCG: Exceeded maximum iteration limit. Error: " << r.norm_inf(fls) << '\n';

	return p;
}