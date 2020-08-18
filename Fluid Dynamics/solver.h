#pragma once

#include "grid.h"
#include "levelset.h"
#include "velocity.h"

struct Stencil {
	int w, h;
	int size;

	double* diag  = nullptr;
	double* plusX = nullptr;
	double* plusY = nullptr;

	Stencil(const int& w, const int& h) : w(w), h(h), size(w * h), diag(new double[size]()),
		plusX(new double[size]()), plusY(new double[size]()) {};

	~Stencil() {
		if (diag  != nullptr) { delete[] diag;  diag  = nullptr; }
		if (plusX != nullptr) { delete[] plusX; plusX = nullptr; }
		if (plusY != nullptr) { delete[] plusY; plusY = nullptr; }
	}

	/* Operator* implements the matrix-vector product, where matrix has the
	   compressed Stencil representation where only diagonal and two bands 
	   around it are actually stored exploiting symmetry. */
	inline Grid operator* (const Grid& g) const {
		/* Here we prepare the output Grid */
		Grid out(w, h, 0.5, 0.5, 0.);

		for (int j = 0; j < h; j++) {
			for (int i = 0; i < w; i++) {
				int idx = j * w + i;

				out[idx] = diag[idx] * g[idx];
				if (i < w - 1) out[idx] += plusX[idx] * g[idx + 1];
				if (j < h - 1) out[idx] += plusY[idx] * g[idx + w];
				if (i > 0) out[idx] += plusX[idx - 1] * g[idx - 1];
				if (j > 0) out[idx] += plusY[idx - w] * g[idx - w];
			}
		}

		return out;
	}
};

class Solver
{
public:
	Solver(const int& w, const int& h, const double& dx, const double& dt);
	~Solver();

	Velocity project(const LevelSet& fls, const LevelSet& sls, Velocity& V);

public:
	void setA(const LevelSet& fls, const LevelSet& sls);
	void setP(const LevelSet& fls, const LevelSet& sls, const double& tau, const double& sigma);
	Grid applyP(const Grid& r, const LevelSet& sls, const LevelSet& fls);
	Grid PCG(const Grid& rhs, const LevelSet& fls, const LevelSet& sls, const double& eps);

	int w, h;

	double dx, dt;
	double rho;

	size_t maxIt = 1000;

	Stencil A;
	Grid P;
};

