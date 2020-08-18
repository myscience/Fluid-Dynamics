#include "particles.h"

Particles::Particles(const int& w, const int& h, const double& r, const LevelSet& ls) : 
	w(w), h(h), r(r), rng(), rpos(0.1, .9)
{
	P.reserve(w * h * avgPerCell);

	/* Here we populate the particles vector in fluid regions */
	for (int j = 0; j < h; j++) {
		for (int i = 0; i < w; i++) {
			if (ls(i, j) != IN) continue;

			for (int c = 0; c < avgPerCell; c++) {
				double x = i + rpos(rng);
				double y = j + rpos(rng);

				if (x >= 0 && x < w && y >= 0 && y < h && 
					ls[Vec(x, y)] < 0) P.push_back(Part(Vec(x, y), r, std::array<double, 3>()));
				else c--;
			}
		}
	}
}

void Particles::init(const Grid& g, const uint8_t key)
{
	for (auto& p : P) p[key] = g[p.p];
}

void Particles::advect(const double& dt, const Velocity& V, const LevelSet& fls, const LevelSet& sls)
{
	#pragma omp parallel for
	for (auto& p : P) {
		p.p = RK3(dt, p.p, V);

		/* Back project particles landed inside solid boundaries */
		if (sls[p.p] < 0) p.p = sls.project(p, 1e-3);
		assert(sls[p.p] > 0);
	}

	//adjust(fls, sls);
}

void Particles::adjust(const LevelSet& fls, const LevelSet& sls)
{
	/* Here we remove all particles that landed outside fluid */
	for (std::vector<Part>::iterator it = P.begin(); it != P.end();) {
		try { 
			if (fls[it->p] < 0) ++it;
			else if (sls[it->p] < 0) it->p = sls.project(*it, 1e-3);
			else it = P.erase(it);
		}

		/* Remove point landing outside domain */
		catch (std::invalid_argument) {
			it = P.erase(it);
		}
	}

	/* Here we re-seed to avoid too-few particles in cell */
	// TODO: Add reseeding
}

void Particles::fit(const Grid& gnew, const Grid& g, const uint8_t key)
{
	Grid Delta = gnew - g;

	for (auto& p : P) {
		p[key] = alpha * gnew[p.p] + (1. - alpha) * (p[key] + Delta[p.p]);
	}
}
