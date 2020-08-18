#pragma once
#define KEY_VX 0
#define KEY_VY 1
#define KEY_RHO 2

#include <vector>
#include <random>

#include "utils.h"
#include "levelset.h"

class Particles
{
public:
	Particles(const int& w, const int& h, const double& r, const LevelSet& ls);

	void init(const Grid& g, const uint8_t key);

	void advect(const double& dt, const Velocity& V, const LevelSet& fls, const LevelSet& sls);
	void fit(const Grid& gnew, const Grid& g, const uint8_t key);

	friend Grid;
	friend LevelSet;

protected:
	void adjust(const LevelSet& fls, const LevelSet& sls);

	std::vector<Part> P;

private:
	int w, h;

	const double r;

	const size_t maxPerCell = 12;
	const size_t minPerCell = 2;
	const size_t avgPerCell = 6;

	const double alpha = 0.01;

	std::default_random_engine rng;
	std::uniform_real_distribution<> rpos;
};

