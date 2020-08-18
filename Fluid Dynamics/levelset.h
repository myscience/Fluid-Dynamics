#pragma once
#include <vector>

#include "grid.h"

#define IN 0
#define OUT 1
#define MIS 2

class Particles;

class LevelSet : public Grid
{
public:
	LevelSet(int& w, int& h);
	LevelSet(std::string filename);
	~LevelSet();

	/* Operator for accessing grid label */
	uint8_t operator() (int i, int j) const;
	uint8_t operator() (Id idx) const;

	double F(const Id& idx) const;
	double Fl(const Id& idx) const;
	double Fr(const Id& idx) const;
	double Ft(const Id& idx) const;
	double Fb(const Id& idx) const;


	/* Operator for masking a LevelSet */
	LevelSet operator& (const LevelSet& mask) const;

	void fit(const Particles& P);
	void redistance(const double& dx);

	Vec project(const Part& p, const double& eps) const;

	friend unsigned char* operator<< (unsigned char* img, const LevelSet& ls);

private:
	double frac3(double a, double b, double c) const;
	double frac4(double a, double b, double c, double d) const;
	double occupancy(double phi0, double phi1, double phi2, double phi3) const;

	void fill_holes(size_t R);

	double frac(const double& a, const double& b) const;
};

