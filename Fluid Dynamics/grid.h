#pragma once
#include <stdexcept>
#include <string>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>
#include <assert.h>

#include "utils.h"

class Velocity;
class Particles;
class LevelSet;

class Grid
{
public:
	Grid(const int& w, const int& h, const double& offx, const double& offy);
	Grid(const int& w, const int& h, const double& offx, const double& offy, const double& v);
	Grid(const std::string filename);
	Grid(const Grid& g);
	~Grid();

	/* Assigment operator */
	Grid& operator= (const Grid& g);

	/* Operator definition for accessing grid elements */
	double& operator[] (int idx);
	double& operator[] (Id idx);
	double  operator[] (Vec v);

	const double& operator[] (int idx) const;
	const double& operator[] (Id idx) const;
	const double  operator[] (Vec v) const;

	const double cerp(const Vec& v) const;

	/* Comparison operations */
	bool operator== (const double& c) const;

	/* Operation on grid values */
	Grid operator+ (const Grid& g) const;
	Grid operator- (const Grid& g) const;
	Grid operator* (const double& c) const;
	double operator ,(const Grid& g) const;
	void operator+= (const Grid& g);
	void operator-= (const Grid& g);

	void fill(const double& v);
	void extrapolate();
	void advect(const double& dt, const Velocity& V);
	void fit(const Particles& P, const uint8_t& key);
	
	/* Grid to file interation: load and save */
	friend std::ostream& operator<< (std::ostream& stream, const Grid& grid);
	friend std::istream& operator>> (std::ifstream& stream, Grid& grid);

	/* Miscellaneous */
	std::vector<Id> neigh(int i, int j);
	std::vector<Id> neigh(Id idx);
	std::vector<Id> nneigh(Id idx);
	Vec locate(const Vec& v) const;
	Vec grad(const Vec& v) const;
	double norm_inf(const LevelSet& fls) const;
	int getw() const;
	int geth() const;
	
protected:	
	double k(const Vec& v);
	double spline(const double& x);

	int w, h;
	int size;

	Vec off;

	double* data = nullptr;
};

