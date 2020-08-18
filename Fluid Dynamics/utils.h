#pragma once
#include <algorithm>
#include <array>
#include <cmath>

#define UNKWN INT_MAX

class Velocity;

struct Id 
{
	int i, j;

	Id() : i(0), j(0) {};
	Id(int i, int j) : i(i), j(j) {};
	Id(const Id& id) : i(id.i), j(id.j) {};
};

struct Vec 
{
	double x, y;

	Vec() : x(0), y(0) {};
	Vec(double _x, double _y) : x(_x), y(_y) {};
	Vec(const Vec& v) : x(v.x), y(v.y) {};

	/* Operation on vectors */
	inline Vec operator+ (const Vec& v) const { return Vec(x + v.x, y + v.y); };
	inline Vec operator- (const Vec& v) const { return Vec(x - v.x, y - v.y); };
	inline Vec operator* (const double& c) const { return Vec(x * c, y * c); };

	inline double norm() const { return std::sqrt(x * x + y * y); };

	inline Vec normalize() const { double n = norm(); return Vec(x / n, y / n); };
};

struct Part
{
	Vec p;

	double r;
	std::array<double, 3> v;

	Part() : p(), r(1.), v({ 0, 0, 0 }) {};
	Part(Vec p, double r, std::array<double, 3> v) : p(p), r(r), v(v) {};
	Part(const Part& P) : p(P.p), r(P.r), v(P.v) {};

	double& operator[] (const int& key) { return v[key]; }
	double operator[] (const int& key) const { return v[key]; }
};

inline double lerp(double a, double b, double x)
{
	return a * (1. - x) + b * x;
}

inline double ucerp(double a, double b, double c, double d, double x)
{
	double xsq = x * x;
	double xcu = x * x * x;

	double vmin = std::min(a, std::min(b, std::min(c, d)));
	double vmax = std::max(a, std::max(b, std::max(c, d)));

	double t = a * (0.0 - 0.5 * x + 1.0 * xsq - 0.5 * xcu) +
			   b * (1.0 + 0.0 * x - 2.5 * xsq + 1.5 * xcu) +
			   c * (0.0 + 0.5 * x + 2.0 * xsq - 1.5 * xcu) +
			   d * (0.0 + 0.0 * x - 0.5 * xsq + 0.5 * xcu);

	return std::min(std::max(t, vmin), vmax);
}

inline double dist(Vec a, Vec b)
{
	double x = a.x - b.x;
	double y = a.y - b.y;

	return std::sqrt(x * x + y * y);
}

inline int sgn(double val) {
	return (double(0) < val) - (val < double(0));
}

Vec RK3(const double& dt, const Vec& v, const Velocity& V);