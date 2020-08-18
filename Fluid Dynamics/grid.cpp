#include "grid.h"
#include "velocity.h"
#include "particles.h"
#include "levelset.h"

Grid::Grid(const int& w, const int& h, const double& offx, const double& offy) : w(w), h(h), size(w * h),
	off(offx, offy), data(new double[size]()) 
{
}

Grid::Grid(const int& w, const int& h, const double& offx, const double& offy, const double& v) : 
	w(w), h(h), size (w * h), off(offx, offy), data(new double[size]()) 
{
	fill(v);
}

Grid::Grid(const std::string filename) : w(0), h(0)
{
	/* Here we open the filestream */
	std::ifstream stream(filename);

	/* First line is assumed to contain integers specifying the domain */
	double offx, offy;
	char c;
	stream >> c >> w >> h >> offx >> offy;

	size = w * h;
	off.x = offx;
	off.y = offy;

	/* Here we allocate memory for data */
	data = new double[size]();

	/* Here we fill the new created memory space */
	for (int i = 0; i < size; i++)
		stream >> data[i];
}

Grid::Grid(const Grid& g) : w(g.w), h(g.h), size(w * h), off(g.off), data(new double[size]())
{
	if (g.data == nullptr) throw std::invalid_argument("Grid copy: invalid data pointer.");

	std::copy(g.data, g.data + size, data);
}

Grid::~Grid() {
	if (data != nullptr) {
		delete[] data;
		data = nullptr;
	}
}

Grid& Grid::operator=(const Grid& g)
{
	if (w != g.w || h != g.h) throw std::invalid_argument("Operator=: Shape mismatch.");

	std::copy(g.data, g.data + g.size, data);

	return *this;
}

double& Grid::operator[](int idx)
{
	if (idx >= 0 && idx < size) return data[idx];
	else throw std::invalid_argument("Operator[]: idx " + std::to_string(idx) +
		" is out of bounds for size " + std::to_string(size));
}

double& Grid::operator[](Id idx)
{
	if (idx.i >= 0 && idx.i < w && idx.j >= 0 && idx.j < h) return data[idx.j * w + idx.i];
	else throw std::invalid_argument("Operator[]: Idx " + std::to_string(idx.i) + " " + std::to_string(idx.j) +
		" is out of bounds for dim: " + std::to_string(w) + " " + std::to_string(h));
}

double Grid::operator[](Vec v)
{
	/* Here a value not on grid is requested. We thus use lerping. */
	if (v.x < 0 || v.x > w || v.y < 0 || v.y > h) throw std::invalid_argument("Operator[]: Vec " +
		std::to_string(v.x) + " " + std::to_string(v.y) + " is out of bounds");

	/* We remove the offset to account for staggering */
	double x = std::min (std::max (v.x - off.x, 0.), w - 1e-10);
	double y = std::min (std::max (v.y - off.y, 0.), h - 1e-10);

	/* Here we extract fractional part */
	int ix = (int) x;
	int iy = (int) y;
	x -= ix;
	y -= iy;

	/* Here we grab the closest known values */
	double x00 = data[iy * w + ix];
	double x10 = data[iy * w + ix + 1];
	double x01 = data[(iy + 1) * w + ix];
	double x11 = data[(iy + 1) * w + ix + 1];

	/* Finally we return lerped value */
	return lerp(lerp(x00, x10, x), lerp(x01, x11, x), y);
}

const double& Grid::operator[](int idx) const
{
	if (idx >= 0 && idx < size) return data[idx];
	else throw std::invalid_argument("Operator[]: idx " + std::to_string(idx) +
		" is out of bounds for size " + std::to_string(size));
}

const double& Grid::operator[](Id idx) const
{
	if (idx.i >= 0 && idx.i < w && idx.j >= 0 && idx.j < h) return data[idx.j * w + idx.i];
	else throw std::invalid_argument("Operator[]: Idx " + std::to_string(idx.i) + " " + std::to_string(idx.j) +
		" is out of bounds for dim: " + std::to_string(w) + " " + std::to_string(h));
}

const double Grid::operator[](Vec v) const
{
	/* Here are value not on grid is requested. We thus use lerping. */
	if (v.x < 0 || v.x > w || v.y < 0 || v.y > h) throw std::invalid_argument("Operator[]: Vec " +
		std::to_string(v.x) + " " + std::to_string(v.y) + " is out of bounds");

	/* We remove the offset to account for staggering */
	double x = std::min(std::max(v.x - off.x, 0.), w - 1. - 1e-10);
	double y = std::min(std::max(v.y - off.y, 0.), h - 1. - 1e-10);

	/* Here we extract fractional part */
	int ix = (int)x;
	int iy = (int)y;
	x -= ix;
	y -= iy;

	/* Here we grab the closest known values */
	double x00 = data[iy * w + ix];
	double x10 = data[iy * w + ix + 1];
	double x01 = data[(iy + 1) * w + ix];
	double x11 = data[(iy + 1) * w + ix + 1];

	/* Finally we return lerped value */
	return lerp(lerp(x00, x10, x), lerp(x01, x11, x), y);
}

const double Grid::cerp(const Vec& v) const
{
	/* Here are value not on grid is requested. We thus use cerping. */
	if (v.x < 0 || v.x > w || v.y < 0 || v.y > h) throw std::invalid_argument("Operator[]: Vec " +
		std::to_string(v.x) + " " + std::to_string(v.y) + " is out of bounds");

	int ix = (int) (v.x - off.x);
	int iy = (int) (v.y - off.y);

	double x = v.x - off.x - ix;
	double y = v.y - off.y - iy;

	/* Here we collect the id for the four neighbours */
	int x0 = std::max(ix - 1, 0), x1 = ix, x2 = std::min(ix + 1, w - 1), x3 = std::min(ix + 2, w - 1);
	int y0 = std::max(iy - 1, 0), y1 = iy, y2 = std::min(iy + 1, h - 1), y3 = std::min(iy + 2, h - 1);

	double q0 = ucerp(data[x0 + y0 * w], data[x1 + y0 * w], data[x2 + y0 * w], data[x3 + y0 * w], x);
	double q1 = ucerp(data[x0 + y1 * w], data[x1 + y1 * w], data[x2 + y1 * w], data[x3 + y1 * w], x);
	double q2 = ucerp(data[x0 + y2 * w], data[x1 + y2 * w], data[x2 + y2 * w], data[x3 + y2 * w], x);
	double q3 = ucerp(data[x0 + y3 * w], data[x1 + y3 * w], data[x2 + y3 * w], data[x3 + y3 * w], x);

	return ucerp(q0, q1, q2, q3, y);
}

bool Grid::operator==(const double& c) const
{
	for (int i = 0; i < size; i++) if (data[i] != c) return false;

	return true;
}

Grid Grid::operator+(const Grid& g) const
{
	if (w != g.w || h != g.h) throw std::invalid_argument("Grid Operator+: Shape mismatch.");

	Grid out(*this);

	for (int i = 0; i < size; i++) out[i] += g[i];

	return out;
}

Grid Grid::operator-(const Grid& g) const
{
	if (w != g.w || h != g.h) throw std::invalid_argument("Grid Operator+: Shape mismatch.");

	Grid out(*this);

	for (int i = 0; i < size; i++) out[i] -= g[i];

	return out;
}

Grid Grid::operator*(const double& c) const
{
	Grid out(*this);

	for (int i = 0; i < size; i++) out[i] *= c;

	return out;
}

double Grid::operator,(const Grid& g) const
{
	if (w != g.w || h != g.h) throw std::invalid_argument("Dot Operator. Shape mismatch.");

	double out = 0.;
	for (int i = 0; i < size; i++) out += data[i] * g[i];

	return out;
}

void Grid::operator+=(const Grid& g)
{
	if (w != g.w || h != g.h) throw std::invalid_argument("Grid Operator+=: Shape mismatch.");

	for (int i = 0; i < size; i++) data[i] += g[i];
	return;
}

void Grid::operator-=(const Grid& g)
{
	if (w != g.w || h != g.h) throw std::invalid_argument("Grid Operator+=: Shape mismatch.");

	for (int i = 0; i < size; i++) data[i] -= g[i];
	return;
}

void Grid::fill(const double& v)
{
	std::fill_n(data, size, v);
}

std::vector<Id> Grid::neigh(int i, int j)
{
	if (i < 0 || i >= w || j < 0 || j >= h) throw std::invalid_argument("In Neigh. Invalid Idx: out of domain");
	
	/* Here we prepare a vector to store up to 4 neighbours */
	std::vector<Id> neigh;
	neigh.reserve(4);

	/* Here we add based on border constraints */
	if (i > 0) neigh.push_back(Id(i - 1, j));
	if (j > 0) neigh.push_back(Id(i, j - 1));
	if (i < w - 1) neigh.push_back(Id(i + 1, j));
	if (j < h - 1) neigh.push_back(Id(i, j + 1));

	return neigh;
}

std::vector<Id> Grid::neigh(Id idx)
{
	return neigh(idx.i, idx.j);
}

std::vector<Id> Grid::nneigh(Id idx)
{
	int i = idx.i;
	int j = idx.j;

	if (i < 0 || i >= w || j < 0 || j >= h) throw std::invalid_argument("In Neigh. Invalid Idx: out of domain");

	/* Here we prepare a vector to store up to 4 neighbours */
	std::vector<Id> neigh;
	neigh.reserve(12);

	/* Here we add based on border constraints */
	if (i > 0) neigh.push_back(Id(i - 1, j));
	if (j > 0) neigh.push_back(Id(i, j - 1));
	if (i < w - 1) neigh.push_back(Id(i + 1, j));
	if (j < h - 1) neigh.push_back(Id(i, j + 1));
	
	if (i > 1) neigh.push_back(Id(i - 2, j));
	if (j > 1) neigh.push_back(Id(i, j - 2));
	if (i < w - 2) neigh.push_back(Id(i + 2, j));
	if (j < h - 2) neigh.push_back(Id(i, j + 2));

	if (i > 0 && j > 0) neigh.push_back(Id(i - 1, j - 1));
	if (i > 0 && j < h - 1) neigh.push_back(Id(i - 1, j + 1));
	if (i < w - 1 && j > 0) neigh.push_back(Id(i + 1, j - 1));
	if (i < w - 1 && j < h - 1) neigh.push_back(Id(i + 1, j + 1));

	return neigh;
}

Vec Grid::locate(const Vec& v) const
{
	if (v.x < 0 || v.x > w - 1e-10 || v.y < 0 || v.y > h - 1e-10)
		throw std::invalid_argument("Locate fail: Vec out of bounds");

	return Vec(std::round(std::max(v.x - off.x, 0.)), std::round(std::max(v.y - off.y, 0.)));
}

Vec Grid::grad(const Vec& v) const
{
	/* Here are value not on grid is requested. We thus use lerping. */
	if (v.x < 0 || v.x > w || v.y < 0 || v.y > h) throw std::invalid_argument("Operator[]: Vec " +
		std::to_string(v.x) + " " + std::to_string(v.y) + " is out of bounds");

	/* We remove the offset to account for staggering */
	double x = std::min(std::max(v.x - off.x, 0.), w - 1. - 1e-10);
	double y = std::min(std::max(v.y - off.y, 0.), h - 1. - 1e-10);

	/* Here we extract fractional part */
	int ix = (int)x;
	int iy = (int)y;
	x -= ix;
	y -= iy;

	double x00 = data[iy * w + ix], x10 = data[iy * w + ix + 1];
	double x01 = data[(iy + 1) * w + ix], x11 = data[(iy + 1) * w + ix + 1];

	/* Here we compute the Gradient */
	double gx0 = (x10 - x00), gx1 = (x11 - x01);
	double gy0 = (x01 - x00), gy1 = (x11 - x10);

	return Vec(lerp(gx0, gx1, y), lerp(gy0, gy1, x));
}

double Grid::norm_inf(const LevelSet& fls) const
{
	if (fls.geth() != h || fls.getw() != w) throw std::invalid_argument("Norm Inf: Shape Mismatch.\n");

	double max = -1e50;

	for (int i = 0; i < size; i++)
		if (fls[i] < 0 && std::fabs(data[i]) > max) max = std::fabs(data[i]);

	return max;
}

int Grid::getw() const
{
	return w;
}

int Grid::geth() const
{
	return h;
}

double Grid::k(const Vec& v)
{
	return spline(v.x) * spline(v.y);
}

double Grid::spline(const double& r)
{
	if (r >= -1.5 && r < -0.5) return 0.5 * (r + 1.5) * (r + 1.5);
	if (r >= -0.5 && r < 0.5) return 1.5 - r * r;
	if (r >= 0.5 && r < 1.5) return 0.5 * (1.5 - r) * (1.5 - r);
	return 0.;
}

void Grid::extrapolate()
{
	/* We use breadth first search for extrapolation */
	/* We use an helper marker array to track updates */
	Grid m(w, h, 0.5, 0.5, 0.);

	/* We use a wavefront vector to collect the next-to-update values */
	std::vector<Id> W;
	W.reserve(w * h / 4);

	/* Initialization of wavefront vector */
	for (int j = 0; j < h; j++) {
		for (int i = 0; i < w; i++) {
			/* Check if extrapolation is needed */
			if (data[j * w + i] != UNKWN) continue;

			/* Set marker array to unknown for this index */
			m[j * w + i] = UNKWN;

			/* Check if wavefront: at least one neigh known */
			bool flag = false;
			for (Id& idx : neigh(i, j)) flag |= data[idx.j * w + idx.i] != UNKWN;

			/* Update wavefront tracking */
			if (flag) { m[j * w + i] = 1; W.push_back(Id(i, j)); }
		}
	}

	/* Breadth first extrapolation */
	size_t t = 0;
	while (t < W.size()) {
		Id idx = W[t];

		double F = 0.;
		size_t N = 0;
		for (Id& nidx : neigh(idx)) {
			/* Compute average of reliable neighbours */
			if (m[nidx] < m[idx]) { F += data[nidx.j * w + nidx.i]; N++; }

			/* Move the wavefront */
			if (m[nidx] == UNKWN) { m[nidx] = m[idx] + 1; W.push_back(nidx); }
		}

		/* Set value to avg of reliable neighbours */
		assert (N != 0);
		data[idx.j * w + idx.i] = F / N;

		t++;
	}
}

void Grid::advect(const double& dt, const Velocity& V)
{
	/* Build temporary grid where to store advection */
	Grid tmp(w, h, off.x, off.y, 0.);

	for (int j = 0; j < h; j++) {
		for (int i = 0; i < w; i++) {
			/* Look where this point originates due to advection */
			Vec p = RK3(-dt, Vec(i, j) + off, V);

			/* Interpolates new value based on previous position */
			tmp[j * w + i] = cerp(p);
		}
	}

	/* Here we copy the advected grid to current data */
	(*this) = tmp;
}

void Grid::fit(const Particles& P, const uint8_t& key)
{
	fill(0.);

	/* Prepare the weight helper grid */
	Grid W(w, h, off.x, off.y, 0.);

	for (const auto& p : P.P) {
		/* Get neighbours around position of particle p */
		Vec v = locate(p.p);
		Id idx = Id((int)v.x, (int)v.y);

		/* Add contribution to current location */
		data[idx.j * w + idx.i] += p[key] * k(p.p - v);
		W[idx] += k(p.p - v);

		/* Add contribution to neighbours */
		for (Id nidx : nneigh(idx)) {
			Vec nv(nidx.i, nidx.j);
			data[nidx.j * w + nidx.i] += p[key] * k(p.p - nv);
			W[nidx] += k(p.p - nv);
		}
	}

	/* Here we devide by the weights */
	for (int j = 0; j < h; j++) {
		for (int i = 0; i < w; i++) {
			int idx = j * w + i;

			/* Update values where fitting toke place */
			if (W[idx] > 0.) data[idx] /= W[idx];
		}
	}
}

std::ostream& operator<<(std::ostream& stream, const Grid& grid)
{
	for (int j = 0; j < grid.h; j++) {
		for (int i = 0; i < grid.w; i++) 
			stream << grid.data[j * grid.w + i] << ' ';
		
		stream << std::endl;
	}
	
	return stream;
}

std::istream& operator>>(std::ifstream& stream, Grid& grid)
{
	/* First line is assumed to contain integers specifying the domain */
	int w, h;
	char c;
	stream >> c >> w >> h;

	if (w != grid.w || h != grid.h) throw std::invalid_argument("Shape mismatch in loading Grid from file");

	/* Here we fill the new created memory space */
	for (int j = 0; j < grid.h; j++)
		for (int i = 0; i < grid.w; i++)
			stream >> grid.data[j * grid.w + i];

	return stream;
}
