#include "levelset.h"
#include "particles.h"

LevelSet::LevelSet(int& w, int& h) : Grid(w, h, 0.5, 0.5) 
{
}

LevelSet::LevelSet(std::string filename) : Grid(filename)
{
}

LevelSet::~LevelSet() {}

uint8_t LevelSet::operator()(int i, int j) const
{
	if (i < 0 || i > w - 1 || j < 0 || j > h - 1) return MIS;
	else return data[j * w + i] < 0 ? IN : OUT;
}
uint8_t LevelSet::operator()(Id idx) const
{
	return (*this)(idx.i, idx.j);
}

LevelSet LevelSet::operator&(const LevelSet& mask) const
{
	if (mask.getw() != w && mask.geth() != h) 
		throw std::invalid_argument("LevelSet: operator& Shape Mismatch.\n");

	LevelSet out(*this);

	for (int i = 0; i < size; i++)
		if (out[i] < 0 && mask[i] < 0) out[i] = -mask[i];

	return out;
}

void LevelSet::fit(const Particles& Part)
{
	/* Here we construct Signed Distance function from set of points */
	/* We start with very high grid values */
	fill(1e50);

	/* We create a custom closest distance tracker grid */
	Grid t(w, h, 0.5, 0.5, -1);

	/* Here we init the distance close to provided points */
	for (size_t e = 0; e < Part.P.size(); e++) {
		Vec p = Part.P[e].p;
		Vec x = locate(p);
		Id idx = Id((int) x.x, (int) x.y);

		double d = dist(x, p) - Part.P[e].r;
		if (d < (*this)[idx]) {
			(*this)[idx] = d;
			t[idx] = (double) e;
		}
	}

	/* Lambda function to ease fast sweep */
	auto update = [&](int i, int j) {
		/* Here we iterate through the neibours */
		for (Id idx : neigh(i, j)) {
			int e = (int) t[idx];
			double d = e != -1 ? dist(Vec(i, j),  Part.P[e].p - off) - Part.P[e].r : 1e10;
			if (e != -1 && d < data[j * w + i]) {
				data[j * w + i] = d;
				t[j * w + i] = e;
			}
		}
	};

	/* Here we propagate the distance information via fast sweep */
	for (int r = 0; r < 2; r++) {
		/* I ascending */
		for (int j = 0; j < h; j++) {
			for (int i = 0; i < w; i++) {
				update(i, j);
			}
		}

		/* J ascending */
		for (int i = 0; i < w; i++) {
			for (int j = 0; j < h; j++) {
				update(i, j);
			}
		}

		/* I descending */
		for (int j = h - 1; j >= 0; j--) {
			for (int i = w - 1; i >= 0; i--) {
				update(i, j);
			}
		}

		/* J descending */
		for (int i = w - 1; i >= 0; i--) {
			for (int j = h - 1; j >= 0; j--) {
				update(i, j);
			}
		}
	}

	/* Here we subtract particles radius to complete signed distance */
	for (int idx = 0; idx < size; idx++) data[idx] -= Part.r;

	/* Here we redistance the levelset and close small holes */
	fill_holes(2);
}

Vec LevelSet::project(const Part& P, const double& eps) const
{
	Vec p(P.p);
	Vec d = grad(p).normalize();

	double phip = (*this)[p];
	
	for (int n = 0; n < 5; n++) {
		double alpha = 1.;

		for (int m = 0; m < 5; m++) {
			Vec q = p - d * alpha * phip;
			
			if (abs((*this)[q]) < abs((*this)[p])) {
				p = q;
				phip = (*this)[q];
				d = grad(q).normalize();

				if (abs(phip) < eps) return p + d * P.r;
			}
			else alpha *= 0.7;
		}
	}

	std::cerr << "Warning: Unaccuratre projection to closest point on surface.\n";
	return p;
}

double LevelSet::F(const Id& idx) const
{
	double phi0 = (*this)[Vec(idx.i, idx.j)];
	double phi1 = (*this)[Vec(idx.i, idx.j + 1)];
	double phi2 = (*this)[Vec(idx.i + 1, idx.j)];
	double phi3 = (*this)[Vec(idx.i + 1, idx.j + 1)];
	
	double volume = 1. - occupancy(phi0, phi1, phi2, phi3);

	return volume < 0.01 ? 0. : volume;
}

double LevelSet::frac(const double& a, const double& b) const
{
	double min = a < b ? a : b;
	double max = a < b ? b : a;

	if (max <= 0.) return 1.;
	if (min >= 0.) return 0.;

	double F = max / (max - min);

	return F < 0.1 ? 0. : F;
}

double LevelSet::Fl(const Id& idx) const
{
	//std::cerr << "Fl called! With: " << idx.i << " " << idx.j << '\n';

	double philu = (*this)[Vec(idx.i, idx.j + 1)];
	double phild = (*this)[Vec(idx.i, idx.j)];

	//std::cerr << "Lerp: " << philu << " " << phild << '\n';

	return frac(philu, phild);
}

double LevelSet::Fr(const Id& idx) const
{
	//std::cerr << "Fr called! With: " << idx.i << " " << idx.j << '\n';

	double phiru = (*this)[Vec(idx.i + 1, idx.j + 1)];
	double phird = (*this)[Vec(idx.i + 1, idx.j)];

	//std::cerr << "Lerp: " << phiru << " " << phird << '\n';

	return frac(phiru, phird);
}

double LevelSet::Fb(const Id& idx) const
{
	//std::cerr << "Fb called! With: " << idx.i << " " << idx.j << '\n';

	double phibl = (*this)[Vec(idx.i, idx.j)];
	double phibr = (*this)[Vec(idx.i + 1, idx.j)];

	//std::cerr << "Lerp: " << phibl << " " << phibr << '\n';

	return frac(phibl, phibr);
}

double LevelSet::Ft(const Id& idx) const
{
	//std::cerr << "Ft called! With: " << idx.i << " " << idx.j << '\n';

	double phitl = (*this)[Vec(idx.i, idx.j + 1)];
	double phitr = (*this)[Vec(idx.i + 1, idx.j + 1)];

	//std::cerr << "Lerp: " << phitl << " " << phitr << '\n';

	return frac(phitl, phitr);
}

double LevelSet::frac3(double out1, double in, double out2) const
{
	return 0.5 * in * in / ((out1 - in) * (out2 - in));
}

double LevelSet::frac4(double out1, double out2, double in1, double in2) const
{
	return 0.5 * (-in1 / (out1 - in1) - in2 / (out2 - in2));
}

double LevelSet::occupancy(double d11, double d12, double d21, double d22) const
{
	double ds[] = { d11, d12, d22, d21 };

	/* Compute mask */
	uint8_t b = 0;
	for (int i = 3; i >= 0; i--) b = (b << 1) | (ds[i] < 0.0 ? 1 : 0);

	switch (b) 
	{
		/* All outside */
		case 0x0: return 0.0;

		/* One inside */
		case 0x1: return frac3(d21, d11, d12);
		case 0x2: return frac3(d11, d12, d22);
		case 0x4: return frac3(d12, d22, d21);
		case 0x8: return frac3(d22, d21, d11);
		
		/* One outside */
		case 0xE: return 1.0 - frac3(-d21, -d11, -d12);
		case 0xD: return 1.0 - frac3(-d11, -d12, -d22);
		case 0xB: return 1.0 - frac3(-d12, -d22, -d21);
		case 0x7: return 1.0 - frac3(-d22, -d21, -d11);
		
		/* Two adjacent inside */
		case 0x3: return frac4(d21, d22, d11, d12);
		case 0x6: return frac4(d11, d21, d12, d22);
		case 0x9: return frac4(d12, d22, d11, d21);
		case 0xC: return frac4(d11, d12, d21, d22);
		
		/* Two opposed inside */
		case 0x5: return frac3(d11, d12, d22) + frac3(d22, d21, d11);
		case 0xA: return frac3(d21, d11, d12) + frac3(d12, d22, d21);
		
		/* All inside */
		case 0xF: return 1.0;
	}

	return 0.0;
}

void LevelSet::redistance(const double& dx)
{
	/* We create a custom closest distance tracker grid */
	Grid t(w, h, 0.5, 0.5, 0);

	/* Here we init the distance close to surface */
	for (int j = 0; j < h; j++) {
		for (int i = 0; i < w; i++) {
			int idx = j * w + i;

			if (i > 0 && sgn(data[idx]) != sgn(data[idx - 1])) { t[idx] = sgn(data[idx]); t[idx - 1] = sgn(data[idx - 1]); }
			if (j > 0 && sgn(data[idx]) != sgn(data[idx - w])) { t[idx] = sgn(data[idx]); t[idx - w] = sgn(data[idx - w]); }
			if (i < w - 1 && sgn(data[idx]) != sgn(data[idx + 1])) { t[idx] = sgn(data[idx]); t[idx + 1] = sgn(data[idx + 1]); }
			if (j < h - 1 && sgn(data[idx]) != sgn(data[idx + w])) { t[idx] = sgn(data[idx]); t[idx + w] = sgn(data[idx + w]); }
		}
	}

	for (int e = 0; e < size; e++) if (t[e] == 0) { t[e] = sgn(data[e]); data[e] = t[e] * 1e50; }

	/* Lambda function to ease fast sweep */
	auto update = [&](int i, int j) {
		/* Here we iterate through the neibours */
		double phix = std::min(i > 0 ? abs(data[j * w + i - 1]) : 1e10, i < w - 1 ? abs(data[j * w + i + 1]) : 1e10);
		double phiy = std::min(j > 0 ? abs(data[(j - 1) * w + i]) : 1e10, j < h - 1 ? abs(data[(j + 1) * w + i]) : 1e10);

		double phi[] = { phix, phiy };
		std::sort(std::begin(phi), std::end(phi));

		/* Redistancing based on Eikonal equation */
		double d = phi[0] + dx;

		if (d > phi[1]) d = 0.5 * (phi[0] + phi[1] + std::sqrt(2 * dx * dx - (phi[1] - phi[0]) * (phi[1] - phi[0])));
		if (d < abs(data[j * w + i])) data[j * w + i] = d;
	};

	/* Here we propagate the distance information via fast sweep */
	for (int r = 0; r < 2; r++) {
		/* I ascending */
		for (int j = 0; j < h; j++) {
			for (int i = 0; i < w; i++) {
				update(i, j);
			}
		}

		/* J ascending */
		for (int i = 0; i < w; i++) {
			for (int j = 0; j < h; j++) {
				update(i, j);
			}
		}

		/* I descending */
		for (int j = h - 1; j >= 0; j--) {
			for (int i = w - 1; i >= 0; i--) {
				update(i, j);
			}
		}

		/* J descending */
		for (int i = w - 1; i >= 0; i--) {
			for (int j = h - 1; j >= 0; j--) {
				update(i, j);
			}
		}
	}

	/* Here we restore the correct sign */
	for (int idx = 0; idx < size; idx++) data[idx] = t[idx] * abs(data[idx]);
}

void LevelSet::fill_holes(size_t R)
{
	for (size_t r = 0; r < R; r++) {
		Grid Phi(*this);

		for (int j = 0; j < h; j++) {
			for (int i = 0; i < w; i++) {
				double phi_avg = 0.;
				int N = 0;

				for (const Id nidx : neigh(Id(i, j))) {
					phi_avg += Phi[nidx];
					N++;
				}

				phi_avg /= N;

				double phi_old = data[j * w + i];
				data[j * w + i] = phi_avg < phi_old ? phi_avg : phi_old;
			}
		}
	}
}


unsigned char* operator<<(unsigned char* img, const LevelSet& ls)
{
	for(int j = 0; j < ls.h; j++) {
		for (int i = 0; i < ls.w; i++) {
			int idx_rgba = 4 * (i + j * ls.w);

			double shade = ls.F(Id(i, j));
			//double shade = 1. - (ls[i + j * ls.w] < 0);

			shade = std::min(std::max(shade, 0.0), 1.0);
			img[idx_rgba + 0] = (int)(shade * 255.0);
			img[idx_rgba + 1] = (int)(shade * 255.0);
			img[idx_rgba + 2] = (int)(shade * 255.0);
			img[idx_rgba + 3] = 0xFF;
		}
	}

	return img;
}
