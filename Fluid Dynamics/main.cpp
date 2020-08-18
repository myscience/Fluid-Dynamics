#include <iostream>
#include <iomanip>
#include <sstream>

#include "grid.h"
#include "levelset.h"
#include "solver.h"
#include "velocity.h"
#include "lodepng.h"
#include "particles.h"

int main() {
	const int w = 126;
	const int h = 126;

	const double dx = 1.;
	const double dt = 0.1;

	const double r = 0.3;

	unsigned char* img = new unsigned char[w * h * 4];

	std::ofstream out("..\\Debug\\debug.txt");

	LevelSet fls(std::string("..\\res\\fluid.txt"));
	LevelSet sls(std::string("..\\res\\solid.txt"));

	fls = fls & sls;

	fls.redistance(dx);
	sls.redistance(dx);

	//out << sls;

	Grid u(std::string("..\\res\\u_test.txt"));
	Grid v(std::string("..\\res\\v_test.txt"));

	Grid g(w, h + 1, 0.5, 0.5, 9.8);

	Velocity V(u, v);

	Solver solver(w, h, dx, dt);
	Particles P(w, h, r, fls);

	/* Here we initialize the particles */
	P.init(V.u, KEY_VX);
	P.init(V.v, KEY_VY);

	for (int t = 0; t < 5; t++) {
		std::cerr << "\rComputing frame: " << t;

		/* Here we construct the fluid levelset from the particles */
		fls.fit(P);

		/* Here we transfer velocity from particles to grid */
		V.u.fit(P, KEY_VX);
		V.v.fit(P, KEY_VY);

		/* Add body force to the velocity field */
		V.v += g * dt;

		/* Construct solid levelset and solid velocity fields */
		fls = fls & sls;

		/* Solve for pressure to get a divergence free velocity fields */
		Velocity Vnew = solver.project(fls, sls, V);

		/* Update particles velocities */
		P.fit(Vnew.u, V.u, KEY_VX);
		P.fit(Vnew.v, V.v, KEY_VY);

		/* Advect the particles through the divergence free velocity field */
		P.advect(dt, Vnew, fls, sls);

		/* Here we save current frame */
		std::ostringstream spath;
		spath << "output\\" << std::setw(3) << std::setfill('0') << t << ".png" << std::setfill(' ');

		const std::string& path = spath.str();

		img << fls;
		lodepng::encode(path.c_str(), img, w, h);
	}
	
	return 0;
}