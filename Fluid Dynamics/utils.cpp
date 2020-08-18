#include "utils.h"
#include "velocity.h"

Vec RK3(const double& dt, const Vec& p, const Velocity& V)
{
    Vec v1(V.u[p], V.v[p]);
    Vec p1 = p + v1 * 0.5 * dt;

    Vec v2(V.u[p1], V.v[p1]);
    Vec p2 = p + v2 * 0.75 * dt;

    Vec v3(V.u[p2], V.v[p2]);

    return p + (v1 * 2 + v2 * 3 + v3 * 4) * (dt / 9.);
}
