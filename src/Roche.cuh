#pragma once

#include "Header.h"


static inline double Roche_R(double phi)
{
    double c = fabs(cos(phi));
    if (c < 1e-12) return R_POLE;                 // точный полюс
    double u = OMEGA_FRAC * c;
    if (u >= 1.0) u = 1.0 - 1e-12;                // защита при ??1
    if (u < 1e-3) return R_POLE * (1.0 + (4.0 / 27.0) * u * u);  // “ейлор у полюса
    return 3.0 * R_POLE * cos((M_PI + acos(u)) / 3.0) / u;
}

//  оординаты вершины (?_i, ?_k)
static inline void vertex_bf(int i, int k, double* x, double* y)
{
    double phi = PHI_EDGE(k);
    double r = R_EDGE_BF(i, phi);
    *x = r * cos(phi);
    *y = r * sin(phi);
}

// ƒлина ?-грани €чейки (i,j): хорда между (?_i,?_j) и (?_i,?_{j+1})
static inline double L_xi_bf(int i, int j)
{
    double x1, y1, x2, y2;
    vertex_bf(i, j, &x1, &y1);
    vertex_bf(i, j + 1, &x2, &y2);
    return hypot(x2 - x1, y2 - y1);
}

#define L_INNER_BF(i,j)  L_xi_bf(i,   j)   // внутренн€€ ?-грань
#define L_OUTER_BF(i,j)  L_xi_bf(i+1, j)   // внешн€€ ?-грань