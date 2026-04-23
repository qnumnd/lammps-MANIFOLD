/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
   ----------------------------------------------------------------------- */

#include "manifold_wrinkle.h"
#include "autodiff.h"

#include <cmath>

using namespace LAMMPS_NS;

using namespace user_manifold;

// params: = L, n_cell, A, t
#define P_L params[0]
#define P_n_cell params[1]
#define P_A params[2]
#define P_t params[3]

/* ---------------------------------------------------------------------- */
real SinXY(real1 x, real2 y, real3 z, double L, double t, double A)
{
    return (2 * M_PI / L * z) -
           (A * M_PI) * sin(2 * M_PI / L * x) * sin(2 * M_PI / L * y) - M_PI;
}

/* ---------------------------------------------------------------------- */
//
manifold_SinXY::manifold_SinXY(LAMMPS *lmp, int /*argc*/, char ** /*argv*/)
    : manifold(lmp)
{
}

double manifold_SinXY::g(const double *x)
{

    return evaluate(SinXY, x, P_L / P_n_cell, P_t, P_A);
}

void manifold_SinXY::n(const double *x, double *n)
{
    gradient(SinXY, n, x, P_L / P_n_cell, P_t, P_A);
}

#undef P_A
#undef P_t
#undef P_L
#undef P_n_cell
