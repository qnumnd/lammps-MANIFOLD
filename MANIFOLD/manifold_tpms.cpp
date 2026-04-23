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

#include "manifold_tpms.h"
#include "autodiff.h"


#include <cmath>

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

using namespace LAMMPS_NS;

using namespace user_manifold;

// clang-format off
real Primitive (real1 x, real2 y, real3 z, double L, double t) {return cos(2*M_PI/L*x) + cos(2*M_PI/L*y) + cos(2*M_PI/L*z) - t; }; // Schwarz Primitive
real Gyroid    (real1 x, real2 y, real3 z, double L, double t) {return sin(2*M_PI/L*x)*cos(2*M_PI/L*y) + sin(2*M_PI/L*y)*cos(2*M_PI/L*z) + sin(2*M_PI/L*z)*cos(2*M_PI/L*x) - t; }; // Schoen Gyroid
real D         (real1 x, real2 y, real3 z, double L, double t) {return cos(2*M_PI/L*x)*cos(2*M_PI/L*y)*cos(2*M_PI/L*z) - sin(2*M_PI/L*x)*sin(2*M_PI/L*y)*sin(2*M_PI/L*z) - t; }; // Schwarz-D
real Diamond   (real1 x, real2 y, real3 z, double L, double t) {return sin(2*M_PI/L*x)*sin(2*M_PI/L*y)*sin(2*M_PI/L*z) + sin(2*M_PI/L*x)*cos(2*M_PI/L*y)*cos(2*M_PI/L*z) + cos(2*M_PI/L*x)*sin(2*M_PI/L*y)*cos(2*M_PI/L*z) + cos(2*M_PI/L*x)*cos(2*M_PI/L*y)*sin(2*M_PI/L*z) - t; }; // Diamond
real IWP       (real1 x, real2 y, real3 z, double L, double t) {return 2*(cos(2*M_PI/L*x)*cos(2*M_PI/L*y) + cos(2*M_PI/L*y)*cos(2*M_PI/L*z) + cos(2*M_PI/L*z)*cos(2*M_PI/L*x)) - (cos(2*M_PI/L*2*x) + cos(2*M_PI/L*2*y) + cos(2*M_PI/L*2*z)) - t; }; // Schoen IWP
real Neovius   (real1 x, real2 y, real3 z, double L, double t) {return 3*(cos(2*M_PI/L*x) + cos(2*M_PI/L*y) + cos(2*M_PI/L*z)) + 4*cos(2*M_PI/L*x)*cos(2*M_PI/L*y)*cos(2*M_PI/L*z) - t; }; // Neovius
real S         (real1 x, real2 y, real3 z, double L, double t) {return cos(2*M_PI/L*2*x)*sin(2*M_PI/L*y)*cos(2*M_PI/L*z) + cos(2*M_PI/L*x)*cos(2*M_PI/L*2*y)*sin(2*M_PI/L*z) + sin(2*M_PI/L*x)*cos(2*M_PI/L*y)*cos(2*M_PI/L*2*z) - t; }; // Fischer-Koch S
real FRD       (real1 x, real2 y, real3 z, double L, double t) {return 4*cos(2*M_PI/L*x)*cos(2*M_PI/L*y)*cos(2*M_PI/L*z) - (cos(2*M_PI/L*2*x)*cos(2*M_PI/L*2*y) + cos(2*M_PI/L*2*y)*cos(2*M_PI/L*2*z) + cos(2*M_PI/L*2*z)*cos(2*M_PI/L*2*x)) - t; }; // Schoen FRD
real PMY       (real1 x, real2 y, real3 z, double L, double t) {return 2*cos(2*M_PI/L*x)*cos(2*M_PI/L*y)*cos(2*M_PI/L*z) + sin(2*M_PI/L*2*x)*sin(2*M_PI/L*y) + sin(2*M_PI/L*x)*sin(2*M_PI/L*2*z) + sin(2*M_PI/L*2*y)*sin(2*M_PI/L*z) - t; }; // PMY
// clang-format on

// params: = L, n_cell, t
#define P_L params[0]
#define P_n_cell params[1]
#define P_t params[2]

/* ---------------------------------------------------------------------- */
//
manifold_Primitive::manifold_Primitive(LAMMPS *lmp, int /*argc*/,
                                       char ** /*argv*/)
    : manifold(lmp)
{
}

double manifold_Primitive::g(const double *x)
{
    return evaluate(Primitive, x, P_L / P_n_cell, P_t);
}

void manifold_Primitive::n(const double *x, double *n)
{
    gradient(Primitive, n, x, P_L / P_n_cell, P_t);
}

/* ---------------------------------------------------------------------- */
//
manifold_Gyroid::manifold_Gyroid(LAMMPS *lmp, int /*argc*/, char ** /*argv*/)
    : manifold(lmp)
{
}

double manifold_Gyroid::g(const double *x)
{
    return evaluate(Gyroid, x, P_L / P_n_cell, P_t);
}

void manifold_Gyroid::n(const double *x, double *n)
{
    gradient(Gyroid, n, x, P_L / P_n_cell, P_t);
}

/* ---------------------------------------------------------------------- */
//
manifold_D::manifold_D(LAMMPS *lmp, int /*argc*/, char ** /*argv*/)
    : manifold(lmp)
{
}

double manifold_D::g(const double *x)
{
    return evaluate(D, x, P_L / P_n_cell, P_t);
}

void manifold_D::n(const double *x, double *n)
{
    gradient(D, n, x, P_L / P_n_cell, P_t);
}

/* ---------------------------------------------------------------------- */
//
manifold_Diamond::manifold_Diamond(LAMMPS *lmp, int /*argc*/, char ** /*argv*/)
    : manifold(lmp)
{
}

double manifold_Diamond::g(const double *x)
{
    return evaluate(Diamond, x, P_L / P_n_cell, P_t);
}

void manifold_Diamond::n(const double *x, double *n)
{
    gradient(Diamond, n, x, P_L / P_n_cell, P_t);
}

/* ---------------------------------------------------------------------- */
//
manifold_IWP::manifold_IWP(LAMMPS *lmp, int /*argc*/, char ** /*argv*/)
    : manifold(lmp)
{
}

double manifold_IWP::g(const double *x)
{
    return evaluate(IWP, x, P_L / P_n_cell, P_t);
}

void manifold_IWP::n(const double *x, double *n)
{
    gradient(IWP, n, x, P_L / P_n_cell, P_t);
}

/* ---------------------------------------------------------------------- */
//
manifold_Neovius::manifold_Neovius(LAMMPS *lmp, int /*argc*/, char ** /*argv*/)
    : manifold(lmp)
{
}

double manifold_Neovius::g(const double *x)
{
    return evaluate(Neovius, x, P_L / P_n_cell, P_t);
}

void manifold_Neovius::n(const double *x, double *n)
{
    gradient(Neovius, n, x, P_L / P_n_cell, P_t);
}

/* ---------------------------------------------------------------------- */
//
manifold_S::manifold_S(LAMMPS *lmp, int /*argc*/, char ** /*argv*/)
    : manifold(lmp)
{
}

double manifold_S::g(const double *x)
{
    return evaluate(S, x, P_L / P_n_cell, P_t);
}

void manifold_S::n(const double *x, double *n)
{
    gradient(S, n, x, P_L / P_n_cell, P_t);
}

/* ---------------------------------------------------------------------- */
//
manifold_FRD::manifold_FRD(LAMMPS *lmp, int /*argc*/, char ** /*argv*/)
    : manifold(lmp)
{
}

double manifold_FRD::g(const double *x)
{
    return evaluate(FRD, x, P_L / P_n_cell, P_t);
}

void manifold_FRD::n(const double *x, double *n)
{
    gradient(FRD, n, x, P_L / P_n_cell, P_t);
}

/* ---------------------------------------------------------------------- */
//
manifold_PMY::manifold_PMY(LAMMPS *lmp, int /*argc*/, char ** /*argv*/)
    : manifold(lmp)
{
}

double manifold_PMY::g(const double *x)
{
    return evaluate(PMY, x, P_L / P_n_cell, P_t);
}

void manifold_PMY::n(const double *x, double *n)
{
    gradient(PMY, n, x, P_L / P_n_cell, P_t);
}

#undef P_t
#undef P_L
#undef P_n_cell
