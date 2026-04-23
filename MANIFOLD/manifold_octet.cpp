// clang-format off
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

#include "manifold_octet.h"

#include <array>
#include <cmath>

using namespace LAMMPS_NS;

using namespace user_manifold;

//// OCTET
// dot product of 2 vector
double ddot(double x[], double y[])
{
    return x[0] * y[0] + x[1] * y[1] + x[2] * y[2];
}

// Length of vector: |vec|
double norm(double x[])
{
    return sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
}

// Calculate distance from point: x0 to straight line defined by: x1, x2
// rad: vector from (projection of x0 on line) to (point x0)
// Let x3 be the projection of x0 on line x1->x2
// vector: x10 = x1->x0;
// rad = x31 + x10
//
// radius: distance from point x0 to line => radius = |r|
// grad: gradient of scalar field: radius => grad = r / |r|
//
double proj2line(const std::array<double, 3> &x1,
                 const std::array<double, 3> &x2, double x0[], double grad[])
{
    // vector x1->x2, x1->x0
    double x12[3] = {x2[0] - x1[0], x2[1] - x1[1], x2[2] - x1[2]};
    double x10[3] = {x0[0] - x1[0], x0[1] - x1[1], x0[2] - x1[2]};

    // length(x12) and length(x13)
    double x12_norm = norm(x12);
    double x10_proj = ddot(x10, x12) / x12_norm;

    double rad[3];
    for (int i = 0; i < 3; i++) {
        rad[i] = -x12[i] / x12_norm * x10_proj + x10[i];
    }
    double radius = norm(rad);
    for (int i = 0; i < 3; i++) {
        if (radius / x12_norm < 1e-6)
            grad[i] = 0;
        else
            grad[i] = rad[i] / radius;
    }

    return radius;
}

// Index of minimum element in array
int argmin(const double arr[], int size)
{
    int minIndex = 0;
    double minValue = arr[0];

    for (int i = 1; i < size; ++i) {
        if (arr[i] < minValue) {
            minValue = arr[i];
            minIndex = i;
        }
    }

    return minIndex;
}

// Define octet value field and gradient
double octet(double *x, double *n, double L)
{
    double r11, r12, r13, r14, r15, r16, r17, r18;
    double g11[3], g12[3], g13[3], g14[3], g15[3], g16[3], g17[3], g18[3];
    double r21, r22, r23, r24, r25, r26, r27, r28;
    double g21[3], g22[3], g23[3], g24[3], g25[3], g26[3], g27[3], g28[3];
    double r31, r32, r33, r34, r35, r36, r37, r38;
    double g31[3], g32[3], g33[3], g34[3], g35[3], g36[3], g37[3], g38[3];

    // xy plane
    r11 = proj2line({0e0, 0, 0}, {L, L, 0}, x, g11);
    r12 = proj2line({0, L, 0}, {L, 0, 0}, x, g12);
    r13 = proj2line({0, 0, L}, {L, L, L}, x, g13);
    r14 = proj2line({0, L, L}, {L, 0, L}, x, g14);

    r15 = proj2line({L / 2, 0, L / 2}, {0, L / 2, L / 2}, x, g15);
    r16 = proj2line({L / 2, 0, L / 2}, {L, L / 2, L / 2}, x, g16);
    r17 = proj2line({L / 2, L, L / 2}, {0, L / 2, L / 2}, x, g17);
    r18 = proj2line({L, L / 2, L / 2}, {L / 2, L, L / 2}, x, g18);

    // yz plane
    r21 = proj2line({0e0, 0, 0}, {0, L, L}, x, g21);
    r22 = proj2line({0, 0, L}, {0, L, 0}, x, g22);
    r23 = proj2line({L, 0, 0}, {L, L, L}, x, g23);
    r24 = proj2line({L, 0, L}, {L, L, 0}, x, g24);

    r25 = proj2line({L / 2, L / 2, 0}, {L / 2, 0, L / 2}, x, g25);
    r26 = proj2line({L / 2, L / 2, 0}, {L / 2, L, L / 2}, x, g26);
    r27 = proj2line({L / 2, L / 2, L}, {L / 2, 0, L / 2}, x, g27);
    r28 = proj2line({L / 2, L, L / 2}, {L / 2, L / 2, L}, x, g28);

    // zx plane
    r31 = proj2line({0e0, 0, 0}, {L, 0, L}, x, g31);
    r32 = proj2line({L, 0, 0}, {0, 0, L}, x, g32);
    r33 = proj2line({0, L, 0}, {L, L, L}, x, g33);
    r34 = proj2line({L, L, 0}, {0, L, L}, x, g34);

    r35 = proj2line({0, L / 2, L / 2}, {L / 2, L / 2, 0}, x, g35);
    r36 = proj2line({0, L / 2, L / 2}, {L / 2, L / 2, L}, x, g36);
    r37 = proj2line({L, L / 2, L / 2}, {L / 2, L / 2, 0}, x, g37);
    r38 = proj2line({L / 2, L / 2, L}, {L, L / 2, L / 2}, x, g38);

    //
    double r_s[] = {
        r11, r12, r13, r14, r15, r16, r17, r18, r21, r22, r23, r24,
        r25, r26, r27, r28, r31, r32, r33, r34, r35, r36, r37, r38,
    };

    double *g_s[] = {
        g11, g12, g13, g14, g15, g16, g17, g18, g21, g22, g23, g24,
        g25, g26, g27, g28, g31, g32, g33, g34, g35, g36, g37, g38,
    };

    int minIdx = argmin(r_s, 24);

    double r = r_s[minIdx];

    n[0] = g_s[minIdx][0];
    n[1] = g_s[minIdx][1];
    n[2] = g_s[minIdx][2];

    return r;
}

/// LMP

manifold_octet::manifold_octet( LAMMPS *lmp, int /*argc*/,
                                      char **/*argv*/ ) : manifold(lmp)
{}


double manifold_octet::g( const double *x )
{
  double L = params[0];
  double t = params[1];

  double x0[] = {x[0], x[1], x[2]};
  double n0[3];
  return octet(x0, n0, L) - t;
  
  // return t -sin(2*M_PI*x[0]/L)*cos(2*M_PI*x[1]/L) - sin(2*M_PI*x[1]/L)*cos(2*M_PI*x[2]/L) - sin(2*M_PI*x[2]/L)*cos(2*M_PI*x[0]/L);
}

void manifold_octet::n( const double *x, double *n )
{
  double L = params[0];
  double t = params[1];

  double x0[] = {x[0], x[1], x[2]};
  double n0[3];

  octet(x0, n0, L);

  n[0] = n0[0];
  n[1] = n0[1];
  n[2] = n0[2];


  // n[0] = -2*M_PI/L*cos(2*M_PI*x[0]/L)*cos(2*M_PI*x[1]/L) + 2*M_PI/L*sin(2*M_PI*x[2]/L)*sin(2*M_PI*x[0]/L);
  // n[1] = -2*M_PI/L*cos(2*M_PI*x[1]/L)*cos(2*M_PI*x[2]/L) + 2*M_PI/L*sin(2*M_PI*x[0]/L)*sin(2*M_PI*x[1]/L);
  // n[2] = -2*M_PI/L*cos(2*M_PI*x[2]/L)*cos(2*M_PI*x[0]/L) + 2*M_PI/L*sin(2*M_PI*x[1]/L)*sin(2*M_PI*x[2]/L);
}


