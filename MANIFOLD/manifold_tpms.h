/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef LMP_MANIFOLD_TPMS_H
#define LMP_MANIFOLD_TPMS_H

#include "manifold.h"

namespace LAMMPS_NS
{

namespace user_manifold
{

class manifold_Primitive : public manifold
{
  public:
    enum { NPARAMS = 3 }; // Number of parameters.
    manifold_Primitive(LAMMPS *lmp, int, char **);
    double g(const double *x) override;
    void n(const double *x, double *n) override;

    static const char *type() { return "Primitive"; }
    const char *id() override { return type(); }
    static int expected_argc() { return NPARAMS; }
    int nparams() override { return NPARAMS; }
};

class manifold_Gyroid : public manifold
{
  public:
    enum { NPARAMS = 3 }; // Number of parameters.
    manifold_Gyroid(LAMMPS *lmp, int, char **);
    double g(const double *x) override;
    void n(const double *x, double *n) override;

    static const char *type() { return "Gyroid"; }
    const char *id() override { return type(); }
    static int expected_argc() { return NPARAMS; }
    int nparams() override { return NPARAMS; }
};

class manifold_D : public manifold
{
  public:
    enum { NPARAMS = 3 }; // Number of parameters.
    manifold_D(LAMMPS *lmp, int, char **);
    double g(const double *x) override;
    void n(const double *x, double *n) override;

    static const char *type() { return "D"; }
    const char *id() override { return type(); }
    static int expected_argc() { return NPARAMS; }
    int nparams() override { return NPARAMS; }
};

class manifold_Diamond : public manifold
{
  public:
    enum { NPARAMS = 3 }; // Number of parameters.
    manifold_Diamond(LAMMPS *lmp, int, char **);
    double g(const double *x) override;
    void n(const double *x, double *n) override;

    static const char *type() { return "Diamond"; }
    const char *id() override { return type(); }
    static int expected_argc() { return NPARAMS; }
    int nparams() override { return NPARAMS; }
};

class manifold_IWP : public manifold
{
  public:
    enum { NPARAMS = 3 }; // Number of parameters.
    manifold_IWP(LAMMPS *lmp, int, char **);
    double g(const double *x) override;
    void n(const double *x, double *n) override;

    static const char *type() { return "IWP"; }
    const char *id() override { return type(); }
    static int expected_argc() { return NPARAMS; }
    int nparams() override { return NPARAMS; }
};

class manifold_Neovius : public manifold
{
  public:
    enum { NPARAMS = 3 }; // Number of parameters.
    manifold_Neovius(LAMMPS *lmp, int, char **);
    double g(const double *x) override;
    void n(const double *x, double *n) override;

    static const char *type() { return "Neovius"; }
    const char *id() override { return type(); }
    static int expected_argc() { return NPARAMS; }
    int nparams() override { return NPARAMS; }
};

class manifold_S : public manifold
{
  public:
    enum { NPARAMS = 3 }; // Number of parameters.
    manifold_S(LAMMPS *lmp, int, char **);
    double g(const double *x) override;
    void n(const double *x, double *n) override;

    static const char *type() { return "S"; }
    const char *id() override { return type(); }
    static int expected_argc() { return NPARAMS; }
    int nparams() override { return NPARAMS; }
};

class manifold_FRD : public manifold
{
  public:
    enum { NPARAMS = 3 }; // Number of parameters.
    manifold_FRD(LAMMPS *lmp, int, char **);
    double g(const double *x) override;
    void n(const double *x, double *n) override;

    static const char *type() { return "FRD"; }
    const char *id() override { return type(); }
    static int expected_argc() { return NPARAMS; }
    int nparams() override { return NPARAMS; }
};

class manifold_PMY : public manifold
{
  public:
    enum { NPARAMS = 3 }; // Number of parameters.
    manifold_PMY(LAMMPS *lmp, int, char **);
    double g(const double *x) override;
    void n(const double *x, double *n) override;

    static const char *type() { return "PMY"; }
    const char *id() override { return type(); }
    static int expected_argc() { return NPARAMS; }
    int nparams() override { return NPARAMS; }
};

} // namespace user_manifold

} // namespace LAMMPS_NS

#endif
