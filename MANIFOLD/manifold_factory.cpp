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
   -----------------------------------------------------------------------

   This file is a part of the MANIFOLD package.

   Copyright (2013-2014) Stefan Paquay, Eindhoven University of Technology.
   License: GNU General Public License.

   See the README file in the top-level LAMMPS directory.

   This file is part of the user-manifold package written by
   Stefan Paquay at the Eindhoven University of Technology.
   This module makes it possible to do MD with particles constrained
   to pretty arbitrary manifolds characterized by some constraint function
   g(x,y,z) = 0 and its normal grad(g). The number of manifolds available
   right now is limited but can be extended straightforwardly by making
   a new class that inherits from manifold and implements all pure virtual
   methods.

   Thanks to Remy Kusters for beta-testing!

------------------------------------------------------------------------- */

#include "manifold_factory.h"

// #include "manifold_octet.h"
#include "manifold_tpms.h"
#include "manifold_wrinkle.h"

#include <cstring>

namespace LAMMPS_NS {
namespace user_manifold {

  template <typename m_type>
  void make_manifold_if(manifold **man_ptr, const char *name, LAMMPS *lmp, int narg, char **arg)
  {
    if (strcmp(m_type::type(), name) == 0) {
      if (*man_ptr == nullptr) { *man_ptr = new m_type(lmp, narg, arg); }
    }
  }

  manifold *create_manifold(const char *mname, LAMMPS *lmp, int narg, char **arg)
  {
    manifold *man = nullptr;

    // make_manifold_if<manifold_octet>(&man, mname, lmp, narg, arg);

    make_manifold_if<manifold_Primitive>(&man, mname, lmp, narg, arg);
    make_manifold_if<manifold_Gyroid>(&man, mname, lmp, narg, arg);
    make_manifold_if<manifold_D>(&man, mname, lmp, narg, arg);
    make_manifold_if<manifold_Diamond>(&man, mname, lmp, narg, arg);
    make_manifold_if<manifold_IWP>(&man, mname, lmp, narg, arg);
    make_manifold_if<manifold_Neovius>(&man, mname, lmp, narg, arg);
    make_manifold_if<manifold_S>(&man, mname, lmp, narg, arg);
    make_manifold_if<manifold_FRD>(&man, mname, lmp, narg, arg);
    make_manifold_if<manifold_PMY>(&man, mname, lmp, narg, arg);

    make_manifold_if<manifold_SinXY>(&man, mname, lmp, narg, arg);


    return man;
  }
}    // namespace user_manifold

}    // namespace LAMMPS_NS
