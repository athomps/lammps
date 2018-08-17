/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(fitsnap,FixFitSNAP)

#else

#ifndef LMP_FIX_FITSNAP_H
#define LMP_FIX_FITSNAP_H

#include "fix.h"

namespace LAMMPS_NS {

class FixFitSNAP : public Fix {
 public:
  FixFitSNAP(class LAMMPS *, int, char **);
  ~FixFitSNAP();
  int setmask();
  void init();
  void post_force();

 private:
  int firsttime;

};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

UNDOCUMENTED

*/
