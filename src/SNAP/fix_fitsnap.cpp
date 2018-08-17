/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include <cstdlib>
#include <cstring>
#include "fix_fitsnap.h"
#include "force.h"
#include "update.h"
#include "modify.h"
#include "compute.h"
#include "input.h"
#include "variable.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

// This fix provides output to fitsnap.py
// 
// Usage: 
//   fix ID group-ID fitsnap Nevery file outfile
//   fix myfs all fitsnap 1 file myfs.out.$i
// 
// use in combination with pair_style hybrid/overlay with dummy snap pair_style
//
// pair_style   hybrid/overlay lj/cut ${rcutfac} &
//              zbl ${zblcutinner} ${zblcutouter} &
//              snap
//
// pair_coeff 1 1 zbl ${zblz} ${zblz}
// pair_coeff * * snap NULL Ta Ta.snapparam Ta
//
// fix fitsnap will write 1+3N+6 rows to the file outfile
// (it may may sense to an append keyword)
// the rows correspond to potential energy, forces, and stress tensor
// the energy row:
//     Nelem1/atoms,Sum(B_1^i,i=1,Nelem1)/atoms,..NelemM/atoms,...,Sum(B_K^i,i=1,NelemM)/atoms, pe
// the force rows:
//     0,Sum(dB_1^i/dxj,i=1,Nelem1),...,Sum(dB_K^i/dxj,i=1,NelemM), fx[j]
// the stress rows:
//     0,Sum(xj*Sum(dB_1^i/dyj,i=1,Nelem1),j=1,atoms,...,0,...,Sum(xj*Sum(dB_K^i/dyj,i=1,NelemM), pxy

/* ---------------------------------------------------------------------- */

FixFitSNAP::FixFitSNAP(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 11) error->all(FLERR,"Illegal fix fitsnap command");

  vector_flag = 1;
  size_vector = 3;
  global_freq = 1;
  extvector = 0;

  nevery = force->inumeric(FLERR,arg[3]);
  if (nevery <= 0) error->all(FLERR,"Illegal fix fitsnap command");

  // process variable arg

  int iarg = 8;
  if (strncmp(arg[iarg],"c_",2) == 0 ||
      strncmp(arg[iarg],"f_",2) == 0 ||
      strncmp(arg[iarg],"v_",2) == 0) {

    int n = strlen(arg[iarg]);

    iarg++;

  } else error->all(FLERR,"Illegal fix fitsnap command");

  // setpoint arg

  double setpoint = force->numeric(FLERR,arg[iarg]);
  iarg++;

  // error check

  if (1) {
    if (1)
      error->all(FLERR,"Compute ID for fix fitsnap does not exist");
  }
  firsttime = 1;
}

/* ---------------------------------------------------------------------- */

FixFitSNAP::~FixFitSNAP()
{
}

/* ---------------------------------------------------------------------- */

int FixFitSNAP::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixFitSNAP::init()
{
}

/* ---------------------------------------------------------------------- */

void FixFitSNAP::post_force()
{
}

