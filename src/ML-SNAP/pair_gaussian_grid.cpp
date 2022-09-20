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

#include "pair_grid.h"
#include "pair_gaussian_grid.h"
#include "sna.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "force.h"
#include "pair.h"
#include "domain.h"
#include "comm.h"
#include "math_const.h"
#include "math_special.h"
#include "memory.h"
#include "error.h"
#include "tokenizer.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using MathConst::MY_2PI;
using MathSpecial::powint;

PairGaussianGrid::PairGaussianGrid(LAMMPS *lmp) :
  PairGrid(lmp),
  radelem(nullptr), sigmaelem(nullptr), prefacelem(nullptr), argfacelem(nullptr)
{}

/* ---------------------------------------------------------------------- */

PairGaussianGrid::~PairGaussianGrid()
{
  memory->destroy(radelem);
  memory->destroy(sigmaelem);
  memory->destroy(prefacelem);
  memory->destroy(argfacelem);
  memory->destroy(cutsq);
  memory->destroy(map);
}

/* ---------------------------------------------------------------------- */

void PairGaussianGrid::init_style()
{
  if (force->newton_pair == 0)
    error->all(FLERR,"Pair style gaussian/grid requires newton pair on");

  ndesc = ndesc_base + atom->ntypes;
}

/* ---------------------------------------------------------------------- */

void PairGaussianGrid::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void PairGaussianGrid::compute(int eflag, int vflag)
{
  double fij[3];
  
  ev_init(eflag,vflag);

  // compute gaussian for each gridpoint

  double** const x = atom->x;
  double **f = atom->f;
  const int* const mask = atom->mask;
  int * const type = atom->type;
  const int ntotal = atom->nlocal + atom->nghost;

  // first generate fingerprint,
  // which allows calculation of beta

  for (int iz = nzlo; iz <= nzhi; iz++)
    for (int iy = nylo; iy <= nyhi; iy++)
      for (int ix = nxlo; ix <= nxhi; ix++)
	for (int icol = ndesc_base; icol < ndesc; icol++)
	  gridlocal[icol][ix][iy][iz] = 0.0;
	
  int igrid = 0;
  for (int iz = nzlo; iz <= nzhi; iz++)
    for (int iy = nylo; iy <= nyhi; iy++)
      for (int ix = nxlo; ix <= nxhi; ix++) {
	double xgrid[3];
	grid2x(ix, iy, iz, xgrid);
	const double xtmp = xgrid[0];
	const double ytmp = xgrid[1];
	const double ztmp = xgrid[2];

	for (int j = 0; j < ntotal; j++) {

	  const double delx = xtmp - x[j][0];
	  const double dely = ytmp - x[j][1];
	  const double delz = ztmp - x[j][2];
	  const double rsq = delx*delx + dely*dely + delz*delz;
	  int jtype = type[j];
          if (rsq < cutsq[jtype][jtype]) {
	    int icol = ndesc_base + jtype - 1; 
	    gridlocal[icol][ix][iy][iz] += prefacelem[jtype] * exp(-rsq * argfacelem[jtype]);
	  }

	}
	igrid++;
      }

  // this is a proxy for a call to the energy model
  // beta is dE/dB^i, the derivative of the total
  // energy w.r.t. to descriptors of grid point i
  
  compute_beta();
  
  // second compute forces using beta
  
  igrid = 0;
  for (int iz = nzlo; iz <= nzhi; iz++)
    for (int iy = nylo; iy <= nyhi; iy++)
      for (int ix = nxlo; ix <= nxhi; ix++) {
	double xgrid[3];
	grid2x(ix, iy, iz, xgrid);
	const double xtmp = xgrid[0];
	const double ytmp = xgrid[1];
	const double ztmp = xgrid[2];

	for (int j = 0; j < ntotal; j++) {

	  const double delx = xtmp - x[j][0];
	  const double dely = ytmp - x[j][1];
	  const double delz = ztmp - x[j][2];
	  const double rsq = delx*delx + dely*dely + delz*delz;
	  int jtype = type[j];
          if (rsq < cutsq[jtype][jtype]) {

	    // f(Rj) = b * exp(-a * rsq) * beta;
	    // df / dRjx = -2 * a * b * exp(-a * rsq) * Rjx * beta;

	    double factmp = -2.0 * argfacelem[jtype] * prefacelem[jtype] * exp(-rsq * argfacelem[jtype]) * beta[igrid][jtype-1];
	  
	    fij[0] = factmp * delx;
	    fij[1] = factmp * dely;
	    fij[2] = factmp * delz;
	    f[j][0] += fij[0];
	    f[j][1] += fij[1];
	    f[j][2] += fij[2];

	    // tally per-atom virial contribution

	    if (vflag)
	      ev_tally_xyz(-1,j,atom->nlocal,force->newton_pair,0.0,0.0,
			   fij[0],fij[1],fij[2],
			   -delx,-dely,-delz);
	  }

	}
	
	// tally energy contribution

	if (eflag) {

	  // evdwl = energy of grid point I, sum over types

	  double evdwl = 0.0;
	  for (int jtype = 1; jtype <= atom->ntypes; jtype++) {
	    int icol = ndesc_base + jtype - 1; 	    
	    evdwl += beta[igrid][jtype-1] * gridlocal[icol][ix][iy][iz];
	  }
	  ev_tally_full(-1,2.0*evdwl,0.0,0.0,0.0,0.0,0.0);

	}
	igrid++;
      }

  if (vflag_fdotr) virial_fdotr_compute();
  
}


/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairGaussianGrid::settings(int narg, char ** arg)
{

  // call base class first

  PairGrid::settings(narg, arg);
    
  // skip over arguments used by base class
  // so that argument positions are identical to
  // regular per-atom compute
  
  arg += nargbase;
  narg -= nargbase;

  int ntypes = atom->ntypes;
  int nargmin = 1+2*ntypes;

  if (narg < nargmin) error->all(FLERR,"Illegal pair gaussian/grid command");

  // process required arguments

  memory->create(radelem,ntypes+1,"pair:gaussian/grid:radelem"); // offset by 1 to match up with types
  memory->create(sigmaelem, ntypes + 1, "gaussian/atom:sigmaelem");
  memory->create(prefacelem, ntypes + 1, "gaussian/atom:prefacelem");
  memory->create(argfacelem, ntypes + 1, "gaussian/atom:argfacelem");

  rcutfac = utils::numeric(FLERR, arg[0], false, lmp);
  for (int i = 0; i < ntypes; i++) radelem[i + 1] = utils::numeric(FLERR, arg[1 + i], false, lmp);
  for (int i = 0; i < ntypes; i++)
    sigmaelem[i + 1] = utils::numeric(FLERR, arg[ntypes + 1 + i], false, lmp);

  // construct cutsq

  double cut;
  cutmax = 0.0;
  memory->create(cutsq, ntypes + 1, ntypes + 1, "gaussian/grid:cutsq");
  for (int i = 1; i <= ntypes; i++) {
    cut = 2.0 * radelem[i] * rcutfac;
    if (cut > cutmax) cutmax = cut;
    cutsq[i][i] = cut * cut;
    for (int j = i + 1; j <= ntypes; j++) {
      cut = (radelem[i] + radelem[j]) * rcutfac;
      cutsq[i][j] = cutsq[j][i] = cut * cut;
    }
  }

  // pre-compute coefficients
  
  for (int i = 0; i < ntypes; i++) {
    prefacelem[i + 1] = 1.0/powint(sigmaelem[i + 1] * sqrt(MY_2PI), 3);
    argfacelem[i + 1] = 1.0/(2.0 * sigmaelem[i + 1] * sigmaelem[i + 1]);
  }

}
  
/* ----------------------------------------------------------------------
   memory usage
------------------------------------------------------------------------- */

double PairGaussianGrid::memory_usage()
{
  int n = atom->ntypes+1;
  int nbytes = (double)n*sizeof(int);            // map

  return nbytes;
}

