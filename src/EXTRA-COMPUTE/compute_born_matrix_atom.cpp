/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
  Contributing Authors : Germain Clavier (TUe), Aidan Thompson (Sandia)
------------------------------------------------------------------------- */

#include "compute_born_matrix_atom.h"

#include "angle.h"
#include "atom.h"
#include "atom_vec.h"
#include "bond.h"
#include "comm.h"
#include "compute.h"
#include "dihedral.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "improper.h"
#include "kspace.h"
#include "memory.h"
#include "modify.h"
#include "molecule.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "pair.h"
#include "universe.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

#define BIG 1000000000
#define SMALL 1e-16

// this table is used to pick the 3d rij vector indices used to
// compute the 6 indices long Voigt stress vector

static int constexpr sigma_albe[6][2] = {
    {0, 0},    // s11
    {1, 1},    // s22
    {2, 2},    // s33
    {1, 2},    // s44
    {0, 2},    // s55
    {0, 1},    // s66
};

// this table is used to pick the correct indices from the Voigt
// stress vector to compute the Cij matrix (21 terms, see doc) contribution

static int constexpr C_albe[21][2] = {
    {0, 0},    // C11
    {1, 1},    // C22
    {2, 2},    // C33
    {3, 3},    // C44
    {4, 4},    // C55
    {5, 5},    // C66
    {0, 1},    // C12
    {0, 2},    // C13
    {0, 3},    // C14
    {0, 4},    // C15
    {0, 5},    // C16
    {1, 2},    // C23
    {1, 3},    // C24
    {1, 4},    // C25
    {1, 5},    // C26
    {2, 3},    // C34
    {2, 4},    // C35
    {2, 5},    // C36
    {3, 4},    // C45
    {3, 5},    // C46
    {4, 5}     // C56
};

// this table is used to pick the 3d rij vector indices used to
// compute the 21 indices long Cij matrix

static int constexpr albemunu[21][4] = {
    {0, 0, 0, 0},    // C11
    {1, 1, 1, 1},    // C22
    {2, 2, 2, 2},    // C33
    {1, 2, 1, 2},    // C44
    {0, 2, 0, 2},    // C55
    {0, 1, 0, 1},    // C66
    {0, 0, 1, 1},    // C12
    {0, 0, 2, 2},    // C13
    {0, 0, 1, 2},    // C14
    {0, 0, 0, 2},    // C15
    {0, 0, 0, 1},    // C16
    {1, 1, 2, 2},    // C23
    {1, 1, 1, 2},    // C24
    {1, 1, 0, 2},    // C25
    {1, 1, 0, 1},    // C26
    {2, 2, 1, 2},    // C34
    {2, 2, 0, 2},    // C35
    {2, 2, 0, 1},    // C36
    {1, 2, 0, 2},    // C45
    {1, 2, 0, 1},    // C46
    {0, 1, 0, 2}     // C56
};

/* ---------------------------------------------------------------------- */

ComputeBornMatrixAtom::ComputeBornMatrixAtom(LAMMPS *lmp, int narg, char **arg) :
    Compute(lmp, narg, arg), id_virial(nullptr), temp_x(nullptr), temp_f(nullptr)
{
  if (narg < 3) error->all(FLERR, "Illegal compute born/matrix command");


  peratom_flag = 1;
  nvalues = 21;
  size_peratom_cols = nvalues;
  
  numflag = 0;
  numdelta = 0.0;

  pairflag = bondflag = angleflag = dihedflag = impflag = 0;
  if (narg == 3) {
    pairflag = bondflag = angleflag = dihedflag = impflag = 1;
  } else {
    int iarg = 3;
    while (iarg < narg) {
      if (strcmp(arg[iarg], "numdiff") == 0) {
        if (iarg + 3 > narg) error->all(FLERR, "Illegal compute born/matrix command");
        numflag = 1;
        numdelta = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
        if (numdelta <= 0.0) error->all(FLERR, "Illegal compute born/matrix command");
        id_virial = utils::strdup(arg[iarg + 2]);
        int icompute = modify->find_compute(id_virial);
        if (icompute < 0) error->all(FLERR, "Could not find compute born/matrix pressure ID");
        compute_virial = modify->compute[icompute];
        if (compute_virial->pressflag == 0)
          error->all(FLERR, "Compute born/matrix pressure ID does not compute pressure");
        iarg += 3;
      } else if (strcmp(arg[iarg], "pair") == 0) {
        pairflag = 1;
      } else if (strcmp(arg[iarg], "bond") == 0) {
        bondflag = 1;
      } else if (strcmp(arg[iarg], "angle") == 0) {
        angleflag = 1;
      } else if (strcmp(arg[iarg], "dihedral") == 0) {
        dihedflag = 1;
      } else if (strcmp(arg[iarg], "improper") == 0) {
        impflag = 1;
      } else {
        error->all(FLERR, "Illegal compute born/matrix command");
      }
      ++iarg;
    }
  }

  if (pairflag) {
    if (numflag)
      error->all(FLERR, "Illegal compute born/matrix command: cannot mix numflag and other flags");
    if (force->pair) {
      if (force->pair->born_matrix_enable == 0)
        error->all(FLERR, "Pair style {} does not support compute born/matrix", force->pair_style);
    } else {
      pairflag = 0;
    }
  }

  if (bondflag) {
    if (numflag)
      error->all(FLERR, "Illegal compute born/matrix command: cannot mix numflag and other flags");
    if (force->bond) {
      if (force->bond->born_matrix_enable == 0)
        error->all(FLERR, "Bond style {} does not support compute born/matrix", force->bond_style);
    } else {
      bondflag = 0;
    }
  }

  if (angleflag) {
    if (numflag)
      error->all(FLERR, "Illegal compute born/matrix command: cannot mix numflag and other flags");
    if (force->angle) {
      if (force->angle->born_matrix_enable == 0)
        error->all(FLERR, "Angle style {} does not support compute born/matrix",
                   force->angle_style);
    } else {
      angleflag = 0;
    }
  }

  if (dihedflag) {
    if (numflag)
      error->all(FLERR, "Illegal compute born/matrix command: cannot mix numflag and other flags");
    if (force->dihedral) {
      if (force->dihedral->born_matrix_enable == 0)
        error->all(FLERR, "Dihedral style {} does not support compute born/matrix",
                   force->dihedral_style);
    } else {
      dihedflag = 0;
    }
  }

  if (impflag) {
    if (numflag)
      error->all(FLERR, "Illegal compute born/matrix command: cannot mix numflag and other flags");
    if (force->improper) {
      if (force->improper->born_matrix_enable == 0)
        error->all(FLERR, "Improper style {} does not support compute born/matrix",
                   force->improper_style);
    } else {
      impflag = 0;
    }
  }

  if (force->kspace) {
    if (!numflag && (comm->me == 0))
      error->all(FLERR, "KSpace contribution not supported by compute born/matrix");
  }

  // Initialize some variables

  values_local = values_global = vector = nullptr;

  // this fix produces a global vector

  memory->create(vector, nvalues, "born_matrix:vector");
  memory->create(values_global, nvalues, "born_matrix:values_global");
  size_vector = nvalues;

  vector_flag = 1;
  extvector = 0;
  maxatom = 0;

  if (!numflag) {
    memory->create(values_local, nvalues, "born_matrix:values_local");
  } else {

    reallocate();

    // set fixed-point to default = center of cell

    fixedpoint[0] = 0.5 * (domain->boxlo[0] + domain->boxhi[0]);
    fixedpoint[1] = 0.5 * (domain->boxlo[1] + domain->boxhi[1]);
    fixedpoint[2] = 0.5 * (domain->boxlo[2] + domain->boxhi[2]);

    // define the cartesian indices for each strain (Voigt order)

    dirlist[0][0] = 0;
    dirlist[0][1] = 0;
    dirlist[1][0] = 1;
    dirlist[1][1] = 1;
    dirlist[2][0] = 2;
    dirlist[2][1] = 2;

    dirlist[3][0] = 1;
    dirlist[3][1] = 2;
    dirlist[4][0] = 0;
    dirlist[4][1] = 2;
    dirlist[5][0] = 0;
    dirlist[5][1] = 1;
  }
}

/* ---------------------------------------------------------------------- */

ComputeBornMatrixAtom::~ComputeBornMatrixAtom()
{
  memory->destroy(born_matrix);
  if (numflag) {
    memory->destroy(temp_x);
    memory->destroy(temp_f);
    delete[] id_virial;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeBornMatrixAtom::init()
{
  if (!numflag) {

    // need an occasional half neighbor list

    neighbor->add_request(this, NeighConst::REQ_FULL | NeighConst::REQ_OCCASIONAL);

  } else {

    // check for virial compute

    int icompute = modify->find_compute(id_virial);
    if (icompute < 0) error->all(FLERR, "Virial compute ID for compute born/matrix does not exist");
    compute_virial = modify->compute[icompute];

    // set up reverse index lookup
    // This table is used for consistency between numdiff and analytical
    // ordering of the terms.

    for (int m = 0; m < nvalues; m++) {
      int a = C_albe[m][0];
      int b = C_albe[m][1];
      revalbe[a][b] = m;
      revalbe[b][a] = m;
    }

    // reorder LAMMPS virial vector to Voigt order

    virialVtoV[0] = 0;
    virialVtoV[1] = 1;
    virialVtoV[2] = 2;
    virialVtoV[3] = 5;
    virialVtoV[4] = 4;
    virialVtoV[5] = 3;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeBornMatrixAtom::init_list(int /* id */, NeighList *ptr)
{
  list = ptr;
}

/* ----------------------------------------------------------------------
   compute output vector
------------------------------------------------------------------------- */

void ComputeBornMatrixAtom::compute_peratom()
{
  invoked_peratom = update->ntimestep;
  if (update->vflag_atom != invoked_peratom)
    error->all(FLERR, "Per-atom virial was not tallied on needed timestep");

  if (numflag) {

    // calculate per-atom Born matrix using stress finite differences

    compute_numdiff();

    // convert from pressure to energy units

    double inv_nktv2p = 1.0 / force->nktv2p;
    // this needs to loop over all atoms
    for (int m = 0; m < nvalues; m++) { born_matrix[0][m] *= inv_nktv2p; }
  }

}

/* ----------------------------------------------------------------------
  compute Born matrix using virial stress finite differences
------------------------------------------------------------------------- */

void ComputeBornMatrixAtom::compute_numdiff()
{
  int vec_index;

  // grow arrays if necessary

  reallocate();

  // store copy of current forces for owned and ghost atoms

  double **x = atom->x;
  double **f = atom->f;

  for (int i = 0; i < nall; i++)
    for (int k = 0; k < 3; k++) {
      temp_x[i][k] = x[i][k];
      temp_f[i][k] = f[i][k];
    }

  // loop over 6 strain directions
  // compute per-atom stress finite difference in each direction

  for (int idir = 0; idir < NDIR_VIRIAL; idir++) {

    // forward

    displace_atoms(nall, idir, 1.0);
    force_clear(nall);
    update_virial();
    for (int jdir = 0; jdir < NDIR_VIRIAL; jdir++) {
      vec_index = revalbe[idir][jdir];
      values_global[vec_index] = compute_virial->vector[virialVtoV[jdir]];
    }
    restore_atoms(nall, idir);

    // backward

    displace_atoms(nall, idir, -1.0);
    force_clear(nall);
    update_virial();
    for (int jdir = 0; jdir < NDIR_VIRIAL; jdir++) {
      vec_index = revalbe[idir][jdir];
      values_global[vec_index] -= compute_virial->vector[virialVtoV[jdir]];
    }
    restore_atoms(nall, idir);
  }

  // apply derivative factor

  double denominator = -0.5 / numdelta;
  for (int m = 0; m < nvalues; m++) values_global[m] *= denominator;

  // recompute virial so all virial and energy contributions are as before
  // also needed for virial stress addon contributions to Born matrix

  update_virial();

  // add on virial terms

  virial_addon();

  // restore original forces for owned and ghost atoms

  for (int i = 0; i < nall; i++)
    for (int k = 0; k < 3; k++) f[i][k] = temp_f[i][k];
}

/* ----------------------------------------------------------------------
   displace position of all owned and ghost atoms
------------------------------------------------------------------------- */

void ComputeBornMatrixAtom::displace_atoms(int nall, int idir, double magnitude)
{
  double **x = atom->x;

  // NOTE: virial_addon() expressions predicated on
  // shear strain fields (l != k) being symmetric here
  int k = dirlist[idir][0];
  int l = dirlist[idir][1];

  // axial strain

  if (l == k)
    for (int i = 0; i < nall; i++)
      x[i][k] = temp_x[i][k] + numdelta * magnitude * (temp_x[i][l] - fixedpoint[l]);

  // symmetric shear strain

  else
    for (int i = 0; i < nall; i++) {
      x[i][k] = temp_x[i][k] + 0.5 * numdelta * magnitude * (temp_x[i][l] - fixedpoint[l]);
      x[i][l] = temp_x[i][l] + 0.5 * numdelta * magnitude * (temp_x[i][k] - fixedpoint[k]);
    }
}

/* ----------------------------------------------------------------------
   restore position of all owned and ghost atoms
------------------------------------------------------------------------- */

void ComputeBornMatrixAtom::restore_atoms(int nall, int idir)
{

  // reset only idir coord

  int k = dirlist[idir][0];
  int l = dirlist[idir][1];
  double **x = atom->x;
  if (l == k)
    for (int i = 0; i < nall; i++) x[i][k] = temp_x[i][k];
  else
    for (int i = 0; i < nall; i++) {
      x[i][l] = temp_x[i][l];
      x[i][k] = temp_x[i][k];
    }
}

/* ----------------------------------------------------------------------
   evaluate potential forces and virial
   same logic as in Verlet
------------------------------------------------------------------------- */

void ComputeBornMatrixAtom::update_virial()
{
  int eflag = 0;

  // this may not be completely general
  // but it works for lj/cut and sw pair styles
  // and compute stress/atom output is unaffected

  int vflag = VIRIAL_PAIR;

  if (force->pair) force->pair->compute(eflag, vflag);

  if (atom->molecular != Atom::ATOMIC) {
    if (force->bond) force->bond->compute(eflag, vflag);
    if (force->angle) force->angle->compute(eflag, vflag);
    if (force->dihedral) force->dihedral->compute(eflag, vflag);
    if (force->improper) force->improper->compute(eflag, vflag);
  }

  if (force->kspace) force->kspace->compute(eflag, vflag);

  compute_virial->compute_vector();
}

/* ----------------------------------------------------------------------
   calculate virial stress addon terms to the Born matrix
   this is based on original code of Dr. Yubao Zhen
   described here: Comp. Phys. Comm. 183 (2012) 261-265
   as well as Yoshimoto et al., PRB, 71 (2005) 184108, Eq 15.and eq A3.
------------------------------------------------------------------------- */

void ComputeBornMatrixAtom::virial_addon()
{
  double *sigv = compute_virial->vector;

  // This way of doing is not very elegant but is correct.
  // The complete Cijkl terms are the sum of symmetric terms
  // computed in compute_numdiff and virial stress terms.
  // The viral terms are not symmetric in the tensor computation.
  // For example:
  // C2323 = s33+D2323; C3232 = s22+D3232 etc...
  // However there are two symmetry breaking when reducing
  // the 4-rank tensor to a 2-rank tensor
  // Cijkl = (Bijkl+Bjikl+Bijlk+Bjilk)/4. = (Bijkl+Bjilk)/2.
  // and when computing only the 21 independant term.

  // these expressions have been verified
  // correct to about 1e-7 compared
  // to the analytic expressions for lj/cut,
  // predicated on shear strain fields being
  // symmetric in displace_atoms()

  values_global[0] += 2.0 * sigv[0];
  values_global[1] += 2.0 * sigv[1];
  values_global[2] += 2.0 * sigv[2];

  values_global[3] += 0.5 * (sigv[1] + sigv[2]);
  values_global[4] += 0.5 * (sigv[0] + sigv[2]);
  values_global[5] += 0.5 * (sigv[0] + sigv[1]);
  values_global[6] += 0.0;
  values_global[7] += 0.0;
  values_global[8] += 0.0;
  values_global[9] += sigv[4];
  values_global[10] += sigv[3];
  values_global[11] += 0.0;
  values_global[12] += sigv[5];
  values_global[13] += 0.0;
  values_global[14] += sigv[3];
  values_global[15] += sigv[5];
  values_global[16] += sigv[4];
  values_global[17] += 0.0;
  values_global[18] += 0.5 * sigv[3];
  values_global[19] += 0.5 * sigv[4];
  values_global[20] += 0.5 * sigv[5];
}

/* ----------------------------------------------------------------------
   clear forces needed
------------------------------------------------------------------------- */

void ComputeBornMatrixAtom::force_clear(int nall)
{
  double **forces = atom->f;
  size_t nbytes = 3 * sizeof(double) * nall;
  if (nbytes) memset(&forces[0][0], 0, nbytes);
}

/* ----------------------------------------------------------------------
   reallocated local per-atoms arrays
------------------------------------------------------------------------- */

void ComputeBornMatrixAtom::reallocate()
{
  int nall = atom->nlocal + atom->nghost;
  if (nall > maxatom) {
    memory->destroy(temp_x);
    memory->destroy(temp_f);
    maxatom = atom->nmax;
    memory->create(temp_x, maxatom, 3, "born/matrix:temp_x");
    memory->create(temp_f, maxatom, 3, "born/matrix:temp_f");
  }
  
  // grow local stress array if necessary
  // needs to be atom->nmax in length

  if (atom->nmax > nmax) {
    memory->destroy(born_matrix);
    nmax = atom->nmax;
    memory->create(born_matrix, nmax, nvalues, "stress/atom:stress");
    array_atom = born_matrix;
  }
  
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double ComputeBornMatrixAtom::memory_usage()
{
  double bytes = 0.0;
  bytes += (double) 2 * maxatom * 3 * sizeof(double);
  return bytes;
}

