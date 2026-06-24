/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.github.io
   Steve Plimpton, sjplimp@gmail.com, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2014) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level SPARTA directory.
------------------------------------------------------------------------- */

#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "surf_collide_outlet_v_gradT.h"
#include "particle.h"
#include "surf.h"
#include "surf_react.h"
#include "input.h"
#include "modify.h"
#include "comm.h"
#include "update.h"
#include "random_mars.h"
#include "random_knuth.h"
#include "math_extra.h"
#include "error.h"
#include "mpi.h"

using namespace SPARTA_NS;

enum{NUMERIC,CUSTOM,VARIABLE,VAREQUAL,VARSURF};   // surf_collide classes

/* ---------------------------------------------------------------------- */

SurfCollideOutletVgradT::SurfCollideOutletVgradT(SPARTA *sparta, int narg, char **arg) :
  SurfCollide(sparta, narg, arg)
{
  if (narg < 4) error->all(FLERR,"Illegal surf_collide outletvgradT command");

  tmode_grad = MODE_T2;
  t1 = t2 = deltaT = gradT = dxghost = 0.0;
  uinf = ut = wt = 0.0;
  min_un = 1.0e-12;
  max_branch = 20;
  tangent_resample = rot_sample = 0;
  warned_safe_u = 0;

  // Legacy compatibility: old local inputs used "uinf T2/T1".
  // Prefer the keyword syntax documented in custom/extensions/outlet_bc/README.md.
  if (narg == 4) {
    uinf = input->numeric(FLERR,arg[2]);
    t1 = 1.0;
    t2 = input->numeric(FLERR,arg[3]);
    if (t2 <= 0.0 || t2 > t1)
      error->all(FLERR,"Illegal surf_collide outletvgradT command");
  } else {
    int iarg = 2;
    while (iarg < narg) {
      if (strcmp(arg[iarg],"t1") == 0) {
        if (iarg+1 >= narg) error->all(FLERR,"Illegal surf_collide outletvgradT command");
        t1 = input->numeric(FLERR,arg[iarg+1]);
        iarg += 2;
      } else if (strcmp(arg[iarg],"t2") == 0) {
        if (iarg+1 >= narg) error->all(FLERR,"Illegal surf_collide outletvgradT command");
        tmode_grad = MODE_T2;
        t2 = input->numeric(FLERR,arg[iarg+1]);
        iarg += 2;
      } else if (strcmp(arg[iarg],"delta") == 0) {
        if (iarg+1 >= narg) error->all(FLERR,"Illegal surf_collide outletvgradT command");
        tmode_grad = MODE_DELTA;
        deltaT = input->numeric(FLERR,arg[iarg+1]);
        iarg += 2;
      } else if (strcmp(arg[iarg],"grad") == 0) {
        if (iarg+1 >= narg) error->all(FLERR,"Illegal surf_collide outletvgradT command");
        tmode_grad = MODE_GRAD;
        gradT = input->numeric(FLERR,arg[iarg+1]);
        iarg += 2;
      } else if (strcmp(arg[iarg],"dx") == 0) {
        if (iarg+1 >= narg) error->all(FLERR,"Illegal surf_collide outletvgradT command");
        dxghost = input->numeric(FLERR,arg[iarg+1]);
        iarg += 2;
      } else if (strcmp(arg[iarg],"uinf") == 0) {
        if (iarg+1 >= narg) error->all(FLERR,"Illegal surf_collide outletvgradT command");
        uinf = input->numeric(FLERR,arg[iarg+1]);
        iarg += 2;
      } else if (strcmp(arg[iarg],"ut") == 0) {
        if (iarg+2 >= narg) error->all(FLERR,"Illegal surf_collide outletvgradT command");
        ut = input->numeric(FLERR,arg[iarg+1]);
        wt = input->numeric(FLERR,arg[iarg+2]);
        iarg += 3;
      } else if (strcmp(arg[iarg],"min_un") == 0) {
        if (iarg+1 >= narg) error->all(FLERR,"Illegal surf_collide outletvgradT command");
        min_un = input->numeric(FLERR,arg[iarg+1]);
        iarg += 2;
      } else if (strcmp(arg[iarg],"max_branch") == 0) {
        if (iarg+1 >= narg) error->all(FLERR,"Illegal surf_collide outletvgradT command");
        max_branch = input->inumeric(FLERR,arg[iarg+1]);
        iarg += 2;
      } else if (strcmp(arg[iarg],"tangent") == 0) {
        if (iarg+1 >= narg) error->all(FLERR,"Illegal surf_collide outletvgradT command");
        if (strcmp(arg[iarg+1],"scale") == 0) tangent_resample = 0;
        else if (strcmp(arg[iarg+1],"resample") == 0) tangent_resample = 1;
        else error->all(FLERR,"Illegal surf_collide outletvgradT command");
        iarg += 2;
      } else if (strcmp(arg[iarg],"rot") == 0) {
        if (iarg+1 >= narg) error->all(FLERR,"Illegal surf_collide outletvgradT command");
        if (strcmp(arg[iarg+1],"scale") == 0) rot_sample = 0;
        else if (strcmp(arg[iarg+1],"sample") == 0) rot_sample = 1;
        else error->all(FLERR,"Illegal surf_collide outletvgradT command");
        iarg += 2;
      } else error->all(FLERR,"Illegal surf_collide outletvgradT command");
    }
  }

  t2 = compute_t2();
  if (t1 <= 0.0 || t2 <= 0.0 || t2 > t1 || min_un < 0.0 || max_branch < 0)
    error->all(FLERR,"Illegal surf_collide outletvgradT command");

  // initialize RNG
  random = new RanKnuth(update->ranmaster->uniform());
  double seed = update->ranmaster->uniform();
  random->reset(seed,comm->me,100);

  allowreact = 0;
  size_vector = 13;
}

/* ---------------------------------------------------------------------- */

SurfCollideOutletVgradT::~SurfCollideOutletVgradT()
{
  if (copy) return;

  delete random;
}

/* ----------------------------------------------------------------------
   initialize cumulative diagnostics
------------------------------------------------------------------------- */

void SurfCollideOutletVgradT::init()
{
  SurfCollide::init();
  nprocessed = ngenerated = sum_branch = max_branch_seen = 0.0;
  n_s2_bad = n_uprime_positive = n_tbad = n_capped = n_min_un = 0.0;
  last_t1 = t1;
  last_t2 = compute_t2();
  last_safe_u = 0.0;
  warned_safe_u = 0;
}

/* ----------------------------------------------------------------------
   compute outside/ghost temperature
------------------------------------------------------------------------- */

double SurfCollideOutletVgradT::compute_t2() const
{
  if (tmode_grad == MODE_T2) return t2;
  if (tmode_grad == MODE_DELTA) return t1 - deltaT;
  return t1 + gradT*dxghost;
}

/* ----------------------------------------------------------------------
   unbiased stochastic integer count
------------------------------------------------------------------------- */

int SurfCollideOutletVgradT::draw_branch_count(double mean)
{
  if (mean <= 0.0) return 0;
  int n = static_cast<int> (mean);
  if (random->uniform() < mean - n) n++;
  if (n > max_branch) {
    n = max_branch;
    n_capped++;
  }
  return n;
}

/* ----------------------------------------------------------------------
   tangential basis around outward normal
------------------------------------------------------------------------- */

void SurfCollideOutletVgradT::tangent_basis(const double *v, const double *n,
                                           double *tangent1, double *tangent2)
{
  const double dot = MathExtra::dot3(v,n);

  tangent1[0] = v[0] - dot*n[0];
  tangent1[1] = v[1] - dot*n[1];
  tangent1[2] = v[2] - dot*n[2];

  if (MathExtra::lensq3(tangent1) == 0.0) {
    if (fabs(n[0]) < 0.9) {
      tangent2[0] = 1.0;
      tangent2[1] = tangent2[2] = 0.0;
    } else {
      tangent2[1] = 1.0;
      tangent2[0] = tangent2[2] = 0.0;
    }
    MathExtra::cross3(n,tangent2,tangent1);
  }

  MathExtra::norm3(tangent1);
  MathExtra::cross3(n,tangent1,tangent2);
}

/* ----------------------------------------------------------------------
   assign mapped velocity and internal energies
------------------------------------------------------------------------- */

void SurfCollideOutletVgradT::set_velocity(Particle::OnePart *p,
                                           const double *n,
                                           const double *tangent1,
                                           const double *tangent2,
                                           double unew, double r,
                                           double vtan1_old,
                                           double vtan2_old,
                                           double erot_old)
{
  double vtan1,vtan2;

  if (tangent_resample) {
    const double mass = particle->species[p->ispecies].mass;
    const double sigma = sqrt(update->boltz*t2/mass);
    const double twopi = 8.0 * atan(1.0);
    const double theta1 = twopi*random->uniform();
    const double theta2 = twopi*random->uniform();
    vtan1 = ut + sigma * sqrt(-2.0*log(random->uniform())) * cos(theta1);
    vtan2 = wt + sigma * sqrt(-2.0*log(random->uniform())) * cos(theta2);
  } else {
    const double sqrt_r = sqrt(r);
    vtan1 = ut + sqrt_r * (vtan1_old - ut);
    vtan2 = wt + sqrt_r * (vtan2_old - wt);
  }

  p->v[0] = unew*n[0] + vtan1*tangent1[0] + vtan2*tangent2[0];
  p->v[1] = unew*n[1] + vtan1*tangent1[1] + vtan2*tangent2[1];
  p->v[2] = unew*n[2] + vtan1*tangent1[2] + vtan2*tangent2[2];

  if (rot_sample) p->erot = particle->erot(p->ispecies,t2,random);
  else p->erot = erot_old * r;
}

/* ----------------------------------------------------------------------
   gradT outlet boundary with branching shifted-Maxwellian mapping
   ip = particle with current x = collision pt, current v = incident v
   isurf = index of surface element
   norm = SPARTA surface normal; box-face normals point inward
------------------------------------------------------------------------- */

Particle::OnePart *SurfCollideOutletVgradT::
collide(Particle::OnePart *&ip, double &,
        int isurf, double *norm, int isr, int &reaction)
{
  nsingle++;
  reaction = 0;
  if (ip == NULL) return NULL;

  t2 = compute_t2();
  last_t1 = t1;
  last_t2 = t2;
  if (t1 <= 0.0 || t2 <= 0.0 || t2 > t1) {
    n_tbad++;
    ip = NULL;
    return NULL;
  }

  double outward[3];
  outward[0] = -norm[0];
  outward[1] = -norm[1];
  outward[2] = -norm[2];

  double *v = ip->v;
  const double un = MathExtra::dot3(v,outward);
  if (un <= 0.0) return NULL;

  nprocessed++;
  if (un < min_un) {
    n_min_un++;
    ip = NULL;
    return NULL;
  }

  const double r = t2/t1;
  double tangent1[3],tangent2[3];
  tangent_basis(v,outward,tangent1,tangent2);
  const double vtan1_old = MathExtra::dot3(v,tangent1);
  const double vtan2_old = MathExtra::dot3(v,tangent2);
  const double erot_old = ip->erot;

  const int iorig = ip - particle->particles;
  const double gasR = update->boltz / particle->species[ip->ispecies].mass;
  const double A = (t2 == t1) ? 0.0 : 3.0 * gasR * t2 * log(t1/t2);
  const double safe_u = sqrt(A);
  last_safe_u = safe_u;
  if (!warned_safe_u && uinf <= safe_u && t2 < t1) {
    error->warning(FLERR,"surf_collide outletvgradT Uinf is below safe cooling diagnostic threshold");
    warned_safe_u = 1;
  }

  const double du = un - uinf;
  const double S2 = r*du*du + A;
  if (S2 < 0.0) {
    n_s2_bad++;
    ip = NULL;
    return NULL;
  }

  const double uprime = uinf - sqrt(S2);
  if (uprime >= 0.0) {
    n_uprime_positive++;
    ip = NULL;
    return NULL;
  }

  const double mean = -uprime/un;
  sum_branch += mean;
  if (mean > max_branch_seen) max_branch_seen = mean;
  int nbranch = draw_branch_count(mean);
  if (nbranch <= 0) {
    ip = NULL;
    return NULL;
  }

  // Reuse original particle as the first incoming branch.
  set_velocity(ip,outward,tangent1,tangent2,uprime,r,vtan1_old,vtan2_old,erot_old);
  ngenerated++;

  for (int m = 1; m < nbranch; m++) {
    Particle::OnePart *particles_old = particle->particles;
    particle->clone_particle(iorig);
    if (particles_old != particle->particles) ip = &particle->particles[iorig];
    Particle::OnePart *jp = &particle->particles[particle->nlocal-1];
    set_velocity(jp,outward,tangent1,tangent2,uprime,r,vtan1_old,vtan2_old,erot_old);
    ngenerated++;
  }

  return NULL;
}

/* ----------------------------------------------------------------------
   vector diagnostics for stats_style sc_ID[N], 1-based from input
------------------------------------------------------------------------- */

double SurfCollideOutletVgradT::compute_vector(int i)
{
  double value = 0.0;
  if (i == 0) value = nprocessed;
  else if (i == 1) value = ngenerated;
  else if (i == 2) value = sum_branch;
  else if (i == 3) value = (nprocessed > 0.0) ? sum_branch/nprocessed : 0.0;
  else if (i == 4) value = max_branch_seen;
  else if (i == 5) value = n_uprime_positive;
  else if (i == 6) value = n_s2_bad;
  else if (i == 7) value = n_tbad;
  else if (i == 8) value = n_capped;
  else if (i == 9) value = n_min_un;
  else if (i == 10) value = last_t1;
  else if (i == 11) value = last_t2;

  if (i == 3) {
    double allprocessed = 0.0;
    double allbranch = 0.0;
    MPI_Allreduce(&nprocessed,&allprocessed,1,MPI_DOUBLE,MPI_SUM,world);
    MPI_Allreduce(&sum_branch,&allbranch,1,MPI_DOUBLE,MPI_SUM,world);
    return (allprocessed > 0.0) ? allbranch/allprocessed : 0.0;
  }

  if (i == 4) {
    double allmax = 0.0;
    MPI_Allreduce(&value,&allmax,1,MPI_DOUBLE,MPI_MAX,world);
    return allmax;
  }
  if (i == 12) value = last_safe_u;

  double allvalue = 0.0;
  MPI_Allreduce(&value,&allvalue,1,MPI_DOUBLE,MPI_SUM,world);
  if (i == 10 || i == 11 || i == 12) allvalue /= comm->nprocs;
  return allvalue;
}
