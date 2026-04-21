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
#include "surf_collide_outlet_vextra.h"
#include "comm.h"
#include "math_const.h"
#include "math_extra.h"
#include "modify.h"
#include "random_mars.h"
#include "random_knuth.h"
#include "surf.h"
#include "update.h"
#include "error.h"

using namespace SPARTA_NS;
using namespace MathConst;

namespace {

static const int EXTRA_SAMPLE_MAX = 256;

/* ----------------------------------------------------------------------
   ratio of the unpaired incoming flux on v in [2u0,0] to the outgoing
   flux on v > 0 for the shifted 1D Gaussian exp(-(v-u0)^2), u0 < 0
------------------------------------------------------------------------- */

double missing_flux_ratio(double u0n)
{
  const double a = -u0n;
  const double sqpi = sqrt(MY_PI);

  if (a <= 0.0) return 0.0;

  const double numerator = a * sqpi * erf(a);

  double denominator;
  if (a < 6.0) {
    denominator = 0.5 * exp(-a*a) - 0.5 * a * sqpi * erfc(a);
  } else {
    const double ainv2 = 1.0 / (a*a);
    denominator = exp(-a*a) * 0.25 * ainv2 * (1.0 - 1.5 * ainv2);
  }

  if (denominator <= 0.0) return 0.0;
  return numerator / denominator;
}

/* ----------------------------------------------------------------------
   sample from exp(-(v-u0)^2) truncated to v in [2u0,0] by inverting the
   bounded erf CDF on the symmetric interval [u0,-u0] for y = v-u0
------------------------------------------------------------------------- */

double sample_truncated_gaussian(double u0n, RanKnuth *random)
{
  const double a = -u0n;
  const double erfa = erf(a);

  if (erfa == 0.0) return u0n;

  const double target = (2.0*random->uniform() - 1.0) * erfa;
  double lo = -a;
  double hi =  a;

  for (int iter = 0; iter < 60; iter++) {
    const double mid = 0.5 * (lo + hi);
    if (erf(mid) < target) lo = mid;
    else hi = mid;
  }

  return u0n + 0.5 * (lo + hi);
}

/* ----------------------------------------------------------------------
   exact sampler for p(v) proportional to (-v) exp(-(v-u0)^2),
   v in [2u0,0], u0 < 0.  The proposal is the truncated Gaussian above,
   then accept with probability (-v)/(-2u0).
------------------------------------------------------------------------- */

double sample_extra_normal(double u0n, RanKnuth *random)
{
  for (int iter = 0; iter < EXTRA_SAMPLE_MAX; iter++) {
    const double v = sample_truncated_gaussian(u0n,random);
    if (random->uniform() <= v / (2.0*u0n)) return v;
  }

  // Defensive fallback for an effectively unreachable RNG path.
  return u0n;
}

/* ----------------------------------------------------------------------
   build tangential unit vectors from the incident velocity.  If the
   incidence is purely normal, choose an arbitrary orthogonal tangent.
------------------------------------------------------------------------- */

void tangential_basis(const double *v, const double *norm, RanKnuth *random,
                      double *tangent1, double *tangent2)
{
  const double dot = MathExtra::dot3(v,norm);

  tangent1[0] = v[0] - dot*norm[0];
  tangent1[1] = v[1] - dot*norm[1];
  tangent1[2] = v[2] - dot*norm[2];

  if (MathExtra::lensq3(tangent1) == 0.0) {
    tangent2[0] = random->uniform();
    tangent2[1] = random->uniform();
    tangent2[2] = random->uniform();
    MathExtra::cross3(norm,tangent2,tangent1);

    if (MathExtra::lensq3(tangent1) == 0.0) {
      if (fabs(norm[0]) < 0.9) {
        tangent2[0] = 1.0;
        tangent2[1] = tangent2[2] = 0.0;
      } else {
        tangent2[1] = 1.0;
        tangent2[0] = tangent2[2] = 0.0;
      }
      MathExtra::cross3(norm,tangent2,tangent1);
    }
  }

  MathExtra::norm3(tangent1);
  MathExtra::cross3(norm,tangent1,tangent2);
}

/* ----------------------------------------------------------------------
   unbiased integer copy count from a mean ratio
------------------------------------------------------------------------- */

int draw_copy_count(double ratio, RanKnuth *random)
{
  int ncopy = static_cast<int> (ratio);
  const double frac = ratio - ncopy;
  if (random->uniform() < frac) ncopy++;
  return ncopy;
}

}

enum{NUMERIC,CUSTOM,VARIABLE,VAREQUAL,VARSURF};   // surf_collide classes

/* ---------------------------------------------------------------------- */

SurfCollideOutletVExtra::SurfCollideOutletVExtra(SPARTA *sparta,
                                                 int narg, char **arg) :
  SurfCollide(sparta, narg, arg)
{
  if (narg != 5) error->all(FLERR,"Illegal surf_collide outletvextra command");

  vx = atof(arg[2]);
  vy = atof(arg[3]);
  vz = atof(arg[4]);

  random = new RanKnuth(update->ranmaster->uniform());
  double seed = update->ranmaster->uniform();
  random->reset(seed,comm->me,100);

  allowreact = 0;
}

/* ---------------------------------------------------------------------- */

SurfCollideOutletVExtra::~SurfCollideOutletVExtra()
{
  if (copy) return;

  delete random;
}

/* ----------------------------------------------------------------------
   outlet boundary with the original mirrored partner rule plus the
   additional unpaired incoming source needed when u0.n < 0.

   The old paired rule maps f(2u0-ux) = f(ux), which covers only incoming
   normal velocities v <= 2u0.  For u0 < 0 the interval v in [2u0,0]
   has no outgoing partner, so its missing incoming flux must be inserted
   explicitly with flux pdf p(v) proportional to (-v) exp(-(v-u0)^2).
------------------------------------------------------------------------- */

Particle::OnePart *SurfCollideOutletVExtra::
collide(Particle::OnePart *&ip, double &,
        int isurf, double *norm, int isr, int &reaction)
{
  nsingle++;
  reaction = 0;

  if (ip == NULL) return NULL;

  const int iorig = ip - particle->particles;

  double tangent1[3], tangent2[3];
  double *v = ip->v;
  double vincident[3];
  vincident[0] = v[0];
  vincident[1] = v[1];
  vincident[2] = v[2];
  const double dot = MathExtra::dot3(v,norm);

  tangent1[0] = v[0] - dot*norm[0];
  tangent1[1] = v[1] - dot*norm[1];
  tangent1[2] = v[2] - dot*norm[2];

  const double tangent_sq = MathExtra::lensq3(tangent1);
  const double u_p = -dot;
  const double v_p = sqrt(tangent_sq);

  double u0dotn = -(vx*norm[0] + vy*norm[1] + vz*norm[2]);
  double u0dottan = 0.0;

  if (tangent_sq > 0.0) {
    MathExtra::snormalize3(v_p,tangent1,tangent1);
    u0dottan = vx*tangent1[0] + vy*tangent1[1] + vz*tangent1[2];
  }

  const double mirrored = 2.0*u0dotn - u_p;

  // Preserve the existing outletv paired-generation logic for the mirrored
  // partner whenever that partner is incoming.
  v[0] = -mirrored*norm[0] + (2.0*u0dottan - v_p)*tangent1[0];
  v[1] = -mirrored*norm[1] + (2.0*u0dottan - v_p)*tangent1[1];
  v[2] = -mirrored*norm[2] + (2.0*u0dottan - v_p)*tangent1[2];

  if (u_p <= 0.0 || mirrored >= 0.0) {
    ip = NULL;
  } else if (random->uniform() >= fabs(mirrored/u_p)) {
    ip = NULL;
  }

  if (u0dotn >= 0.0) return NULL;

  const double ratio = missing_flux_ratio(u0dotn);
  if (ratio <= 0.0) return NULL;

  const int ncopy = draw_copy_count(ratio,random);
  if (ncopy <= 0) return NULL;

  tangential_basis(vincident,norm,random,tangent1,tangent2);

  // The extra u0.n < 0 source only fills the missing normal interval
  // v in [2u0,0].  Its tangential components are taken from the outgoing
  // particle associated with this generation event and are not resampled.
  const double vtan1 = MathExtra::dot3(vincident,tangent1);
  const double vtan2 = MathExtra::dot3(vincident,tangent2);

  for (int m = 0; m < ncopy; m++) {
    particle->clone_particle(iorig);
    Particle::OnePart *jp = &particle->particles[particle->nlocal-1];

    const double vnormal = sample_extra_normal(u0dotn,random);

    jp->v[0] = -vnormal*norm[0] + vtan1*tangent1[0] + vtan2*tangent2[0];
    jp->v[1] = -vnormal*norm[1] + vtan1*tangent1[1] + vtan2*tangent2[1];
    jp->v[2] = -vnormal*norm[2] + vtan1*tangent1[2] + vtan2*tangent2[2];
  }

  return NULL;
}
