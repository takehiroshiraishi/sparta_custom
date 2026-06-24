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

#ifdef SURF_COLLIDE_CLASS

SurfCollideStyle(outletvgradT,SurfCollideOutletVgradT)

#else

#ifndef SPARTA_SURF_COLLIDE_OUTLET_V_GRADT_H
#define SPARTA_SURF_COLLIDE_OUTLET_V_GRADT_H

#include "surf_collide.h"
#include "particle.h"

namespace SPARTA_NS {

class SurfCollideOutletVgradT : public SurfCollide {
 public:
  SurfCollideOutletVgradT(class SPARTA *, int, char **);
  SurfCollideOutletVgradT(class SPARTA *sparta) : SurfCollide(sparta) {} // needed for Kokkos
  virtual ~SurfCollideOutletVgradT();
  void init();
  Particle::OnePart *collide(Particle::OnePart *&, double &,
                             int, double *, int, int &);
  double compute_vector(int);

 protected:
  enum{MODE_T2,MODE_DELTA,MODE_GRAD};

  int tmode_grad;
  double t1,t2,deltaT,gradT,dxghost;
  double uinf,ut,wt;
  double min_un;
  int max_branch;
  int tangent_resample,rot_sample;
  int warned_safe_u;

  double nprocessed,ngenerated,sum_branch,max_branch_seen;
  double n_s2_bad,n_uprime_positive,n_tbad,n_capped,n_min_un;
  double last_t1,last_t2,last_safe_u;

  class RanKnuth *random;

  double compute_t2() const;
  int draw_branch_count(double);
  void tangent_basis(const double *, const double *, double *, double *);
  void set_velocity(Particle::OnePart *, const double *, const double *,
                    const double *, double, double, double, double, double);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPARTA to see the offending line.

*/
