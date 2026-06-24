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
#include "surf_collide_outlet_v_gradT.h"
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

#include <iostream>

using namespace SPARTA_NS;

enum{NUMERIC,CUSTOM,VARIABLE,VAREQUAL,VARSURF};   // surf_collide classes

/* ---------------------------------------------------------------------- */

SurfCollideOutletVgradT::SurfCollideOutletVgradT(SPARTA *sparta, int narg, char **arg) :
  SurfCollide(sparta, narg, arg)
{
  if (narg != 5 ) error->all(FLERR,"Illegal surf_collide outletv command");

  // macroscopic velocity
  uinf = atof(arg[2]); // velocity in normal direction
  C = atof(arg[3]); // C = T2/T1

  // initialize RNG
  random = new RanKnuth(update->ranmaster->uniform());
  double seed = update->ranmaster->uniform();
  random->reset(seed,comm->me,100);

  allowreact = 0;
}

/* ---------------------------------------------------------------------- */

SurfCollideOutletVgradT::~SurfCollideOutletVgradT()
{
  if (copy) return;

  delete random;
}

/* ----------------------------------------------------------------------
   general outlet boundary condition with prescribed velocity
   ip = particle with current x = collision pt, current v = incident v
   isurf = index of surface element
   norm = surface normal unit vector
------------------------------------------------------------------------- */

Particle::OnePart *SurfCollideOutletVgradT::
collide(Particle::OnePart *&ip, double &,
        int isurf, double *norm, int isr, int &reaction)
{
  nsingle++;

  double tangent1[3];
  double *v = ip->v; // particle velocity
  double dot = MathExtra::dot3(v,norm);

  tangent1[0] = v[0] - dot*norm[0];
  tangent1[1] = v[1] - dot*norm[1];
  tangent1[2] = v[2] - dot*norm[2];

  double u_p = -dot; // > 0 if particle is moving into surface
  double v_p = MathExtra::len3(tangent1);

  MathExtra::norm3(tangent1); // normalize vector

  double new_un = ((1+sqrt(C))*uinf - sqrt(C)*u_p);

  v[0] = new_un*-norm[0] + sqrt(C)*tangent1[0];
  v[1] = new_un*-norm[1] + sqrt(C)*tangent1[1];
  v[2] = new_un*-norm[2] + sqrt(C)*tangent1[2];

  if( new_un >= 0){
    //remove particle
    ip = NULL;
  }else if( random->uniform() * u_p * C >= abs(new_un)){
    //remove particle
    ip = NULL;
  }
  
  return NULL;
}
