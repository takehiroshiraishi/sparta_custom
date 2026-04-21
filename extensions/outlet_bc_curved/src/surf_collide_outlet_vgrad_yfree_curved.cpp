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
#include "surf_collide_outlet_vgrad_yfree_curved.h"
#include "surf.h"
#include "input.h"
#include "modify.h"
#include "comm.h"
#include "update.h"
#include "random_mars.h"
#include "random_knuth.h"
#include "math_extra.h"
#include "surf_collide_curved_utils.h"
#include "error.h"

#include <iostream>

using namespace SPARTA_NS;

enum{NUMERIC,CUSTOM,VARIABLE,VAREQUAL,VARSURF};   // surf_collide classes

/* ---------------------------------------------------------------------- */

SurfCollideOutletVgradYfreeCurved::SurfCollideOutletVgradYfreeCurved(
  SPARTA *sparta, int narg, char **arg) :
  SurfCollide(sparta, narg, arg)
{
  if (narg != 8)
    error->all(FLERR,"Illegal surf_collide outletvgradyfree/curved command");

  // macroscopic velocity
  directionv = atoi(arg[2]);
  directiongrad = atoi(arg[3]);
  grad = atof(arg[4]);
  refvel = atof(arg[5]);
  width = atof(arg[6]);
  maxwidth = atof(arg[7]);
  margin = (maxwidth - width) * 0.5;

  // initialize RNG
  random = new RanKnuth(update->ranmaster->uniform());
  double seed = update->ranmaster->uniform();
  random->reset(seed,comm->me,100);

  allowreact = 0;
}

/* ---------------------------------------------------------------------- */

SurfCollideOutletVgradYfreeCurved::~SurfCollideOutletVgradYfreeCurved()
{
  if (copy) return;

  delete random;
}

/* ----------------------------------------------------------------------
   particle collision with surface with optional chemistry
   ip = particle with current x = collision pt, current v = incident v
   isurf = index of surface element
   norm = surface normal unit vector
   isr = index of reaction model if >= 0, -1 for no chemistry
   ip = reset to NULL if destroyed by chemistry
   return jp = new particle if created by chemistry
   return reaction = index of reaction (1 to N) that took place, 0 = no reaction
   resets particle(s) to post-collision outward velocity
------------------------------------------------------------------------- */

Particle::OnePart *SurfCollideOutletVgradYfreeCurved::
collide(Particle::OnePart *&ip, double &,
        int isurf, double *norm, int isr, int &reaction)
{
  nsingle++;
  double utarget[3];
  double speed;
  const double coord =
    SurfCollideCurvedUtils::tangential_coordinate(ip->x,norm,directiongrad);
  const double normal_in = -MathExtra::dot3(ip->v,norm);

  if (coord < margin) speed = refvel;
  else if (coord > maxwidth - margin) speed = width * grad + refvel;
  else speed = (coord - margin) * grad + refvel;

  SurfCollideCurvedUtils::velocity_from_direction(norm,directionv,speed,utarget);
  SurfCollideCurvedUtils::shift_normal_component(ip->v,utarget,norm,ip->v);
  if (SurfCollideCurvedUtils::remove_after_outlet_update(ip,norm,normal_in,
                                                         random))
    ip = NULL;
  
  return NULL;
}
