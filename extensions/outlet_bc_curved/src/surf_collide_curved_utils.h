/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.github.io
------------------------------------------------------------------------- */

#ifndef SPARTA_SURF_COLLIDE_CURVED_UTILS_H
#define SPARTA_SURF_COLLIDE_CURVED_UTILS_H

#include "math.h"
#include "math_extra.h"
#include "particle.h"
#include "random_knuth.h"

namespace SPARTA_NS {
namespace SurfCollideCurvedUtils {

inline void axis_vector(int dir, double *axis)
{
  axis[0] = axis[1] = axis[2] = 0.0;
  axis[dir] = 1.0;
}

inline void fallback_tangent(const double *norm, double *tangent)
{
  double axis[3];

  if (fabs(norm[0]) <= fabs(norm[1]) && fabs(norm[0]) <= fabs(norm[2])) {
    axis[0] = 1.0; axis[1] = 0.0; axis[2] = 0.0;
  } else if (fabs(norm[1]) <= fabs(norm[2])) {
    axis[0] = 0.0; axis[1] = 1.0; axis[2] = 0.0;
  } else {
    axis[0] = 0.0; axis[1] = 0.0; axis[2] = 1.0;
  }

  MathExtra::cross3(axis,norm,tangent);
  MathExtra::norm3(tangent);
}

inline void tangent_from_direction(const double *norm, int dir, double *tangent)
{
  double axis[3];
  axis_vector(dir,axis);

  const double dot = MathExtra::dot3(axis,norm);
  tangent[0] = axis[0] - dot*norm[0];
  tangent[1] = axis[1] - dot*norm[1];
  tangent[2] = axis[2] - dot*norm[2];

  if (MathExtra::lensq3(tangent) < 1.0e-12) fallback_tangent(norm,tangent);
  else MathExtra::norm3(tangent);
}

inline double tangential_coordinate(const double *x, const double *norm, int dir)
{
  double tangent[3];
  tangent_from_direction(norm,dir,tangent);
  return MathExtra::dot3(x,tangent);
}

inline void velocity_from_direction(const double *norm, int dir, double speed,
                                    double *velocity)
{
  double tangent[3];
  tangent_from_direction(norm,dir,tangent);
  velocity[0] = speed*tangent[0];
  velocity[1] = speed*tangent[1];
  velocity[2] = speed*tangent[2];
}

inline void mirror_about_velocity(const double *vin, const double *utarget,
                                  double *vout)
{
  vout[0] = 2.0*utarget[0] - vin[0];
  vout[1] = 2.0*utarget[1] - vin[1];
  vout[2] = 2.0*utarget[2] - vin[2];
}

inline void shift_normal_component(const double *vin, const double *utarget,
                                   const double *norm, double *vout)
{
  const double delta =
    2.0 * (MathExtra::dot3(utarget,norm) - MathExtra::dot3(vin,norm));

  vout[0] = vin[0] + delta*norm[0];
  vout[1] = vin[1] + delta*norm[1];
  vout[2] = vin[2] + delta*norm[2];
}

inline int remove_after_outlet_update(Particle::OnePart *ip, const double *norm,
                                      double normal_in, class RanKnuth *random)
{
  const double normal_out = MathExtra::dot3(ip->v,norm);
  if (normal_out <= 0.0) return 1;
  if (normal_in <= 0.0) return 1;
  if (random->uniform() >= normal_out/normal_in) return 1;
  return 0;
}

}
}

#endif
