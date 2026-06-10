/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(drop/conduction,FixDropConduction)

#else

#ifndef SPARTA_FIX_DROP_CONDUCTION_H
#define SPARTA_FIX_DROP_CONDUCTION_H

#include "fix.h"

namespace SPARTA_NS {

enum{TRANSIENT,STEADY,STEADY2D};

class FixDropConduction : public Fix {
 public:
  FixDropConduction(class SPARTA *, int, char **);
  virtual ~FixDropConduction();
  int setmask();
  void init();
  void end_of_step();
  double memory_usage();

 private:
  int groupbit;
  int ifix, source_index;
  int tindex, nindex, created_custom, created_density_custom;
  int nbins, firstflag, mode, surface_bins;
  int nx2d, ny2d, ncell2d, maxiter2d;
  char *id_source;
  char *id_custom;
  char *id_density_custom;

  double twall, latent, conductivity, liquid_rho, liquid_cp, relaxation;
  double ylo, yhi, dy, dtcond;
  double grid2d, tolerance2d, xhi2d;

  class Fix *source_fix;

  double *temperature;
  double *heat_local;
  double *heat_global;
  double *volume;
  double *width;
  double *yedge;
  double *lower;
  double *diag;
  double *upper;
  double *rhs;
  double *cp;
  double *dp;
  double *temperature2d;
  double *source2d;
  double *source2d_local;
  int *mask2d;

  void build_geometry();
  void build_geometry_2d();
  int find_bin(double) const;
  int cell_index_2d(int, int) const;
  int find_cell_2d(double, double) const;
  double width_at_y(double) const;
  void initialize_temperature();
  void accumulate_heat();
  void accumulate_heat_2d();
  void solve_temperature();
  void solve_temperature_2d();
  void write_surface_state();
  void write_surface_state_2d();
  double saturation_pressure(double) const;
  double source_value(int) const;
};

}

#endif
#endif
