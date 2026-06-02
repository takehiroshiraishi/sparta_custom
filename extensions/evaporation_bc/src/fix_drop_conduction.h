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
  int tindex, created_custom;
  int nbins, firstflag;
  char *id_source;
  char *id_custom;

  double twall, latent, conductivity, liquid_rho, liquid_cp;
  double ylo, yhi, dy, dtcond;

  class Fix *source_fix;

  double *temperature;
  double *heat_local;
  double *heat_global;
  double *volume;
  double *width;
  double *lower;
  double *diag;
  double *upper;
  double *rhs;
  double *cp;
  double *dp;

  void build_geometry();
  void initialize_temperature();
  void accumulate_heat();
  void solve_temperature();
  void write_surface_temperature();
  double source_value(int) const;
};

}

#endif
#endif
