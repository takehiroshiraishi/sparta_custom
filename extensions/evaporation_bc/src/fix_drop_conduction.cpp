/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
------------------------------------------------------------------------- */

#include "fix_drop_conduction.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "input.h"
#include "memory.h"
#include "modify.h"
#include "surf.h"
#include "update.h"

#include "math.h"
#include "mpi.h"
#include "stdlib.h"
#include "string.h"

using namespace SPARTA_NS;

enum{INT,DOUBLE};

/* ---------------------------------------------------------------------- */

FixDropConduction::FixDropConduction(SPARTA *sparta, int narg, char **arg) :
  Fix(sparta,narg,arg)
{
  if (narg != 12 && narg != 13) error->all(FLERR,"Illegal fix drop/conduction command");
  if (surf->implicit)
    error->all(FLERR,"Cannot use fix drop/conduction with implicit surfs");
  if (domain->dimension != 2)
    error->all(FLERR,"Fix drop/conduction currently requires 2d simulation");

  int igroup = surf->find_group(arg[2]);
  if (igroup < 0) error->all(FLERR,"Fix drop/conduction group ID does not exist");
  groupbit = surf->bitmask[igroup];

  nevery = input->inumeric(FLERR,arg[3]);
  if (nevery <= 0) error->all(FLERR,"Illegal fix drop/conduction command");

  if (strncmp(arg[4],"f_",2) != 0)
    error->all(FLERR,"Fix drop/conduction source must be a per-surf fix");
  int n = strlen(arg[4]);
  id_source = new char[n];
  strcpy(id_source,&arg[4][2]);

  char *ptr = strchr(id_source,'[');
  if (ptr) {
    if (id_source[strlen(id_source)-1] != ']')
      error->all(FLERR,"Invalid source in fix drop/conduction command");
    source_index = atoi(ptr+1);
    *ptr = '\0';
  } else source_index = 0;

  n = strlen(arg[5]) + 1;
  id_custom = new char[n];
  strcpy(id_custom,arg[5]);

  id_density_custom = NULL;
  int ioffset = 0;
  if (narg == 13) {
    n = strlen(arg[6]) + 1;
    id_density_custom = new char[n];
    strcpy(id_density_custom,arg[6]);
    ioffset = 1;
  }

  twall = input->numeric(FLERR,arg[6+ioffset]);
  latent = input->numeric(FLERR,arg[7+ioffset]);
  conductivity = input->numeric(FLERR,arg[8+ioffset]);
  liquid_rho = input->numeric(FLERR,arg[9+ioffset]);
  liquid_cp = input->numeric(FLERR,arg[10+ioffset]);
  nbins = input->inumeric(FLERR,arg[11+ioffset]);

  if (twall <= 0.0 || latent <= 0.0 || conductivity <= 0.0 ||
      liquid_rho <= 0.0 || liquid_cp <= 0.0 || nbins <= 0)
    error->all(FLERR,"Illegal fix drop/conduction command");

  ifix = -1;
  source_fix = NULL;
  tindex = surf->find_custom(id_custom);
  if (tindex < 0) {
    tindex = surf->add_custom(id_custom,DOUBLE,0);
    created_custom = 1;
  } else {
    if (surf->etype[tindex] != DOUBLE || surf->esize[tindex] != 0)
      error->all(FLERR,"Fix drop/conduction custom temperature is not a double vector");
    created_custom = 0;
  }
  if (id_density_custom) {
    nindex = surf->find_custom(id_density_custom);
    if (nindex < 0) {
      nindex = surf->add_custom(id_density_custom,DOUBLE,0);
      created_density_custom = 1;
    } else {
      if (surf->etype[nindex] != DOUBLE || surf->esize[nindex] != 0)
        error->all(FLERR,"Fix drop/conduction custom density is not a double vector");
      created_density_custom = 0;
    }
  } else {
    nindex = -1;
    created_density_custom = 0;
  }
  firstflag = 1;

  temperature = heat_local = heat_global = volume = width = NULL;
  lower = diag = upper = rhs = cp = dp = NULL;
}

/* ---------------------------------------------------------------------- */

FixDropConduction::~FixDropConduction()
{
  delete [] id_source;
  delete [] id_custom;
  delete [] id_density_custom;
  memory->destroy(temperature);
  memory->destroy(heat_local);
  memory->destroy(heat_global);
  memory->destroy(volume);
  memory->destroy(width);
  memory->destroy(lower);
  memory->destroy(diag);
  memory->destroy(upper);
  memory->destroy(rhs);
  memory->destroy(cp);
  memory->destroy(dp);
  if (created_custom && tindex >= 0) surf->remove_custom(tindex);
  if (created_density_custom && nindex >= 0) surf->remove_custom(nindex);
}

/* ---------------------------------------------------------------------- */

int FixDropConduction::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixDropConduction::init()
{
  ifix = modify->find_fix(id_source);
  if (ifix < 0) error->all(FLERR,"Could not find fix drop/conduction source fix ID");
  source_fix = modify->fix[ifix];
  if (source_fix->per_surf_flag == 0)
    error->all(FLERR,"Fix drop/conduction source does not compute per-surf info");
  if (source_index == 0 && source_fix->size_per_surf_cols != 0)
    error->all(FLERR,"Fix drop/conduction source is not a per-surf vector");
  if (source_index > 0 && source_fix->size_per_surf_cols == 0)
    error->all(FLERR,"Fix drop/conduction source is not a per-surf array");
  if (source_index > 0 && source_index > source_fix->size_per_surf_cols)
    error->all(FLERR,"Fix drop/conduction source array is accessed out-of-range");
  if (nevery % source_fix->per_surf_freq)
    error->all(FLERR,"Fix drop/conduction source not computed at compatible times");

  if (firstflag) {
    firstflag = 0;
    memory->create(temperature,nbins,"drop/conduction:temperature");
    memory->create(heat_local,nbins,"drop/conduction:heat_local");
    memory->create(heat_global,nbins,"drop/conduction:heat_global");
    memory->create(volume,nbins,"drop/conduction:volume");
    memory->create(width,nbins+1,"drop/conduction:width");
    memory->create(lower,nbins,"drop/conduction:lower");
    memory->create(diag,nbins,"drop/conduction:diag");
    memory->create(upper,nbins,"drop/conduction:upper");
    memory->create(rhs,nbins,"drop/conduction:rhs");
    memory->create(cp,nbins,"drop/conduction:cp");
    memory->create(dp,nbins,"drop/conduction:dp");
    build_geometry();
    initialize_temperature();
    write_surface_state();
  }

  dtcond = nevery * update->dt;
}

/* ---------------------------------------------------------------------- */

void FixDropConduction::build_geometry()
{
  double local_ylo = 1.0e100;
  double local_yhi = -1.0e100;

  Surf::Line *lines = surf->lines;
  int nlocal = surf->nlocal;
  for (int i = 0; i < nlocal; i++) {
    if (!(lines[i].mask & groupbit)) continue;
    if (lines[i].p1[1] < local_ylo) local_ylo = lines[i].p1[1];
    if (lines[i].p2[1] < local_ylo) local_ylo = lines[i].p2[1];
    if (lines[i].p1[1] > local_yhi) local_yhi = lines[i].p1[1];
    if (lines[i].p2[1] > local_yhi) local_yhi = lines[i].p2[1];
  }

  MPI_Allreduce(&local_ylo,&ylo,1,MPI_DOUBLE,MPI_MIN,world);
  MPI_Allreduce(&local_yhi,&yhi,1,MPI_DOUBLE,MPI_MAX,world);
  if (yhi <= ylo) error->all(FLERR,"Fix drop/conduction could not determine droplet y extent");
  dy = (yhi - ylo) / nbins;

  double *width_local;
  memory->create(width_local,nbins+1,"drop/conduction:width_local");
  for (int i = 0; i <= nbins; i++) width_local[i] = 0.0;

  for (int ibin = 0; ibin <= nbins; ibin++) {
    double y = ylo + ibin*dy;
    double xmax = 0.0;
    for (int i = 0; i < nlocal; i++) {
      if (!(lines[i].mask & groupbit)) continue;
      double y1 = lines[i].p1[1];
      double y2 = lines[i].p2[1];
      double ymin = y1 < y2 ? y1 : y2;
      double ymax = y1 > y2 ? y1 : y2;
      if (y < ymin - 1.0e-14 || y > ymax + 1.0e-14) continue;
      double x;
      if (fabs(y2-y1) < 1.0e-30) x = 0.5*(lines[i].p1[0] + lines[i].p2[0]);
      else {
        double frac = (y-y1)/(y2-y1);
        if (frac < 0.0) frac = 0.0;
        if (frac > 1.0) frac = 1.0;
        x = lines[i].p1[0] + frac*(lines[i].p2[0]-lines[i].p1[0]);
      }
      if (x > xmax) xmax = x;
    }
    width_local[ibin] = xmax;
  }

  MPI_Allreduce(width_local,width,nbins+1,MPI_DOUBLE,MPI_MAX,world);
  memory->destroy(width_local);

  for (int i = 0; i < nbins; i++) {
    volume[i] = 0.5*(width[i] + width[i+1])*dy;
    if (volume[i] <= 0.0) volume[i] = 1.0e-300;
  }
}

/* ---------------------------------------------------------------------- */

void FixDropConduction::initialize_temperature()
{
  for (int i = 0; i < nbins; i++) temperature[i] = twall;
}

/* ---------------------------------------------------------------------- */

void FixDropConduction::end_of_step()
{
  if (update->ntimestep % nevery) return;

  accumulate_heat();
  solve_temperature();
  write_surface_state();
}

/* ---------------------------------------------------------------------- */

void FixDropConduction::accumulate_heat()
{
  for (int i = 0; i < nbins; i++) heat_local[i] = 0.0;

  Surf::Line *lines;
  if (surf->distributed) lines = surf->mylines;
  else lines = surf->lines;

  int nown = surf->nown;
  for (int i = 0; i < nown; i++) {
    int m = surf->distributed ? i : comm->me + i*comm->nprocs;
    if (!(lines[m].mask & groupbit)) continue;
    double ymid = 0.5*(lines[m].p1[1] + lines[m].p2[1]);
    int ibin = static_cast<int>((ymid - ylo) / dy);
    if (ibin < 0) ibin = 0;
    if (ibin >= nbins) ibin = nbins - 1;
    double area = surf->line_size(&lines[m]);
    // SPARTA's surface mflux is positive for net evaporation out of the surface
    // and negative for net condensation into the surface.  Heat added to the
    // liquid is therefore -mflux*latent: evaporation cools, condensation heats.
    heat_local[ibin] += -source_value(i) * area * latent;
  }

  MPI_Allreduce(heat_local,heat_global,nbins,MPI_DOUBLE,MPI_SUM,world);
}

/* ---------------------------------------------------------------------- */

void FixDropConduction::solve_temperature()
{
  for (int i = 0; i < nbins; i++) {
    double capdt = liquid_rho * liquid_cp * volume[i] / dtcond;
    double gdown = conductivity * width[i] / dy;
    double gup = (i == nbins-1) ? 0.0 : conductivity * width[i+1] / dy;

    lower[i] = (i == 0) ? 0.0 : -gdown;
    upper[i] = (i == nbins-1) ? 0.0 : -gup;
    diag[i] = capdt + gdown + gup;
    rhs[i] = capdt * temperature[i] + heat_global[i];
    if (i == 0) rhs[i] += gdown * twall;
  }

  cp[0] = upper[0] / diag[0];
  dp[0] = rhs[0] / diag[0];
  for (int i = 1; i < nbins; i++) {
    double denom = diag[i] - lower[i]*cp[i-1];
    if (fabs(denom) < 1.0e-300)
      error->all(FLERR,"Fix drop/conduction singular tridiagonal solve");
    cp[i] = upper[i] / denom;
    dp[i] = (rhs[i] - lower[i]*dp[i-1]) / denom;
  }

  temperature[nbins-1] = dp[nbins-1];
  for (int i = nbins-2; i >= 0; i--) temperature[i] = dp[i] - cp[i]*temperature[i+1];

  for (int i = 0; i < nbins; i++)
    if (temperature[i] <= 0.0 || !isfinite(temperature[i]))
      error->all(FLERR,"Fix drop/conduction produced invalid temperature");
}

/* ---------------------------------------------------------------------- */

void FixDropConduction::write_surface_state()
{
  double *tcustom = surf->edvec[surf->ewhich[tindex]];
  double *ncustom = NULL;
  if (nindex >= 0) ncustom = surf->edvec[surf->ewhich[nindex]];
  Surf::Line *lines;
  if (surf->distributed) lines = surf->mylines;
  else lines = surf->lines;

  int nown = surf->nown;
  for (int i = 0; i < nown; i++) {
    int m = surf->distributed ? i : comm->me + i*comm->nprocs;
    if (!(lines[m].mask & groupbit)) continue;
    double ymid = 0.5*(lines[m].p1[1] + lines[m].p2[1]);
    int ibin = static_cast<int>((ymid - ylo) / dy);
    if (ibin < 0) ibin = 0;
    if (ibin >= nbins) ibin = nbins - 1;
    double t = temperature[ibin];
    tcustom[i] = t;
    if (ncustom) ncustom[i] = saturation_pressure(t) / (update->boltz * t);
  }
  surf->estatus[tindex] = 0;
  if (nindex >= 0) surf->estatus[nindex] = 0;
}

/* ---------------------------------------------------------------------- */

double FixDropConduction::saturation_pressure(double temp) const
{
  // Antoine equation for water, pressure in Pa, temperature in K.
  // This covers the near-room-temperature range used by these cases.
  const double A = 8.07131;
  const double B = 1730.63;
  const double C = 233.426;
  const double mmhg_to_pa = 133.32236842105263;
  double temp_c = temp - 273.15;
  return pow(10.0,A - B/(C + temp_c)) * mmhg_to_pa;
}

/* ---------------------------------------------------------------------- */

double FixDropConduction::source_value(int i) const
{
  if (source_index == 0) return source_fix->vector_surf[i];
  return source_fix->array_surf[i][source_index-1];
}

/* ---------------------------------------------------------------------- */

double FixDropConduction::memory_usage()
{
  return (9.0*nbins + nbins + 1.0) * sizeof(double);
}
