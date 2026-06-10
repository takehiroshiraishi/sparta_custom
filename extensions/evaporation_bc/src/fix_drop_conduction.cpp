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

#include <algorithm>
#include <vector>

using namespace SPARTA_NS;

enum{INT,DOUBLE};

/* ---------------------------------------------------------------------- */

FixDropConduction::FixDropConduction(SPARTA *sparta, int narg, char **arg) :
  Fix(sparta,narg,arg)
{
  if (narg < 12) error->all(FLERR,"Illegal fix drop/conduction command");
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
  int remaining = narg - 12;
  if (remaining > 0 && strcmp(arg[12],"mode") && strcmp(arg[12],"relax")) {
    n = strlen(arg[6]) + 1;
    id_density_custom = new char[n];
    strcpy(id_density_custom,arg[6]);
    ioffset = 1;
  }

  int required = 12 + ioffset;
  if (narg < required) error->all(FLERR,"Illegal fix drop/conduction command");

  twall = input->numeric(FLERR,arg[6+ioffset]);
  latent = input->numeric(FLERR,arg[7+ioffset]);
  conductivity = input->numeric(FLERR,arg[8+ioffset]);
  liquid_rho = input->numeric(FLERR,arg[9+ioffset]);
  liquid_cp = input->numeric(FLERR,arg[10+ioffset]);
  if (strcmp(arg[11+ioffset],"surface") == 0) {
    surface_bins = 1;
    nbins = 0;
  } else {
    surface_bins = 0;
    nbins = input->inumeric(FLERR,arg[11+ioffset]);
  }

  if (twall <= 0.0 || latent <= 0.0 || conductivity <= 0.0 ||
      liquid_rho <= 0.0 || liquid_cp <= 0.0 || (!surface_bins && nbins <= 0))
    error->all(FLERR,"Illegal fix drop/conduction command");

  mode = TRANSIENT;
  relaxation = 1.0;
  grid2d = 0.0;
  tolerance2d = 1.0e-8;
  maxiter2d = 10000;
  int iarg = required;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"mode") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix drop/conduction command");
      if (strcmp(arg[iarg+1],"transient") == 0) mode = TRANSIENT;
      else if (strcmp(arg[iarg+1],"steady") == 0) mode = STEADY;
      else if (strcmp(arg[iarg+1],"steady2d") == 0) mode = STEADY2D;
      else error->all(FLERR,"Illegal fix drop/conduction command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"relax") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix drop/conduction command");
      relaxation = input->numeric(FLERR,arg[iarg+1]);
      if (relaxation <= 0.0 || relaxation > 1.0)
        error->all(FLERR,"Illegal fix drop/conduction command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"grid") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix drop/conduction command");
      grid2d = input->numeric(FLERR,arg[iarg+1]);
      if (grid2d <= 0.0) error->all(FLERR,"Illegal fix drop/conduction command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"tol") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix drop/conduction command");
      tolerance2d = input->numeric(FLERR,arg[iarg+1]);
      if (tolerance2d <= 0.0) error->all(FLERR,"Illegal fix drop/conduction command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"maxiter") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix drop/conduction command");
      maxiter2d = input->inumeric(FLERR,arg[iarg+1]);
      if (maxiter2d <= 0) error->all(FLERR,"Illegal fix drop/conduction command");
      iarg += 2;
    } else error->all(FLERR,"Illegal fix drop/conduction command");
  }

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

  temperature = heat_local = heat_global = volume = width = yedge = NULL;
  lower = diag = upper = rhs = cp = dp = NULL;
  temperature2d = source2d = source2d_local = NULL;
  mask2d = NULL;
  nx2d = ny2d = ncell2d = 0;
  xhi2d = 0.0;
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
  memory->destroy(yedge);
  memory->destroy(lower);
  memory->destroy(diag);
  memory->destroy(upper);
  memory->destroy(rhs);
  memory->destroy(cp);
  memory->destroy(dp);
  memory->destroy(temperature2d);
  memory->destroy(source2d);
  memory->destroy(source2d_local);
  memory->destroy(mask2d);
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
    build_geometry();
    if (mode == STEADY2D) build_geometry_2d();
    memory->create(temperature,nbins,"drop/conduction:temperature");
    memory->create(heat_local,nbins,"drop/conduction:heat_local");
    memory->create(heat_global,nbins,"drop/conduction:heat_global");
    memory->create(lower,nbins,"drop/conduction:lower");
    memory->create(diag,nbins,"drop/conduction:diag");
    memory->create(upper,nbins,"drop/conduction:upper");
    memory->create(rhs,nbins,"drop/conduction:rhs");
    memory->create(cp,nbins,"drop/conduction:cp");
    memory->create(dp,nbins,"drop/conduction:dp");
    if (mode == STEADY2D) {
      memory->create(temperature2d,ncell2d,"drop/conduction:temperature2d");
      memory->create(source2d,ncell2d,"drop/conduction:source2d");
      memory->create(source2d_local,ncell2d,"drop/conduction:source2d_local");
    }
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
  std::vector<double> yvalues_local;

  Surf::Line *lines = surf->lines;
  int nlocal = surf->nlocal;
  for (int i = 0; i < nlocal; i++) {
    if (!(lines[i].mask & groupbit)) continue;
    if (lines[i].p1[1] < local_ylo) local_ylo = lines[i].p1[1];
    if (lines[i].p2[1] < local_ylo) local_ylo = lines[i].p2[1];
    if (lines[i].p1[1] > local_yhi) local_yhi = lines[i].p1[1];
    if (lines[i].p2[1] > local_yhi) local_yhi = lines[i].p2[1];
    if (surface_bins) {
      yvalues_local.push_back(lines[i].p1[1]);
      yvalues_local.push_back(lines[i].p2[1]);
    }
  }

  MPI_Allreduce(&local_ylo,&ylo,1,MPI_DOUBLE,MPI_MIN,world);
  MPI_Allreduce(&local_yhi,&yhi,1,MPI_DOUBLE,MPI_MAX,world);
  if (yhi <= ylo) error->all(FLERR,"Fix drop/conduction could not determine droplet y extent");

  if (surface_bins) {
    int nlocal_values = static_cast<int>(yvalues_local.size());
    std::vector<int> counts(comm->nprocs);
    MPI_Allgather(&nlocal_values,1,MPI_INT,counts.data(),1,MPI_INT,world);
    std::vector<int> displs(comm->nprocs);
    int ntotal_values = 0;
    for (int i = 0; i < comm->nprocs; i++) {
      displs[i] = ntotal_values;
      ntotal_values += counts[i];
    }
    std::vector<double> yvalues(ntotal_values);
    MPI_Allgatherv(
      yvalues_local.data(),nlocal_values,MPI_DOUBLE,
      yvalues.data(),counts.data(),displs.data(),MPI_DOUBLE,world);

    std::sort(yvalues.begin(),yvalues.end());
    std::vector<double> unique_y;
    for (double y : yvalues) {
      if (unique_y.empty() || fabs(y - unique_y.back()) > 1.0e-14)
        unique_y.push_back(y);
    }
    if (unique_y.size() < 2)
      error->all(FLERR,"Fix drop/conduction could not determine surface-based bins");
    nbins = static_cast<int>(unique_y.size()) - 1;
    memory->create(yedge,nbins+1,"drop/conduction:yedge");
    for (int i = 0; i <= nbins; i++) yedge[i] = unique_y[i];
  } else {
    memory->create(yedge,nbins+1,"drop/conduction:yedge");
    for (int i = 0; i <= nbins; i++) yedge[i] = ylo + (yhi-ylo) * i / nbins;
  }
  dy = (yhi - ylo) / nbins;

  memory->create(volume,nbins,"drop/conduction:volume");
  memory->create(width,nbins+1,"drop/conduction:width");

  double *width_local;
  memory->create(width_local,nbins+1,"drop/conduction:width_local");
  for (int i = 0; i <= nbins; i++) width_local[i] = 0.0;

  for (int ibin = 0; ibin <= nbins; ibin++) {
    double y = yedge[ibin];
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
    double dyi = yedge[i+1] - yedge[i];
    volume[i] = 0.5*(width[i] + width[i+1])*dyi;
    if (volume[i] <= 0.0) volume[i] = 1.0e-300;
  }
}

/* ---------------------------------------------------------------------- */

double FixDropConduction::width_at_y(double y) const
{
  if (y <= yedge[0]) return width[0];
  if (y >= yedge[nbins]) return width[nbins];
  double *upper_edge = std::upper_bound(yedge,yedge+nbins+1,y);
  int ibin = static_cast<int>(upper_edge - yedge) - 1;
  if (ibin < 0) ibin = 0;
  if (ibin >= nbins) ibin = nbins - 1;
  double denom = yedge[ibin+1] - yedge[ibin];
  if (denom <= 0.0) return width[ibin];
  double frac = (y - yedge[ibin]) / denom;
  return width[ibin] + frac*(width[ibin+1] - width[ibin]);
}

/* ---------------------------------------------------------------------- */

int FixDropConduction::cell_index_2d(int ix, int iy) const
{
  return iy*nx2d + ix;
}

/* ---------------------------------------------------------------------- */

void FixDropConduction::build_geometry_2d()
{
  if (grid2d <= 0.0) grid2d = dy;
  if (grid2d <= 0.0) error->all(FLERR,"Illegal fix drop/conduction command");

  xhi2d = 0.0;
  for (int i = 0; i <= nbins; i++)
    if (width[i] > xhi2d) xhi2d = width[i];
  if (xhi2d <= 0.0) error->all(FLERR,"Fix drop/conduction could not determine 2d liquid width");

  nx2d = static_cast<int>(ceil(xhi2d/grid2d));
  ny2d = static_cast<int>(ceil((yhi-ylo)/grid2d));
  if (nx2d <= 0 || ny2d <= 0) error->all(FLERR,"Illegal fix drop/conduction command");
  ncell2d = nx2d * ny2d;

  memory->create(mask2d,ncell2d,"drop/conduction:mask2d");
  int ninside = 0;
  for (int iy = 0; iy < ny2d; iy++) {
    double y = ylo + (iy + 0.5)*grid2d;
    double local_width = width_at_y(y);
    for (int ix = 0; ix < nx2d; ix++) {
      double x = (ix + 0.5)*grid2d;
      int index = cell_index_2d(ix,iy);
      mask2d[index] = (x <= local_width + 1.0e-14) ? 1 : 0;
      if (mask2d[index]) ninside++;
    }
  }

  if (ninside <= 0) error->all(FLERR,"Fix drop/conduction found no 2d liquid cells");
}

/* ---------------------------------------------------------------------- */

int FixDropConduction::find_cell_2d(double x, double y) const
{
  int ix = static_cast<int>(floor(x/grid2d));
  int iy = static_cast<int>(floor((y-ylo)/grid2d));
  if (ix < 0) ix = 0;
  if (ix >= nx2d) ix = nx2d - 1;
  if (iy < 0) iy = 0;
  if (iy >= ny2d) iy = ny2d - 1;

  int index = cell_index_2d(ix,iy);
  if (mask2d[index]) return index;

  int best = -1;
  double best_dist2 = 1.0e100;
  for (int radius = 1; radius <= nx2d + ny2d; radius++) {
    for (int jy = iy-radius; jy <= iy+radius; jy++) {
      if (jy < 0 || jy >= ny2d) continue;
      for (int jx = ix-radius; jx <= ix+radius; jx++) {
        if (jx < 0 || jx >= nx2d) continue;
        if (abs(jx-ix) != radius && abs(jy-iy) != radius) continue;
        int candidate = cell_index_2d(jx,jy);
        if (!mask2d[candidate]) continue;
        double cx = (jx + 0.5)*grid2d;
        double cy = ylo + (jy + 0.5)*grid2d;
        double dist2 = (cx-x)*(cx-x) + (cy-y)*(cy-y);
        if (dist2 < best_dist2) {
          best = candidate;
          best_dist2 = dist2;
        }
      }
    }
    if (best >= 0) return best;
  }

  return -1;
}

/* ---------------------------------------------------------------------- */

void FixDropConduction::initialize_temperature()
{
  for (int i = 0; i < nbins; i++) temperature[i] = twall;
  if (mode == STEADY2D)
    for (int i = 0; i < ncell2d; i++) temperature2d[i] = twall;
}

/* ---------------------------------------------------------------------- */

int FixDropConduction::find_bin(double y) const
{
  if (y <= yedge[0]) return 0;
  if (y >= yedge[nbins]) return nbins - 1;
  double *upper = std::upper_bound(yedge,yedge+nbins+1,y);
  int ibin = static_cast<int>(upper - yedge) - 1;
  if (ibin < 0) ibin = 0;
  if (ibin >= nbins) ibin = nbins - 1;
  return ibin;
}

/* ---------------------------------------------------------------------- */

void FixDropConduction::end_of_step()
{
  if (update->ntimestep % nevery) return;

  if (mode == STEADY2D) {
    accumulate_heat_2d();
    solve_temperature_2d();
    write_surface_state_2d();
  } else {
    accumulate_heat();
    solve_temperature();
    write_surface_state();
  }
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
    int ibin = find_bin(ymid);
    double area = surf->line_size(&lines[m]);
    // SPARTA compute surf mflux is positive for net mass entering the surface
    // in this surface-collision/emission setup.  Positive condensation heats
    // the liquid; negative evaporation cools it.
    heat_local[ibin] += source_value(i) * area * latent;
  }

  MPI_Allreduce(heat_local,heat_global,nbins,MPI_DOUBLE,MPI_SUM,world);
}

/* ---------------------------------------------------------------------- */

void FixDropConduction::accumulate_heat_2d()
{
  for (int i = 0; i < ncell2d; i++) source2d_local[i] = 0.0;

  Surf::Line *lines;
  if (surf->distributed) lines = surf->mylines;
  else lines = surf->lines;

  int nown = surf->nown;
  for (int i = 0; i < nown; i++) {
    int m = surf->distributed ? i : comm->me + i*comm->nprocs;
    if (!(lines[m].mask & groupbit)) continue;

    double xmid = 0.5*(lines[m].p1[0] + lines[m].p2[0]);
    double ymid = 0.5*(lines[m].p1[1] + lines[m].p2[1]);
    double local_width = width_at_y(ymid);

    // Deposit the boundary heat flux to the nearest liquid cell just inside
    // the curved surface.  This is a conservative FV source in W per depth.
    double xin = xmid;
    double yin = ymid;
    if (xmid >= 0.5*local_width) xin = xmid - 0.25*grid2d;
    else yin = ymid - 0.25*grid2d;
    if (xin < 0.0) xin = 0.0;
    if (yin < ylo) yin = ylo;
    if (yin > yhi) yin = yhi;

    int cell = find_cell_2d(xin,yin);
    if (cell < 0) continue;
    double area = surf->line_size(&lines[m]);
    source2d_local[cell] += source_value(i) * area * latent;
  }

  MPI_Allreduce(source2d_local,source2d,ncell2d,MPI_DOUBLE,MPI_SUM,world);
}

/* ---------------------------------------------------------------------- */

void FixDropConduction::solve_temperature()
{
  for (int i = 0; i < nbins; i++) {
    double capdt = 0.0;
    if (mode == TRANSIENT) capdt = liquid_rho * liquid_cp * volume[i] / dtcond;
    double center = 0.5*(yedge[i] + yedge[i+1]);
    double down_distance = center - yedge[i];
    if (i > 0) down_distance = center - 0.5*(yedge[i-1] + yedge[i]);
    double gdown = conductivity * width[i] / down_distance;
    double gup = 0.0;
    if (i < nbins-1) {
      double up_distance = 0.5*(yedge[i+1] + yedge[i+2]) - center;
      gup = conductivity * width[i+1] / up_distance;
    }

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

  for (int i = nbins-2; i >= 0; i--) dp[i] -= cp[i]*dp[i+1];

  for (int i = 0; i < nbins; i++)
    temperature[i] = (1.0 - relaxation)*temperature[i] + relaxation*dp[i];

  for (int i = 0; i < nbins; i++)
    if (temperature[i] <= 0.0 || !isfinite(temperature[i]))
      error->all(FLERR,"Fix drop/conduction produced invalid temperature");
}

/* ---------------------------------------------------------------------- */

void FixDropConduction::solve_temperature_2d()
{
  const double conductance = conductivity;

  for (int iter = 0; iter < maxiter2d; iter++) {
    double max_delta = 0.0;

    for (int iy = 0; iy < ny2d; iy++) {
      for (int ix = 0; ix < nx2d; ix++) {
        int index = cell_index_2d(ix,iy);
        if (!mask2d[index]) continue;

        double sum_g = 0.0;
        double sum_gt = source2d[index];

        int left = (ix > 0) ? cell_index_2d(ix-1,iy) : -1;
        if (left >= 0 && mask2d[left]) {
          sum_g += conductance;
          sum_gt += conductance * temperature2d[left];
        }

        int right = (ix < nx2d-1) ? cell_index_2d(ix+1,iy) : -1;
        if (right >= 0 && mask2d[right]) {
          sum_g += conductance;
          sum_gt += conductance * temperature2d[right];
        }

        if (iy == 0) {
          // Fixed wall/base temperature one half-cell below the first liquid
          // cell center.  For square cells this gives 2*k conductance.
          sum_g += 2.0*conductance;
          sum_gt += 2.0*conductance * twall;
        } else {
          int down = cell_index_2d(ix,iy-1);
          if (mask2d[down]) {
            sum_g += conductance;
            sum_gt += conductance * temperature2d[down];
          }
        }

        int up = (iy < ny2d-1) ? cell_index_2d(ix,iy+1) : -1;
        if (up >= 0 && mask2d[up]) {
          sum_g += conductance;
          sum_gt += conductance * temperature2d[up];
        }

        if (sum_g <= 0.0) continue;
        double solved = sum_gt / sum_g;
        double updated = (1.0 - relaxation)*temperature2d[index] + relaxation*solved;
        double delta = fabs(updated - temperature2d[index]);
        if (delta > max_delta) max_delta = delta;
        temperature2d[index] = updated;
      }
    }

    double global_delta;
    MPI_Allreduce(&max_delta,&global_delta,1,MPI_DOUBLE,MPI_MAX,world);
    if (global_delta < tolerance2d) break;
  }

  for (int i = 0; i < ncell2d; i++)
    if (mask2d[i] && (temperature2d[i] <= 0.0 || !isfinite(temperature2d[i])))
      error->all(FLERR,"Fix drop/conduction produced invalid 2d temperature");
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
    int ibin = find_bin(ymid);
    double t = temperature[ibin];
    tcustom[i] = t;
    if (ncustom) ncustom[i] = saturation_pressure(t) / (update->boltz * t);
  }
  surf->estatus[tindex] = 0;
  if (nindex >= 0) surf->estatus[nindex] = 0;
}

/* ---------------------------------------------------------------------- */

void FixDropConduction::write_surface_state_2d()
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
    double xmid = 0.5*(lines[m].p1[0] + lines[m].p2[0]);
    double ymid = 0.5*(lines[m].p1[1] + lines[m].p2[1]);
    double local_width = width_at_y(ymid);

    double xin = xmid;
    double yin = ymid;
    if (xmid >= 0.5*local_width) xin = xmid - 0.25*grid2d;
    else yin = ymid - 0.25*grid2d;
    if (xin < 0.0) xin = 0.0;
    if (yin < ylo) yin = ylo;
    if (yin > yhi) yin = yhi;

    int cell = find_cell_2d(xin,yin);
    double t = twall;
    if (cell >= 0) t = temperature2d[cell];
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
  return (9.0*nbins + nbins + 1.0 + 3.0*ncell2d) * sizeof(double) +
    ncell2d * sizeof(int);
}
