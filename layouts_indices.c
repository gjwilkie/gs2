// This can be made by just replacing nakx by ntheta0
//   from AstroGK's layouts_indices.c
#include <string.h>
#include "layouts_type.h"

int ik_idx_g_y     (struct g_layout_type*,int*);
int ik_idx_g_xy    (struct g_layout_type*,int*);
int ik_idx_g_mxy   (struct g_layout_type*,int*);
int ik_idx_g_my    (struct g_layout_type*,int*);
int ik_idx_g_xmy   (struct g_layout_type*,int*);
int imu_idx_g_m    (struct g_layout_type*,int*);
int imu_idx_g_ym   (struct g_layout_type*,int*);
int imu_idx_g_xm   (struct g_layout_type*,int*);
int imu_idx_g_xym  (struct g_layout_type*,int*);
int imu_idx_g_yxm  (struct g_layout_type*,int*);
int is_idx_g_yxms  (struct g_layout_type*,int*);
int is_idx_g_xyms  (struct g_layout_type*,int*);
int is_idx_g_xmys  (struct g_layout_type*,int*);
int is_idx_g_ymxs  (struct g_layout_type*,int*);
int is_idx_g_myxs  (struct g_layout_type*,int*);
int is_idx_g_mxys  (struct g_layout_type*,int*);
int idx_g_yxms     (struct g_layout_type*,int*);
int idx_g_xyms     (struct g_layout_type*,int*,int*,int*,int*,int*);
int idx_g_xmys     (struct g_layout_type*,int*,int*,int*,int*,int*);
int idx_g_ymxs     (struct g_layout_type*,int*,int*,int*,int*,int*);
int idx_g_myxs     (struct g_layout_type*,int*,int*,int*,int*,int*);
int idx_g_mxys     (struct g_layout_type*,int*,int*,int*,int*,int*);

int (*pik_g)  (struct g_layout_type*,int*);
int (*pit_g)  (struct g_layout_type*,int*);
int (*pimu_g) (struct g_layout_type*,int*);
int (*pis_g)  (struct g_layout_type*,int*);
int (*pidx_g) (struct g_layout_type*,int*,int*,int*,int*,int*);

# ifdef NO_UNDERSCORE
int init_indices_glo_c (char *layout)
# elif UNDERSCORE != 2
int init_indices_glo_c_ (char *layout)
# else
int init_indices_glo_c__ (char *layout)
# endif
{
  int ierr;
  ierr = 0;
  if (!strcmp(layout,"yxms")) {
    pik_g = ik_idx_g_y;
    pimu_g = imu_idx_g_yxm;
    pis_g = is_idx_g_yxms;
    pidx_g = idx_g_yxms;
  } else if (!strcmp(layout,"xyms")) {
    pik_g = ik_idx_g_xy;
    pimu_g = imu_idx_g_xym;
    pis_g = is_idx_g_xyms;
    pidx_g = idx_g_xyms;
  } else if (!strcmp(layout,"xmys")) {
    pik_g = ik_idx_g_xmy;
    pimu_g = imu_idx_g_xm;
    pis_g = is_idx_g_xmys;
    pidx_g = idx_g_xmys;
  } else if (!strcmp(layout,"ymxs")) {
    pik_g = ik_idx_g_y;
    pimu_g = imu_idx_g_y;
    pis_g = is_idx_g_ymxs;
    pidx_g = idx_g_ymxs;
  } else if (!strcmp(layout,"mxys")) {
    pik_g = ik_idx_g_mxy;
    pimu_g = imu_idx_g_m;
    pis_g = is_idx_g_mxys;
    pidx_g = idx_g_mxys;
  } else if (!strcmp(layout,"myxs")) {
    pik_g = ik_idx_g_my;
    pimu_g = imu_idx_g_m;
    pis_g = is_idx_g_myxs;
    pidx_g = idx_g_myxs;
  } else {
    ierr = 1;
  }
  return ierr;
}

# ifdef NO_UNDERSCORE
int ik_idx_g_c (struct g_layout_type *lo, int *i)
# elif UNDERSCORE != 2
int ik_idx_g_c_ (struct g_layout_type *lo, int *i)
# else
int ik_idx_g_c__ (struct g_layout_type *lo, int *i)
# endif
{ return (*pik_g)(lo,i); }

# ifdef NO_UNDERSCORE
int imu_idx_g_c (struct g_layout_type *lo, int *i)
# elif UNDERSCORE != 2
int imu_idx_g_c_ (struct g_layout_type *lo, int *i)
# else
int imu_idx_g_c__ (struct g_layout_type *lo, int *i)
# endif
{ return (*pimu_g)(lo,i); }

# ifdef NO_UNDERSCORE
int is_idx_g_c (struct g_layout_type *lo, int *i)
# elif UNDERSCORE != 2
int is_idx_g_c_ (struct g_layout_type *lo, int *i)
# else
int is_idx_g_c__ (struct g_layout_type *lo, int *i)
# endif
{ return (*pis_g)(lo,i); }

# ifdef NO_UNDERSCORE
int idx_g_c (struct g_layout_type *lo, int *ik, int *it, int *il, int *ie, int *is)
# elif UNDERSCORE != 2
int idx_g_c_ (struct g_layout_type *lo, int *ik, int *it, int *il, int *ie, int *is)
# else
int idx_g_c__ (struct g_layout_type *lo, int *ik, int *it, int *il, int *ie, int *is)
# endif
{ return (*pidx_g)(lo,ik,it,imu,is); }

/* Internal functions */

/* dist_fn layout functions */

int ik_idx_g_y (struct g_layout_type *lo, int *i)
{ return 1 + (*i - lo->llim_world) % lo->naky; }
int ik_idx_g_xy (struct g_layout_type *lo, int *i)
{ return 1 + (*i - lo->llim_world)/lo->ntheta0 % lo->naky; }
int ik_idx_g_mxy (struct g_layout_type *lo, int *i)
{ return 1 + (*i - lo->llim_world)/lo->nmu/lo->ntheta0 % lo->naky; }
int ik_idx_g_my (struct g_layout_type *lo, int *i)
{ return 1 + (*i - lo->llim_world)/lo->nmu % lo->naky; }
int ik_idx_g_xmy (struct g_layout_type *lo, int *i)
{ return 1 + (*i - lo->llim_world)/lo->ntheta0/lo->nmu % lo->naky; }

int imu_idx_g_m (struct g_layout_type *lo, int *i)
{ return 1 + (*i - lo->llim_world) % lo->nmu; }
int imu_idx_g_ym (struct g_layout_type *lo, int *i)
{ return 1 + (*i - lo->llim_world)/lo->naky % lo->nmu; }
int imu_idx_g_xm (struct g_layout_type *lo, int *i)
{ return 1 + (*i - lo->llim_world)/lo->ntheta0 % lo->nmu; }
int imu_idx_g_xym (struct g_layout_type *lo, int *i)
{ return 1 + (*i - lo->llim_world)/lo->ntheta0/lo->naky % lo->nmu; }
int imu_idx_g_yxm (struct g_layout_type *lo, int *i)
{ return 1 + (*i - lo->llim_world)/lo->naky/lo->ntheta0 % lo->nmu; }

int is_idx_g_yxms (struct g_layout_type *lo, int *i)
{ return 1 + (*i - lo->llim_world)/lo->naky/lo->ntheta0/lo->nmu % lo->nspec; }
int is_idx_g_xyms (struct g_layout_type *lo, int *i)
{ return 1 + (*i - lo->llim_world)/lo->ntheta0/lo->naky/lo->nmu % lo->nspec; }
int is_idx_g_xmys (struct g_layout_type *lo, int *i)
{ return 1 + (*i - lo->llim_world)/lo->ntheta0/lo->nmu/lo->naky % lo->nspec; }
int is_idx_g_ymxs (struct g_layout_type *lo, int *i)
{ return 1 + (*i - lo->llim_world)/lo->naky/lo->nmu/lo->ntheta0 % lo->nspec; }
int is_idx_g_myxs (struct g_layout_type *lo, int *i)
{ return 1 + (*i - lo->llim_world)/lo->nmu/lo->naky/lo->ntheta0 % lo->nspec; }
int is_idx_g_mxys (struct g_layout_type *lo, int *i)
{ return 1 + (*i - lo->llim_world)/lo->nmu/lo->ntheta0/lo->naky % lo->nspec; }

int idx_g_yxms (struct g_layout_type *lo, int *ik, int *it, int *imu, int *is)
{ return *ik-1 + lo->naky*(*it-1 + lo->ntheta0*(*imu-1 + lo->nmu*(*is-1))); }
int idx_g_xyms (struct g_layout_type *lo, int *ik, int *it, int *imu, int *is)
{ return *it-1 + lo->ntheta0*(*ik-1 + lo->naky*(*imu-1 + lo->nmu*(*is-1))); }
int idx_g_xmys (struct g_layout_type *lo, int *ik, int *it, int *imu, int *is)
{ return *it-1 + lo->ntheta0*(*imu-1 + lo->nmu*(*ik-1 + lo->naky*(*is-1))); }
int idx_g_ymxs (struct g_layout_type *lo, int *ik, int *it, int *imu, int *is)
{ return *ik-1 + lo->naky*(*imu-1 + lo->nmu*(*it-1 + lo->ntheta0*(*is-1))); }
int idx_g_myxs (struct g_layout_type *lo, int *ik, int *it, int *imu, int *is)
{ return *imu-1 + lo->nmu*(*ik-1 + lo->naky*(*it-1 + lo->ntheta0*(*is-1))); }
int idx_g_mxys (struct g_layout_type *lo, int *ik, int *it, int *imu, int *is)
{ return *imu-1 + lo->nmu*(*it-1 + lo->ntheta0*(*ik-1 + lo->naky*(*is-1))); }
