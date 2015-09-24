// This can be made by just replacing nakx by ntheta0
//   from AstroGK's layouts_indices.c
#include <string.h>
#include "layouts_type.h"

int ik_idx_g_y     (struct g_layout_type*,int*);
int ik_idx_g_lexy  (struct g_layout_type*,int*);
int ik_idx_g_lxy   (struct g_layout_type*,int*);
int ik_idx_g_ly    (struct g_layout_type*,int*);
int it_idx_g_yx    (struct g_layout_type*,int*);
int it_idx_g_lex   (struct g_layout_type*,int*);
int it_idx_g_lx    (struct g_layout_type*,int*);
int it_idx_g_lyx   (struct g_layout_type*,int*);
int il_idx_g_yxel  (struct g_layout_type*,int*);
int il_idx_g_yxl   (struct g_layout_type*,int*);
int il_idx_g_l     (struct g_layout_type*,int*);
int ie_idx_g_yxe   (struct g_layout_type*,int*);
int ie_idx_g_yxle  (struct g_layout_type*,int*);
int ie_idx_g_le    (struct g_layout_type*,int*);
int ie_idx_g_lxye  (struct g_layout_type*,int*);
int ie_idx_g_lyxe  (struct g_layout_type*,int*);
int is_idx_g_yxels (struct g_layout_type*,int*);
int is_idx_g_yxles (struct g_layout_type*,int*);
int is_idx_g_lexys (struct g_layout_type*,int*);
int is_idx_g_lxyes (struct g_layout_type*,int*);
int is_idx_g_lyxes (struct g_layout_type*,int*);
int idx_g_yxels    (struct g_layout_type*,int*,int*,int*,int*,int*);
int idx_g_yxles    (struct g_layout_type*,int*,int*,int*,int*,int*);
int idx_g_lexys    (struct g_layout_type*,int*,int*,int*,int*,int*);
int idx_g_lxyes    (struct g_layout_type*,int*,int*,int*,int*,int*);
int idx_g_lyxes    (struct g_layout_type*,int*,int*,int*,int*,int*);

int ik_idx_lz_y    (struct lz_layout_type*,int*);
int ik_idx_lz_exy  (struct lz_layout_type*,int*);
int ik_idx_lz_xy   (struct lz_layout_type*,int*);
int it_idx_lz_yx   (struct lz_layout_type*,int*);
int it_idx_lz_ex   (struct lz_layout_type*,int*);
int it_idx_lz_x    (struct lz_layout_type*,int*);
int is_idx_lz_yxes (struct lz_layout_type*,int*);
int is_idx_lz_exys (struct lz_layout_type*,int*);
int is_idx_lz_xyes (struct lz_layout_type*,int*);
int ie_idx_lz_yxe  (struct lz_layout_type*,int*);
int ie_idx_lz_e    (struct lz_layout_type*,int*);
int ie_idx_lz_xye  (struct lz_layout_type*,int*);
int idx_lz_yxes    (struct lz_layout_type*,int*,int*,int*,int*,int*);
int idx_lz_exys    (struct lz_layout_type*,int*,int*,int*,int*,int*);
int idx_lz_xyes    (struct lz_layout_type*,int*,int*,int*,int*,int*);

int ik_idx_e_y     (struct e_layout_type*,int*);
int ik_idx_e_lxy   (struct e_layout_type*,int*);
int ik_idx_e_ly    (struct e_layout_type*,int*);
int it_idx_e_yx    (struct e_layout_type*,int*);
int it_idx_e_lx    (struct e_layout_type*,int*);
int it_idx_e_lyx   (struct e_layout_type*,int*);
int il_idx_e_yxl   (struct e_layout_type*,int*);
int il_idx_e_l     (struct e_layout_type*,int*);
int is_idx_e_yxls  (struct e_layout_type*,int*);
int is_idx_e_lxys  (struct e_layout_type*,int*);
int is_idx_e_lyxs  (struct e_layout_type*,int*);
int idx_e_yxls     (struct e_layout_type*,int*,int*,int*,int*,int*,int*);
int idx_e_lxys     (struct e_layout_type*,int*,int*,int*,int*,int*,int*);
int idx_e_lyxs     (struct e_layout_type*,int*,int*,int*,int*,int*,int*);

int ik_idx_le_y    (struct le_layout_type*,int*);
int ik_idx_le_xy   (struct le_layout_type*,int*);
int it_idx_le_x    (struct le_layout_type*,int*);
int it_idx_le_yx   (struct le_layout_type*,int*);
int idx_le_yxs     (struct le_layout_type*,int*,int*,int*,int*);
int idx_le_xys     (struct le_layout_type*,int*,int*,int*,int*);

int (*pik_g)  (struct g_layout_type*,int*);
int (*pit_g)  (struct g_layout_type*,int*);
int (*pil_g)  (struct g_layout_type*,int*);
int (*pie_g)  (struct g_layout_type*,int*);
int (*pis_g)  (struct g_layout_type*,int*);
int (*pidx_g) (struct g_layout_type*,int*,int*,int*,int*,int*);

int (*pik_lz)  (struct lz_layout_type*,int*);
int (*pit_lz)  (struct lz_layout_type*,int*);
int (*pie_lz)  (struct lz_layout_type*,int*);
int (*pis_lz)  (struct lz_layout_type*,int*);
int (*pidx_lz) (struct lz_layout_type*,int*,int*,int*,int*,int*);

int (*pik_e)  (struct e_layout_type*,int*);
int (*pit_e)  (struct e_layout_type*,int*);
int (*pil_e)  (struct e_layout_type*,int*);
int (*pis_e)  (struct e_layout_type*,int*);
int (*pidx_e) (struct e_layout_type*,int*,int*,int*,int*,int*,int*);

int (*pik_le)  (struct le_layout_type*,int*);
int (*pit_le)  (struct le_layout_type*,int*);
int (*pis_le)  (struct le_layout_type*,int*);
int (*pidx_le) (struct le_layout_type*,int*,int*,int*,int*);

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
  if (!strcmp(layout,"yxels")) {
    pik_g = ik_idx_g_y;
    pit_g = it_idx_g_yx;
    pil_g = il_idx_g_yxel;
    pie_g = ie_idx_g_yxe;
    pis_g = is_idx_g_yxels;
    pidx_g = idx_g_yxels;
  } else if (!strcmp(layout,"yxles")) {
    pik_g = ik_idx_g_y;
    pit_g = it_idx_g_yx;
    pil_g = il_idx_g_yxl;
    pie_g = ie_idx_g_yxle;
    pis_g = is_idx_g_yxles;
    pidx_g = idx_g_yxles;
  } else if (!strcmp(layout,"lexys")) {
    pik_g = ik_idx_g_lexy;
    pit_g = it_idx_g_lex;
    pil_g = il_idx_g_l;
    pie_g = ie_idx_g_le;
    pis_g = is_idx_g_lexys;
    pidx_g = idx_g_lexys;
  } else if (!strcmp(layout,"lxyes")) {
    pik_g = ik_idx_g_lxy;
    pit_g = it_idx_g_lx;
    pil_g = il_idx_g_l;
    pie_g = ie_idx_g_lxye;
    pis_g = is_idx_g_lxyes;
    pidx_g = idx_g_lxyes;
  } else if (!strcmp(layout,"lyxes")) {
    pik_g = ik_idx_g_ly;
    pit_g = it_idx_g_lyx;
    pil_g = il_idx_g_l;
    pie_g = ie_idx_g_lyxe;
    pis_g = is_idx_g_lyxes;
    pidx_g = idx_g_lyxes;
  } else {
    ierr = 1;
  }
  return ierr;
}

# ifdef NO_UNDERSCORE
int init_indices_lzlo_c (char *layout)
# elif UNDERSCORE != 2
int init_indices_lzlo_c_ (char *layout)
# else
int init_indices_lzlo_c__ (char *layout)
# endif
{
  int ierr;
  ierr = 0;
  if (!strcmp(layout,"yxels")) {
    pik_lz = ik_idx_lz_y;
    pit_lz = it_idx_lz_yx;
    pie_lz = ie_idx_lz_yxe;
    pis_lz = is_idx_lz_yxes;
    pidx_lz = idx_lz_yxes;
  } else if (!strcmp(layout,"yxles")) {
    pik_lz = ik_idx_lz_y;
    pit_lz = it_idx_lz_yx;
    pie_lz = ie_idx_lz_yxe;
    pis_lz = is_idx_lz_yxes;
    pidx_lz = idx_lz_yxes;
  } else if (!strcmp(layout,"lexys")) {
    pik_lz = ik_idx_lz_exy;
    pit_lz = it_idx_lz_ex;
    pie_lz = ie_idx_lz_e;
    pis_lz = is_idx_lz_exys;
    pidx_lz = idx_lz_exys;
  } else if (!strcmp(layout,"lxyes")) {
    pik_lz = ik_idx_lz_xy;
    pit_lz = it_idx_lz_x;
    pie_lz = ie_idx_lz_xye;
    pis_lz = is_idx_lz_xyes;
    pidx_lz = idx_lz_xyes;
  } else if (!strcmp(layout,"lyxes")) {
    pik_lz = ik_idx_lz_y;
    pit_lz = it_idx_lz_yx;
    pie_lz = ie_idx_lz_yxe;
    pis_lz = is_idx_lz_yxes;
    pidx_lz = idx_lz_yxes;
  } else {
    ierr = 1;
  }
  return ierr;
}

# ifdef NO_UNDERSCORE
int init_indices_elo_c (char *layout)
# elif UNDERSCORE != 2
int init_indices_elo_c_ (char *layout)
# else
int init_indices_elo_c__ (char *layout)
# endif
{
  int ierr;
  ierr = 0;
  if (!strcmp(layout,"yxels")) {
    pik_e = ik_idx_e_y;
    pit_e = it_idx_e_yx;
    pil_e = il_idx_e_yxl;
    pis_e = is_idx_e_yxls;
    pidx_e = idx_e_yxls;
  } else if (!strcmp(layout,"yxles")) {
    pik_e = ik_idx_e_y;
    pit_e = it_idx_e_yx;
    pil_e = il_idx_e_yxl;
    pis_e = is_idx_e_yxls;
    pidx_e = idx_e_yxls;
  } else if (!strcmp(layout,"lexys")) {
    pik_e = ik_idx_e_lxy;
    pit_e = it_idx_e_lx;
    pil_e = il_idx_e_l;
    pis_e = is_idx_e_lxys;
    pidx_e = idx_e_lxys;
  } else if (!strcmp(layout,"lxyes")) {
    pik_e = ik_idx_e_lxy;
    pit_e = it_idx_e_lx;
    pil_e = il_idx_e_l;
    pis_e = is_idx_e_lxys;
    pidx_e = idx_e_lxys;
  } else if (!strcmp(layout,"lyxes")) {
    pik_e = ik_idx_e_ly;
    pit_e = it_idx_e_lyx;
    pil_e = il_idx_e_l;
    pis_e = is_idx_e_lyxs;
    pidx_e = idx_e_lyxs;
  } else {
    ierr = 1;
  }
  return ierr;
}

# ifdef NO_UNDERSCORE
int init_indices_lelo_c (char *layout)
# elif UNDERSCORE != 2
int init_indices_lelo_c_ (char *layout)
# else
int init_indices_lelo_c__ (char *layout)
# endif
{
  int ierr;
  ierr = 0;
  if (!strcmp(layout,"yxels")) {
    pik_le = ik_idx_le_y;
    pit_le = it_idx_le_yx;
    pidx_le = idx_le_yxs;
  } else if (!strcmp(layout,"yxles")) {
    pik_le = ik_idx_le_y;
    pit_le = it_idx_le_yx;
    pidx_le = idx_le_yxs;
  } else if (!strcmp(layout,"lexys")) {
    pik_le = ik_idx_le_xy;
    pit_le = it_idx_le_x;
    pidx_le = idx_le_xys;
  } else if (!strcmp(layout,"lxyes")) {
    pik_le = ik_idx_le_xy;
    pit_le = it_idx_le_x;
    pidx_le = idx_le_xys;
  } else if (!strcmp(layout,"lyxes")) {
    pik_le = ik_idx_le_y;
    pit_le = it_idx_le_yx;
    pidx_le = idx_le_yxs;
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
int it_idx_g_c (struct g_layout_type *lo, int *i)
# elif UNDERSCORE != 2
int it_idx_g_c_ (struct g_layout_type *lo, int *i)
# else
int it_idx_g_c__ (struct g_layout_type *lo, int *i)
# endif
{ return (*pit_g)(lo,i); }

# ifdef NO_UNDERSCORE
int il_idx_g_c (struct g_layout_type *lo, int *i)
# elif UNDERSCORE != 2
int il_idx_g_c_ (struct g_layout_type *lo, int *i)
# else
int il_idx_g_c__ (struct g_layout_type *lo, int *i)
# endif
{ return (*pil_g)(lo,i); }

# ifdef NO_UNDERSCORE
int ie_idx_g_c (struct g_layout_type *lo, int *i)
# elif UNDERSCORE != 2
int ie_idx_g_c_ (struct g_layout_type *lo, int *i)
# else
int ie_idx_g_c__ (struct g_layout_type *lo, int *i)
# endif
{ return (*pie_g)(lo,i); }

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
{ return (*pidx_g)(lo,ik,it,il,ie,is); }

# ifdef NO_UNDERSCORE
int ik_idx_lz_c (struct lz_layout_type *lo, int *i)
# elif UNDERSCORE != 2
int ik_idx_lz_c_ (struct lz_layout_type *lo, int *i)
# else
int ik_idx_lz_c__ (struct lz_layout_type *lo, int *i)
# endif
{ return (*pik_lz)(lo,i); }

# ifdef NO_UNDERSCORE
int it_idx_lz_c (struct lz_layout_type *lo, int *i)
# elif UNDERSCORE != 2
int it_idx_lz_c_ (struct lz_layout_type *lo, int *i)
# else
int it_idx_lz_c__ (struct lz_layout_type *lo, int *i)
# endif
{ return (*pit_lz)(lo,i); }

# ifdef NO_UNDERSCORE
int ie_idx_lz_c (struct lz_layout_type *lo, int *i)
# elif UNDERSCORE != 2
int ie_idx_lz_c_ (struct lz_layout_type *lo, int *i)
# else
int ie_idx_lz_c__ (struct lz_layout_type *lo, int *i)
# endif
{ return (*pie_lz)(lo,i); }

# ifdef NO_UNDERSCORE
int is_idx_lz_c (struct lz_layout_type *lo, int *i)
# elif UNDERSCORE != 2
int is_idx_lz_c_ (struct lz_layout_type *lo, int *i)
# else
int is_idx_lz_c__ (struct lz_layout_type *lo, int *i)
# endif
{ return (*pis_lz)(lo,i); }

# ifdef NO_UNDERSCORE
int idx_lz_c (struct lz_layout_type *lo, int *ig, int *ik, int *it, int *ie, int *is)
# elif UNDERSCORE != 2
int idx_lz_c_ (struct lz_layout_type *lo, int *ig, int *ik, int *it, int *ie, int *is)
# else
int idx_lz_c__ (struct lz_layout_type *lo, int *ig, int *ik, int *it, int *ie, int *is)
# endif
{ return (*pidx_lz)(lo,ig,ik,it,ie,is); }

# ifdef NO_UNDERSCORE
int ik_idx_e_c (struct e_layout_type *lo, int *i)
# elif UNDERSCORE != 2
int ik_idx_e_c_ (struct e_layout_type *lo, int *i)
# else
int ik_idx_e_c__ (struct e_layout_type *lo, int *i)
# endif
{ return (*pik_e)(lo,i); }

# ifdef NO_UNDERSCORE
int it_idx_e_c (struct e_layout_type *lo, int *i)
# elif UNDERSCORE != 2
int it_idx_e_c_ (struct e_layout_type *lo, int *i)
# else
int it_idx_e_c__ (struct e_layout_type *lo, int *i)
# endif
{ return (*pit_e)(lo,i); }

# ifdef NO_UNDERSCORE
int il_idx_e_c (struct e_layout_type *lo, int *i)
# elif UNDERSCORE != 2
int il_idx_e_c_ (struct e_layout_type *lo, int *i)
# else
int il_idx_e_c__ (struct e_layout_type *lo, int *i)
# endif
{ return (*pil_e)(lo,i); }

# ifdef NO_UNDERSCORE
int is_idx_e_c (struct e_layout_type *lo, int *i)
# elif UNDERSCORE != 2
int is_idx_e_c_ (struct e_layout_type *lo, int *i)
# else
int is_idx_e_c__ (struct e_layout_type *lo, int *i)
# endif
{ return (*pis_e)(lo,i); }

# ifdef NO_UNDERSCORE
int idx_e_c (struct e_layout_type *lo, int *ig, int *isign, int *ik, int *it, int *il, int *is)
# elif UNDERSCORE != 2
int idx_e_c_ (struct e_layout_type *lo, int *ig, int *isign, int *ik, int *it, int *il, int *is)
# else
int idx_e_c__ (struct e_layout_type *lo, int *ig, int *isign, int *ik, int *it, int *il, int *is)
# endif
{ return (*pidx_e)(lo,ig,isign,ik,it,il,is); }

# ifdef NO_UNDERSCORE
int ik_idx_le_c (struct le_layout_type *lo, int *i)
# elif UNDERSCORE != 2
int ik_idx_le_c_ (struct le_layout_type *lo, int *i)
# else
int ik_idx_le_c__ (struct le_layout_type *lo, int *i)
# endif
{ return (*pik_le)(lo,i); }

# ifdef NO_UNDERSCORE
int it_idx_le_c (struct le_layout_type *lo, int *i)
# elif UNDERSCORE != 2
int it_idx_le_c_ (struct le_layout_type *lo, int *i)
# else
int it_idx_le_c__ (struct le_layout_type *lo, int *i)
# endif
{ return (*pit_le)(lo,i); }

# ifdef NO_UNDERSCORE
int idx_le_c (struct le_layout_type *lo, int *ig, int *ik, int *it, int *is)
# elif UNDERSCORE != 2
int idx_le_c_ (struct le_layout_type *lo, int *ig, int *ik, int *it, int *is)
# else
int idx_le_c__ (struct le_layout_type *lo, int *ig, int *ik, int *it, int *is)
# endif
{ return (*pidx_le)(lo,ig,ik,it,is); }


/* Internal functions */

/* dist_fn layout functions */
int ik_idx_g_y (struct g_layout_type *lo, int *i)
{ return 1 + (*i - lo->llim_world) % lo->naky; }
int ik_idx_g_lexy (struct g_layout_type *lo, int *i)
{ return 1 + (*i - lo->llim_world)/lo->nlambda/lo->negrid/lo->ntheta0 % lo->naky; }
int ik_idx_g_lxy (struct g_layout_type *lo, int *i)
{ return 1 + (*i - lo->llim_world)/lo->nlambda/lo->ntheta0 % lo->naky; }
int ik_idx_g_ly (struct g_layout_type *lo, int *i)
{ return 1 + (*i - lo->llim_world)/lo->nlambda % lo->naky; }

int it_idx_g_yx (struct g_layout_type *lo, int *i)
{ return 1 + (*i - lo->llim_world)/lo->naky % lo->ntheta0; }
int it_idx_g_lex (struct g_layout_type *lo, int *i)
{ return 1 + (*i - lo->llim_world)/lo->nlambda/lo->negrid % lo->ntheta0; }
int it_idx_g_lx (struct g_layout_type *lo, int *i)
{ return 1 + (*i - lo->llim_world)/lo->nlambda % lo->ntheta0; }
int it_idx_g_lyx (struct g_layout_type *lo, int *i)
{ return 1 + (*i - lo->llim_world)/lo->nlambda/lo->naky % lo->ntheta0; }

int il_idx_g_yxel (struct g_layout_type *lo, int *i)
{ return 1 + (*i - lo->llim_world)/lo->naky/lo->ntheta0/lo->negrid % lo->nlambda; }
int il_idx_g_yxl (struct g_layout_type *lo, int *i)
{ return 1 + (*i - lo->llim_world)/lo->naky/lo->ntheta0 % lo->nlambda; }
int il_idx_g_l (struct g_layout_type *lo, int *i)
{ return 1 + (*i - lo->llim_world) % lo->nlambda; }

int ie_idx_g_yxe (struct g_layout_type *lo, int *i)
{ return 1 + (*i - lo->llim_world)/lo->naky/lo->ntheta0 % lo->negrid; }
int ie_idx_g_yxle (struct g_layout_type *lo, int *i)
{ return 1 + (*i - lo->llim_world)/lo->naky/lo->ntheta0/lo->nlambda % lo->negrid; }
int ie_idx_g_le (struct g_layout_type *lo, int *i)
{ return 1 + (*i - lo->llim_world)/lo->nlambda % lo->negrid; }
int ie_idx_g_lxye (struct g_layout_type *lo, int *i)
{ return 1 + (*i - lo->llim_world)/lo->nlambda/lo->ntheta0/lo->naky % lo->negrid; }
int ie_idx_g_lyxe (struct g_layout_type *lo, int *i)
{ return 1 + (*i - lo->llim_world)/lo->nlambda/lo->naky/lo->ntheta0 % lo->negrid; }

int is_idx_g_yxels (struct g_layout_type *lo, int *i)
{ return 1 + (*i - lo->llim_world)/lo->naky/lo->ntheta0/lo->negrid/lo->nlambda % lo->nspec; }
int is_idx_g_yxles (struct g_layout_type *lo, int *i)
{ return 1 + (*i - lo->llim_world)/lo->naky/lo->ntheta0/lo->nlambda/lo->negrid % lo->nspec; }
int is_idx_g_lexys (struct g_layout_type *lo, int *i)
{ return 1 + (*i - lo->llim_world)/lo->nlambda/lo->negrid/lo->ntheta0/lo->naky % lo->nspec; }
int is_idx_g_lxyes (struct g_layout_type *lo, int *i)
{ return 1 + (*i - lo->llim_world)/lo->nlambda/lo->ntheta0/lo->naky/lo->negrid % lo->nspec; }
int is_idx_g_lyxes (struct g_layout_type *lo, int *i)
{ return 1 + (*i - lo->llim_world)/lo->nlambda/lo->naky/lo->ntheta0/lo->negrid % lo->nspec; }

int idx_g_yxels (struct g_layout_type *lo, int *ik, int *it, int *il, int *ie, int *is)
{ return *ik-1 + lo->naky*(*it-1 + lo->ntheta0*(*ie-1 + lo->negrid*(*il-1 + lo->nlambda*(*is-1)))); }
int idx_g_yxles (struct g_layout_type *lo, int *ik, int *it, int *il, int *ie, int *is)
{ return *ik-1 + lo->naky*(*it-1 + lo->ntheta0*(*il-1 + lo->nlambda*(*ie-1 + lo->negrid*(*is-1)))); }
int idx_g_lexys (struct g_layout_type *lo, int *ik, int *it, int *il, int *ie, int *is)
{ return *il-1 + lo->nlambda*(*ie-1 + lo->negrid*(*it-1 + lo->ntheta0*(*ik-1 + lo->naky*(*is-1)))); }
int idx_g_lxyes (struct g_layout_type *lo, int *ik, int *it, int *il, int *ie, int *is)
{ return *il-1 + lo->nlambda*(*it-1 + lo->ntheta0*(*ik-1 + lo->naky*(*ie-1 + lo->negrid*(*is-1)))); }
int idx_g_lyxes (struct g_layout_type *lo, int *ik, int *it, int *il, int *ie, int *is)
{ return *il-1 + lo->nlambda*(*ik-1 + lo->naky*(*it-1 + lo->ntheta0*(*ie-1 + lo->negrid*(*is-1)))); }

/* Lorentz layout functions */
int ik_idx_lz_y (struct lz_layout_type *lo, int *i)
{ return 1 + (*i - lo->llim_world)/(2*lo->ntgrid + 1) % lo->naky; }
int ik_idx_lz_exy (struct lz_layout_type *lo, int *i)
{ return 1 + (*i - lo->llim_world)/(2*lo->ntgrid + 1)/lo->negrid/lo->ntheta0 % lo->naky; }
int ik_idx_lz_xy (struct lz_layout_type *lo, int *i)
{ return 1 + (*i - lo->llim_world)/(2*lo->ntgrid + 1)/lo->ntheta0 % lo->naky; }

int it_idx_lz_yx (struct lz_layout_type *lo, int *i)
{ return 1 + (*i - lo->llim_world)/(2*lo->ntgrid + 1)/lo->naky % lo->ntheta0; }
int it_idx_lz_ex (struct lz_layout_type *lo, int *i)
{ return 1 + (*i - lo->llim_world)/(2*lo->ntgrid + 1)/lo->negrid % lo->ntheta0; }
int it_idx_lz_x (struct lz_layout_type *lo, int *i)
{ return 1 + (*i - lo->llim_world)/(2*lo->ntgrid + 1) % lo->ntheta0; }

int ie_idx_lz_yxe (struct lz_layout_type *lo, int *i)
{ return 1 + (*i - lo->llim_world)/(2*lo->ntgrid + 1)/lo->naky/lo->ntheta0 % lo->negrid; }
int ie_idx_lz_e (struct lz_layout_type *lo, int *i)
{ return 1 + (*i - lo->llim_world)/(2*lo->ntgrid + 1) % lo->negrid; }
int ie_idx_lz_xye (struct lz_layout_type *lo, int *i)
{ return 1 + (*i - lo->llim_world)/(2*lo->ntgrid + 1)/lo->ntheta0/lo->naky % lo->negrid; }

int is_idx_lz_yxes (struct lz_layout_type *lo, int *i)
{ return 1 + (*i - lo->llim_world)/(2*lo->ntgrid + 1)/lo->naky/lo->ntheta0/lo->negrid % lo->nspec; }
int is_idx_lz_exys (struct lz_layout_type *lo, int *i)
{ return 1 + (*i - lo->llim_world)/(2*lo->ntgrid + 1)/lo->negrid/lo->ntheta0/lo->naky % lo->nspec; }
int is_idx_lz_xyes (struct lz_layout_type *lo, int *i)
{ return 1 + (*i - lo->llim_world)/(2*lo->ntgrid + 1)/lo->ntheta0/lo->naky/lo->negrid % lo->nspec; }

int idx_lz_yxes (struct lz_layout_type *lo, int *ig, int *ik, int *it, int *ie, int *is)
{ return *ig+lo->ntgrid + (2*lo->ntgrid+1)*(*ik-1 + lo->naky*(*it-1 + lo->ntheta0*(*ie-1 + lo->negrid*(*is-1)))); }
int idx_lz_exys (struct lz_layout_type *lo, int *ig, int *ik, int *it, int *ie, int *is)
{ return *ig+lo->ntgrid + (2*lo->ntgrid+1)*(*ie-1 + lo->negrid*(*it-1 + lo->ntheta0*(*ik-1 + lo->naky*(*is-1)))); }
int idx_lz_xyes (struct lz_layout_type *lo, int *ig, int *ik, int *it, int *ie, int *is)
{ return *ig+lo->ntgrid + (2*lo->ntgrid+1)*(*it-1 + lo->ntheta0*(*ik-1 + lo->naky*(*ie-1 + lo->negrid*(*is-1)))); }

/* energy layout functions */
int ik_idx_e_y (struct e_layout_type *lo, int *i)
{ return 1 + (*i - lo->llim_world)/(2*lo->ntgrid + 1)/lo->nsign % lo->naky; }
int ik_idx_e_lxy (struct e_layout_type *lo, int *i)
{ return 1 + (*i - lo->llim_world)/(2*lo->ntgrid + 1)/lo->nsign/lo->nlambda/lo->ntheta0 % lo->naky; }
int ik_idx_e_ly (struct e_layout_type *lo, int *i)
{ return 1 + (*i - lo->llim_world)/(2*lo->ntgrid + 1)/lo->nsign/lo->nlambda % lo->naky; }

int it_idx_e_yx (struct e_layout_type *lo, int *i)
{ return 1 + (*i - lo->llim_world)/(2*lo->ntgrid + 1)/lo->nsign/lo->naky % lo->ntheta0; }
int it_idx_e_lx (struct e_layout_type *lo, int *i)
{ return 1 + (*i - lo->llim_world)/(2*lo->ntgrid + 1)/lo->nsign/lo->nlambda % lo->ntheta0; }
int it_idx_e_lyx (struct e_layout_type *lo, int *i)
{ return 1 + (*i - lo->llim_world)/(2*lo->ntgrid + 1)/lo->nsign/lo->nlambda/lo->naky % lo->ntheta0; }

int il_idx_e_yxl (struct e_layout_type *lo, int *i)
{ return 1 + (*i - lo->llim_world)/(2*lo->ntgrid + 1)/lo->nsign/lo->naky/lo->ntheta0 % lo->nlambda; }
int il_idx_e_l (struct e_layout_type *lo, int *i)
{ return 1 + (*i - lo->llim_world)/(2*lo->ntgrid + 1)/lo->nsign % lo->nlambda; }

int is_idx_e_yxls (struct e_layout_type *lo, int *i)
{ return 1 + (*i - lo->llim_world)/(2*lo->ntgrid + 1)/lo->nsign/lo->naky/lo->ntheta0/lo->nlambda % lo->nspec; }
int is_idx_e_lxys (struct e_layout_type *lo, int *i)
{ return 1 + (*i - lo->llim_world)/(2*lo->ntgrid + 1)/lo->nsign/lo->nlambda/lo->ntheta0/lo->naky % lo->nspec; }
int is_idx_e_lyxs (struct e_layout_type *lo, int *i)
{ return 1 + (*i - lo->llim_world)/(2*lo->ntgrid + 1)/lo->nsign/lo->nlambda/lo->naky/lo->ntheta0 % lo->nspec; }

int idx_e_yxls (struct e_layout_type *lo, int *ig, int *isign, int *ik, int *it, int *il, int *is)
{ return *ig+lo->ntgrid + (2*lo->ntgrid+1)*(*isign-1 + lo->nsign*(*ik-1 + lo->naky*(*it-1 + lo->ntheta0*(*il-1 + lo->nlambda*(*is-1))))); }
int idx_e_lxys (struct e_layout_type *lo, int *ig, int *isign, int *ik, int *it, int *il, int *is)
{ return *ig+lo->ntgrid + (2*lo->ntgrid+1)*(*isign-1 + lo->nsign*(*il-1 + lo->nlambda*(*it-1 + lo->ntheta0*(*ik-1 + lo->naky*(*is-1))))); }
int idx_e_lyxs (struct e_layout_type *lo, int *ig, int *isign, int *ik, int *it, int *il, int *is)
{ return *ig+lo->ntgrid + (2*lo->ntgrid+1)*(*isign-1 + lo->nsign*(*il-1 + lo->nlambda*(*ik-1 + lo->naky*(*it-1 + lo->ntheta0*(*is-1))))); }

/* lorentz-energy layout functions */
int ik_idx_le_y (struct le_layout_type *lo, int *i)
{ return 1 + (*i - lo->llim_world)/(2*lo->ntgrid + 1) % lo->naky; }
int ik_idx_le_xy (struct le_layout_type *lo, int *i)
{ return 1 + (*i - lo->llim_world)/(2*lo->ntgrid + 1)/lo->ntheta0 % lo->naky; }

int it_idx_le_yx (struct le_layout_type *lo, int *i)
{ return 1 + (*i - lo->llim_world)/(2*lo->ntgrid + 1)/lo->naky % lo->ntheta0; }
int it_idx_le_x (struct le_layout_type *lo, int *i)
{ return 1 + (*i - lo->llim_world)/(2*lo->ntgrid + 1) % lo->ntheta0; }

int idx_le_yxs (struct le_layout_type *lo, int *ig, int *ik, int *it, int *is)
{ return *ig+lo->ntgrid + (2*lo->ntgrid+1)*(*ik-1 + lo->naky*(*it-1 + lo->ntheta0*(*is-1))); }
int idx_le_xys (struct le_layout_type *lo, int *ig, int *ik, int *it, int *is)
{ return *ig+lo->ntgrid + (2*lo->ntgrid+1)*(*it-1 + lo->ntheta0*(*ik-1 + lo->naky*(*is-1))); }
