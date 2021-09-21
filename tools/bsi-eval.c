/* This file is part of BSI, a library for Biot-Savart Integral evaluation
 *
 * Copyright (C) 2021 Michael Carley
 *
 * BSI is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version. BSI is distributed in the
 * hope that it will be useful, but WITHOUT ANY WARRANTY; without even
 * the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 * PURPOSE.  See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with BSI.  If not, see <https://www.gnu.org/licenses/>.
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /*HAVE_CONFIG_H*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <glib.h>

#include <bsi.h>

#include <wbfmm.h>

#include <gqr.h>

#include "quadrature.h"
#include "tetquad.h"

#include "duffy.h"

#define REAL gdouble
REAL orient3d(REAL *pa, REAL *pb, REAL *pc, REAL *pd);

#define SIGN(_x) ((_x) < 0.0 ? -1 : ((_x) > 0.0 ? 1 : 0))

#include "bsi-private.h"

/* #define BSI_UNCORRECTED */

gchar *progname ;

GTimer *timer ;

/* #define BSI_INTERPOLATE_RBF */

#define BSI_DATA_SIZE            24
#define BSI_DATA_VERTICES         0
#define BSI_DATA_QUAD_NODES       1
#define BSI_DATA_QUAD_RULE        2
#define BSI_DATA_QUAD_NUMBER      3
#define BSI_DATA_VECTOR_SIZE      4
#define BSI_DATA_TET_COUNT        5
#define BSI_DATA_NODE_STRIDE      6
#define BSI_DATA_NODE_DATA        7
#define BSI_DATA_VELOCITY         8
#define BSI_DATA_VELOCITY_STRIDE  9
#define BSI_DATA_INTERPOLATOR    10
#define BSI_DATA_WORKSPACE       11
#define BSI_DATA_INTERP_ORDER    12
#define BSI_DATA_VERTEX_DATA     13
#define BSI_DATA_VOLUME          14
#define BSI_DATA_QUADRATURE_PHI  15
#define BSI_DATA_QUADRATURE_TH   16
#define BSI_DATA_QUADRATURE_R    17

gint bsi_eval_velocity(gdouble *x, gint xstr, gdouble *f, gint fstr,
		       gint nx, gdouble *xe, gdouble *v) ;

#define bsi_vector_cross(BSI_C,BSI_A,BSI_B)				\
  ((BSI_C)[0] = (BSI_A)[1]*(BSI_B)[2] - (BSI_A)[2]*(BSI_B)[1],		\
   (BSI_C)[1] = (BSI_A)[2]*(BSI_B)[0] - (BSI_A)[0]*(BSI_B)[2],		\
   (BSI_C)[2] = (BSI_A)[0]*(BSI_B)[1] - (BSI_A)[1]*(BSI_B)[0])

static void invert3x3(gdouble *Ai, gdouble *A)

{
  gdouble det ;

  det = 
    A[0]*(A[8]*A[4]-A[7]*A[5]) - 
    A[3]*(A[8]*A[1]-A[7]*A[2]) +
    A[6]*(A[5]*A[1]-A[4]*A[2]) ;
  det = 1.0/det ;

  Ai[0] =  det*(A[8]*A[4] - A[7]*A[5]) ;
  Ai[1] = -det*(A[8]*A[1] - A[7]*A[2]) ;
  Ai[2] =  det*(A[5]*A[1] - A[4]*A[2]) ;

  Ai[3] = -det*(A[8]*A[3] - A[6]*A[5]) ;
  Ai[4] =  det*(A[8]*A[0] - A[6]*A[2]) ;
  Ai[5] = -det*(A[5]*A[0] - A[3]*A[2]) ;

  Ai[6] =  det*(A[7]*A[3] - A[6]*A[4]) ;
  Ai[7] = -det*(A[7]*A[0] - A[6]*A[1]) ;
  Ai[8] =  det*(A[4]*A[0] - A[3]*A[1]) ;

  return ;
}

static gint read_nodes(FILE *f, gdouble **nodes, gdouble **data,
		       gint *nnodes, gint *nd)

{
  gint i, j ;

  fscanf(f, "%d", nnodes) ;
  fscanf(f, "%d", nd) ;

  *nodes = (gdouble *)g_malloc0((*nnodes)*3*sizeof(gdouble)) ;
  *data  = (gdouble *)g_malloc0((*nnodes)*(*nd)*sizeof(gdouble)) ;

  for ( i = 0 ; i < *nnodes ; i ++ ) {
    fscanf(f, "%lg", &((*nodes)[3*i+0])) ;
    fscanf(f, "%lg", &((*nodes)[3*i+1])) ;
    fscanf(f, "%lg", &((*nodes)[3*i+2])) ;
    for ( j = 0 ; j < (*nd) ; j ++ )
      fscanf(f, "%lg", &((*data)[i*(*nd)+j])) ;
  }
  
  return 0 ;
}

static gint tetrahedron_bs_func(gdouble *y, gdouble *s, gdouble R,
				gdouble wt, gdouble *u, gint nu,
				gpointer data[])

{
  /* bsi_polynomial_t *p = data[0] ; */
  gdouble rw[3], w[3] ;

  wt *= 0.25*M_1_PI/(R*R) ;
  
  /*vector s is the unit vector from the evaluation point to the
    quadrature point (Rs = y-x)*/
  
  /* bsi_polynomial_evaluate(p, y, w) ;   */

  bsi_source_func_ring_gaussian(y, w, 3, NULL) ;
  tq_vector_cross(rw, w, s) ;

  /* u[0] -= rw[0]/R/R*0.25*M_1_PI*wt ; */
  /* u[1] -= rw[1]/R/R*0.25*M_1_PI*wt ; */
  /* u[2] -= rw[2]/R/R*0.25*M_1_PI*wt ; */
  u[0] -= rw[0]*wt ;
  u[1] -= rw[1]*wt ;
  u[2] -= rw[2]*wt ;

  return 0 ;
}

static gint tetrahedron_duffy_func(gdouble *y, gdouble *x, gdouble wt,
				   gdouble *u, gint nu, gpointer data)

{
  gdouble s[3], w[3], sw[3], r ;

  bsi_source_func_ring_gaussian(y, w, 3, NULL) ;

  duffy_vector3(s, x, y) ;
  r = sqrt(duffy_length3sqr(s)) ;
  s[0] /= r ; s[1] /= r ; s[2] /= r ; 
  duffy_cross3(sw, w, s) ;

  u[0] -= sw[0]*0.25*M_1_PI*wt ;
  u[1] -= sw[1]*0.25*M_1_PI*wt ;
  u[2] -= sw[2]*0.25*M_1_PI*wt ;
  
  return 0 ;
}

static void quadrature_nodes_tetgen(gint tet[4], gpointer *data)

{
  gdouble *x = data[BSI_DATA_VERTICES] ;
  gdouble *xq =  data[BSI_DATA_QUAD_NODES] ;
  gdouble *q =   data[BSI_DATA_QUAD_RULE] ;
  gint nq =     *((gint *)data[BSI_DATA_QUAD_NUMBER]) ;
  gint nf =     *((gint *)data[BSI_DATA_VECTOR_SIZE]) ;
  /* gint iorder = *((gint *)data[BSI_DATA_INTERP_ORDER]) ; */
  gint fstr =   *((gint *)data[BSI_DATA_NODE_STRIDE]) ;
  gint *ntf =    data[BSI_DATA_TET_COUNT] ;
  gdouble *f =   data[BSI_DATA_NODE_DATA] ;
  gdouble *u =   data[BSI_DATA_VELOCITY] ;
  /* gdouble *fn =  data[BSI_DATA_VERTEX_DATA] ; */
#ifdef BSI_INTERPOLATE_RBF
  bsi_interpolator_t *rbf = data[BSI_DATA_INTERPOLATOR] ;
  gdouble *work = data[BSI_DATA_WORKSPACE] ;
#else /*BSI_INTERPOLATE_RBF*/
  bsi_polynomial_t *pb = data[BSI_DATA_INTERPOLATOR] ;
  /* bsi_polynomial_workspace_t *wpb = data[BSI_DATA_WORKSPACE]    ; */
#endif /*BSI_INTERPOLATE_RBF*/
  gqr_rule_t *qp = data[BSI_DATA_QUADRATURE_PHI] ;
  gqr_rule_t *qt = data[BSI_DATA_QUADRATURE_TH] ;
  gqr_rule_t *qr = data[BSI_DATA_QUADRATURE_R] ;
  gdouble L[4], J, r[3], *xquad, *fquad ;
  gint i, j, str, rot[] = {0, 1, 2, 3, 0, 1, 2} ;

  str = 3+nf ;

  J = fabs(orient3d(&(x[3*tet[0]]), &(x[3*tet[1]]), &(x[3*tet[2]]),
		    &(x[3*tet[3]])))/6.0 ;
  
  if ( J < 1e-12 ) return ;
  
#ifdef BSI_INTERPOLATE_RBF
  /* bsi_rbf_parameter_fit_cell(c, x, fn, nf, nf, rbf, work) ; */
  bsi_rbf_parameter_fit(x, idx[0], v, fn, nf, nf, rbf, work) ;
#else /*BSI_INTERPOLATE_RBF*/
  /* bsi_polynomial_cell_interpolant(c, x, fn, nf, nf, 5, pb, */
  /* 				  iorder, wpb) ; */
#endif /*BSI_INTERPOLATE_RBF*/  
  
  for ( i = 0 ; i < nq ; i ++ ) {
    L[1] = q[4*i+0] ; L[2] = q[4*i+1] ; L[3] = q[4*i+2] ;
    L[0] = 1.0 - L[1] - L[2] - L[3] ;
    xquad = &(xq[(*ntf)*str]) ;
    fquad = &(f[(*ntf)*fstr]) ;
    xquad[0] = L[0]*x[3*tet[0]+0] + L[1]*x[3*tet[1]+0] +
      L[2]*x[3*tet[2]+0] + L[3]*x[3*tet[3]+0] ; 
    xquad[1] = L[0]*x[3*tet[0]+1] + L[1]*x[3*tet[1]+1] +
      L[2]*x[3*tet[2]+1] + L[3]*x[3*tet[3]+1] ; 
    xquad[2] = L[0]*x[3*tet[0]+2] + L[1]*x[3*tet[1]+2] +
      L[2]*x[3*tet[2]+2] + L[3]*x[3*tet[3]+2] ; 

    fquad[0] = fquad[1] = fquad[2] = 0.0 ;

#ifdef BSI_INTERPOLATE_RBF
    bsi_rbf_evaluate(rbf, x, xquad, nf, fquad) ;
#else /*BSI_INTERPOLATE_RBF*/
    /* bsi_polynomial_evaluate(pb, xquad, fquad) ; */
#endif /*BSI_INTERPOLATE_RBF*/

    bsi_source_func_ring_gaussian(xquad, fquad, nf, NULL) ;
    
    fquad[0] *= J*q[4*i+3] ;
    fquad[1] *= J*q[4*i+3] ;
    fquad[2] *= J*q[4*i+3] ;

#ifndef BSI_UNCORRECTED
    for ( j = 0 ; j < 4 ; j ++ ) {
      gdouble rw[3], R ;

      tq_vector_init(r, xquad, &(x[3*tet[j]])) ;
      
      R = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]) ;
      R *= R*R*4.0*M_PI ;
    
      bsi_vector_cross(rw, r, fquad) ;
      u[tet[j]*3+0] += rw[0]/R ;
      u[tet[j]*3+1] += rw[1]/R ;
      u[tet[j]*3+2] += rw[2]/R ;
    }
#endif /*BSI_UNCORRECTED*/
    (*ntf) ++ ;
  }
#ifdef BSI_UNCORRECTED
  return ;
#endif /*BSI_UNCORRECTED*/

  /*local volume interpolation and integration*/
  for ( i = 0 ; i < 4 ; i ++ ) {
    gpointer fdata[2] ;
    fdata[0] = pb ;
    /* duffy_tet_quad(&(x[3*tet[rot[i+0]]]), &(x[3*tet[rot[i+1]]]), */
    /* 		   &(x[3*tet[rot[i+2]]]), &(x[3*tet[rot[i+3]]]), */
    /* 		   2.0, 2.0, */
    /* 		   &(gqr_rule_abscissa(qp,0)), 1, &(gqr_rule_weight(qp,0)), 1, */
    /* 		   gqr_rule_length(qp), */
    /* 		   &(gqr_rule_abscissa(qt,0)), 1, &(gqr_rule_weight(qt,0)), 1, */
    /* 		   gqr_rule_length(qt), */
    /* 		   &(gqr_rule_abscissa(qr,0)), 1, &(gqr_rule_weight(qr,0)), 1, */
    /* 		   gqr_rule_length(qr), */
    /* 		   (duffy_func_t)tetrahedron_duffy_func, NULL, */
    /* 		   &(u[tet[i]*3]), 3) ; */
    tq_tet_quad(&(x[3*tet[rot[i+0]]]), &(x[3*tet[rot[i+1]]]),
    		&(x[3*tet[rot[i+2]]]), &(x[3*tet[rot[i+3]]]),
    		&(gqr_rule_abscissa(qp,0)), 1, &(gqr_rule_weight(qp,0)), 1,
    		gqr_rule_length(qp),
    		&(gqr_rule_abscissa(qt,0)), 1, &(gqr_rule_weight(qt,0)), 1,
    		gqr_rule_length(qt),
    		&(gqr_rule_abscissa(qr,0)), 1, &(gqr_rule_weight(qr,0)), 1,
    		gqr_rule_length(qr),
    		(tq_tetquad_func_t)tetrahedron_bs_func, fdata,
    		&(u[tet[i]*3]), 3) ;
  }
  
  return ;
}

gint bsi_eval_velocity(gdouble *x, gint xstr, gdouble *f, gint fstr,
		       gint nx, gdouble *xe, gdouble *v)

{
  gint i ;
  gdouble R, r[3], rw[3] ;
  
  for ( i = 0 ; i < nx ; i ++ ) {
    r[0] = xe[0] - x[i*xstr+0] ;
    r[1] = xe[1] - x[i*xstr+1] ;
    r[2] = xe[2] - x[i*xstr+2] ;
    R = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]) ;
    R *= R*R*4.0*M_PI ;
    
    bsi_vector_cross(rw, r, &(f[i*fstr])) ;

    v[0] -= rw[0]/R ;
    v[1] -= rw[1]/R ;
    v[2] -= rw[2]/R ;
  }
  
  return 0 ;
}

static gint points_origin_width(gdouble *x, gint nx, gdouble *xmin,
				gdouble *xmax, gdouble *D, gboolean init)
{
  gint i ;

  if ( init ) {
    xmin[0] = xmin[1] = xmin[2] =  G_MAXDOUBLE ;
    xmax[0] = xmax[1] = xmax[2] = -G_MAXDOUBLE ;
  }

  for ( i = 0 ; i < nx ; i ++ ) {
    xmin[0] = MIN(xmin[0], x[3*i+0]) ;
    xmin[1] = MIN(xmin[1], x[3*i+1]) ;
    xmin[2] = MIN(xmin[2], x[3*i+2]) ;
    xmax[0] = MAX(xmax[0], x[3*i+0]) ;
    xmax[1] = MAX(xmax[1], x[3*i+1]) ;
    xmax[2] = MAX(xmax[2], x[3*i+2]) ;
  }    

  *D = xmax[0] - xmin[0] ;
  *D = MAX(*D, xmax[1] - xmin[1]) ;
  *D = MAX(*D, xmax[2] - xmin[2]) ;

  return 0 ;
}

static gboolean point_in_tet(gdouble *x, gint tet[4], gdouble *p)

{
  gint i0, i1, i2, i3, i, s ;
  gdouble min, max ;

  i0 = tet[0] ; i1 = tet[1] ; i2 = tet[2] ; i3 = tet[3] ; 

  /*bounding box test here to eliminate most predicate calls*/
  for ( i = 0 ; i < 3 ; i ++ ) {
    min = MIN(x[3*i0+i], MIN(x[3*i1+i], MIN(x[3*i2+i], x[3*i3+i]))) ;
    max = MAX(x[3*i0+i], MAX(x[3*i1+i], MAX(x[3*i2+i], x[3*i3+i]))) ;
    if ( p[i] < min ) return FALSE ;
    if ( p[i] > max ) return FALSE ;
  }
  
  /*add a circumsphere test here to speed things up a bit?*/
  if ( (s = SIGN(orient3d(&(x[3*i0]), &(x[3*i1]), &(x[3*i2]), p))) == 0 )
    return TRUE ;  
  if ( s != SIGN(orient3d(&(x[3*i0]), &(x[3*i1]), &(x[3*i2]), &(x[3*i3]))) )
    return FALSE ;

  if ( (s = SIGN(orient3d(&(x[3*i1]), &(x[3*i2]), &(x[3*i3]), p))) == 0 )
    return TRUE ;  
  if ( s != SIGN(orient3d(&(x[3*i1]), &(x[3*i2]), &(x[3*i3]), &(x[3*i0]))) )
    return FALSE ;

  if ( (s = SIGN(orient3d(&(x[3*i2]), &(x[3*i3]), &(x[3*i0]), p))) == 0 )
    return TRUE ;  
  if ( s != SIGN(orient3d(&(x[3*i2]), &(x[3*i3]), &(x[3*i0]), &(x[3*i1]))) )
    return FALSE ;

  if ( (s = SIGN(orient3d(&(x[3*i3]), &(x[3*i0]), &(x[3*i1]), p))) == 0 )
    return TRUE ;  
  if ( s != SIGN(orient3d(&(x[3*i3]), &(x[3*i0]), &(x[3*i1]), &(x[3*i2]))) )
    return FALSE ;
  
  return TRUE ;
}

static gint tet_barycentric(gdouble *x1, gdouble *x2, gdouble *x3, gdouble *x4,
			    gdouble *p, gdouble *L)

{
  gdouble A[9], Ai[9], r[3] ;

  A[0] = x1[0] - x4[0] ; 
  A[1] = x2[0] - x4[0] ; 
  A[2] = x3[0] - x4[0] ; 
  A[3] = x1[1] - x4[1] ;  
  A[4] = x2[1] - x4[1] ;  
  A[5] = x3[1] - x4[1] ;  
  A[6] = x1[2] - x4[2] ;  
  A[7] = x2[2] - x4[2] ;  
  A[8] = x3[2] - x4[2] ;  

  invert3x3(Ai, A) ;
  r[0] = p[0] - x4[0] ;
  r[1] = p[1] - x4[1] ;
  r[2] = p[2] - x4[2] ;

  L[0] = Ai[0]*r[0] + Ai[1]*r[1] + Ai[2]*r[2] ;
  L[1] = Ai[3]*r[0] + Ai[4]*r[1] + Ai[5]*r[2] ;
  L[2] = Ai[6]*r[0] + Ai[7]*r[1] + Ai[8]*r[2] ;

  L[3] = 1.0 - L[0] - L[1] - L[2] ;
  
  return 0 ;
}

static gint mesh_interp_linear(gdouble *x, gint *tet, gint ntet, gdouble *f,
			       gint nf, gint fstr, gdouble *xi, gdouble *fi)

{
  gint i, j ;
  gdouble L[4] ;
  
  memset(fi, 0, nf*sizeof(gdouble)) ;
  for ( i = 0 ; i < ntet ; i ++ ) {
    if ( point_in_tet(x, &(tet[4*i]), xi) ) break ;
  }
  if ( i == ntet ) return 0 ;

  tet_barycentric(&(x[3*tet[4*i+0]]), &(x[3*tet[4*i+1]]),
		  &(x[3*tet[4*i+2]]), &(x[3*tet[4*i+3]]),
		  xi, L) ;

  for ( j = 0 ; j < nf ; j ++ ) {
    fi[j] =
      L[0]*f[fstr*tet[4*i+0]+j] + 
      L[1]*f[fstr*tet[4*i+1]+j] + 
      L[2]*f[fstr*tet[4*i+2]+j] + 
      L[3]*f[fstr*tet[4*i+3]+j] ;
  }
  
  return 0 ;
}

gint main(gint argc, gchar **argv)

{
  FILE *input, *output ;
  gdouble *f, *q, *xq, *u, xg[3], ui[3] ;
  gdouble *x ;
  gint nx, nf, order, nq, ntet, ntf, str, i, j, ngx, ngz ;
  gint nnmax, iorder ;
  gint ngp, ngt, ngr ;
#ifdef BSI_INTERPOLATE_RBF
  bsi_interpolator_t *rbf ;
  gdouble work[4096] ;
#else /*BSI_INTERPOLATE_RBF*/
  bsi_polynomial_t *pb ;
  bsi_polynomial_workspace_t *wpb ;
#endif /*BSI_INTERPOLATE_RBF*/
  gpointer data[BSI_DATA_SIZE] ;
  gqr_rule_t *qp, *qt, *qr ;
  wbfmm_tree_t *tree ;
  wbfmm_target_list_t *targets ;
  wbfmm_shift_operators_t *shifts ;
  gdouble xtree[3], xtmax[3], D, del, *fmmwork ;
  gint pstr, torder[32], order_r, order_s, order_max, depth ;
  gint sizew, level, nthreads ;
  tetwrap_tetgenio_t *in, *out ;
  gchar *outfile = "velocity.dat" ;
  gboolean sort_sources ;
  
  progname = g_strdup(g_path_get_basename(argv[0])) ;
  timer = g_timer_new() ;

  nthreads = 2 ;
  
  input = stdin ;
  output = stdout ;

  sort_sources = TRUE ;
  /*interpolation parameters*/
  nnmax = 10 ; iorder = 2 ;
  /* mop_logging_init(stderr, "", G_LOG_LEVEL_DEBUG, NULL) ;   */

  /*global quadrature*/
  order = 3 ;

  /*local quadratures*/
  /* ngp = 16 ; ngt = 16 ; ngr = 16 ; */
  /* ngp = 12 ; ngt = 12 ; ngr = 12 ; */
  /* ngp = 8 ; ngt = 8 ; ngr = 8 ; */
  ngp = 4 ; ngt = 4 ; ngr = 4 ;
  /* ngp = 2 ; ngt = 2 ; ngr = 2 ; */
  /* ngp = 1 ; ngt = 2 ; ngr = 2 ; */
  
  /*tree parameters*/
  del = 1e-2 ;
  order_r = 4 ; order_s = 4 ;
  depth = 6 ;
  tree = NULL ;
  
  fprintf(stderr, "%s: reading nodes; %lg\n",
	  progname, g_timer_elapsed(timer, NULL)) ;
  /* read_nodes(input, &x, &f, &nx, &nf) ; */
  read_nodes(input, &x, &f, &nx, &nf) ;

  fprintf(stderr, "%s: tetrahedralizing %d nodes; %lg\n",
	  progname, nx, g_timer_elapsed(timer, NULL)) ;
  /* v = bsi_make_delaunay_volume(x, nx) ; */
  in = tetwrap_new(x, nx) ;
  out = tetwrap_new(NULL, 0) ;
  tetwrap_tetrahedralize("", in, out, NULL, NULL) ;
  
  bsi_quadrature_js_select(order, &q, &nq) ;
  fprintf(stderr, "%s: %dth order quadrature, %d nodes; %lg\n",
	  progname, order, nq, g_timer_elapsed(timer, NULL)) ;

  ntet = out->numberoftetrahedra ;
  
  xq = (gdouble *)g_malloc0(ntet*(3+nf)*nq*sizeof(gdouble)) ;
  u  = (gdouble *)g_malloc0(3*nx*sizeof(gdouble)) ;

#ifdef BSI_INTERPOLATE_RBF
  rbf = bsi_interpolator_new(BSI_RBF_CP_C6, iorder, 3, nnmax) ;
#else /*BSI_INTERPOLATE_RBF*/
  pb  = bsi_polynomial_new(iorder, nf) ;
  wpb = bsi_polynomial_workspace_new(iorder, nf, nnmax) ;  
#endif /*BSI_INTERPOLATE_RBF*/
  
  str = 3+nf ;

  ntf = ntet*nq ;

  /*initialize the FMM tree*/
  
#if 1
  fprintf(stderr, "%s: building FMM tree; %lg\n",
	  progname, g_timer_elapsed(timer, NULL)) ;
  wbfmm_points_origin_width(xq, str, ntf, xtree, xtmax, &D, TRUE) ;
  /* bsi_points_origin_width(x, nx, xtree, xtmax, &D, FALSE) ; */
  points_origin_width(x, nx, xtree, xtmax, &D, FALSE) ;
  xtree[0] -= del ; xtree[1] -= del ; xtree[2] -= del ;
  D += 2.0*del ;

  pstr = str*sizeof(gdouble) ;
  tree = wbfmm_tree_new(xtree, D, ntf) ;
  torder[2*depth+0] = order_s ; 
  torder[2*depth+1] = order_r ; 
  order_max = MAX(order_s, order_r) ;
  for ( i = depth-1 ; i > 0 ; i -- ) {
    torder[2*i+0] = torder[2*(i+1)+0] + 2 ;
    torder[2*i+1] = torder[2*(i+1)+1] + 2 ;
    order_max = MAX(order_max, torder[2*i+0]) ;
    order_max = MAX(order_max, torder[2*i+1]) ;
  }
  
#endif
  
  sizew = wbfmm_element_number_rotation(2*order_max) ;
  sizew = MAX(sizew, (order_max+1)*(order_max+1)*nf*16) ;
  fmmwork = (gdouble *)g_malloc0(2*MAX(1,nthreads)*sizew*sizeof(gdouble)) ;

  ntf = 0 ;

  qp = gqr_rule_alloc(ngp) ;
  qr = gqr_rule_alloc(ngr) ;
  qt = gqr_rule_alloc(ngt) ;
  
  gqr_rule_select(qp, GQR_GAUSS_LEGENDRE, ngp, NULL) ;
  gqr_rule_select(qr, GQR_GAUSS_LEGENDRE, ngr, NULL) ;
  gqr_rule_select(qt, GQR_GAUSS_LEGENDRE, ngt, NULL) ;

  data[BSI_DATA_VERTICES]     = x ;
  data[BSI_DATA_QUAD_NODES]   = xq ;
  data[BSI_DATA_QUAD_RULE]    = q ;
  data[BSI_DATA_QUAD_NUMBER]  = &nq ;
  data[BSI_DATA_VECTOR_SIZE]  = &nf ;
  data[BSI_DATA_INTERP_ORDER] = &iorder ;
  data[BSI_DATA_NODE_STRIDE]  = &str ;
  data[BSI_DATA_NODE_DATA]    = &(xq[3]) ;
  data[BSI_DATA_VELOCITY]     = u ;
#ifdef BSI_INTERPOLATE_RBF
  data[BSI_DATA_INTERPOLATOR] = rbf ;
  data[BSI_DATA_WORKSPACE]    = fmmwork ;
  data[BSI_DATA_VOLUME]       = v ;
#else /*BSI_INTERPOLATE_RBF*/
  data[BSI_DATA_INTERPOLATOR] = pb ;
  data[BSI_DATA_WORKSPACE]    = wpb ;
#endif /*BSI_INTERPOLATE_RBF*/
  data[BSI_DATA_VERTEX_DATA]  = f ;
  data[BSI_DATA_TET_COUNT]    = &ntf ;
  data[BSI_DATA_QUADRATURE_PHI] = qp ;
  data[BSI_DATA_QUADRATURE_TH]  = qt ; 
  data[BSI_DATA_QUADRATURE_R]   = qr ;

  fprintf(stderr, "%s: generating quadrature nodes; %lg\n",
	  progname, g_timer_elapsed(timer, NULL)) ;
  for ( i = 0 ; i < ntet ; i ++ ) {
    quadrature_nodes_tetgen(&(out->tetrahedronlist[4*i]), data) ;
  }
  
  fprintf(stderr, "ntf == %d (%d)\n", ntf, ntet*nq) ;
  
#if 1
  fprintf(stderr, "%s: initializing shift rotation operators; %lg\n",
  	  progname, g_timer_elapsed(timer, NULL)) ;
  wbfmm_shift_angle_table_init() ;
  shifts = wbfmm_shift_operators_new(order_max, TRUE, fmmwork) ;
  fprintf(stderr, "%s: shift rotation operators initialized; %lg\n",
  	  progname, g_timer_elapsed(timer, NULL)) ;

  fprintf(stderr, "%s: initializing coaxial translation coefficients; %lg\n",
  	  progname, g_timer_elapsed(timer, NULL)) ;
  wbfmm_laplace_coaxial_translate_init(order_max+2) ;
  fprintf(stderr, "%s: coaxial translation coefficients initialized; %lg\n",
  	  progname, g_timer_elapsed(timer, NULL)) ;

  /*add source points to the tree and refine to allocate sources to
    leaf boxes*/
  if ( sort_sources )
    wbfmm_tree_sort_points(tree, xq, pstr, ntf) ;
  wbfmm_tree_add_points(tree, (gpointer)xq, pstr, NULL, 0, ntf, sort_sources) ;
  for ( i = 0 ; i < depth ; i ++ ) wbfmm_tree_refine(tree) ;

  wbfmm_tree_problem(tree) = WBFMM_PROBLEM_LAPLACE ;
  wbfmm_tree_source_size(tree) = nf ;
  for ( i = 1 ; i <= depth ; i ++ ) {
    wbfmm_tree_laplace_coefficient_init(tree, i,
					torder[2*i+1], torder[2*i+0]) ;
  }
  fprintf(stderr, "%s: building target list; %lg\n",
	  progname, g_timer_elapsed(timer, NULL)) ;
  targets = wbfmm_target_list_new(tree, nx) ;
  wbfmm_target_list_add_points(targets, x, 3*sizeof(gdouble), nx) ;

  fprintf(stderr, "%s: initializing leaf expansions; %lg\n",
	  progname, g_timer_elapsed(timer, NULL)) ;  
  wbfmm_tree_laplace_leaf_expansions(tree,
				     &(xq[3]), str, NULL, 0, NULL, 0,
				     TRUE, fmmwork) ;  
  fprintf(stderr, "%s: leaf expansions initialized; %lg\n",
	  progname, g_timer_elapsed(timer, NULL)) ;

  fprintf(stderr, "%s: upward pass; %lg\n",
	  progname, g_timer_elapsed(timer, NULL)) ;
  for ( level = depth ; level >= 3 ; level -- ) {
    wbfmm_laplace_upward_pass(tree, shifts, level, fmmwork) ;
  }  
  fprintf(stderr, "%s: upward pass completed; %lg\n",
	  progname, g_timer_elapsed(timer, NULL)) ;

  fprintf(stderr, "%s: downward pass; %lg\n",
	  progname, g_timer_elapsed(timer, NULL)) ;
  for ( level = 2 ; level <= depth ; level ++ ) {
    wbfmm_laplace_downward_pass(tree, shifts, level, fmmwork, nthreads) ;
  }
  fprintf(stderr, "%s: downward pass completed; %lg\n",
	  progname, g_timer_elapsed(timer, NULL)) ;
#endif
  
  if ( tree == NULL ) {
    fprintf(stderr, "%s: starting point source summation; %lg\n",
	    progname, g_timer_elapsed(timer, NULL)) ;
    for ( i = 0 ; i < nx ; i ++ ) {
      bsi_eval_velocity(xq, str, &(xq[3]), str, ntf,
			&(x[3*i]), &(u[3*i])) ;
    }
    fprintf(stderr, "%s: point source summation completed; %lg\n",
	    progname, g_timer_elapsed(timer, NULL)) ;
  } else {
    fprintf(stderr, "%s: starting FMM point source summation; %lg\n",
	    progname, g_timer_elapsed(timer, NULL)) ;

    for ( i = 0 ; i < nx ; i ++ ) {
      guint64 box ;
      /* box = wbfmm_point_box(tree, tree->depth, &(x[3*i])) ; */
      box = wbfmm_target_point_box(targets,i) ;
      wbfmm_tree_laplace_box_local_curl(tree, tree->depth, box,
					&(x[3*i]),
					&(u[3*i]), 3,
					&(xq[3]), str,
					NULL, 0, NULL, 0,
					TRUE, fmmwork) ;
    }
    
    fprintf(stderr, "%s: FMM point source summation completed; %lg\n",
	    progname, g_timer_elapsed(timer, NULL)) ;
  }
  
  ngx = 0 ; ngz = 1025 ;

  if ( outfile != NULL )
    output = fopen(outfile, "w") ;
  else
    output = stdout ;
  
  for ( i = 0 ; i <= ngx ; i ++ ) {
    xg[0] = 1.0 ; /* -2 + (4.0)*i/ngx ; */
    xg[1] = 0.0 ;
    for ( j = 0 ; j <= ngz ; j ++ ) {
      xg[2] = -1 + (2.0)*j/ngz ;
      ui[0] = ui[1] = ui[2] = 0.0 ;

      mesh_interp_linear(x, out->tetrahedronlist, ntet, u, 3, 3, xg, ui) ;
      
      fprintf(output,
  	      "%e %e %e %e %e %e\n",
  	      xg[0], xg[1], xg[2],
  	      ui[0], ui[1], ui[2]) ;
    }
  }
  
  if ( output != stdout ) fclose(output) ;
  
  return 0 ;
}
