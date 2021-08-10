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

static gint read_nodes(FILE *f, GtsVertex **nodes, gdouble **data,
		       gint *nnodes, gint *nd)

{
  gint i, j ;
  gdouble x, y, z ;
  GtsVertex *v ;

  fscanf(f, "%d", nnodes) ;
  fscanf(f, "%d", nd) ;

  *nodes = (GtsVertex *)g_malloc0((*nnodes)*sizeof(GtsVertex)) ;
  *data  = (gdouble   *)g_malloc0((*nnodes)*(*nd)*sizeof(gdouble)) ;

  for ( i = 0 ; i < *nnodes ; i ++ ) {
    fscanf(f, "%lg %lg %lg", &x, &y, &z) ;
    v = gts_vertex_new(gts_vertex_class(), x, y, z) ;
    memcpy(&((*nodes)[i]), v, sizeof(GtsVertex)) ;
    gts_object_destroy(GTS_OBJECT(v)) ;
    for ( j = 0 ; j < (*nd) ; j ++ )
      fscanf(f, "%lg", &((*data)[i*(*nd)+j])) ;
  }
  
  return 0 ;
}

gboolean cell_match(GtsVertex *v[])

{
  gdouble xc[] = {1, 0, 0,
		  1.01894, 0.0170455, -0.0104167,
		  1.00568, 0.0246212, -0.00347222,
		  0.996212, 0.0132576, -0.0173611} ;
  gdouble tol = 1e-3 ;
  gint i, j ;
  gboolean match ;
  
  for ( i = 0 ; i < 4 ; i ++ ) {
    match = FALSE ;
    for ( j = 0 ; j < 4 ; j ++ ) {
      if ( fabs(GTS_POINT(v[i])->x - xc[3*j+0]) < tol &&
	   fabs(GTS_POINT(v[i])->y - xc[3*j+1]) < tol &&
	   fabs(GTS_POINT(v[i])->z - xc[3*j+2]) < tol ) {
	match = TRUE ;
      }	
    }

    if ( match == FALSE ) return FALSE ;
  }

  return TRUE ;
}

gint tetrahedron_bs_func(gdouble *y, gdouble *s, gdouble R,
			 gdouble wt, gdouble *u, gint nu,
			 gpointer data[])
/* bsi_polynomial_t *p) */

{
  bsi_polynomial_t *p = data[0] ;
  gdouble rw[3], w[3] ;

  /* g_assert(wt > 0.0) ; */


  /*vector s is the unit vector from the evaluation point to the
    quadrature point (Rs = y-x)*/
  
  bsi_polynomial_evaluate(p, y, w) ;  

#if 1
  tq_vector_cross(rw, w, s) ;
  /* tq_vector_cross(rw, s, w) ; */

  u[0] -= rw[0]/R/R*0.25*M_1_PI*wt ;
  u[1] -= rw[1]/R/R*0.25*M_1_PI*wt ;
  u[2] -= rw[2]/R/R*0.25*M_1_PI*wt ;

  return 0 ;
#else
  gdouble *x = data[1] ;
  gdouble r[3], R3 ;
  tq_vector_init(r, y, x) ;
  bsi_vector_cross(rw, r, w) ;
  
  R3 = tq_vector_length2(r) ;
  R3 *= sqrt(R3)*4.0*M_PI ;
  
  u[0] -= rw[0]/R3*wt ;
  u[1] -= rw[1]/R3*wt ;
  u[2] -= rw[2]/R3*wt ;
#endif
  
  return 0 ;
}

static void quadrature_nodes(GtvCell *c, gpointer *data)

{
  GtsVertex *x = data[BSI_DATA_VERTICES] ;
  gdouble *xq =  data[BSI_DATA_QUAD_NODES] ;
  gdouble *q =   data[BSI_DATA_QUAD_RULE] ;
  gint nq =     *((gint *)data[BSI_DATA_QUAD_NUMBER]) ;
  gint nf =     *((gint *)data[BSI_DATA_VECTOR_SIZE]) ;
  gint iorder = *((gint *)data[BSI_DATA_INTERP_ORDER]) ;
  gint fstr =   *((gint *)data[BSI_DATA_NODE_STRIDE]) ;
  gint *ntf =    data[BSI_DATA_TET_COUNT] ;
  gdouble *f =   data[BSI_DATA_NODE_DATA] ;
  gdouble *u =   data[BSI_DATA_VELOCITY] ;
  gdouble *fn =  data[BSI_DATA_VERTEX_DATA] ;
#ifdef BSI_INTERPOLATE_RBF
  bsi_interpolator_t *rbf = data[BSI_DATA_INTERPOLATOR] ;
  gdouble *work = data[BSI_DATA_WORKSPACE] ;
  GtvVolume *v = data[BSI_DATA_VOLUME] ;
#else /*BSI_INTERPOLATE_RBF*/
  bsi_polynomial_t *pb = data[BSI_DATA_INTERPOLATOR] ;
  bsi_polynomial_workspace_t *wpb = data[BSI_DATA_WORKSPACE]    ;
#endif /*BSI_INTERPOLATE_RBF*/
  gqr_rule_t *qp = data[BSI_DATA_QUADRATURE_PHI] ;
  gqr_rule_t *qt = data[BSI_DATA_QUADRATURE_TH] ;
  gqr_rule_t *qr = data[BSI_DATA_QUADRATURE_R] ;
  GtsVertex *vt[4] ;
  gdouble L[4], J, r[3], *xquad, *fquad ;
  gint i, j, str, idx[4], rot[] = {0, 1, 2, 3, 0, 1, 2} ;

  str = 3+nf ;
  J = fabs(gtv_tetrahedron_volume(GTV_TETRAHEDRON(c))) ;

  if ( J < 1e-12 ) return ;
  
  gtv_tetrahedron_vertices(GTV_TETRAHEDRON(c),
			   &(vt[0]), &(vt[1]), &(vt[2]), &(vt[3])) ;
  /*for local corrections on cell vertices*/
  bsi_cell_indices(c, x, &(idx[0]), &(idx[1]), &(idx[2]), &(idx[3])) ;
  
#ifdef BSI_INTERPOLATE_RBF
  /* bsi_rbf_parameter_fit_cell(c, x, fn, nf, nf, rbf, work) ; */
  bsi_rbf_parameter_fit(x, idx[0], v, fn, nf, nf, rbf, work) ;
#else /*BSI_INTERPOLATE_RBF*/
  bsi_polynomial_cell_interpolant(c, x, fn, nf, nf, 5, pb,
				  iorder, wpb) ;
#endif /*BSI_INTERPOLATE_RBF*/  

  /* if ( cell_match(vt) ) { */
  /*   for ( i = 0 ; i < 4 ; i ++ )  */
  /*     fprintf(stdout, "%lg %lg %lg\n",  */
  /* 	      GTS_POINT(vt[i])->x, GTS_POINT(vt[i])->y, GTS_POINT(vt[i])->z) ; */

  /*   for ( i = 0 ; i < bsi_poly_coefficient_3d_offset(pb->order+1) ; i ++ ) { */
  /*     fprintf(stdout, "%lg %lg %lg\n", */
  /* 	      pb->c[3*i+0], pb->c[3*i+1], pb->c[3*i+2]) ; */
  /*   } */
    
  /*   exit(0) ; */
  /* } */
  
  for ( i = 0 ; i < nq ; i ++ ) {
    L[1] = q[4*i+0] ; L[2] = q[4*i+1] ; L[3] = q[4*i+2] ;
    L[0] = 1.0 - L[1] - L[2] - L[3] ;
    xquad = &(xq[(*ntf)*str]) ;
    fquad = &(f[(*ntf)*fstr]) ;
    xquad[0] = 
      L[0]*GTS_POINT(vt[0])->x + L[1]*GTS_POINT(vt[1])->x + 
      L[2]*GTS_POINT(vt[2])->x + L[3]*GTS_POINT(vt[3])->x ;
    xquad[1] = 
      L[0]*GTS_POINT(vt[0])->y + L[1]*GTS_POINT(vt[1])->y + 
      L[2]*GTS_POINT(vt[2])->y + L[3]*GTS_POINT(vt[3])->y ;
    xquad[2] = 
      L[0]*GTS_POINT(vt[0])->z + L[1]*GTS_POINT(vt[1])->z + 
      L[2]*GTS_POINT(vt[2])->z + L[3]*GTS_POINT(vt[3])->z ;

    fquad[0] = fquad[1] = fquad[2] = 0.0 ;

#ifdef BSI_INTERPOLATE_RBF
    bsi_rbf_evaluate(rbf, x, xquad, nf, fquad) ;
#else /*BSI_INTERPOLATE_RBF*/
    bsi_polynomial_evaluate(pb, xquad, fquad) ;
#endif /*BSI_INTERPOLATE_RBF*/

    /* bsi_source_func_ring_gaussian(xquad, fquad, nf, NULL) ; */
    
    fquad[0] *= J*q[4*i+3] ;
    fquad[1] *= J*q[4*i+3] ;
    fquad[2] *= J*q[4*i+3] ;

#ifndef BSI_UNCORRECTED
    for ( j = 0 ; j < 4 ; j ++ ) {
    /* for ( j = 0 ; j < 0 ; j ++ ) { */
      gdouble rw[3], R ;

      r[0] = GTS_POINT(vt[j])->x - xquad[0] ;
      r[1] = GTS_POINT(vt[j])->y - xquad[1] ;
      r[2] = GTS_POINT(vt[j])->z - xquad[2] ;

      R = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]) ;
      R *= R*R*4.0*M_PI ;
    
      bsi_vector_cross(rw, r, fquad) ;
      u[idx[j]*3+0] += rw[0]/R ;
      u[idx[j]*3+1] += rw[1]/R ;
      u[idx[j]*3+2] += rw[2]/R ;
    }
#endif /*BSI_UNCORRECTED*/
    /* fprintf(stdout, */
    /* 	    "%1.16e %1.16e %1.16e %1.16e %1.16e %1.16e\n", */
    /* 	    xquad[0], xquad[1], xquad[2],  */
    /* 	    fquad[0], fquad[1], fquad[2] */
    /* 	    ) ; */
    
    (*ntf) ++ ;
  }
#ifdef BSI_UNCORRECTED
  return ;
#endif /*BSI_UNCORRECTED*/

  /*local volume interpolation and integration*/
  for ( i = 0 ; i < 4 ; i ++ ) {
    gpointer fdata[2] ;
    fdata[0] = pb ;
    fdata[1] = &(GTS_POINT(vt[rot[i+0]])->x) ;
    tq_tet_quad(&(GTS_POINT(vt[rot[i+0]])->x),
		&(GTS_POINT(vt[rot[i+1]])->x),
		&(GTS_POINT(vt[rot[i+2]])->x),
		&(GTS_POINT(vt[rot[i+3]])->x),
		&(gqr_rule_abscissa(qp,0)), 1, &(gqr_rule_weight(qp,0)), 1,
		gqr_rule_length(qp),
		&(gqr_rule_abscissa(qt,0)), 1, &(gqr_rule_weight(qt,0)), 1,
		gqr_rule_length(qt),
		&(gqr_rule_abscissa(qr,0)), 1, &(gqr_rule_weight(qr,0)), 1,
		gqr_rule_length(qr),
		tetrahedron_bs_func, fdata,
		/* pb, */
		&(u[idx[i]*3]), 3) ;
  /* tq_tet_quad(&(x[3*rot[i+0]]), &(x[3*rot[i+3]]), */
  /* 		&(x[3*rot[i+2]]), &(x[3*rot[i+1]]),  */
  /* 		&(gqr_rule_abscissa(qp,0)), 1, &(gqr_rule_weight(qp,0)), 1, */
  /* 		ngp, */
  /* 		&(gqr_rule_abscissa(qt,0)), 1, &(gqr_rule_weight(qt,0)), 1, */
  /* 		ngt, */
  /* 		&(gqr_rule_abscissa(qr,0)), 1, &(gqr_rule_weight(qr,0)), 1, */
  /* 		ngr, */
  /* 		volume_test_func, NULL, */
  /* 		&(q[i*nq]), nq) ; */
  }

#if 0
  /*local volume integration using interpolant*/
  for ( id = 0 ; id < 4 ; id ++ ) {
    /* gdouble V ; */
    /* V = gtv_tetrahedron_volume(GTV_TETRAHEDRON(c)) ; */
    
    bsi_tet_integrate_powers(&(GTS_POINT(vt[rot[id+0]])->x),
			     &(GTS_POINT(vt[rot[id+1]])->x),
  			     &(GTS_POINT(vt[rot[id+2]])->x),
			     &(GTS_POINT(vt[rot[id+3]])->x),
  			     iorder+1, 175, 1e-12, 12, P) ;

    /* fprintf(stderr, "%lg %lg\n", V, P[0]) ; */
    
    bsi_polynomial_shift_origin(pb, &(GTS_POINT(vt[rot[id]])->x), wpb) ;
    for ( n = 0 ; n <= iorder ; n ++ ) {
      for ( i = 0 ; i <= n ; i ++ ) {
	for ( j = 0 ; j <= n-i ; j ++ ) {
	  k = n - i - j ;
	  ic = bsi_poly_coefficient_3d_index(i,j  ,k) ;
	  ip = bsi_poly_coefficient_3d_index(i  ,j+1,k  ) ;
	  u[3*idx[rot[id]]+0] -= P[ip]*pb->c[3*ic+2] ;
	  /* u[3*idx[rot[id]]+0] += P[ip]*pb->c[3*ic+2] ; */
	  ip = bsi_poly_coefficient_3d_index(i  ,j  ,k+1) ;
	  u[3*idx[rot[id]]+0] += P[ip]*pb->c[3*ic+1] ;
	  /* u[3*idx[rot[id]]+0] -= P[ip]*pb->c[3*ic+1] ; */

	  ip = bsi_poly_coefficient_3d_index(i  ,j  ,k+1) ;
	  u[3*idx[rot[id]]+1] -= P[ip]*pb->c[3*ic+0] ;
	  /* u[3*idx[rot[id]]+1] += P[ip]*pb->c[3*ic+0] ; */
	  ip = bsi_poly_coefficient_3d_index(i+1,j  ,k  ) ;
	  u[3*idx[rot[id]]+1] += P[ip]*pb->c[3*ic+2] ;
	  /* u[3*idx[rot[id]]+1] -= P[ip]*pb->c[3*ic+2] ; */
	  
	  ip = bsi_poly_coefficient_3d_index(i+1,j  ,k  ) ;
	  u[3*idx[rot[id]]+2] -= P[ip]*pb->c[3*ic+1] ;
	  /* u[3*idx[rot[id]]+2] += P[ip]*pb->c[3*ic+1] ; */
	  ip = bsi_poly_coefficient_3d_index(i  ,j+1,k  ) ;
	  u[3*idx[rot[id]]+2] += P[ip]*pb->c[3*ic+0] ;
	  /* u[3*idx[rot[id]]+2] -= P[ip]*pb->c[3*ic+0] ; */
	}
      }
    }
  }
  
#endif
  /* exit(0) ; */
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

gint main(gint argc, gchar **argv)

{
  FILE *input, *output ;
  GtvVolume *v ;
  GtsVertex *x ;
  GtsPoint *xi ;
  gdouble *f, *q, *xq, *u, xg[3], ui[3] ;
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
  wbfmm_shift_operators_t *shifts ;
  gdouble xtree[3], xtmax[3], D, del, *fmmwork ;
  gint pstr, torder[32], order_r, order_s, order_max, depth ;
  gint sizew, level, nthreads ;
  
  progname = g_strdup(g_path_get_basename(argv[0])) ;
  timer = g_timer_new() ;

  nthreads = 2 ;
  
  input = stdin ;
  output = stdout ;

  /*interpolation parameters*/
  nnmax = 10 ; iorder = 3 ;
  /* mop_logging_init(stderr, "", G_LOG_LEVEL_DEBUG, NULL) ;   */
  
  /*local quadratures*/
  /* ngp = 16 ; ngt = 16 ; ngr = 16 ; */
  /* ngp = 12 ; ngt = 12 ; ngr = 12 ; */
  ngp = 4 ; ngt = 4 ; ngr = 4 ;
  
  /*tree parameters*/
  del = 1e-2 ;
  order_r = 8 ; order_s = 8 ;
  depth = 6 ;
  tree = NULL ;
  
  fprintf(stderr, "%s: reading nodes; %lg\n",
	  progname, g_timer_elapsed(timer, NULL)) ;
  read_nodes(input, &x, &f, &nx, &nf) ;

  fprintf(stderr, "%s: tetrahedralizing %d nodes; %lg\n",
	  progname, nx, g_timer_elapsed(timer, NULL)) ;
  v = bsi_make_delaunay_volume(x, nx) ;

  order = 4 ;
  bsi_quadrature_js_select(order, &q, &nq) ;
  fprintf(stderr, "%s: %dth order quadrature, %d nodes; %lg\n",
	  progname, order, nq, g_timer_elapsed(timer, NULL)) ;

  ntet = gtv_volume_cell_number(v) ;

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
  bsi_points_origin_width(x, nx, xtree, xtmax, &D, FALSE) ;
  xtree[0] -= del ; xtree[1] -= del ; xtree[2] -= del ;
  D += 2.0*del ;

  pstr = str*sizeof(gdouble) ;
  tree = wbfmm_tree_new(xtree, D, ntf) ;
  torder[2*depth+0] = order_s ; 
  torder[2*depth+1] = order_r ; 
  order_max = MAX(order_s, order_r) ;
  for ( i = depth-1 ; i > 0 ; i -- ) {
    torder[2*i+0] = torder[2*(i+1)+0] + 4 ;
    torder[2*i+1] = torder[2*(i+1)+1] + 4 ;
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
  
  gtv_volume_foreach_cell(v, (GtsFunc)quadrature_nodes, data) ;

  /* return 0 ; */
  
  /* g_assert(ntf == ntet*nq) ; */

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
  wbfmm_tree_add_points(tree, (gpointer)xq, pstr, NULL, 0, ntf) ;
  for ( i = 0 ; i < depth ; i ++ ) wbfmm_tree_refine(tree) ;

  wbfmm_tree_problem(tree) = WBFMM_PROBLEM_LAPLACE ;
  wbfmm_tree_source_size(tree) = nf ;
  for ( i = 1 ; i <= depth ; i ++ ) {
    wbfmm_tree_laplace_coefficient_init(tree, i,
					torder[2*i+1], torder[2*i+0]) ;
  }

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
			&(GTS_POINT(&(x[i]))->x), &(u[3*i])) ;    
    }
    fprintf(stderr, "%s: point source summation completed; %lg\n",
	    progname, g_timer_elapsed(timer, NULL)) ;
  } else {
    fprintf(stderr, "%s: starting FMM point source summation; %lg\n",
	    progname, g_timer_elapsed(timer, NULL)) ;

    for ( i = 0 ; i < nx ; i ++ ) {
      guint64 box ;
      box = wbfmm_point_box(tree, tree->depth, &(GTS_POINT(&(x[i]))->x)) ;
      wbfmm_tree_laplace_box_local_curl(tree, tree->depth, box,
					&(GTS_POINT(&(x[i]))->x),
					&(u[3*i]), 3,
					&(xq[3]), str,
					NULL, 0, NULL, 0,
					TRUE, fmmwork) ;
    }
    
    fprintf(stderr, "%s: FMM point source summation completed; %lg\n",
	    progname, g_timer_elapsed(timer, NULL)) ;
  }

  /* fprintf(stderr, "%lg %lg %lg %1.16e %1.16e %1.16e\n", */
  /* 	  GTS_POINT(&(x[0]))->x, GTS_POINT(&(x[0]))->y, GTS_POINT(&(x[0]))->z, */
  /* 	  u[0], u[1], u[2]) ; */

  /* return 0 ; */
  
  ngx = 0 ; ngz = 1025 ;
  xi = gts_point_new(gts_point_class(), 0, 0, 0) ;
  for ( i = 0 ; i <= ngx ; i ++ ) {
    xg[0] = 1.0 ; /* -2 + (4.0)*i/ngx ; */
    xg[1] = 0.0 ;
    for ( j = 0 ; j <= ngz ; j ++ ) {
      xg[2] = -1 + (2.0)*j/ngz ;
      ui[0] = ui[1] = ui[2] = 0.0 ;

      xi->x = xg[0] ; xi->y = xg[1] ; xi->z = xg[2] ; 
      bsi_mesh_interp_linear(x, v, u, 3, 3, xi, ui) ;
    
      fprintf(output,
  	      "%e %e %e %e %e %e\n",
  	      xg[0], xg[1], xg[2],
  	      ui[0], ui[1], ui[2]) ;
    }
  }

  return 0 ;
}
