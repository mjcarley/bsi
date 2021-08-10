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
#include <string.h>
#include <math.h>

#include <glib.h>

#include <bsi.h>

#include <mop.h>

#include "quadrature.h"

gchar *progname ;

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

static void quadrature_nodes(GtvCell *c, gpointer *data)

{
  /* GtsVertex *x = data[0] ; */
  gdouble *xq = data[1] ;
  gint *ntf = data[2] ;
  gdouble *q = data[3] ;
  gint nq = *((gint *)data[4]) ;
  gint nf = *((gint *)data[5]) ;
  gdouble *V = data[6] ;
  GtsVertex *v1, *v2, *v3, *v0 ;
  gdouble L[4], J ;
  gint i, str ;

  str = 3+nf ;
  gtv_tetrahedron_vertices(GTV_TETRAHEDRON(c), &v0, &v1, &v2, &v3) ;

  J = fabs(gtv_tetrahedron_volume(GTV_TETRAHEDRON(c))) ;
  
  for ( i = 0 ; i < nq ; i ++ ) {
    L[1] = q[4*i+0] ; L[2] = q[4*i+1] ; L[3] = q[4*i+2] ;
    L[0] = 1.0 - L[1] - L[2] - L[3] ;
    xq[(*ntf)*str+0] =
      L[0]*GTS_POINT(v0)->x + L[1]*GTS_POINT(v1)->x + 
      L[2]*GTS_POINT(v2)->x + L[3]*GTS_POINT(v3)->x ;
    xq[(*ntf)*str+1] =
      L[0]*GTS_POINT(v0)->y + L[1]*GTS_POINT(v1)->y + 
      L[2]*GTS_POINT(v2)->y + L[3]*GTS_POINT(v3)->y ;
    xq[(*ntf)*str+2] =
      L[0]*GTS_POINT(v0)->z + L[1]*GTS_POINT(v1)->z + 
      L[2]*GTS_POINT(v2)->z + L[3]*GTS_POINT(v3)->z ;

    xq[(*ntf)*str+3] = J*q[4*i+3] ;
    
    (*ntf) ++ ;
    
    *V += q[4*i+3]*J ;
  }
  
  return ;
}

gint init_node_source(gdouble *x, gint xstr, gdouble *f, gint nf,
		      gint fstr, gint nx,
		      bsi_source_func_t sfunc, gpointer data)

{
  gint i, j ;
  gdouble w ;

  for ( i = 0 ; i < nx ; i ++ ) {
    w = f[i*fstr] ;
    sfunc(&(x[i*xstr]), &(f[i*fstr]), nf, data) ;
    for ( j = 0 ; j < nf ; j ++ ) f[i*fstr+j] *= w ;
  }
  
  return 0 ;
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

gint find_nearest_vertex(GtsVertex *x, gint nx, gdouble *xs)

{
  gint i, n ;
  gdouble R2, R2min ;

  R2min = G_MAXDOUBLE ; n = 0 ;

  for ( i = 0 ; i < nx ; i ++ ) {
    R2 =
      (GTS_POINT(&(x[i]))->x - xs[0])*(GTS_POINT(&(x[i]))->x - xs[0]) +
      (GTS_POINT(&(x[i]))->y - xs[1])*(GTS_POINT(&(x[i]))->y - xs[1]) +
      (GTS_POINT(&(x[i]))->z - xs[2])*(GTS_POINT(&(x[i]))->z - xs[2]) ;
    if ( R2 < R2min ) { n = i ; R2min = R2 ; }
  }
  
  return n ;
}

gint bsi_source_func_quadratic(gdouble *x, gdouble *f, gint nf,
			       gpointer data)

{
  f[0] = x[0]*x[0] ;
  f[1] = x[1]*x[1] ;
  f[2] = x[2]*x[2] ;

  return 0 ;
}

gint polynomial_test(gint N)

{
  gdouble p[1024], ps[1024], f[8], fs[8], x[3], xs[3], xc[3] ;
  gint i, np ;

  np = 3 ;
  x[0] = -0.9 ; x[1] = 0.3 ; x[2] = 1.7 ;
  
  for ( i = 0 ; i < bsi_poly_coefficient_3d_offset(N+1)*np ; i ++ )
    p[i] = g_random_double() ;
  
  bsi_poly_evaluate_3d(p, np, N, x, f) ;

  xc[0] = 0.1 ; xc[1] = -0.4 ; xc[2] = 3.1 ;
  bsi_poly_shift_3d(p, np, N, xc, ps) ;

  xs[0] = x[0] - xc[0] ;
  xs[1] = x[1] - xc[1] ;
  xs[2] = x[2] - xc[2] ;
  bsi_poly_evaluate_3d(ps, np, N, xs, fs) ;

  for ( i = 0 ; i < np ; i ++ ) 
    fprintf(stderr, "%lg %lg (%lg)\n", f[i], fs[i], fabs(f[i]-fs[i])) ;
  
  return 0 ;
}

static gint cell_neighbour_test(GtsVertex *x, gint i,
				GtvVolume *v)

{
  GtvCell *c ;
  GSList *cells ;
  gint nbrs[256], nnbrs, nnmax, j ;

  nnmax = 256 ;
  cells = gtv_vertex_cells(&(x[i]), v) ;
  c = cells->data ;
  
  bsi_cell_neighbours(c, x, 4, nbrs, &nnbrs, nnmax) ;

  for ( j = 0 ; j < nnbrs ; j ++ ) {
    fprintf(stdout, "%lg %lg %lg\n",
	    GTS_POINT(&(x[nbrs[j]]))->x,
	    GTS_POINT(&(x[nbrs[j]]))->y,
	    GTS_POINT(&(x[nbrs[j]]))->z) ;
  }
  
  return 0 ;
}

gint compare_distance(gconstpointer a, gconstpointer b,
		      gpointer data[])

{
  gint i = *((gint *)a) ;
  gint j = *((gint *)b) ;
  GtsVertex *x = data[0] ;
  gdouble *xc = data[1] ;
  gdouble ri, rj ;

  ri =
    (xc[0] - GTS_POINT(&(x[i]))->x)*(xc[0] - GTS_POINT(&(x[i]))->x) +
    (xc[1] - GTS_POINT(&(x[i]))->y)*(xc[1] - GTS_POINT(&(x[i]))->y) +
    (xc[2] - GTS_POINT(&(x[i]))->z)*(xc[2] - GTS_POINT(&(x[i]))->z) ;
  rj =
    (xc[0] - GTS_POINT(&(x[j]))->x)*(xc[0] - GTS_POINT(&(x[j]))->x) +
    (xc[1] - GTS_POINT(&(x[j]))->y)*(xc[1] - GTS_POINT(&(x[j]))->y) +
    (xc[2] - GTS_POINT(&(x[j]))->z)*(xc[2] - GTS_POINT(&(x[j]))->z) ;

  if ( ri < rj ) return -1 ;
  if ( ri > rj ) return  1 ;
  
  return 0 ;
}

static gint merge_interpolant(mop_polynomial_t *p, gdouble *c, gint nc,
			      gdouble *pm)

{
  gint n, i, j, k, idx ;
  gdouble pc ;
  
  n = mop_polynomial_order(p) ;

  memset(pm, 0, (n+1)*(n+2)*(n+3)/6*sizeof(gdouble)) ;

  for ( i = 0 ; i < mop_polynomial_term_number(p) ; i ++ ) {
    idx = bsi_poly_coefficient_3d_index(mop_polynomial_monomial_power(p,i,0),
					mop_polynomial_monomial_power(p,i,1),
					mop_polynomial_monomial_power(p,i,2)) ;
    for ( j = 0 ; j < mop_polynomial_term_number(p) ; j ++ ) {
      for ( k = 0 ; k < nc ; k ++ ) {
	pm[idx*nc+k] += c[j*nc+k]*mop_polynomial_coefficient(p,j,i) ;
      }
    }
  }
  
  return 0 ;
}

static gint cell_interpolation_test(GtsVertex *x, gint i,
				    GtvVolume *v,
				    gdouble *f, gint nf, gint fstr)

{
  GtvCell *c ;
  GtsVertex *v1, *v2, *v3, *v4 ;
  GSList *cells ;
  gint nbrs[256], nnbrs, nnmax, m, j, k, order, idx ;
  mop_polynomial_t *p ;
  mop_polynomial_workspace_t *wp ;
  gdouble work[2048], *xi, *fi, *ci, fc[3], xc[3], P[2048], fr[3], xs[3] ;
  gdouble intp[2048], intps[2048], fm[3], fms[3], L[4], xe[3] ;
  gpointer udata[2] ;
  bsi_polynomial_t *pb ;
  bsi_polynomial_workspace_t *wpb ;
  
  nnmax = 64 ; order = 4 ;
  cells = gtv_vertex_cells(&(x[i]), v) ;
  c = cells->data ;

  gtv_tetrahedron_vertices(GTV_TETRAHEDRON(c), &v1, &v2, &v3, &v4) ;
  
  L[0] = 0.1 ; L[1] = 0.3 ; L[2] = 0.05 ; L[3] = 1.0 - L[0] - L[1] - L[2] ;
  xe[0] = L[0]*GTS_POINT(v1)->x + L[1]*GTS_POINT(v2)->x +
    L[2]*GTS_POINT(v3)->x + L[3]*GTS_POINT(v4)->x ;
  xe[1] = L[0]*GTS_POINT(v1)->y + L[1]*GTS_POINT(v2)->y +
    L[2]*GTS_POINT(v3)->y + L[3]*GTS_POINT(v4)->y ;
  xe[2] = L[0]*GTS_POINT(v1)->z + L[1]*GTS_POINT(v2)->z +
    L[2]*GTS_POINT(v3)->z + L[3]*GTS_POINT(v4)->z ;

  bsi_source_func_ring_gaussian(xe, fr, 3, NULL) ;
  
  pb  = bsi_polynomial_new(order, nf) ;
  wpb = bsi_polynomial_workspace_new(order, nf, nnmax) ;

  /*generate local interpolant centred on cell centroid*/
  bsi_polynomial_cell_interpolant(c, x, f, nf, fstr, 5, pb, order, wpb) ;
  bsi_polynomial_evaluate(pb, xe, fc) ;

  /*shift polynomial origin to v2*/
  bsi_polynomial_shift_origin(pb, &(GTS_POINT(v2)->x), wpb) ;
  bsi_polynomial_evaluate(pb, xe, fm) ;

  fprintf(stderr, "%lg %lg %lg\n", fc[0], fc[1], fc[2]) ;
  fprintf(stderr, "%lg %lg %lg\n", fr[0], fr[1], fr[2]) ;
  fprintf(stderr, "%lg %lg %lg\n", fm[0], fm[1], fm[2]) ;
  fprintf(stderr, "(%lg %lg %lg)\n",
	  fabs(fm[0]-fc[0]), fabs(fm[1]-fc[1]), fabs(fm[2]-fc[2])) ;

  return 0 ;
  
  /*quadrature on tet*/
  bsi_tet_integrate_powers(&(GTS_POINT(v1)->x), &(GTS_POINT(v2)->x),
			   &(GTS_POINT(v3)->x), &(GTS_POINT(v4)->x),
			   order,
			   25, 1e-6, 4,
			   work) ;
  /* bsi_tet_integrate_powers(&(GTS_POINT(v2)->x), &(GTS_POINT(v3)->x), */
  /* 			   &(GTS_POINT(v4)->x), &(GTS_POINT(v1)->x), */
  /* 			   order, */
  /* 			   25, 1e-6, 4, */
  /* 			   work) ; */
  /* bsi_tet_integrate_powers(&(GTS_POINT(v3)->x), &(GTS_POINT(v4)->x), */
  /* 			   &(GTS_POINT(v1)->x), &(GTS_POINT(v2)->x), */
  /* 			   order, */
  /* 			   25, 1e-6, 4, */
  /* 			   work) ; */
  /* bsi_tet_integrate_powers(&(GTS_POINT(v4)->x), &(GTS_POINT(v1)->x), */
  /* 			   &(GTS_POINT(v2)->x), &(GTS_POINT(v3)->x), */
  /* 			   order, */
  /* 			   25, 1e-6, 4, */
  /* 			   work) ; */

  for ( m = 0 ; m <= order ; m ++ ) {
    for ( j = 0 ; j <= m ; j ++ ) {
      for ( k = 0 ; k <= m-j ; k ++ ) {
	idx = bsi_poly_coefficient_3d_index(j, k, (m-j-k)) ;
	fprintf(stderr, "%d %d %d %d %lg\n",
		idx, j, k, m-j-k, work[idx]) ;
      }
    }
  }    
  
  return 0 ;
}

gint main(gint argc, gchar **argv)

{
  FILE *input, *output ;
  GtvVolume *v ;
  GtsVertex *x ;
  gdouble *f, *q, *xq, V, u[3], xg[3], fg[3], fi[3], *work, emax[3], fmax[3] ;
  gdouble xt[3] ;
  gint nx, nf, order, nq, ntet, ntf, str, i, j, k, ngx, ngz ;
  bsi_source_func_t sfunc ;
  bsi_interpolator_t *interp ;
  gpointer data[16] ;
  
  progname = g_strdup(g_path_get_basename(argv[0])) ;

  input = stdin ;
  output = stdout ;

  /* polynomial_test(4) ; */
  
  /* return 0 ; */

  sfunc = bsi_source_func_ring_gaussian ;
  /* sfunc = bsi_source_func_quadratic ; */
  
  fprintf(stderr, "reading nodes\n") ;
  read_nodes(input, &x, &f, &nx, &nf) ;

  /* for ( i = 0 ; i < nx ; i ++ ) { */
  /*   sfunc(&(GTS_POINT(&(x[i]))->x), &(f[i*nf]), 3, NULL) ; */
  /* } */
  
  fprintf(stderr, "tetrahedralizing nodes\n") ;
  v = bsi_make_delaunay_volume(x, nx) ;
  
  order = 2 ;
  bsi_quadrature_js_select(order, &q, &nq) ;
  fprintf(stderr, "%s: %dth order quadrature, %d nodes\n",
	  progname, order, nq) ;

  ntet = gtv_volume_cell_number(v) ;

  xq = (gdouble *)g_malloc0(ntet*(3+nf)*nq*sizeof(gdouble)) ;

  ntf = 0 ;
  data[0] = x ;
  data[1] = xq ;
  data[2] = &ntf ;
  data[3] = q ;
  data[4] = &nq ;
  data[5] = &nf ;
  data[6] = &V ; V = 0.0 ;
  /* gtv_volume_foreach_cell(v, (GtsFunc)quadrature_nodes, data) ; */

  str = 3+nf ;
  init_node_source(xq, str, &(xq[3]), nf, str, ntf, sfunc, NULL) ;

  interp = bsi_interpolator_new(BSI_RBF_CP_C6, 2, 3,128) ;

  xg[0] = 0.6 ; xg[1] = 0.0 ; xg[2] = 0.0 ;
  /* xg[0] = M_SQRT1_2 + 0.2 ; xg[1] = M_SQRT1_2 ; xg[2] = 0.0 ; */
  i = find_nearest_vertex(x, nx, xg) ;
  fprintf(stderr, "node %d: %lg %lg %lg\n",
	  i,
	  GTS_POINT(&(x[i]))->x, GTS_POINT(&(x[i]))->y,
	  GTS_POINT(&(x[i]))->z) ;

  cell_interpolation_test(x, i, v, f, nf, nf) ;
  
  return 0 ;
  
  cell_neighbour_test(x, i, v) ;

  return 0 ;
  
  work = (gdouble *)g_malloc((interp->nnmax+interp->np)*
			     (interp->nnmax+interp->np)*sizeof(gdouble)) ;

  /* bsi_polynomial_fit(x, i, v, f, 3, 3, interp, work) ; */

  return 0 ;
  
  bsi_rbf_parameter_fit(x, i, v, f, 3, 3, interp, work) ;

  emax[0] = emax[1] = emax[2] = 
    fmax[0] = fmax[1] = fmax[2] = 0.0 ;

  xt[0] = xt[1] = xt[2] = 1.0/sqrt(3.0) ;
  
  for ( j = 0 ; j < 33 ; j ++ ) {
    xg[0] = GTS_POINT(&(x[i]))->x + ((gdouble)j/32-0.5)*xt[0]*0.01 ;
    xg[1] = GTS_POINT(&(x[i]))->y + ((gdouble)j/32-0.5)*xt[1]*0.01 ;
    xg[2] = GTS_POINT(&(x[i]))->z + ((gdouble)j/32-0.5)*xt[2]*0.01 ;

    bsi_rbf_evaluate(interp, x, xg, 3, fi) ;
    fprintf(stdout, "%lg %lg %lg", xg[0], xg[1], xg[2]) ;
    fprintf(stdout, " %lg %lg %lg", fi[0], fi[1], fi[2]) ;

    sfunc(xg, fg, 3, NULL) ;
    fprintf(stdout, " %lg %lg %lg\n", fg[0], fg[1], fg[2]) ;

    for ( k = 0 ; k < 3 ; k ++ ) {
      emax[k] = MAX(fabs(fi[k]-fg[k]), emax[k]) ;
      fmax[k] = MAX(fabs(fg[k]), fmax[k]) ;
    }
    
  }

  fprintf(stderr, "errors:") ;
  for ( i = 0 ; i < 3 ; i ++ ) fprintf(stderr, " %lg", emax[i]) ;
  fprintf(stderr, "\n") ;
  
  return 0 ;
  
  /* for ( i = 0 ; i < nx ; i ++ ) { */
  /*   u[0] = u[1] = u[2] = 0.0 ; */
  /*   bsi_eval_velocity(xq, str, &(xq[3]), str, ntf, &(x[i].p.x), u) ; */
    
  /*   fprintf(output, */
  /* 	    "%e %e %e %e %e %e\n", */
  /* 	    x[i].p.x, x[i].p.y, x[i].p.z, */
  /* 	    u[0], u[1], u[2]) ; */
  /* } */

  ngx = 0 ; ngz = 257 ;
  for ( i = 0 ; i <= ngx ; i ++ ) {
    xg[0] = 1.0 ; /* -2 + (4.0)*i/ngx ; */
    xg[1] = 0.0 ;
    for ( j = 0 ; j <= ngz ; j ++ ) {
      xg[2] = -1 + (2.0)*j/ngz ;
      u[0] = u[1] = u[2] = 0.0 ;
      bsi_eval_velocity(xq, str, &(xq[3]), str, ntf, xg, u) ;
    
      fprintf(output,
	      "%e %e %e %e %e %e\n",
	      xg[0], xg[1], xg[2],
	      u[0], u[1], u[2]) ;
    }
  }
  
  fprintf(stderr, "volume: %lg %lg (%lg)\n",
	  V, gtv_volume_volume(v), fabs(V - gtv_volume_volume(v))) ;

  g_assert(ntf == ntet*nq) ;
  
  /* fprintf(stderr, "writing volume\n") ; */
  /* gtv_volume_write(v, output) ; */

  return 0 ;
}
