/* This file is part of Tetquad, a set of functions for evaluating
 * integrals on tetrahedra.
 *
 * Copyright (C) 2021 Michael Carley
 *
 * Tetquad is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version. Tetquad is distributed in the
 * hope that it will be useful, but WITHOUT ANY WARRANTY; without even
 * the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 * PURPOSE.  See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Tetquad.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <glib.h>

#include "tetquad.h"

#define EPSILON 1e-15

gint index_sort3(gdouble *x, gint *idx)

{
  idx[0] = idx[1] = idx[2] = -1 ; 
  if ( x[0] <= x[1] && x[0] <= x[2] ) {
    idx[0] = 0 ;
    if ( x[1] <= x[2] ) {
      idx[1] = 1 ; idx[2] = 2 ;
    } else {
      idx[1] = 2 ; idx[2] = 1 ;
    }
  }

  if ( x[1] <= x[2] && x[1] <= x[0] ) {
    idx[0] = 1 ;
    if ( x[2] <= x[0] ) {
      idx[1] = 2 ; idx[2] = 0 ;
    } else {
      idx[1] = 0 ; idx[2] = 2 ;
    }
  }

  if ( x[2] <= x[0] && x[2] <= x[1] ) {
    idx[0] = 2 ;
    if ( x[0] <= x[1] ) {
      idx[1] = 0 ; idx[2] = 1 ;
    } else {
      idx[1] = 1 ; idx[2] = 0 ;
    }
  }

  if (idx[0] != -1 && idx[1] != -1 && idx[2] != -1) return 0 ;

  return -1 ;
  
  /* g_error("%s: %lg %lg %lg\n", __FUNCTION__, x[0], x[1], x[2]) ; */
  
  return 0 ;
}

gint multiply_matrices3x3(gdouble *A, gdouble *B, gdouble *C,
			  gboolean tc)

/*A = B*C or A = B*C'*/
  
{
  if ( !tc ) {
    A[0] = B[3*0+0]*C[3*0+0] + B[3*0+1]*C[3*1+0] + B[3*0+2]*C[3*2+0] ; 
    A[1] = B[3*0+0]*C[3*0+1] + B[3*0+1]*C[3*1+1] + B[3*0+2]*C[3*2+1] ; 
    A[2] = B[3*0+0]*C[3*0+2] + B[3*0+1]*C[3*1+2] + B[3*0+2]*C[3*2+2] ; 
    
    A[3] = B[3*1+0]*C[3*0+0] + B[3*1+1]*C[3*1+0] + B[3*1+2]*C[3*2+0] ; 
    A[4] = B[3*1+0]*C[3*0+1] + B[3*1+1]*C[3*1+1] + B[3*1+2]*C[3*2+1] ; 
    A[5] = B[3*1+0]*C[3*0+2] + B[3*1+1]*C[3*1+2] + B[3*1+2]*C[3*2+2] ; 
    
    A[6] = B[3*2+0]*C[3*0+0] + B[3*2+1]*C[3*1+0] + B[3*2+2]*C[3*2+0] ; 
    A[7] = B[3*2+0]*C[3*0+1] + B[3*2+1]*C[3*1+1] + B[3*2+2]*C[3*2+1] ; 
    A[8] = B[3*2+0]*C[3*0+2] + B[3*2+1]*C[3*1+2] + B[3*2+2]*C[3*2+2] ;

    return 0 ;
  }

  /*multiplication with C transposed*/
  A[0] = B[3*0+0]*C[3*0+0] + B[3*0+1]*C[3*0+1] + B[3*0+2]*C[3*0+2] ; 
  A[1] = B[3*0+0]*C[3*1+0] + B[3*0+1]*C[3*1+1] + B[3*0+2]*C[3*1+2] ; 
  A[2] = B[3*0+0]*C[3*2+0] + B[3*0+1]*C[3*2+1] + B[3*0+2]*C[3*2+2] ; 
    
  A[3] = B[3*1+0]*C[3*0+0] + B[3*1+1]*C[3*0+1] + B[3*1+2]*C[3*0+2] ; 
  A[4] = B[3*1+0]*C[3*1+0] + B[3*1+1]*C[3*1+1] + B[3*1+2]*C[3*1+2] ; 
  A[5] = B[3*1+0]*C[3*2+0] + B[3*1+1]*C[3*2+1] + B[3*1+2]*C[3*2+2] ; 
    
  A[6] = B[3*2+0]*C[3*0+0] + B[3*2+1]*C[3*0+1] + B[3*2+2]*C[3*0+2] ; 
  A[7] = B[3*2+0]*C[3*1+0] + B[3*2+1]*C[3*1+1] + B[3*2+2]*C[3*1+2] ; 
  A[8] = B[3*2+0]*C[3*2+0] + B[3*2+1]*C[3*2+1] + B[3*2+2]*C[3*2+2] ;

  return 0 ;
}

gint multiply_matrix_vector3(gdouble *y, gdouble *A, gdouble *x,
				    gboolean tr)

{
  if ( !tr ) {
    y[0] = A[0]*x[0] + A[1]*x[1] + A[2]*x[2] ; 
    y[1] = A[3]*x[0] + A[4]*x[1] + A[5]*x[2] ; 
    y[2] = A[6]*x[0] + A[7]*x[1] + A[8]*x[2] ;

    return 0 ;
  }

  y[0] = A[0]*x[0] + A[3]*x[1] + A[6]*x[2] ; 
  y[1] = A[1]*x[0] + A[4]*x[1] + A[7]*x[2] ; 
  y[2] = A[2]*x[0] + A[5]*x[1] + A[8]*x[2] ;

  return 0 ;
}

gint rotation_x(gdouble th, gdouble *A)

{
  gdouble R[9], At[9] ;

  memcpy(At, A, 9*sizeof(gdouble)) ;
  
  R[0] = 1.0 ; R[1] = 0.0 ; R[2] = 0.0 ;
  R[3] = 0.0 ; R[4] = cos(th) ; R[5] = -sin(th) ;
  R[6] = 0.0 ; R[7] = sin(th) ; R[8] =  cos(th) ;

  multiply_matrices3x3(A, R, At, FALSE) ;
  
  return 0 ;
}

gint rotation_y(gdouble th, gdouble *A)

{
  gdouble R[9], At[9] ;

  memcpy(At, A, 9*sizeof(gdouble)) ;

  R[0] =  cos(th) ; R[1] = 0.0 ; R[2] = sin(th) ;
  R[3] =      0.0 ; R[4] = 1.0 ; R[5] =     0.0 ;
  R[6] = -sin(th) ; R[7] = 0.0 ; R[8] = cos(th) ;
  
  multiply_matrices3x3(A, R, At, FALSE) ;

  return 0 ;
}

gint rotation_z(gdouble th, gdouble *A)

{
  gdouble R[9], At[9] ;

  memcpy(At, A, 9*sizeof(gdouble)) ;

  R[0] =  cos(th) ; R[1] = -sin(th) ; R[2] = 0.0 ;
  R[3] =  sin(th) ; R[4] =  cos(th) ; R[5] = 0.0 ;
  R[6] = 0.0 ; R[7] = 0.0 ; R[8] = 1.0 ;
  
  multiply_matrices3x3(A, R, At, FALSE) ;

  return 0 ;
}

/* gdouble turning_point(gdouble *x1, gdouble *x2) */

/* /\* */
/*  * find the point on a line through two points where \phi is minimum */
/*  * or maximum, d\cos\phi/d u = 0, \cos\phi = x_{3}/|x| */
/*  *\/ */
  
/* { */
/*   gdouble u, v[3], d[3] ; */

/*   tq_vector_init(d, x1, x2) ; */

/*   v[0] = d[2]*x1[0] - x1[2]*d[0] ; */
/*   v[1] = d[2]*x1[1] - x1[2]*d[1] ; */
/*   v[2] = d[2]*x1[2] - x1[2]*d[2] ; */

/*   u = -tq_vector_scalar(x1,v)/tq_vector_scalar(d,v) ; */
  
/*   return u ; */
/* } */

static gboolean intersect_triangle(gdouble orig[3], gdouble dir[3],
				   gdouble vert0[3], gdouble vert1[3], 
				   gdouble vert2[3], gdouble *t)
				   /* gdouble *u, gdouble *v) */
/*
 * Taken from Moller and Trumbore, Journal of Graphics Tools, 1997,
 * 2(1):21--28, http://jgt.akpeters.com/papers/MollerTrumbore97/
 */

{
  gdouble edge1[3], edge2[3], tvec[3], pvec[3], qvec[3]; 
  gdouble det, inv_det;
  gdouble u, v ;
  
  /* find vectors for two edges sharing vert0 */ 
  tq_vector_init(edge1, vert1, vert0);
  tq_vector_init(edge2, vert2, vert0);

  /* begin calculating determinant - also used to calculate U parameter */ 
  tq_vector_cross(pvec, dir, edge2);
  
  /* if determinant is near zero, ray lies in plane of triangle */ 
  det = tq_vector_scalar(edge1, pvec);

  /*this is the in-plane check*/
  if (det > -EPSILON && det < EPSILON) return FALSE; 

  inv_det = 1.0 / det ;
  /* calculate distance from vert0 to ray origin */ 
  tq_vector_init(tvec, orig, vert0) ;

  /* calculate U parameter and test bounds */ 
  /* *u = tq_vector_scalar(tvec, pvec) * inv_det;  */
  /* if (*u < 0.0 || *u > 1.0) return FALSE ; */
  u = tq_vector_scalar(tvec, pvec) * inv_det; 
  if (u < 0.0 || u > 1.0) return FALSE ;

  /* prepare to test V parameter */ 
  tq_vector_cross(qvec, tvec, edge1);

  /* calculate V parameter and test bounds */ 
  /* *v = tq_vector_scalar(dir, qvec) * inv_det;  */
  /* if (*v < 0.0 || *u + *v > 1.0) return FALSE  ; */
  v = tq_vector_scalar(dir, qvec) * inv_det; 
  if (v < 0.0 || u + v > 1.0) return FALSE  ;

  /* calculate t, ray intersects triangle */ 
  *t = tq_vector_scalar(edge2, qvec) * inv_det; 

  return TRUE ;
}

gdouble segment_intersect_2d(gdouble *x1, gdouble *x2, gdouble C, gdouble S)

{
  gdouble d[2], u ;

  d[0] = x2[0] - x1[0] ; d[1] = x2[1] - x1[1] ;

  u = (x1[0]*S - x1[1]*C)/(d[1]*C - d[0]*S) ;
  
  return u ;
}

gint segment_interp(gdouble *x1, gdouble *x2, gdouble u, gdouble *x)

{
  x[0] = x1[0] + u*(x2[0] - x1[0]) ;
  x[1] = x1[1] + u*(x2[1] - x1[1]) ;
  x[2] = x1[2] + u*(x2[2] - x1[2]) ;
  
  return 0 ;
}

gint in_plane_points(gdouble *x, gint *idx)

{
  gint n, j ;
  
  /*find ordering of points, with in-plane points first*/
  n = 0 ; j = -1 ;

  if ( x[3*0+2] < 1e-12 ) {
    idx[n] = 0 ; n ++ ;
  } else {
    j = 0 ;
  }
  
  if ( x[3*1+2] < 1e-12 ) {
    idx[n] = 1 ; n ++ ;
  } else {
    j = 1 ;
  }

  if ( x[3*2+2] < 1e-12 ) {
    idx[n] = 2 ; n ++ ;
  } else {
    j = 2 ;
  }

  idx[2] = j ;
  
  g_assert(n == 2) ;
  g_assert(j != -1) ;

  if ( x[3*idx[0]+1] < 1e-12 ) return 0 ;

  /*swap first two rows*/
  n = idx[0] ; idx[0] = idx[1] ; idx[1] = n ;
  
  return 0 ;
}

gint order_points(gdouble *x, gint *idx)

{
  gint j ;
  gdouble th[3] ;

  for ( j = 0 ; j < 3 ; j ++ ) {
    th[j] = atan2(x[3*j+1], x[3*j+0]) ;
  }

  index_sort3(th, idx) ;

  /* j = idx[1] ; idx[1] = idx[2] ; idx[2] = j ; */
  
  return 0 ;
}

static gint subtet_quad_th(gdouble *x, gdouble *xs,
			   gdouble th0, gdouble th1,
			   gint i0, gint i1,
			   gint j0, gint j1,
			   gdouble *qph, gint qpstr, gdouble *wph, gint wpstr,
			   gint nph,
			   gdouble *qth, gint qtstr, gdouble *wth, gint wtstr,
			   gint nth,
			   gdouble *qr, gint qrstr, gdouble *wr, gint wrstr,
			   gint nr,
			   tq_tetquad_func_t qfunc, gpointer qdata,
			   gdouble *q, gint nq, gdouble *A)

{
  gint i, j, k ;
  gdouble th, dth, thb, u, Cth, Sth, xi[3] ;
  gdouble ph0, ph1, phb, dph, ph, Cph, Sph ;
  gdouble r2, r, dr, rb ;
  gdouble s[3], orig[3] = {0,0,0}, wt, yg[3], sg[3] ;
  
  thb = 0.5*(th1 + th0) ; dth = 0.5*(th1 - th0) ;
  dth = fabs(dth) ;

  /* fprintf(stderr, "dth: %lg\n", dth) ; */
  
  if  ( dth < 1e-6 ) return 0 ;
  
  for ( i = 0 ; i < nth ; i ++ ) {
    th = thb + dth*qth[i*qtstr] ;
    Cth = cos(th) ; Sth = sin(th) ;
    u = segment_intersect_2d(&(x[3*i0]), &(x[3*i1]), Cth, Sth) ;
    if ( !BETWEEN(u,0.0,1.0) ) return -1 ;
      /* g_error("%s: u out of bounds: %lg (%lg)", __FUNCTION__, u, u-1.0) ; */
    segment_interp(&(x[3*i0]), &(x[3*i1]), u, xi) ;
    ph0 = acos(xi[2]/sqrt(tq_vector_scalar(xi,xi))) ;

    u = segment_intersect_2d(&(x[3*j0]), &(x[3*j1]), Cth, Sth) ;
    /* g_assert(u >= 0 && u <= 1) ; */
    if ( !BETWEEN(u,0.0,1.0) ) return -1 ;
      /* g_error("%s: u out of bounds: %lg (%lg)", __FUNCTION__, u, u-1.0) ; */
    segment_interp(&(x[3*j0]), &(x[3*j1]), u, xi) ;
    ph1 = acos(xi[2]/sqrt(tq_vector_scalar(xi,xi))) ;
    phb = 0.5*(ph1 + ph0) ; dph = 0.5*(ph1 - ph0) ;
    dph = fabs(dph) ;

    for ( j = 0 ; j < nph ; j ++ ) {
      ph = phb + dph*qph[j*qpstr] ;
      Cph = cos(ph) ; Sph = sin(ph) ; 
      s[0] = Sph*Cth ;
      s[1] = Sph*Sth ;
      s[2] = Cph ;
      if ( !intersect_triangle(orig, s,
			       &(x[3*0]), &(x[3*1]), &(x[3*2]),
			       &r2) ) {
	return 3 ;
	/* g_error("%s: intersection failed", __FUNCTION__) ; */
      }
      r2 = fabs(r2) ;
      rb = 0.5*r2 ; dr = 0.5*r2 ;
      for ( k = 0 ; k < nr ; k ++ ) {
	r = rb + dr*qr[k*qrstr] ;
	/*global coordinate system*/
	multiply_matrix_vector3(sg, A,  s, TRUE) ;
	yg[0] = r*sg[0] + xs[0] ; 
	yg[1] = r*sg[1] + xs[1] ; 
	yg[2] = r*sg[2] + xs[2] ; 
	wt = r*r*Sph*wth[i*wtstr]*wph[j*wpstr]*wr[k*wrstr]*dph*dth*dr ;
	qfunc(yg, sg, r, wt, q, nq, qdata) ;
      }
    }    
  }
  
  return 0 ;
}

gint zero_points(gdouble *x, gdouble tol, gint n)

{
  gint i ;

  for ( i = 0 ; i < n ; i ++ ) if ( fabs(x[i]) < tol ) x[i] = 0.0 ;
  
  return 0 ;
}

gint tet_quad_origin_th(gdouble *x, gdouble *xs,
			       gdouble *qph, gint qpstr, gdouble *wph,
			       gint wpstr, gint nph,
			       gdouble *qth, gint qtstr, gdouble *wth,
			       gint wtstr, gint nth,
			       gdouble *qr, gint qrstr, gdouble *wr,
			       gint wrstr, gint nr,
			       tq_tetquad_func_t qfunc, gpointer qdata,
			       gdouble *q, gint nq)

{
  gint i, perm[3], bounds[8] ;
  gdouble th[3], xr[9], thlim[4], A[9] ;
  
  tq_transform_matrix(x, A) ;
  
  /*tetrahedron nodes in rotated coordinate system*/
  multiply_matrices3x3(xr, x, A, TRUE) ;  

  zero_points(xr, 1e-15, 9) ;
  
  order_points(xr, perm) ;
  
  /* fprintf(stderr, "xr = [%lg %lg %lg\n%lg %lg %lg\n%lg %lg %lg] ;\n", */
  /* 	  xr[3*perm[0]+0], xr[3*perm[0]+1], xr[3*perm[0]+2], */
  /* 	  xr[3*perm[1]+0], xr[3*perm[1]+1], xr[3*perm[1]+2], */
  /* 	  xr[3*perm[2]+0], xr[3*perm[2]+1], xr[3*perm[2]+2]) ; */
  
  /*find \theta limits and corresponding line segments*/
  for ( i = 0 ; i < 3 ; i ++ ) {
    th[i] = atan2(xr[3*i+1], xr[3*i+0]) ;
  }

  if ( BETWEEN(th[perm[2]], th[perm[0]], th[perm[1]]) ) {
    /* fprintf(stderr, "between\n") ; */
    thlim[0] = th[perm[0]] ; thlim[1] = th[perm[2]] ;
    bounds[0] = perm[0] ; bounds[1] = perm[1] ; 
    bounds[2] = perm[0] ; bounds[3] = perm[2] ; 

    thlim[2] = th[perm[2]] ; thlim[3] = th[perm[1]] ;
    bounds[4] = perm[0] ; bounds[5] = perm[1] ; 
    bounds[6] = perm[1] ; bounds[7] = perm[2] ; 
  } else {
    /* fprintf(stderr, "outside\n") ; */
    thlim[0] = th[perm[0]] ; thlim[1] = th[perm[1]] ;
    bounds[0] = perm[0] ; bounds[1] = perm[1] ; 
    bounds[2] = perm[0] ; bounds[3] = perm[2] ; 

    thlim[2] = th[perm[1]] ; thlim[3] = th[perm[2]] ;
    bounds[4] = perm[1] ; bounds[5] = perm[2] ; 
    bounds[6] = perm[0] ; bounds[7] = perm[2] ; 
  }

  /* fprintf(stderr, "subtet 1\n") ; */
  /* bounds[0] = 0 ; bounds[1] = 1 ; */
  /* bounds[2] = 2 ; bounds[3] = 1 ; */
  
  i = subtet_quad_th(xr, xs, thlim[0], thlim[1],
		     bounds[0], bounds[1], bounds[2], bounds[3],
		     qph, qpstr, wph, wpstr, nph,
		     qth, qtstr, wth, wtstr, nth,
		     qr,  qrstr, wr,  wrstr, nr,
		     qfunc, qdata,
		     q, nq, A) ;
  if ( i != 0 ) return 1 ;
  /* return 0 ; */
  /* fprintf(stderr, "subtet 2\n") ; */
  i = subtet_quad_th(xr, xs, thlim[2], thlim[3],
		     bounds[4], bounds[5], bounds[6], bounds[7],
		     qph, qpstr, wph, wpstr, nph,
		     qth, qtstr, wth, wtstr, nth,
		     qr,  qrstr, wr,  wrstr, nr,
		     qfunc, qdata,
		     q, nq, A) ;
  if ( i != 0 ) return 2 ;
  
  return 0 ;
}

gdouble select_rotation(gdouble *xr, gdouble *th, gdouble *dih,
			gint j, gint k)

{
  gdouble thr, xb[2], dth ;

  /*check quadrants*/
  if ( xr[3*j+2] >= 0 && xr[3*k+2] >= 0 ) {
    /* g_assert_not_reached() ;  */
    thr = -th[j] + 0.5*M_PI ;
    g_assert(th[k] + thr > -0.5*M_PI) ;
    return thr ;
  }

  if ( xr[3*j+2] >= 0 && xr[3*k+2] < 0 ) {
    /* g_assert_not_reached() ;  */
    /* thr = -0.5*M_PI - th[k] ; */
    thr = -0.5*M_PI - th[j] ;
    /* g_assert(th[j] + thr < 0.5*M_PI) ; */
    /* if ( th[j] + thr > 0.5*M_PI ) { */
    dth = th[k] + thr ;
    if ( dth < -0.5*M_PI ) dth += 2.0*M_PI ;
    if ( dth >  0.5*M_PI+1e-12 ||
	 dth < -0.5*M_PI-1e-12 ) {
    
    /* if ( th[k] + thr >  0.5*M_PI+1e-12 || */
    /* 	 th[k] + thr < -0.5*M_PI-1e-12 ) { */
      thr += M_PI ;
    }
    return thr ;
  }

  if ( xr[3*j+2] < 0 && xr[3*k+2] >= 0 ) {
    /* g_assert_not_reached() ;  */
    thr = -th[j] - 0.5*M_PI ;
    if ( th[k] + thr >  0.5*M_PI ||
	 th[k] + thr < -0.5*M_PI ) {
    /* if ( th[k] + thr > 0.5*M_PI ) { */
      thr += M_PI ;
    }
    
    return thr ;
  }

  if ( xr[3*j+2] < 0 && xr[3*k+2] < 0 ) {
    /* g_assert_not_reached() ;  */
    thr = -th[j] - 0.5*M_PI ;
    /*check node k*/
    g_assert ( th[k] + thr <= 0.5*M_PI ) ;
    return thr ;
  }

  g_assert_not_reached() ;
  
  xb[0] = 0.5*(xr[3*j+1] + xr[3*k+1]) ;
  xb[1] = 0.5*(xr[3*j+2] + xr[3*k+2]) ;

  thr = -atan2(xb[1], xb[0]) ;

  return thr ;
}
  
gint matrix_define(gdouble *x, gdouble *A, gint *p)

{
  gdouble th, thmin, xr[9], dih[3], ph[3], phr ;
  gint i, j, k, idih[3] ;

  tq_tet_dihedral_angles(&(x[3*0]), &(x[3*1]), &(x[3*2]),
  			 &(dih[0]), &(dih[1]), &(dih[2]),
  			 NULL, NULL, NULL) ;

  if ( index_sort3(dih, idih) != 0 ) {
    fprintf(stderr, "dihedral error:") ;
    fprintf(stderr,
	    "%lg %lg %lg\n"
	    "%lg %lg %lg\n"
	    "%lg %lg %lg\n",
	    x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8]) ;
  }
  /* fprintf(stderr, "dih: %d %d %d\n", idih[0], idih[1], idih[2]) ; */
  /*align the edge with largest dihedral with the x axis*/
  i = idih[2] ; j = idih[0] ; k = idih[1] ;
  /*rotate to bring node to z==0*/
  rotation_y(atan2(x[3*i+2], x[3*i+0]), A) ;
  
  multiply_matrices3x3(xr, x, A, TRUE) ;
  /* return 0 ; */

  /*rotate to bring node to y==0*/
  rotation_z(-atan2(xr[3*i+1], xr[3*i+0]), A) ;
  
  multiply_matrices3x3(xr, x, A, TRUE) ;

  if ( xr[3*0+1] > -1e-12 &&
       xr[3*1+1] > -1e-12 &&
       xr[3*2+1] > -1e-12 ) return 0 ;

  /*check for points with negative y*/
  ph[0] = atan2(xr[3*0+2], xr[3*0+1]) ;
  ph[1] = atan2(xr[3*1+2], xr[3*1+1]) ;
  ph[2] = atan2(xr[3*2+2], xr[3*2+1]) ;

  /* fprintf(stderr, "negative y (%lg, %lg, %lg)\n", ph[0], ph[1], ph[2]) ; */

  if ( xr[3*j+1] < -1e-12 && xr[3*k+1] < -1e-12 ) {
    /*both values of y negative, easier to make a 180 or a reflection*/
    /* fprintf(stderr, "flipping\n") ; */
    /* g_assert_not_reached() ; /\*unchecked code*\/ */
    rotation_x(M_PI, A) ;

    return 0 ;
  }

  if ( xr[3*k+1] < -1e-12 ) {
    phr = select_rotation(xr, ph, dih, k, j) ;
  } else {
    phr = select_rotation(xr, ph, dih, j, k) ;
  }
  
  /* if ( xr[3*k+1] < -1e-12 ) { */
  /*   phr = -ph[k] + 0.5*dih[k] ; */
  /*   if ( xr[3*k+2] <= 0.0 ) phr -= 0.5*M_PI ; */
  /* } */
  /* if ( xr[3*j+1] < -1e-12 ) { */
  /*   phr = -ph[j] + 0.5*dih[j] ; */
  /*   if ( xr[3*j+2] <= 0.0 ) phr -= 0.5*M_PI ; */
  /* } */
  
  /*rotate to bring node to x==0*/
  rotation_x(phr, A) ;
  
  multiply_matrices3x3(xr, x, A, TRUE) ;

  return 0 ;
  
  /* i = p[0] ; j = p[1] ; k = p[2] ; */
  /* /\* fprintf(stderr, "%d %d %d\n", i, j, k) ; *\/ */
  /* /\*rotate to bring node to z==0*\/ */
  /* rotation_y(atan2(x[3*i+2], x[3*i+0]), A) ; */
  
  /* multiply_matrices3x3(xr, x, A, TRUE) ; */
  /* /\* return 0 ; *\/ */

  /* /\*rotate to bring node to y==0*\/ */
  /* rotation_z(-atan2(xr[3*i+1], xr[3*i+0]), A) ; */
  
  /* multiply_matrices3x3(xr, x, A, TRUE) ; */
  /* /\*rotate to bring one more node to z==0*\/ */
  /* rotation_x(-atan2(xr[3*j+2], xr[3*j+1]), A) ; */
  
  /* multiply_matrices3x3(xr, x, A, TRUE) ; */

  /* /\*check for negative theta*\/ */
  /* thmin = 0.0 ; */
  /* for ( i = 0 ; i < 3 ; i ++ ) { */
  /*   th = atan2(xr[3*i+1], xr[3*i+0]) ; */
  /*   if ( th < 0.0 ) thmin = MIN(thmin, th) ; */
  /* } */

  /* if ( thmin < 0.0 ) { */
  /*   rotation_z(-thmin, A) ; */
  
  /*   multiply_matrices3x3(xr, x, A, TRUE) ; */
  /* } */
  
  /* /\* return 0 ; *\/ */
  
  /* /\*check all y are positive (to tolerance)*\/ */
  /* if ( xr[3*j+1] > -1e-12 && xr[3*k+1] > -1e-12 ) return 0 ; */

  /* if ( xr[3*j+1] < -1e-12 && xr[3*k+1] < -1e-12 ) { */
  /*   /\*flip in y*\/ */
  /*   fprintf(stderr, "flip\n") ; */
  /*   A[3] *= -1 ; A[4] *= -1 ; A[5] *= -1 ;     */
  /*   return 0 ; */
  /* } */

  /* /\*one y is negative, one is positive*\/ */
  /* if ( xr[3*j+1] < -1e-12 ) { */
  /*   /\* th = atan2(xr[3*j+2], xr[3*j+0]) ; *\/ */
  /*   th = atan2(xr[3*j+1], xr[3*j+2]) ; */
  /* } else { */
  /*   /\* th = atan2(xr[3*k+2], xr[3*k+0]) ; *\/ */
  /*   th = atan2(xr[3*k+1], xr[3*k+2]) ; */
  /* } */
  
  /* rotation_x(th, A) ; */

  /* multiply_matrices3x3(xr, x, A, TRUE) ; */

  /* thmin = 0.0 ; */
  /* for ( i = 0 ; i < 3 ; i ++ ) { */
  /*   th = atan2(xr[3*i+1], xr[3*i+0]) ; */
  /*   if ( th < 0.0 ) thmin = MIN(thmin, th) ; */
  /* } */

  /* if ( thmin < 0.0 ) { */
  /*   rotation_z(-thmin, A) ; */
  
  /*   multiply_matrices3x3(xr, x, A, TRUE) ; */
  /* } */

  
  /* if ( xr[3*j+1] < -1e-12 ) { */
  /*   th = atan2(xr[3*j+2], xr[3*j+0]) ; */
  /* } else { */
  /*   th = atan2(xr[3*k+2], xr[3*k+0]) ; */
  /* } */
  
  /* rotation_x(-th, A) ; */

  /* multiply_matrices3x3(xr, x, A, TRUE) ; */
  
  /* if ( xr[3*j+1] < -1e-12 ) { */
  /*   th = atan2(xr[3*j+1], xr[3*j+0]) ; */
  /* } else { */
  /*   th = atan2(xr[3*k+1], xr[3*k+0]) ; */
  /* } */
  
  /* rotation_z(-th, A) ; */

  /* multiply_matrices3x3(xr, x, A, TRUE) ; */

  if ( xr[3*0+1] < -1e-12 ) return 1 ;
  if ( xr[3*1+1] < -1e-12 ) return 2 ;
  if ( xr[3*2+1] < -1e-12 ) return 3 ;
  
  return 0 ;
}

gint tq_transform_matrix(gdouble *x, gdouble *A)

{
  gint p[] = {0, 1, 2, 1, 2, 0, 2, 0, 1}, i ;

  for ( i = 0 ; i < 3 ; i ++ ) {
    /* fprintf(stderr, "%d\n", i) ; */

    A[0] = 1.0 ; A[1] = 0.0 ; A[2] = 0.0 ;
    A[3] = 0.0 ; A[4] = 1.0 ; A[5] = 0.0 ;
    A[6] = 0.0 ; A[7] = 0.0 ; A[8] = 1.0 ;

    /* if ( matrix_define(x, A, p) == 0 ) return 0 ; */

    /* if ( matrix_define(x, A, &(p[3*(2-i)])) == 0 ) { */
    /*   fprintf(stderr, "perm: %d\n", 2-i) ; */
    /*   return 0 ; */
    /* } */
    if ( matrix_define(x, A, &(p[3*i])) == 0 ) return 0 ;
  }

    g_assert_not_reached() ;
    
    /* g_assert(xr[3*0+1] > -1e-12) ; */
    /* g_assert(xr[3*1+1] > -1e-12) ; */
    /* g_assert(xr[3*2+1] > -1e-12) ; */
  
  return 0 ;
}

gint tq_tet_quad(gdouble *x1, gdouble *x2, gdouble *x3, gdouble *x4,
		 gdouble *qph, gint qpstr, gdouble *wph, gint wpstr,
		 gint nph,
		 gdouble *qth, gint qtstr, gdouble *wth, gint wtstr,
		 gint nth,
		 gdouble *qr, gint qrstr, gdouble *wr, gint wrstr,
		 gint nr,
		 tq_tetquad_func_t qfunc, gpointer qdata,		     
		 gdouble *q, gint nq)

/*
 * integrate function qfunc over tetrahedron x1,2,3,4 with x1 taken as
 * origin ("singular point")
 *
 * (qph,wph), etc are Gauss-Legendre quadrature rules
 *
 * integrals are added to array q
 */
  
{
  gdouble x[9] ;
  gint i ;
  
  tq_vector_init(&(x[3*0]), x1, x2) ;
  tq_vector_init(&(x[3*1]), x1, x3) ;
  tq_vector_init(&(x[3*2]), x1, x4) ;

  i = tet_quad_origin_th(x, x1,
			 qph, qpstr, wph, wpstr, nph,
			 qth, qtstr, wth, wtstr, nth,
			 qr , qtstr, wr , wtstr, nr ,
			 qfunc, qdata, q, nq) ;

  if ( i != 0 ) {
    fprintf(stderr,
	    "%1.16e %1.16e %1.16e\n"
	    "%1.16e %1.16e %1.16e\n"
	    "%1.16e %1.16e %1.16e\n"
	    "%1.16e %1.16e %1.16e\n",
	    x1[0], x1[1], x1[2],
	    x2[0], x2[1], x2[2],
	    x3[0], x3[1], x3[2],
	    x4[0], x4[1], x4[2]) ;
    g_error("%s: integration failure", __FUNCTION__) ;
    
  }
  
  return 0 ;
}

gint tq_tet_dihedral_angles(gdouble *xa, gdouble *xb, gdouble *xc,
			    gdouble *oa, gdouble *ob, gdouble *oc,
			    gdouble *ab, gdouble *bc, gdouble *ca)

/* 
 * From 
 * 
 * https://math.stackexchange.com/questions/49330/the-dihedral-angles-of-a-tetrahedron-in-terms-of-its-edge-lengths 
*/
  
{
  gdouble W, X, Y, Z, x[9], a2, b2, c2, d2, e2, f2, H2, J2, K2 ;

  a2 = tq_vector_length2(xa) ;
  b2 = tq_vector_length2(xb) ;
  c2 = tq_vector_length2(xc) ;
  d2 = tq_vector_distance2(xb,xc) ;
  e2 = tq_vector_distance2(xc,xa) ;
  f2 = tq_vector_distance2(xa,xb) ;

  /*four face areas*/
  tq_vector_cross(x, xb, xc) ;
  X = 0.5*sqrt(tq_vector_length2(x)) ;
  
  tq_vector_cross(x, xc, xa) ;
  Y = 0.5*sqrt(tq_vector_length2(x)) ;

  tq_vector_cross(x, xa, xb) ;
  Z = 0.5*sqrt(tq_vector_length2(x)) ;

  tq_vector_init(&(x[3]), xa, xb) ;
  tq_vector_init(&(x[6]), xb, xc) ;
  tq_vector_cross(x, &(x[3]), &(x[6])) ;
  W = 0.5*sqrt(tq_vector_length2(x)) ;

  H2 = (4.0*a2*d2 - ((b2+e2) - (c2+f2))*((b2+e2) - (c2+f2)))/16.0 ;
  J2 = (4.0*b2*e2 - ((c2+f2) - (a2+d2))*((c2+f2) - (a2+d2)))/16.0 ;
  K2 = (4.0*c2*f2 - ((a2+d2) - (b2+e2))*((a2+d2) - (b2+e2)))/16.0 ;

  *oa = (Y*Y + Z*Z - H2)/(2.0*Y*Z) ;
  *oa = acos(*oa) ;

  *ob = (Z*Z + X*X - J2)/(2.0*Z*X) ;
  *ob = acos(*ob) ;

  *oc = (X*X + Y*Y - K2)/(2.0*X*Y) ;
  *oc = acos(*oc) ;

  if ( ab == NULL ) return 0 ;
  
  *bc = (W*W + X*X - H2)/(2.0*W*X) ;
  *bc = acos(*bc) ;

  *ca = (W*W + Y*Y - J2)/(2.0*W*Y) ;
  *ca = acos(*ca) ;

  *ab = (W*W + Z*Z - K2)/(2.0*W*Z) ;
  *ab = acos(*ab) ;

  return 0 ;
}
