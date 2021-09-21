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

#include <blaswrap.h>

#include "tetquad.h"

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /*HAVE_CONFIG_H*/

#ifdef HAVE_AVX_INSTRUCTIONS
#include <immintrin.h>
#endif /*HAVE_AVX_INSTRUCTIONS*/

/*tolerance for in-plane check in ray-tracing calculation*/
#define EPSILON 1e-15
/*tolerance for zeroing entries in transformed coordinates*/
#define TOL_ZERO 1e-12
/*
 * tolerance for integrals in theta: do not integrate over intervals
 * smaller than this (treat interval as so small the integral is zero)
 */
#define TOL_THETA 1e-6

/* #define TRANSFORM_PHI */

#ifndef _GNU_SOURCE
#define sincos(_x,_S,_C)					\
  do { (*(_S)) = sin((_x)) ; (*(_C)) = cos((_x)) ; } while (0)
#endif

static gint index_sort3(gdouble *x, gint *idx)

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

static inline void multiply_matrices3x3(gdouble *A, gdouble *B, gdouble *C,
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

    return ;
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

  return ;
}

static inline void multiply_matrix_vector3(gdouble *y, gdouble *A, gdouble *x,
					   gboolean tr)

{
  if ( !tr ) {
    y[0] = A[0]*x[0] + A[1]*x[1] + A[2]*x[2] ; 
    y[1] = A[3]*x[0] + A[4]*x[1] + A[5]*x[2] ; 
    y[2] = A[6]*x[0] + A[7]*x[1] + A[8]*x[2] ;

    return ;
  }

  y[0] = A[0]*x[0] + A[3]*x[1] + A[6]*x[2] ; 
  y[1] = A[1]*x[0] + A[4]*x[1] + A[7]*x[2] ; 
  y[2] = A[2]*x[0] + A[5]*x[1] + A[8]*x[2] ;

  return ;
}

static inline void multiply_matrix4_vector3(gdouble *y, gdouble *A, gdouble *x)


{
  y[0] = A[0]*x[0] + A[1]*x[1] + A[2]*x[2] ; 
  y[1] = A[3]*x[0] + A[4]*x[1] + A[5]*x[2] ; 
  y[2] = A[6]*x[0] + A[7]*x[1] + A[8]*x[2] ;
  y[3] = A[9]*x[0] + A[10]*x[1] + A[11]*x[2] ;
  
  return ;
}

static void rotation_x(gdouble th, gdouble *A)

{
  gdouble R[9], At[9], C, S ;

  memcpy(At, A, 9*sizeof(gdouble)) ;

  sincos(th, &S, &C) ;
  R[0] = 1.0 ; R[1] = 0.0 ; R[2] = 0.0 ;
  R[3] = 0.0 ; R[4] = C ; R[5] = -S ;
  R[6] = 0.0 ; R[7] = S ; R[8] =  C ;

  multiply_matrices3x3(A, R, At, FALSE) ;
  
  return ;
}

static void rotation_y(gdouble th, gdouble *A)

{
  gdouble R[9], At[9], S, C ;

  memcpy(At, A, 9*sizeof(gdouble)) ;

  sincos(th, &S, &C) ;

  R[0] =  C   ; R[1] = 0.0 ; R[2] = S   ;
  R[3] =  0.0 ; R[4] = 1.0 ; R[5] = 0.0 ;
  R[6] = -S   ; R[7] = 0.0 ; R[8] = C   ;
  
  multiply_matrices3x3(A, R, At, FALSE) ;

  return ;
}

static void rotation_z(gdouble th, gdouble *A)

{
  gdouble R[9], At[9], S, C ;

  memcpy(At, A, 9*sizeof(gdouble)) ;
  sincos(th, &S, &C) ;

  R[0] = C   ; R[1] = -S   ; R[2] = 0.0 ;
  R[3] = S   ; R[4] =  C   ; R[5] = 0.0 ;
  R[6] = 0.0 ; R[7] =  0.0 ; R[8] = 1.0 ;
  
  multiply_matrices3x3(A, R, At, FALSE) ;

  return ;
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

#ifndef HAVE_AVX_INSTRUCTIONS
/* static gboolean intersect_triangle(gdouble dir[3], gdouble vert0[3], */
/* 				   gdouble edge1[3], gdouble edge2[3], */
/* 				   gdouble qvec[3], */
/* 				   gdouble *t) */

/* /\* */
/*  * Taken from Moller and Trumbore, Journal of Graphics Tools, 1997, */
/*  * 2(1):21--28, http://jgt.akpeters.com/papers/MollerTrumbore97/ */
/*  * */
/*  * adapted to return only distance to plane and to take precomputed */
/*  * edge vectors */
/*  *\/ */

/* { */
/*   gdouble *tvec, pvec[3], det, inv_det ; */
/* #ifdef RAY_TRACE_CHECK_BOUNDS */
/*   gdouble u, v ; */
/* #endif /\*RAY_TRACE_CHECK_BOUNDS*\/ */
  
/*   /\* begin calculating determinant - also used to calculate U parameter *\/  */
/*   tq_vector_cross(pvec, dir, edge2); */
  
/*   /\* if determinant is near zero, ray lies in plane of triangle *\/  */
/*   det = tq_vector_scalar(edge1, pvec); */

/*   /\*this is the in-plane check*\/ */
/*   if (det > -EPSILON && det < EPSILON) return FALSE;  */

/*   inv_det = 1.0 / det ; */
/*   /\* calculate distance from vert0 to ray origin *\/  */
/*   /\* tq_vector_init(tvec, orig, vert0) ; *\/ */
/*   tvec = vert0 ; */

/* #ifdef RAY_TRACE_CHECK_BOUNDS   */
/*   /\* calculate U parameter and test bounds *\/  */
/*   /\* u = tq_vector_scalar(tvec, pvec) * inv_det;  *\/ */
/*   /\* if (u < 0.0 || u > 1.0) return FALSE ; *\/ */
/* #endif /\*RAY_TRACE_CHECK_BOUNDS  *\/ */
  
/*   /\* prepare to test V parameter *\/  */
/*   tq_vector_cross(qvec, tvec, edge1); */

/* #ifdef RAY_TRACE_CHECK_BOUNDS     */
/*   /\* calculate V parameter and test bounds *\/  */
/*   v = tq_vector_scalar(dir, qvec) * inv_det; */
/*   if (v < 0.0 || u + v > 1.0) return FALSE  ; */
/* #endif /\*RAY_TRACE_CHECK_BOUNDS*\/ */
  
/*   /\* calculate t, ray intersects triangle *\/  */
/*   *t = tq_vector_scalar(edge2, qvec) * inv_det;  */

/*   return TRUE ; */
/* } */

static gdouble segment_intersect_2d(gdouble *x1, gdouble *d,
				    gdouble C, gdouble S)

/*
 * calculate intersection point in two dimensions of ray through
 * origin with direction (C,S) and line segment x1--(x1+d)
 */
  
{
  gdouble u ;

  u = -(x1[0]*S - x1[1]*C)/(d[0]*S - d[1]*C) ;
  
  return u ;
}

static void segment_interp(gdouble *x1, gdouble *d, gdouble u, gdouble *x)

/*
 * interpolate on line segment x1--(x1+d), 0 <= u <= 1
 */
  
{
  x[0] = x1[0] + u*d[0] ;
  x[1] = x1[1] + u*d[1] ;
  x[2] = x1[2] + u*d[2] ;
  
  return ;
}

#endif /*HAVE_AVX_INSTRUCTIONS*/


/* static void in_plane_points(gdouble *x, gint *idx) */

/* { */
/*   gint n, j ; */
  
/*   /\*find ordering of points, with in-plane points first*\/ */
/*   n = 0 ; j = -1 ; */

/*   if ( x[3*0+2] < TOL_ZERO ) { */
/*     idx[n] = 0 ; n ++ ; */
/*   } else { */
/*     j = 0 ; */
/*   } */
  
/*   if ( x[3*1+2] < TOL_ZERO ) { */
/*     idx[n] = 1 ; n ++ ; */
/*   } else { */
/*     j = 1 ; */
/*   } */

/*   if ( x[3*2+2] < TOL_ZERO ) { */
/*     idx[n] = 2 ; n ++ ; */
/*   } else { */
/*     j = 2 ; */
/*   } */

/*   idx[2] = j ; */
  
/*   g_assert(n == 2) ; */
/*   g_assert(j != -1) ; */

/*   if ( x[3*idx[0]+1] < TOL_ZERO ) return 0 ; */

/*   /\*swap first two rows*\/ */
/*   n = idx[0] ; idx[0] = idx[1] ; idx[1] = n ; */
  
/*   return ; */
/* } */

static void order_points(gdouble *x, gint *idx)

{
  gint j ;
  gdouble th[3] ;

  for ( j = 0 ; j < 3 ; j ++ ) {
    th[j] = atan2(x[3*j+1], x[3*j+0]) ;
  }

  index_sort3(th, idx) ;
  
  return ;
}

#ifndef HAVE_AVX_INSTRUCTIONS
static gint subtet_quad_th(gdouble *x, gdouble *xs,
			   gdouble th0, gdouble th1,
			   gint i0, gint i1,
			   gint j0, gint j1,
			   gdouble *qph, gint qpstr,
			   gdouble *wph, gint wpstr,
			   gint nph,
			   gdouble *qth, gint qtstr,
			   gdouble *wth, gint wtstr,
			   gint nth,
			   gdouble *qr, gint qrstr,
			   gdouble *wr, gint wrstr,
			   gint nr,
			   tq_tetquad_func_t qfunc, gpointer qdata,
			   gdouble *q, gint nq, gdouble *A)

{
  gint i, j, k ;
  gdouble th, dth, thb, u, Cth, Sth, xi[3], xj[3] ;
  gdouble ph0, ph1, phb, dph, ph, Cph, Sph ;
  gdouble r2, r, dr, rb ;
  gdouble s[3], wt, yg[3], sg[3] ;
  gdouble dxi[3], dxj[3] ;
  gdouble dx, dz, xp ;
  
  thb = 0.5*(th1 + th0) ; dth = 0.5*(th1 - th0) ;
  dth = fabs(dth) ;
  
  if  ( dth < TOL_THETA ) return 0 ;

  /*intermediate quantities for geometric calculations*/
  tq_vector_init(dxi, &(x[3*i0]), &(x[3*i1])) ;
  tq_vector_init(dxj, &(x[3*j0]), &(x[3*j1])) ;
  
  for ( i = 0 ; i < nth ; i ++ ) {
    th = thb + dth*qth[i*qtstr] ;
    Cth = cos(th) ; Sth = sqrt(1.0 - Cth*Cth) ;
    u = segment_intersect_2d(&(x[3*i0]), dxi, Cth, Sth) ;
    segment_interp(&(x[3*i0]), dxi, u, xi) ;
    ph0 = acos(xi[2]/sqrt(tq_vector_scalar(xi,xi))) ;

    u = segment_intersect_2d(&(x[3*j0]), dxj, Cth, Sth) ;
    segment_interp(&(x[3*j0]), dxj, u, xj) ;
    ph1 = acos(xj[2]/sqrt(tq_vector_scalar(xj,xj))) ;

    phb = 0.5*(ph1 + ph0) ; dph = 0.5*(ph1 - ph0) ;
    dph = fabs(dph) ;

    dx = xj[0] - xi[0] ; dz = xj[2] - xi[2] ; 

    xp = sqrt(tq_vector_scalar(xj,xj)*tq_vector_scalar(xi,xi))*
      sin(2.0*dph)*Cth ;
    
    for ( j = 0 ; j < nph ; j ++ ) {
      ph = phb + dph*qph[j*qpstr] ;
      Cph = cos(ph) ; Sph = sqrt(1.0 - Cph*Cph) ;
      s[0] = Sph*Cth ;
      s[1] = Sph*Sth ;
      s[2] = Cph ;
      r2 = xp/(dz*s[0] - dx*s[2]) ;
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
	/* wt = r*r*wth[i*wtstr]*wph[j*wpstr]*wr[k*wrstr]*dtp*dth*dr ; */
	qfunc(yg, sg, r, wt, q, nq, qdata) ;
      }
    }    
  }

  return 0 ;
}
#else /*HAVE_AVX_INSTRUCTIONS*/

/* static inline gboolean triangle_det4(gdouble dir[12], */
/* 				     __m256d rdvx, __m256d rdvy, __m256d rdvz,  */
/* 				     __m256d re2q, */
/* 				     gdouble *t) */

/* { */
/*   __m256d op1, acc, rdx, rdy, rdz ; */

/*   rdx = _mm256_load_pd(&(dir[0])) ; */
/*   rdy = _mm256_load_pd(&(dir[4])) ; */
/*   rdz = _mm256_load_pd(&(dir[8])) ; */
  
/*   acc = _mm256_mul_pd(rdx, rdvx) ; */
/*   op1 = _mm256_mul_pd(rdy, rdvy) ; */
/*   acc = _mm256_add_pd(acc, op1) ; */
/*   op1 = _mm256_mul_pd(rdz, rdvz) ; */
/*   acc = _mm256_add_pd(acc, op1) ; */
/*   acc = _mm256_div_pd(re2q, acc) ; */
  
/*   _mm256_store_pd(t, acc) ; */

/*   return TRUE ; */
/* } */

static __m256d calc_cos4(gdouble thb, gdouble dth, gdouble *q)

{
  gint i ;
  __m256d rCth ;
  __attribute__ ((aligned (32))) gdouble Cth[4] ;

  /*this can be done with AVX 512 when I get to a machine that lets
    me*/
  
  for ( i = 0 ; i < 4 ; i ++ ) {
    Cth[i] = cos(thb + dth*q[i]) ;
  }

  rCth = _mm256_load_pd(Cth) ;

  return rCth ;
}

static __m256d calc_sin4(__m256d rCth)

{
  __m256d rSth, op1 ;
  const __m256d op2 = _mm256_set1_pd(1.0) ;

  op1 = _mm256_mul_pd(rCth,rCth) ;
  op1 = _mm256_sub_pd(op2,op1) ;
  rSth = _mm256_sqrt_pd(op1) ;
  
  return rSth ;
}

static void segment_intersect_2d4(gdouble *x1, gdouble *d,
				  __m256d rCth, __m256d rSth, 
				  gdouble u[4])

/*
 * calculate intersection point in two dimensions of ray through
 * origin with direction (C,S) and line segment x1--(x1+d)
 *
 * this could probably be rewritten to return u in a __m256d
 */
  
{
  __m256d op1, op2, op3 ;

  op1 = _mm256_set1_pd(x1[0]) ;
  op1 = _mm256_mul_pd(op1, rSth) ;
  op2 = _mm256_set1_pd(d[0]) ;
  op2 = _mm256_mul_pd(op2, rSth) ;
  op3 = _mm256_set1_pd(d[1]) ;
  op3 = _mm256_mul_pd(op3, rCth) ;
  /*denominator*/
  op3 = _mm256_sub_pd(op3, op2) ;
  /*numerator*/
  op2 = _mm256_set1_pd(x1[1]) ;
  op2 = _mm256_mul_pd(op2, rCth) ;
  op2 = _mm256_sub_pd(op1, op2) ;
  op2 = _mm256_div_pd(op2, op3) ;
  _mm256_store_pd(u, op2) ;

  return ;
}

static void segment_interp4(gdouble *x1, gdouble *d, gdouble u[4],
			    gdouble tp[4], gdouble r[4])

/*
 * interpolate on line segment x1--(x1+d), 0 <= u <= 1
 */
  
{
  gint i ;
  gdouble x[3], z[4] ;
  __m256d op1, op2 ;
  
  for ( i = 0 ; i < 4 ; i ++ ) {
    x[0] = x1[0] + u[i]*d[0] ;
    x[1] = x1[1] + u[i]*d[1] ;
    x[2] = x1[2] + u[i]*d[2] ;
    z[i] = x[2]  ;
    r[i] = tq_vector_scalar(x,x) ;
  }

  op1 = _mm256_load_pd(r) ;
  op1 = _mm256_sqrt_pd(op1) ;
  _mm256_store_pd(r, op1) ;
  op2 = _mm256_load_pd(z) ;
  op1 = _mm256_div_pd(op2, op1) ;
  _mm256_store_pd(tp, op1) ;

  return ;
}

static void intersect4(__m256d rxp, __m256d ra, __m256d rb,
		       gdouble S[4], gdouble C[4], gdouble r2[4])

{
  __m256d rS, rC ;

  rS = _mm256_load_pd(S) ;
  rS = _mm256_mul_pd(rS, ra) ;
  rC = _mm256_load_pd(C) ;  
  rC = _mm256_mul_pd(rC, rb) ;
  rS = _mm256_sub_pd(rS, rC) ;
  rC = _mm256_div_pd(rxp, rS) ;
  _mm256_store_pd(r2, rC) ;

  return ;
}

static gint subtet_quad_th_avx(gdouble *x, gdouble *xs,
			       gdouble th0, gdouble th1,
			       gint i0, gint i1,
			       gint j0, gint j1,
			       gdouble *qph, gint qpstr,
			       gdouble *wph, gint wpstr,
			       gint nph,
			       gdouble *qth, gint qtstr,
			       gdouble *wth, gint wtstr,
			       gint nth,
			       gdouble *qr, gint qrstr,
			       gdouble *wr, gint wrstr,
			       gint nr,
			       tq_tetquad_func_t qfunc, gpointer qdata,
			       gdouble *q, gint nq, gdouble *A)

{
  gint i, j, k, m, n ;
#ifdef TRANSFORM_PHI
  gdouble dtp, tpb ;
#else /*TRANSFORM_PHI*/
  gdouble ph0[4], ph1[4], phb, dph, ph ;
#endif /*TRANSFORM_PHI*/
  gdouble dth, thb ;
  gdouble Cph[4], Sph[4] ;
  __attribute__ ((aligned (32))) gdouble r2[4], u[4], Cth[4], Sth[4], s[12] ;
  __attribute__ ((aligned (32))) gdouble tp0[4], tp1[4] ;
  gdouble r, rp1[4], rp0[4] ;
  gdouble wt, wj, yg[3], sg[12], dxi[3], dxj[3] ;
  /* gdouble edge1[3], edge2[3], qvec[3], e2q, dv[3] ; */
  /* __m256d rdvx, rdvy, rdvz, re2q ; */
  __m256d ra, rb, rxp, rCth, rSth ;
  gint i3 = 3, i4 = 4 ;
  gdouble d1 = 1.0, d0 = 0.0 ;	
  
  g_assert(qpstr == 1) ; g_assert(wpstr == 1) ;
  g_assert(qtstr == 1) ; g_assert(wtstr == 1) ;
  g_assert(qrstr == 1) ; g_assert(wrstr == 1) ;

  g_assert(nth % 4 == 0) ;
  g_assert(nph % 4 == 0) ;
  
  thb = 0.5*(th1 + th0) ; dth = 0.5*(th1 - th0) ;
  dth = fabs(dth) ;
  
  if  ( dth < TOL_THETA ) return 0 ;

  /*intermediate quantities for geometric calculations*/
  tq_vector_init(dxi, &(x[3*i0]), &(x[3*i1])) ;
  tq_vector_init(dxj, &(x[3*j0]), &(x[3*j1])) ;
  /* tq_vector_init(edge1, &(x[3*1]), &(x[3*0])) ; */
  /* tq_vector_init(edge2, &(x[3*2]), &(x[3*0])) ; */
  /* tq_vector_cross(qvec, &(x[3*0]), edge1) ; */
  /* e2q = tq_vector_scalar(edge2, qvec) ;   */

  /* re2q = _mm256_set1_pd(e2q) ; */
  
  /*
   * precompute terms used in vector triple product
   * a.(b x c) = b.dv (with b == direction vector)
   */
  /* dv[0] = edge1[2]*edge2[1] - edge1[1]*edge2[2] ;  */
  /* dv[1] = edge1[0]*edge2[2] - edge1[2]*edge2[0] ;  */
  /* dv[2] = edge1[1]*edge2[0] - edge1[0]*edge2[1] ;  */
  /* rdvx = _mm256_set1_pd(dv[0]) ; */
  /* rdvy = _mm256_set1_pd(dv[1]) ; */
  /* rdvz = _mm256_set1_pd(dv[2]) ; */

  for ( i = 0 ; i < nth ; i += 4 ) {
    rCth = calc_cos4(thb, dth, &(qth[i])) ;
    rSth = calc_sin4(rCth) ;
    segment_intersect_2d4(&(x[3*i0]), dxi, rCth, rSth, u) ;
    segment_interp4(&(x[3*i0]), dxi, u, tp0, rp0) ;
#ifndef TRANSFORM_PHI
    for ( m = 0 ; m < 4 ; m ++ ) ph0[m] = acos(tp0[m]) ;
#endif /*TRANSFORM_PHI*/
    
    segment_intersect_2d4(&(x[3*j0]), dxj, rCth, rSth, u) ;
    segment_interp4(&(x[3*j0]), dxj, u, tp1, rp1) ;
#ifndef TRANSFORM_PHI
    for ( m = 0 ; m < 4 ; m ++ ) ph1[m] = acos(tp1[m]) ;
#endif /*TRANSFORM_PHI*/

    _mm256_store_pd(Sth, rSth) ;
    _mm256_store_pd(Cth, rCth) ;
  
    for ( n = 0 ; n < 4 ; n ++ ) {
#ifdef TRANSFORM_PHI
      tpb = 0.5*(tp1[n] + tp0[n]) ; dtp = 0.5*(tp1[n] - tp0[n]) ;
      dtp = fabs(dtp) ;
#else /*TRANSFORM_PHI*/
      phb = 0.5*(ph1[n] + ph0[n]) ; dph = 0.5*(ph1[n] - ph0[n]) ;
      dph = fabs(dph) ;
#endif /*TRANSFORM_PHI*/

      r = rp0[n]*rp1[n]*(tp1[n]*sqrt(1.0-tp0[n]*tp0[n]) -
			 tp0[n]*sqrt(1.0-tp1[n]*tp1[n])) ;
      rxp = _mm256_set1_pd(r) ;
      r = rp1[n]*tp1[n] - rp0[n]*tp0[n] ;
      ra = _mm256_set1_pd(r) ;
      r = rp1[n]*sqrt(1.0-tp1[n]*tp1[n]) - rp0[n]*sqrt(1.0-tp0[n]*tp0[n]) ;
      rb = _mm256_set1_pd(r) ;
      
      for ( j = 0 ; j < nph ; j += 4 ) {
	for ( m = 0 ; m < 4 ; m ++ ) {
#ifdef TRANSFORM_PHI
	  Cph[m] = tpb + dtp*qph[j+m] ;
#else /*TRANSFORM_PHI*/
	  ph = phb + dph*qph[j+m] ; Cph[m] = cos(ph) ;
#endif /*TRANSFORM_PHI*/
	  Sph[m] = sqrt(1.0 - Cph[m]*Cph[m]) ;
	  s[4*0+m] = Sph[m]*Cth[n] ;
	  s[4*1+m] = Sph[m]*Sth[n] ;
	  s[4*2+m] = Cph[m] ;
	}

	intersect4(rxp, ra, rb, Sph, Cph, r2) ;
	
	blaswrap_dgemm(TRUE, FALSE, i4, i3, i3, d1, s, i4, A, i3, d0, sg, i3) ;
      
	for ( m = 0 ; m < 4 ; m ++ ) {
	  /* r2[m] = fabs(r2[m]) ; */
#ifdef TRANSFORM_PHI
	  wj = wth[i+n]*wph[j+m]*dtp*dth*0.5*r2[m] ;
#else /*TRANSFORM_PHI*/
	  wj = Sph[m]*wth[i+n]*wph[j+m]*dph*dth*0.5*r2[m] ;
#endif /*TRANSFORM_PHI*/
	  
	  for ( k = 0 ; k < nr ; k ++ ) {
	    r = 0.5*r2[m]*(1.0 + qr[k]) ;
	    /*global coordinate system*/
	    yg[0] = r*sg[3*m+0] + xs[0] ; 
	    yg[1] = r*sg[3*m+1] + xs[1] ; 
	    yg[2] = r*sg[3*m+2] + xs[2] ; 
	    wt = r*r*wr[k]*wj ;
	    qfunc(yg, &(sg[3*m]), r, wt, q, nq, qdata) ;
	  }
	}      
      }
    }
  }
  
  return 0 ;
}
#endif /*HAVE_AVX_INSTRUCTIONS*/

static void zero_points(gdouble *x, gdouble tol, gint n)

{
  gint i ;

  for ( i = 0 ; i < n ; i ++ ) if ( fabs(x[i]) < tol ) x[i] = 0.0 ;
  
  return ;
}

static gint tet_quad_origin_th(gdouble *x, gdouble *xs,
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

  zero_points(xr, TOL_ZERO, 9) ;
  
  order_points(xr, perm) ;
  
  /*find \theta limits and corresponding line segments*/
  for ( i = 0 ; i < 3 ; i ++ ) {
    th[i] = atan2(xr[3*i+1], xr[3*i+0]) ;
  }

  if ( th[0] < 0 || th[1] < 0 || th[2] < 0 ) {
    /* fprintf(stderr, */
    /* 	    "%1.32e %1.32e %1.32e\n" */
    /* 	    "%1.32e %1.32e %1.32e\n" */
    /* 	    "%1.32e %1.32e %1.32e\n" */
    /* 	    "%1.32e %1.32e %1.32e\n", */
    /* 	    x[3*0+0], x[3*0+1], x[3*0+2],  */
    /* 	    x[3*1+0], x[3*1+1], x[3*1+2],  */
    /* 	    x[3*2+0], x[3*2+1], x[3*2+2],  */
    /* 	    x[3*3+0], x[3*3+1], x[3*3+2]) ; */
    /* g_error("%s: integration failure", __FUNCTION__) ; */
    return 3 ;
  }
  
  /* g_assert(th[0] >= 0) ; */
  /* g_assert(th[1] >= 0) ; */
  /* g_assert(th[2] >= 0) ; */
  
  if ( BETWEEN(th[perm[2]], th[perm[0]], th[perm[1]]) ) {
    thlim[0] = th[perm[0]] ; thlim[1] = th[perm[2]] ;
    bounds[0] = perm[0] ; bounds[1] = perm[1] ; 
    bounds[2] = perm[0] ; bounds[3] = perm[2] ; 

    thlim[2] = th[perm[2]] ; thlim[3] = th[perm[1]] ;
    bounds[4] = perm[0] ; bounds[5] = perm[1] ; 
    bounds[6] = perm[1] ; bounds[7] = perm[2] ; 
  } else {
    thlim[0] = th[perm[0]] ; thlim[1] = th[perm[1]] ;
    bounds[0] = perm[0] ; bounds[1] = perm[1] ; 
    bounds[2] = perm[0] ; bounds[3] = perm[2] ; 

    thlim[2] = th[perm[1]] ; thlim[3] = th[perm[2]] ;
    bounds[4] = perm[1] ; bounds[5] = perm[2] ; 
    bounds[6] = perm[0] ; bounds[7] = perm[2] ; 
  }

  /* thlim[0] = 0 ; */
  /* thlim[1] = MAX(th[0], MAX(th[1], th[2])) ; */
  /* thlim[2] = 0 ; */
  /* thlim[3] = 0 ; */
  
#ifdef HAVE_AVX_INSTRUCTIONS
  
  i = subtet_quad_th_avx(xr, xs, thlim[0], thlim[1],
			 bounds[0], bounds[1], bounds[2], bounds[3],
			 qph, qpstr, wph, wpstr, nph,
			 qth, qtstr, wth, wtstr, nth,
			 qr,  qrstr, wr,  wrstr, nr,
			 qfunc, qdata,
			 q, nq, A) ;
  if ( i != 0 ) return 1 ;

  i = subtet_quad_th_avx(xr, xs, thlim[2], thlim[3],
  			 bounds[4], bounds[5], bounds[6], bounds[7],
  			 qph, qpstr, wph, wpstr, nph,
  			 qth, qtstr, wth, wtstr, nth,
  			 qr,  qrstr, wr,  wrstr, nr,
  			 qfunc, qdata,
  			 q, nq, A) ;
  if ( i != 0 ) return 2 ;
#else /*HAVE_AVX_INSTRUCTIONS*/
  i = subtet_quad_th(xr, xs, thlim[0], thlim[1],
  		     bounds[0], bounds[1], bounds[2], bounds[3],
  		     qph, qpstr, wph, wpstr, nph,
  		     qth, qtstr, wth, wtstr, nth,
  		     qr,  qrstr, wr,  wrstr, nr,
  		     qfunc, qdata,
  		     q, nq, A) ;
  if ( i != 0 ) return 1 ;

  i = subtet_quad_th(xr, xs, thlim[2], thlim[3],
  		     bounds[4], bounds[5], bounds[6], bounds[7],
  		     qph, qpstr, wph, wpstr, nph,
  		     qth, qtstr, wth, wtstr, nth,
  		     qr,  qrstr, wr,  wrstr, nr,
  		     qfunc, qdata,
  		     q, nq, A) ;
  if ( i != 0 ) return 2 ;
#endif /*HAVE_AVX_INSTRUCTIONS*/
  
  return 0 ;
}

static gdouble select_rotation(gdouble *xr, gdouble *th, gint j, gint k)

{
  gdouble thr, xb[2], dth ;

  /*check quadrants*/
  if ( xr[3*j+2] >= 0 && xr[3*k+2] >= 0 ) {
    thr = -th[j] + 0.5*M_PI ;
    g_assert(th[k] + thr > -0.5*M_PI) ;
    return thr ;
  }

  if ( xr[3*j+2] >= 0 && xr[3*k+2] < 0 ) {
    thr = -0.5*M_PI - th[j] ;
    dth = th[k] + thr ;
    if ( dth < -0.5*M_PI ) dth += 2.0*M_PI ;
    if ( dth >  0.5*M_PI+TOL_ZERO ||
	 dth < -0.5*M_PI-TOL_ZERO ) {

      thr += M_PI ;
    }
    return thr ;
  }

  if ( xr[3*j+2] < 0 && xr[3*k+2] >= 0 ) {
    thr = -th[j] - 0.5*M_PI ;
    if ( th[k] + thr >  0.5*M_PI ||
	 th[k] + thr < -0.5*M_PI ) {
      thr += M_PI ;
    }
    
    return thr ;
  }

  if ( xr[3*j+2] < 0 && xr[3*k+2] < 0 ) {
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

static void cart2sph(gdouble *x, gdouble *r, gdouble *th, gdouble *ph)

{
  *r = sqrt(tq_vector_length2(x)) ;
  *th = atan2(x[1], x[0]) ;
  *ph = acos(x[2]/(*r)) ;
	     
  return ;
}

static gint matrix_define(gdouble *x, gdouble *A, gint *p)

{
  gdouble xr[9], ph[3], phr, dih[3] ;
  gint i, j, k, idih[3] ;

  tq_tet_dihedral_angles(&(x[3*0]), &(x[3*1]), &(x[3*2]),
  			 &(dih[0]), &(dih[1]), &(dih[2]),
  			 NULL, NULL, NULL) ;

  /* idih[0] = 0 ; idih[1] = 1 ; idih[2] = 2 ;  */
  idih[0] = p[0] ; idih[1] = p[1] ; idih[2] = p[2] ; 
  
  /* if ( index_sort3(dih, idih) != 0 ) { */
  /*   fprintf(stderr, "dihedral error:") ; */
  /*   fprintf(stderr, */
  /* 	    "%lg %lg %lg\n" */
  /* 	    "%lg %lg %lg\n" */
  /* 	    "%lg %lg %lg\n", */
  /* 	    x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8]) ; */
  /* } */

  /*align the edge with largest dihedral with the x axis*/
  i = idih[2] ; j = idih[0] ; k = idih[1] ;
  /*rotate to bring node to z==0*/
  rotation_y(atan2(x[3*i+2], x[3*i+0]), A) ;
  
  multiply_matrices3x3(xr, x, A, TRUE) ;

  /*rotate to bring node to y==0*/
  rotation_z(-atan2(xr[3*i+1], xr[3*i+0]), A) ;
  
  multiply_matrices3x3(xr, x, A, TRUE) ;

/* #if 0 */
  /*rotation in x to bring nodes 3 and 4 to same azimuth*/
  gdouble r3, th3, ph3, r4, th4, ph4 ;
  cart2sph(&(xr[3*j]), &r3, &th3, &ph3) ;
  cart2sph(&(xr[3*k]), &r4, &th4, &ph4) ;

  /* fprintf(stderr, "th3=%lg; ph3=%lg\n", th3, ph3) ; */
  /* fprintf(stderr, "th4=%lg; ph4=%lg\n", th4, ph4) ; */

  if ( fabs(th3 - 0.5*M_PI) < 1e-9 || fabs(th4 - 0.5*M_PI) < 1e-9 ) {
    rotation_x(-0.5*M_PI, A) ;
    multiply_matrices3x3(xr, x, A, TRUE) ;
    cart2sph(&(xr[3*j]), &r3, &th3, &ph3) ;
    cart2sph(&(xr[3*k]), &r4, &th4, &ph4) ;

    if ( fabs(th3-th4) < 1e-9 ) return 0 ;
  }

  if ( fabs(th3 + 0.5*M_PI) < 1e-9 || fabs(th4 + 0.5*M_PI) < 1e-9 ) {
    rotation_x( 0.5*M_PI, A) ;
    multiply_matrices3x3(xr, x, A, TRUE) ;
    cart2sph(&(xr[3*j]), &r3, &th3, &ph3) ;
    cart2sph(&(xr[3*k]), &r4, &th4, &ph4) ;

    if ( fabs(th3-th4) < 1e-9 ) return 0 ;
  }

  if ( fabs(ph3) < 1e-9 || fabs(ph4) < 1e-9 ) {
    if ( th3 > -1e-9 && th4 > -1e-9 ) return 0 ;
    rotation_x(M_PI, A) ;
    return 0 ;
  }
  
  /* phr = (sin(th3)*cos(th4) - cos(th3)*sin(th4))*sin(ph3)*sin(ph4) ; */
  /* phr /= cos(ph3)*sin(ph4)*cos(th4) - sin(ph3)*cos(ph4)*cos(th3) ; */
  /* phr = atan(phr) ; */
  phr = atan2((sin(th3)*cos(th4) - cos(th3)*sin(th4))*sin(ph3)*sin(ph4),
  	      cos(ph3)*sin(ph4)*cos(th4) - sin(ph3)*cos(ph4)*cos(th3)) ;

  /* rr = cos(phr)*sin(ph3)*sin(th3) - sin(phr)*cos(ph3) ; */
  /* rr /= sin(ph3)*cos(th3) ; */
  /* if ( rr < 0.0 ) phr += M_PI ; */

  r3 = atan2(cos(phr)*sin(ph3)*sin(th3) - sin(phr)*cos(ph3),
	     sin(ph3)*cos(th3)) ;
  /* if ( rr < 0 ) phr += M_PI ; */
  
  /* rr = cos(phr)*sin(ph4)*sin(th4) - sin(phr)*cos(ph4) ; */
  /* rr /= sin(ph4)*cos(th4) ; */
  
  r4 = atan2(cos(phr)*sin(ph4)*sin(th4) - sin(phr)*cos(ph4),
	     sin(ph4)*cos(th4)) ;

  if ( fabs(r3 - r4) > 1e-9 ) return 1 ;

  if ( r3 < 0 || r4 < 0 ) phr += M_PI ;
  
  /* if ( rr < 0 ) phr += M_PI ; */
  rotation_x(phr, A) ;
  
  multiply_matrices3x3(xr, x, A, TRUE) ;
  
  return 0 ;
/* #endif */
  
  if ( xr[3*0+1] > -TOL_ZERO &&
       xr[3*1+1] > -TOL_ZERO &&
       xr[3*2+1] > -TOL_ZERO ) return 0 ;

  /* fprintf(stderr, "negative y (%lg, %lg, %lg)\n", ph[0], ph[1], ph[2]) ; */

  if ( xr[3*j+1] < -TOL_ZERO && xr[3*k+1] < -TOL_ZERO ) {
    /*both values of y negative, easier to make a 180 or a reflection*/
    /* fprintf(stderr, "flipping\n") ; */
    rotation_x(M_PI, A) ;

    return 0 ;
  }

  /* check for points with negative y */
  ph[0] = atan2(xr[3*0+2], xr[3*0+1]) ;
  ph[1] = atan2(xr[3*1+2], xr[3*1+1]) ;
  ph[2] = atan2(xr[3*2+2], xr[3*2+1]) ;

  if ( xr[3*k+1] < -TOL_ZERO ) {
    phr = select_rotation(xr, ph, k, j) ;
  } else {
    phr = select_rotation(xr, ph, j, k) ;
  }
  
  /*rotate to bring node to x==0*/
  rotation_x(phr, A) ;
  
  multiply_matrices3x3(xr, x, A, TRUE) ;

  return 0 ;
}

static gint transform_nodes(gdouble *x, gdouble y[9], gdouble *A)

{
  gdouble p[3], z, n[3], s[3], t[3], Cth, Sth ;
  gdouble pd[3], zd, nd[3], sd[3], td[3] ;
  gdouble x1, x2, x3, len, ds, dt, d23 ;
  gint i ;

  A[0] = 1.0 ; A[1] = 0.0 ; A[2] = 0.0 ; 
  A[3] = 0.0 ; A[4] = 1.0 ; A[5] = 0.0 ; 
  A[6] = 0.0 ; A[7] = 0.0 ; A[8] = 1.0 ;
  /*check that tet is not already in required form*/
  if ( fabs(x[3*0+1]) < TOL_ZERO && fabs(x[3*0+2]) < TOL_ZERO ) {
    if ( fabs(atan2(x[3*1+1], x[3*1+0]) - atan2(x[3*2+1], x[3*2+0])) < 1e-9 ) {
      memcpy(y, x, 9*sizeof(gdouble)) ;
      
      return 0 ;
    }
    if ( fabs(atan2(x[3*1+2], x[3*1+1]) - atan2(x[3*2+2], x[3*2+1])) < 1e-9 ) {
      rotation_x(0.5*M_PI, A) ;      
      multiply_matrices3x3(y, x, A, TRUE) ;

      return 0 ;
    }
  }
  
  x1 = sqrt(tq_vector_length2(&(x[0]))) ;
  x2 = sqrt(tq_vector_length2(&(x[3]))) ;
  x3 = sqrt(tq_vector_length2(&(x[6]))) ;

  /*normal to plane 023*/
  tq_vector_cross(n, &(x[3]), &(x[6])) ;
  len = sqrt(tq_vector_length2(n)) ;
  n[0] /= len ; n[1] /= len ; n[2] /= len ;
  /*projection of point 1 on plane 023*/
  z = tq_vector_scalar(&(x[0]), n) ;
  p[0] = x[0] - z*n[0] ; p[1] = x[1] - z*n[1] ; p[2] = x[2] - z*n[2] ;

  /*check for p lying too close to origin*/
  /* if ( (len = tq_vector_length2(p)) < 1e-12 ) return 1 ; */
  len = tq_vector_length2(p) ;
  /*vectors for in-plane coordinate system*/
  if ( len > 1e-12 ) {
    len = sqrt(len) ;
    s[0] = p[0]/len ; s[1] = p[1]/len ; s[2] = p[2]/len ;
  } else {
    if ( fabs(n[0]) > 1e-12 || fabs(n[1]) > 1e-12 ) {
      s[0] = -n[1] ; s[1] = n[0] ; s[2] = 0.0 ;
      len = sqrt(tq_vector_length2(s)) ;
      s[0] /= len ; s[1] /= len ; s[2] /= len ; 
    } else {
      s[0] = 0.0 ; s[1] = 0.0 ; s[2] = 1.0 ;
    }
  }
  tq_vector_cross(t, n, s) ;
  /*angle between edge 1 and plane 023*/
  Cth = tq_vector_scalar(&(x[0]), s)/x1 ;
  Sth = sqrt(1.0 - Cth*Cth) ;

  /*vectors for rotated coordinate system*/
  nd[0] = -Sth ; nd[1] = Cth ; nd[2] = 0.0 ;
  zd = -x1*Sth ;
  pd[0] = x1 - zd*nd[0] ; pd[1] = -zd*nd[1] ; pd[2] = -zd*nd[2] ;
  /*this should be the same as the length of p*/
  len = sqrt(tq_vector_length2(pd)) ;
  if ( len > 1e-6 ) {
    sd[0] = pd[0]/len ; sd[1] = pd[1]/len ; sd[2] = pd[2]/len ;
  } else {
    if ( fabs(nd[0]) > 1e-12 || fabs(nd[1]) > 1e-12 ) {
      sd[0] = -nd[1] ; sd[1] = nd[0] ; sd[2] = 0.0 ;
      len = sqrt(tq_vector_length2(sd)) ;
      sd[0] /= len ; sd[1] /= len ; sd[2] /= len ; 
    } else {
      sd[0] = 0.0 ; sd[1] = 0.0 ; sd[2] = 1.0 ;
    }
    /* sd[0] = -nd[1] ; sd[1] = nd[0] ; sd[2] = 0.0 ; */
  }
  tq_vector_cross(td, nd, sd) ;

  /*coordinates of rotated tet*/
  y[3*0+0] = x1 ; y[3*0+1] = 0.0 ; y[3*0+2] = 0.0 ;
  for ( i = 1 ; i < 3 ; i ++ ) {
    ds = tq_vector_diff_scalar(s, p, &(x[3*i])) ;
    dt = tq_vector_diff_scalar(t, p, &(x[3*i])) ;
    y[3*i+0] = pd[0] + ds*sd[0] + dt*td[0] ; 
    y[3*i+1] = pd[1] + ds*sd[1] + dt*td[1] ; 
    y[3*i+2] = pd[2] + ds*sd[2] + dt*td[2] ; 
  }

  /*solve for rotation matrix*/
  A[3*0+0] = x[3*0+0]/x1 ; 
  A[3*0+1] = x[3*0+1]/x1 ; 
  A[3*0+2] = x[3*0+2]/x1 ; 

  d23 = y[3*1+2]*y[3*2+1] - y[3*1+1]*y[3*2+2] ;
  for ( i = 0 ; i < 3 ; i ++ ) {
    A[3*2+i] =
      (x[3*1+i] - A[3*0+i]*y[3*1+0])*y[3*2+1] -
      (x[3*2+i] - A[3*0+i]*y[3*2+0])*y[3*1+1] ;
    A[3*2+i] /= d23 ;
    if ( fabs(y[3*2+1]) > 1e-9 ) 
      A[3*1+i] = (x[3*2+i] - A[3*0+i]*y[3*2+0] - A[3*2+i]*y[3*2+2])/y[3*2+1] ;
    else
      A[3*1+i] = (x[3*1+i] - A[3*0+i]*y[3*1+0] - A[3*2+i]*y[3*1+2])/y[3*1+1] ;
  }

#if 0
  /*internal check on rotation*/
  gdouble xc[9], err = 0.0 ;
  /*transpose A for rotation from global to transformed coordinates*/
  multiply_matrices3x3(xc, x, A, TRUE) ;

  for ( i = 0 ; i < 9 ; i ++ ) err = MAX(err, fabs(xc[i]-y[i])) ;

  fprintf(stderr, "forward rotation error: %lg\n", err) ;
#endif
  
  return 0 ;
}

gint tq_transform_matrix(gdouble *x, gdouble *A)

{
  gdouble y[9], xt[9] ;
  gint p[] = {0, 1, 2, 1, 2, 0, 2, 1, 0, 1, 0, 2, 0, 2, 1, 2, 1, 0} ;
  gint i, j ;

  for ( i = 0 ; i < 6 ; i ++ ) {
    memcpy(&(xt[3*0]), &(x[3*p[3*i+0]]), 3*sizeof(gdouble)) ;
    memcpy(&(xt[3*1]), &(x[3*p[3*i+1]]), 3*sizeof(gdouble)) ;
    memcpy(&(xt[3*2]), &(x[3*p[3*i+2]]), 3*sizeof(gdouble)) ;
    j = transform_nodes(xt, y, A) ;

    multiply_matrices3x3(y, xt, A, TRUE) ;

    zero_points(y, TOL_ZERO, 9) ;

    if ( j == 0 ) {
      /*check for sign of y*/
      if ( y[3*1+1] >= 0.0 && y[3*2+1] >= 0.0 ) return 0 ;

      /*check for both values negative, easiest to flip about xz*/
      if ( (y[3*1+1] <= 0.0 && y[3*2+1] <= 0.0 ) ) {
	/* rotation_x(M_PI, A) ; */
	A[3*1+0] *= -1 ; A[3*1+1] *= -1 ; A[3*1+2] *= -1 ;
	/* A[3*0+1] *= -1 ; A[3*1+1] *= -1 ; A[3*2+1] *= -1 ; */
	multiply_matrices3x3(y, xt, A, TRUE) ;
	
	return 0 ;
      }
    }
  }

  fprintf(stderr,
	  "%1.32e %1.32e %1.32e\n"
	  "%1.32e %1.32e %1.32e\n"
	  "%1.32e %1.32e %1.32e\n",
	  xt[3*0+0], xt[3*0+1], xt[3*0+2], 
	  xt[3*1+0], xt[3*1+1], xt[3*1+2], 
	  xt[3*2+0], xt[3*2+1], xt[3*2+2]) ;
  
  g_assert_not_reached() ;
  
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
	    "%1.32e %1.32e %1.32e\n"
	    "%1.32e %1.32e %1.32e\n"
	    "%1.32e %1.32e %1.32e\n"
	    "%1.32e %1.32e %1.32e\n",
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
