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

#include <sqt.h>

#include <blaswrap.h>

static gint tri_quad_monomials(gdouble s, gdouble t, gdouble w,
			       gdouble *y, gdouble *n,
			       gdouble *K, gint nk,
			       gdouble *quad, gint nc,
			       gint init,
			       gpointer data)
{
  gdouble R3, pw[48] ;
  gint N = *((gint *)data) ;
  gint i, j, k, idx, m ;

  g_assert(N < 16) ;
  
  R3 = y[0]*y[0] + y[1]*y[1] + y[2]*y[2] ;
  R3 *= sqrt(R3)*4.0*M_PI ;
  R3 = w/R3 ;
  
  pw[0] = pw[1] = pw[2] = 1.0 ;
  for ( i = 1 ; i <= N ; i ++ ) {
    pw[3*i+0] = pw[3*(i-1)+0]*y[0] ;
    pw[3*i+1] = pw[3*(i-1)+1]*y[1] ;
    pw[3*i+2] = pw[3*(i-1)+2]*y[2] ;
  }

  /* quad[0] = 0.0 ; */
  quad[0] += w ;
  for ( m = 1 ; m <= N ; m ++ ) {
    for ( i = 0 ; i <= m ; i ++ ) {
      for ( j = 0 ; j <= m-i ; j ++ ) {
	k = m - i - j ;
	idx = bsi_poly_coefficient_3d_index(i, j, k) ;
	quad[idx] += pw[3*i+0]*pw[3*j+1]*pw[3*k+2]*R3 ;
      }
    }
  }

  return 0 ;
}

gint bsi_tet_integrate_powers(gdouble *x1, gdouble *x2, gdouble *x3,
			      gdouble *x4,
			      gint N,
			      gint nq, gdouble tol, gint dmax,
			      gdouble *quad)

/*
  volume integration of (x-x_1)^n/R^3 n = 1, ..., N

  quad[0] is zeroed to help the indexing, extract the integral for
  (x^i)(y^j)(z^k) with

  bsi_poly_coefficient_3d_index(i, j, k) ;
*/
  
{
  gdouble *q, xt[9], nn[3], s ;
  gint order, xstr, i, j, k, n, idx, np ;
  sqt_quadrature_func_t func = (sqt_quadrature_func_t)tri_quad_monomials ;
  
  sqt_quadrature_select(nq, &q, &order) ;

  /*opposite face in shifted coordinates*/
  xt[3*0+0] = x2[0]-x1[0] ; xt[3*0+1] = x2[1]-x1[1] ; xt[3*0+2] = x2[2]-x1[2] ;
  xt[3*1+0] = x3[0]-x1[0] ; xt[3*1+1] = x3[1]-x1[1] ; xt[3*1+2] = x3[2]-x1[2] ;
  xt[3*2+0] = x4[0]-x1[0] ; xt[3*2+1] = x4[1]-x1[1] ; xt[3*2+2] = x4[2]-x1[2] ;

  /*initialize quadrature*/
  np = (N+1)*(N+2)*(N+3)/6 ;
  memset(quad, 0, np*sizeof(gdouble)) ;

  xstr = 3 ;
  sqt_adaptive_quad_tri(xt, xstr, 3, q, nq, func, quad, np, tol, dmax, &N) ;

  /*distance to plane of xt*/
  xt[3*1+0] -= xt[3*0+0] ; xt[3*1+1] -= xt[3*0+1] ; xt[3*1+2] -= xt[3*0+2] ;
  xt[3*2+0] -= xt[3*0+0] ; xt[3*2+1] -= xt[3*0+1] ; xt[3*2+2] -= xt[3*0+2] ;

  gts_vector_cross(nn, &(xt[3]), &(xt[6])) ;
  gts_vector_normalize(nn) ;
  s = gts_vector_scalar(xt, nn) ;
  s = fabs(s) ;
  
  /* return 0 ; */
  
  /* /\*use sign of quadratic terms to check sign of result*\/ */
  /* if ( quad[bsi_poly_coefficient_3d_index(2,0,0)] > 0.0 ) { */
  /*   s = -s ; */
  /* } */

  quad[0] *= s/3 ;
  
  for ( n = 1 ; n <= N ; n ++ ) {
    for ( i = 0 ; i <= n ; i ++ ) {
      for ( j = 0 ; j <= n-i ; j ++ ) {
	k = n - i -j ;
	idx = bsi_poly_coefficient_3d_index(i, j, k) ;
	quad[idx] *= s/n ;
      }
    }
  }

  g_assert(quad[bsi_poly_coefficient_3d_index(2,0,0)] >= 0.0) ;
  g_assert(quad[bsi_poly_coefficient_3d_index(0,2,0)] >= 0.0) ;
  g_assert(quad[bsi_poly_coefficient_3d_index(0,0,2)] >= 0.0) ;
  
  return 0 ;
}
