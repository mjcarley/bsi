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

#include <blaswrap.h>

#include "rbf.h"

extern gint dgesv_(gint *n, gint *nrhs, gdouble *A, gint *lda,
		   gint *ip, gdouble *b, gint *ldb, gint *info) ;

gdouble rbf_basis_func_r1(gdouble r)

{
  return r ;
}

gdouble rbf_basis_func_r3(gdouble r)

{
  return r*r*r ;
}

gdouble rbf_basis_func_tps(gdouble r)

{
  if ( r == 0.0 ) return 0.0 ;

  return r*r*log(r) ;
}

gdouble rbf_basis_func_q(gdouble r)

{
  return 1.0 + r*r ;
}

gdouble rbf_basis_func_mq(gdouble r)

{
  return sqrt(1.0 + r*r) ;
}

gdouble rbf_basis_func_imq(gdouble r)

{
  return 1.0/sqrt(1.0 + r*r) ;
}

gdouble rbf_basis_func_iq(gdouble r)

{
  return 1.0/(1.0 + r*r) ;
}

gdouble rbf_basis_func_g(gdouble r)

{
  return exp(-r*r) ;
}

gdouble rbf_basis_func_cp_c0(gdouble r)

{
  if ( r > 1.0 ) return 0.0 ;
  
  return 1.0 - r*r ;
}

gdouble rbf_basis_func_cp_c2(gdouble r)

{
  if ( r > 1.0 ) return 0.0 ;

  gdouble rm12 = 1.0 - r ;

  rm12 *= rm12 ;
  
  return rm12*rm12*(4.0*r + 1.0) ;
  
  /* return rm1*rm1*rm1*rm1*(4.0*r + 1.0) ; */
}

gdouble rbf_basis_func_cp_c4(gdouble r)

{
  if ( r > 1.0 ) return 0.0 ;

  gdouble rm12 = 1.0 - r ;

  rm12 *= rm12 ;
  
  return rm12*rm12*rm12*(r*(35.0/3.0*r + 6.0) + 1.0) ;
  
  /* return rm1*rm1*rm1*rm1*rm1*rm1*(r*(35.0/3.0*r + 6.0) + 1.0) ; */
}

gdouble rbf_basis_func_cp_c6(gdouble r)

{
  if ( r > 1.0 ) return 0.0 ;
  gdouble rm14 = 1.0 - r ;

  rm14 *= rm14 ;
  rm14 *= rm14 ;
  
  return rm14*rm14*(r*(r*(32.0*r + 25.0) + 8.0) + 1.0) ;
  /* return rm1*rm1*rm1*rm1*rm1*rm1*rm1*rm1*(r*(r*(32.0*r + 25.0) + 8.0) + 1.0) ; */
}

gdouble rbf_basis_func_ctps_c0(gdouble r)

{
  if ( r > 1.0 ) return 0.0 ;
  gdouble rm1 = 1.0 - r ;

  return rm1*rm1*rm1*rm1*rm1 ;
}

gdouble rbf_basis_func_ctps_c1(gdouble r)

{
  if ( r > 1.0 ) return 0.0 ;
  if ( r == 0.0 ) return 1.0 ;
  
  return 1.0 + r*(r*(80.0/3.0 - 40.0*r + 15.0*r*r -
		     8.0/3.0*r*r*r + 20.0*log(r))) ;
}

gdouble rbf_basis_func_ctps_c2a(gdouble r)

{
  if ( r > 1.0 ) return 0.0 ;
  if ( r == 0.0 ) return 1.0 ;

  return 1.0 + r*(-30.0*r + r*r*(-10.0 + r*(45.0 - 6.0*r) - 60.0*log(r))) ;

  /* return 1.0 - 30.0*r*r - 10.0*r*r*r + 45.0*r*r*r*r - 6.0*r*r*r*r*r - */
  /*    60.0*r*r*r*log(r) ; */
}

gdouble rbf_basis_func_ctps_c2b(gdouble r)

{
  if ( r > 1.0 ) return 0.0 ;
  if ( r == 0.0 ) return 1.0 ;

  return 1.0 + r*r*(-20.0 + r*(80.0 + r*(-45.0 - 16.0*r + 60.0*log(r)))) ;
  /* return 1.0 - 20.0*r*r + 80*r*r*r - 45.0*r*r*r*r - 16.0*r*r*r*r*r + */
  /*   60.0*r*r*r*r*log(r) ; */
}

static gint monomials_3d(gdouble *A, gdouble *x, gdouble *x0, gint m)

{
  if ( m < 0 ) return 0 ;

  gdouble r[3] ;

  r[0] = x[0] - x0[0] ; r[1] = x[1] - x0[1] ; r[2] = x[2] - x0[2] ; 

  switch ( m ) {
  default: g_assert_not_reached() ; break ;
  case 2:
    A[4] = r[0]*r[0] ; A[5] = r[1]*r[1] ; A[6] = r[2]*r[2] ;
    A[7] = r[0]*r[1] ; A[8] = r[1]*r[2] ; A[9] = r[2]*r[0] ; 
  case 1:
    A[1] = r[0] ; A[2] = r[1] ; A[3] = r[2] ;
  case 0:
    A[0] = 1.0 ;
    break ;
  }
  
  return 0 ;
}

gint rbf_parameter_fit_3d(gdouble *xs, gint xstr, gint ns,
			  gdouble *fs, gint fstr, gint nf,
			  gdouble *x0,
			  rbf_basis_func_t basis,
			  gdouble R, gint m,
			  gdouble *c,
			  gdouble *work)
{
  gdouble *A, r ;
  gint i, j, n, info, ip[4096], np ;
  
  g_assert(xstr > 2) ;
  g_assert(fstr >= nf) ;

  g_assert(m < 3) ;

  A = work ;
  np = 0 ;
  if ( m >= 0 ) np = (m+1)*(m+2)*(m+3)/6 ;
  /*size and initialize matrix*/
  n = ns + np ;

  memset(c, 0, nf*n*sizeof(gdouble)) ;
  
  for ( i = 0 ; i < ns ; i ++ ) {
    for ( j = i ; j < ns ; j ++ ) {
      r = sqrt((xs[i*xstr+0] - xs[j*xstr+0])*(xs[i*xstr+0] - xs[j*xstr+0]) +
	       (xs[i*xstr+1] - xs[j*xstr+1])*(xs[i*xstr+1] - xs[j*xstr+1]) +
	       (xs[i*xstr+2] - xs[j*xstr+2])*(xs[i*xstr+2] - xs[j*xstr+2])) ;
      /*this is FORTRAN indexing because we are making a direct call
	to a LAPACK matrix solver (though A is symmetric ...)*/
      A[i*n+j] = A[j*n+i] = basis(r/R) ;
    }
    monomials_3d(&(A[i*n+ns]), &(xs[i*xstr+0]), x0, m) ;
    for ( j = 0 ; j < np ; j ++ ) A[(ns+j)*n+i] = A[i*n+ns+j] ;
  }

  for ( j = 0 ; j < nf ; j ++ ) {
    /*FORTRAN indexing here too*/
    gint one = 1 ;
    blaswrap_dcopy(ns, &(fs[j]), fstr, &(c[j*n]), one) ;
  }

  dgesv_(&n, &nf, A, &n, ip, c, &n, &info) ;

  g_assert(info == 0) ;
  
  return 0 ;
}

gint rbf_evaluate_3d(gdouble *xs, gint xstr, gint ns,
		     gdouble *x0,
		     rbf_basis_func_t basis, gdouble R, gint m,
		     gdouble *c, gint nf,
		     gdouble *x, gdouble *f)

{
  gdouble r, rbf, p[20] ;
  gint i, j, n, np ;

  g_assert(m < 3) ;
  np = 0 ;
  if ( m >= 0 ) np = (m+1)*(m+2)*(m+3)/6 ;
  n = ns + np ;

  monomials_3d(p, x, x0, m) ;
  for ( j = 0 ; j < nf ; j ++ ) {
    f[j] = 0.0 ;
    for ( i = 0 ; i < np ; i ++ ) f[j] += p[i]*c[j*n+ns+i] ;
  }

  for ( i = 0 ; i < ns ; i ++ ) {
    gint one = 1 ;
    r = sqrt((xs[i*xstr+0] - x[0])*(xs[i*xstr+0] - x[0]) +
	     (xs[i*xstr+1] - x[1])*(xs[i*xstr+1] - x[1]) +
	     (xs[i*xstr+2] - x[2])*(xs[i*xstr+2] - x[2])) ;
    rbf = basis(r/R) ;
    blaswrap_daxpy(nf, rbf, &(c[i]), n, f, one) ;
  }

  return 0 ;
}
