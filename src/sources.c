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

#include <math.h>
#include <bsi.h>

gint bsi_source_func_ring_gaussian(gdouble *x, gdouble *f, gint nf,
				   gpointer data)

{
  gdouble s, r, r2, r0, th, om ;

  if ( data == NULL ) {
    s = 0.3 ; r0 = 1.0 ;
  } else {
    gdouble *d = data ;
    s = d[0] ; r0 = d[1] ;
  }

  r = sqrt(x[0]*x[0] + x[1]*x[1]) ;
  th = atan2(x[1], x[0]) ;

  r2 = (r-r0)*(r-r0) + x[2]*x[2] ;

  om = exp(-r2/s/s) ;

  f[0] = -om*sin(th) ; 
  f[1] =  om*cos(th) ; 
  f[2] =  0.0 ;
  
  return 0 ;
}
