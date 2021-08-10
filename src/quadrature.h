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

/*
 * High order quadrature rules for tetrahedra, from the supplementary
 * material of
 * 
 * Jan Jaskowiec and N. Sukumar, High-order cubature rules for
 * tetrahedra, Int J Numer Methods Eng. 121:2418-2436, 2020, 
 * http://dx.doi.org/10.1002/nme.6313
 */

#ifndef _QUADRATURE_H_INCLUDED_
#define _QUADRATURE_H_INCLUDED_

#include <glib.h>

gint bsi_quadrature_js_select(gint p, gdouble **q, gint *n) ;

#endif /*_QUADRATURE_H_INCLUDED_*/
