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
 * Radial Basis Function interpolation, based on the Matlab code at
 * https://github.com/CCM-Deusto/RBF-interp
 */

#ifndef _RBF_H_INCLUDED_
#define _RBF_H_INCLUDED_

typedef gdouble (*rbf_basis_func_t)(gdouble r) ;

gdouble rbf_basis_func_r1(gdouble r) ;
gdouble rbf_basis_func_r3(gdouble r) ;
gdouble rbf_basis_func_tps(gdouble r) ;
gdouble rbf_basis_func_q(gdouble r) ;
gdouble rbf_basis_func_mq(gdouble r) ;
gdouble rbf_basis_func_imq(gdouble r) ;
gdouble rbf_basis_func_iq(gdouble r) ;
gdouble rbf_basis_func_g(gdouble r) ;
gdouble rbf_basis_func_cp_c0(gdouble r) ;
gdouble rbf_basis_func_cp_c2(gdouble r) ;
gdouble rbf_basis_func_cp_c4(gdouble r) ;
gdouble rbf_basis_func_cp_c6(gdouble r) ;
gdouble rbf_basis_func_ctps_c0(gdouble r) ;
gdouble rbf_basis_func_ctps_c1(gdouble r) ;
gdouble rbf_basis_func_ctps_c2a(gdouble r) ;
gdouble rbf_basis_func_ctps_c2b(gdouble r) ;

gint rbf_parameter_fit_3d(gdouble *xs, gint xstr, gint ns,
			  gdouble *fs, gint fstr, gint nf,
			  gdouble *x0,
			  rbf_basis_func_t basis,
			  gdouble R, gint m,
			  gdouble *c,
			  gdouble *work) ;
gint rbf_evaluate_3d(gdouble *xs, gint xstr, gint ns,
		     gdouble *x0,
		     rbf_basis_func_t basis, gdouble R, gint m,
		     gdouble *c, gint nf,
		     gdouble *x, gdouble *f) ;

#endif /*_RBF_H_INCLUDED_*/
