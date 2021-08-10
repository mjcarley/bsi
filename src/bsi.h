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

#ifndef BSI_H_INCLUDED
#define BSI_H_INCLUDED

#include <glib.h>

#include <gtv.h>

#include <mop.h>

typedef gint (*bsi_source_func_t)(gdouble *x, gdouble *f, gint nf,
				  gpointer data) ;

typedef struct _bsi_point_t bsi_point_t ;

typedef struct _bsi_interpolator_t bsi_interpolator_t ;

struct _bsi_interpolator_t {
  gdouble *c, R ;
  gint m, np, nc, nf, nn, nnmax, *idx ;
  gpointer rbf ;
  GtsVertex *xc ;
} ;

typedef struct _bsi_polynomial_t bsi_polynomial_t ;

struct _bsi_polynomial_t {
  gint order, omax, nf, nfmax ;
  gdouble *c, xc[3] ;
} ;

typedef struct _bsi_polynomial_workspace_t bsi_polynomial_workspace_t ;

struct _bsi_polynomial_workspace_t {
  gint order, omax, nnbrs, nnmax, nf, nfmax, *nbrs ;
  gdouble *work ;
  mop_polynomial_t *mp ;
  mop_polynomial_workspace_t *wp ;
} ;

typedef enum
  {
   BSI_POLY,
   BSI_RBF_R1,
   BSI_RBF_R3,
   BSI_RBF_TPS,
   BSI_RBF_Q,
   BSI_RBF_MQ,
   BSI_RBF_IQ,
   BSI_RBF_IMQ,
   BSI_RBF_G,
   BSI_RBF_CP_C0,
   BSI_RBF_CP_C2,
   BSI_RBF_CP_C4,
   BSI_RBF_CP_C6
  } bsi_interpolation_t ;

#define bsi_poly_coefficient_3d_offset(_n)	\
  ((_n)*((_n)+1)*((_n)+2)/6)

#define bsi_poly_coefficient_3d_index(_i,_j,_k)		\
  (bsi_poly_coefficient_3d_offset(((_i)+(_j)+(_k))) +	\
   ((_i)*(2*((_i)+(_j)+(_k))+3-(_i))/2+(_j)))

gint bsi_distribution_refine(GtvVolume *v, GtsVertex *x, gdouble *f,
			     gint nf, gint nx, gint nxmax,
			     bsi_source_func_t sfunc, gpointer sdata,
			     gdouble tol) ;
GtvVolume *bsi_make_delaunay_volume(GtsVertex *x, gint nx) ;

gint bsi_source_func_ring_gaussian(gdouble *x, gdouble *f, gint nf,
				   gpointer data) ;

gint bsi_points_origin_width(GtsVertex *x, gint nx,
			     gdouble *xmin, gdouble *xmax,
			     gdouble *D, gboolean init) ;

gint bsi_node_neighbours(GtsVertex *x, gint i,
			 GtvVolume *v,
			 gint *nbrs, gint *nnbrs, gint nnmax) ;
gint bsi_cell_neighbours(GtvCell *c, GtsVertex *x, gint d,
			 gint *nbrs, gint *nnbrs, gint nnmax) ;
gint bsi_cell_indices(GtvCell *c, GtsVertex *x,
		      gint *i1, gint *i2, gint *i3, gint *i4) ;

bsi_interpolator_t *bsi_interpolator_new(bsi_interpolation_t interp,
					 gint m, gint nf, gint nnmax) ;
gint bsi_rbf_parameter_fit(GtsVertex *x, gint i, GtvVolume *v,
			   gdouble *f, gint nf, gint fstr, 
			   bsi_interpolator_t *interp, gdouble *work) ;
gint bsi_rbf_parameter_fit_cell(GtvCell *c, GtsVertex *x,
				gdouble *f, gint nf, gint fstr, 
				bsi_interpolator_t *interp, gdouble *work) ;
gint bsi_rbf_evaluate(bsi_interpolator_t *interp,
		      GtsVertex *xi, gdouble *x,
		      gint nf, gdouble *f) ;

gint bsi_poly_evaluate_3d(gdouble *p, gint np, gint N, gdouble *x,
			  gdouble *f) ;
gint bsi_poly_evaluate_derivative_3d(gdouble *p, gint np, gint N,
				     gint dx, gint dy, gint dz,
				     gdouble *x, gdouble *f) ;
gint bsi_poly_shift_3d(gdouble *p, gint np, gint N, gdouble *x, gdouble *ps) ;
/* gint bsi_monomials_3d(gint N, gdouble *x, gdouble *m, gint mstr) ; */
/* gint bsi_polynomial_fit(GtsVertex *x, gint i, GtvVolume *v, */
/* 			gdouble *f, gint nf, gint fstr,  */
/* 			bsi_interpolator_t *interp, gdouble *work) ; */

bsi_polynomial_t *bsi_polynomial_new(gint order, gint nf) ;
bsi_polynomial_workspace_t *bsi_polynomial_workspace_new(gint order,
							 gint nf,
							 gint nnbrs) ;
gint bsi_polynomial_cell_interpolant(GtvCell *c, GtsVertex *x,
				     gdouble *f, gint nf, gint fstr,
				     gint depth,
				     bsi_polynomial_t *p, gint order,
				     bsi_polynomial_workspace_t *w) ;
gint bsi_polynomial_evaluate(bsi_polynomial_t *p, gdouble *x, gdouble *f) ;
gint bsi_polynomial_shift_origin(bsi_polynomial_t *p, gdouble *x,
				 bsi_polynomial_workspace_t *w) ;
gint bsi_mesh_interp_linear(GtsVertex *x, GtvVolume *v,
			    gdouble *f, gint nf, gint fstr,
			    GtsPoint *xi, gdouble *fi) ;


gint bsi_tet_integrate_powers(gdouble *x1, gdouble *x2, gdouble *x3,
			      gdouble *x4,
			      gint N,
			      /* gint *p, gint np, */
			      gint nq, gdouble tol, gint dmax,
			      gdouble *quad) ;
/* gint bsi_tet_integrate_powers(gdouble *x1, gdouble *x2, gdouble *x3, */
/* 			      gdouble *x4, */
/* 			      gint *p, gint np, */
/* 			      gint nq, gdouble tol, gint dmax, */
/* 			      gdouble *quad) ; */

#endif /*BSI_H_INCLUDED*/
