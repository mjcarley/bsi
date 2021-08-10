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

#include <blaswrap.h>

#include "rbf.h"

extern gint dgesv_(gint *n, gint *nrhs, gdouble *A, gint *lda,
		   gint *ip, gdouble *b, gint *ldb, gint *info) ;

extern gint dgels_(gchar *trans, gint *m, gint *n, gint *nrhs,
		   gdouble *A, gint *lda, gdouble *B, gint *ldb,
		   gdouble *work, gint *lwork, gint *info) ;

static gint compare_distance(gconstpointer a, gconstpointer b,
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

bsi_interpolator_t *bsi_interpolator_new(bsi_interpolation_t interp,
					 gint m, gint nf, gint nnmax)

{
  bsi_interpolator_t *i ;

  i = (bsi_interpolator_t *)g_malloc0(sizeof(bsi_interpolator_t)) ;

  i->np = 0 ;
  i->nn = 0 ;
  i->nnmax = nnmax ;
  i->m = m ;
  i->nf = nf ;
  
  if ( m > 0 ) i->np = (m+1)*(m+2)*(m+3)/6 ;

  i->c = (gdouble *)g_malloc0((i->np+i->nnmax)*(i->nf)*sizeof(gdouble)) ;
  i->idx = (gint *)g_malloc0(nnmax*sizeof(gint)) ;

  switch ( interp ) {
  default: g_assert_not_reached() ; break ;
  case BSI_POLY: i->rbf = NULL ; break ;
  case BSI_RBF_R1: i->rbf = rbf_basis_func_r1 ; break ;
  case BSI_RBF_R3: i->rbf = rbf_basis_func_r3 ; break ;
  case BSI_RBF_TPS: i->rbf = rbf_basis_func_tps ; break ;
  case BSI_RBF_Q: i->rbf = rbf_basis_func_q ; break ;
  case BSI_RBF_MQ: i->rbf = rbf_basis_func_mq ; break ;
  case BSI_RBF_IQ: i->rbf = rbf_basis_func_iq ; break ;
  case BSI_RBF_IMQ: i->rbf = rbf_basis_func_imq ; break ;
  case BSI_RBF_G: i->rbf = rbf_basis_func_g ; break ;
  case BSI_RBF_CP_C0: i->rbf = rbf_basis_func_cp_c0 ; break ;
  case BSI_RBF_CP_C2: i->rbf = rbf_basis_func_cp_c2 ; break ;
  case BSI_RBF_CP_C4: i->rbf = rbf_basis_func_cp_c4 ; break ;
  case BSI_RBF_CP_C6: i->rbf = rbf_basis_func_cp_c6 ; break ;
  }
  
  return i ;
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

gint bsi_rbf_parameter_fit(GtsVertex *x, gint i, GtvVolume *v,
			   gdouble *f, gint nf, gint fstr, 
			   bsi_interpolator_t *interp, gdouble *work)

{
  rbf_basis_func_t basis ;
  gint n, info, ip[4096], j, k ;
  gdouble *A, r, R ;
  gpointer udata[2] ;
  
  g_assert(fstr >= nf) ;
  g_assert(nf <= interp->nf) ;
  g_assert(interp->m < 3) ;

  R = 2.0 ;
  interp->R = R ;
  
  basis = interp->rbf ;
  
  bsi_node_neighbours(x, i, v, interp->idx, &(interp->nn), interp->nnmax) ;
  /* interp->idx[0] = i ; */
  /* bsi_node_neighbours(x, i, v, interp->idx, &(interp->nn), interp->nnmax) ; */

  interp->xc = &(x[i]) ;
  udata[0] = x ; udata[1] = interp->xc ;
  g_qsort_with_data(&(interp->idx[0]), interp->nn, sizeof(gint),
		    compare_distance, udata) ;
  
  A = work ;

  n = interp->nn + interp->np ;

  memset(interp->c, 0, nf*n*sizeof(gdouble)) ;
  
  for ( j = 0 ; j < interp->nn ; j ++ ) {
    for ( k = j ; k < interp->nn ; k ++ ) {
      r = gts_point_distance(GTS_POINT(&(x[interp->idx[j]])),
			     GTS_POINT(&(x[interp->idx[k]]))) ;
      /*this is FORTRAN indexing because we are making a direct call
	to a LAPACK matrix solver (though A is symmetric ...)*/
      A[j*n+k] = A[k*n+j] = basis(r/R) ;
    }
    monomials_3d(&(A[j*n+interp->nn]),
		 &(GTS_POINT(&(x[interp->idx[j]]))->x),
		 &(GTS_POINT(interp->xc)->x), interp->m) ;
    for ( k = 0 ; k < interp->np ; k ++ )
      A[(interp->nn+k)*n+j] = A[j*n+interp->nn+k] ;
  }

  for ( j = 0 ; j < interp->nn ; j ++ ) {
    for ( k = 0 ; k < nf ; k ++ ) {
      interp->c[k*n+j] = f[interp->idx[j]*fstr+k] ;
    }
  }

  dgesv_(&n, &nf, A, &n, ip, interp->c, &n, &info) ;

  if ( info != 0 )
    g_error("%s: dgesv failed with error code %d", __FUNCTION__, info) ;
  
  return 0 ;
}

gint bsi_rbf_parameter_fit_cell(GtvCell *c, GtsVertex *x,
				gdouble *f, gint nf, gint fstr, 
				bsi_interpolator_t *interp, gdouble *work)

{
  rbf_basis_func_t basis ;
  gint n, info, ip[4096], j, k ;
  gdouble *A, r, R ;
  gpointer udata[2] ;
  
  g_assert(fstr >= nf) ;
  g_assert(nf <= interp->nf) ;
  g_assert(interp->m < 3) ;

  R = 1.0 ;
  interp->R = R ;
  
  basis = interp->rbf ;
  
  /* bsi_node_neighbours(x, i, v, interp->idx, &(interp->nn), interp->nnmax) ; */
  bsi_cell_neighbours(c, x, 4, interp->idx, &(interp->nn), interp->nnmax) ;
  A = work ;

  g_assert(interp->nn <= interp->nnmax) ;
  
  n = interp->nn + interp->np ;

  memset(interp->c, 0, nf*n*sizeof(gdouble)) ;

  /* g_assert(n*n < 4096) ; */
  
  interp->xc = &(x[interp->idx[0]]) ;

  udata[0] = x ; udata[1] = interp->xc ;
  g_qsort_with_data(&(interp->idx[0]), interp->nn, sizeof(gint),
		    compare_distance, udata) ;
  
  for ( j = 0 ; j < interp->nn ; j ++ ) {
    for ( k = j ; k < interp->nn ; k ++ ) {
      r = gts_point_distance(GTS_POINT(&(x[interp->idx[j]])),
			     GTS_POINT(&(x[interp->idx[k]]))) ;
      /*this is FORTRAN indexing because we are making a direct call
	to a LAPACK matrix solver (though A is symmetric ...)*/
      A[j*n+k] = A[k*n+j] = basis(r/R) ;
    }
    monomials_3d(&(A[j*n+interp->nn]),
		 &(GTS_POINT(&(x[interp->idx[j]]))->x),
		 &(GTS_POINT(interp->xc)->x), interp->m) ;
    for ( k = 0 ; k < interp->np ; k ++ )
      A[(interp->nn+k)*n+j] = A[j*n+interp->nn+k] ;
  }

  for ( j = 0 ; j < interp->nn ; j ++ ) {
    for ( k = 0 ; k < nf ; k ++ ) {
      interp->c[k*n+j] = f[interp->idx[j]*fstr+k] ;
    }
  }

  dgesv_(&n, &nf, A, &n, ip, interp->c, &n, &info) ;

  if ( info != 0 )
    g_error("%s: dgesv failed with error code %d", __FUNCTION__, info) ;
  
  return 0 ;
}

gint bsi_rbf_evaluate(bsi_interpolator_t *interp,
		      GtsVertex *xi,
		      gdouble *x,
		      gint nf, gdouble *f)

{
  rbf_basis_func_t basis ;
  gint n, i, j ;
  gdouble r, R, p[20], rbf ;
  
  g_assert(nf <= interp->nf) ;

  R = interp->R ;
  
  basis = interp->rbf ;

  n = interp->nn + interp->np ;

  monomials_3d(p, x, &(GTS_POINT(interp->xc)->x), interp->m) ;
  for ( j = 0 ; j < nf ; j ++ ) {
    f[j] = 0.0 ;
    for ( i = 0 ; i < interp->np ; i ++ )
      f[j] += p[i]*(interp->c[j*n+interp->nn+i]) ;
  }

  for ( i = 0 ; i < interp->nn ; i ++ ) {
    gint one = 1 ;
    r = sqrt(((GTS_POINT(&(xi[interp->idx[i]]))->x)-x[0])*
	     ((GTS_POINT(&(xi[interp->idx[i]]))->x)-x[0]) +
	     ((GTS_POINT(&(xi[interp->idx[i]]))->y)-x[1])*
	     ((GTS_POINT(&(xi[interp->idx[i]]))->y)-x[1]) +
	     ((GTS_POINT(&(xi[interp->idx[i]]))->z)-x[2])*
	     ((GTS_POINT(&(xi[interp->idx[i]]))->z)-x[2])) ;
    rbf = basis(r/R) ;
    blaswrap_daxpy(nf, rbf, &(interp->c[i]), n, f, one) ;
  }
  
  return 0 ;
}

gint bsi_poly_evaluate_3d(gdouble *p, gint np, gint N, gdouble *x,
			  gdouble *f)

{
  gdouble pw[64] ;
  gint n, i, j, k, m ;

  for ( m = 0 ; m < np ; m ++ ) f[m] = p[np*0+m] ;
  /* f = p[0] ; */

  pw[0] = pw[1] = pw[2] = 1.0 ;

  for ( n = 1 ; n <= N ; n ++ ) {

    pw[3*n+0] = pw[3*(n-1)+0]*x[0] ;
    pw[3*n+1] = pw[3*(n-1)+1]*x[1] ;
    pw[3*n+2] = pw[3*(n-1)+2]*x[2] ;
    
    for ( i = 0 ; i <= n ; i ++ ) {
      for ( j = 0 ; j <= n-i ; j ++ ) {
	k = n - i - j ;
	for ( m = 0 ; m < np ; m ++ ) {
	  f[m] += p[np*bsi_poly_coefficient_3d_index(i,j,k)+m]*
	    pw[3*i+0]*pw[3*j+1]*pw[3*k+2] ;
	}
      }
    }
  }
  
  return 0 ;
}

gint bsi_poly_evaluate_derivative_3d(gdouble *p, gint np,
				     gint N,
				     gint dx, gint dy, gint dz,
				     gdouble *x, gdouble *f)


/*
 * this is not strictly the derivative but the derivative scaled on
 * dx!dy!dz! (since it is mostly used in finding coefficients of
 * Taylor series)
*/
  
{
  gdouble pw[64] ;
  gint n, i, j, k, m ;

  if ( dx == 0 && dy == 0 && dz == 0 )
    return bsi_poly_evaluate_3d(p, np, N, x, f) ;
  
  for ( m = 0 ; m < np ; m ++ ) f[m] = 0.0 ;
  
  pw[0] = pw[1] = pw[2] = 1.0 ;

  for ( n = 1 ; n <= N ; n ++ ) {

    pw[3*n+0] = pw[3*(n-1)+0]*x[0] ;
    pw[3*n+1] = pw[3*(n-1)+1]*x[1] ;
    pw[3*n+2] = pw[3*(n-1)+2]*x[2] ;
    
    for ( i = dx ; i <= n ; i ++ ) {
      for ( j = dy ; j <= n-i ; j ++ ) {
	k = n - i - j ;
	if ( k >= dz ) {
	  for ( m = 0 ; m < np ; m ++ ) {
	    f[m] += p[np*bsi_poly_coefficient_3d_index(i,j,k)+m]*
	      pw[3*(i-dx)+0]*pw[3*(j-dy)+1]*pw[3*(k-dz)+2]*
	      tgamma(i+1)/tgamma(i-dx+1)/tgamma(dx+1)*
	      tgamma(j+1)/tgamma(j-dy+1)/tgamma(dy+1)*
	      tgamma(k+1)/tgamma(k-dz+1)/tgamma(dz+1) ;
	  }
	}
      }
    }
  }
  
  return 0 ;
}

gint bsi_poly_shift_3d(gdouble *p, gint np, gint N, gdouble *x, gdouble *ps)

{
  gint n, i, j, k, idx ;
  
  memset(ps, 0, bsi_poly_coefficient_3d_offset(N+1)*np*sizeof(gdouble)) ;

  bsi_poly_evaluate_3d(p, np, N, x, &(ps[np*0])) ;

  for ( n = 1 ; n <= N ; n ++ ) {
    for ( i = 0 ; i <= n ; i ++ ) {
      for ( j = 0 ; j <= n-i ; j ++ ) {
	k = n - i - j ;
	idx = bsi_poly_coefficient_3d_index(i,j,k) ;
	bsi_poly_evaluate_derivative_3d(p, np, N, i, j, k, x, &(ps[np*idx])) ;
      }      
    }
  }
  
  return 0 ;
}

gint bsi_monomials_3d(gint N, gdouble *x, gdouble *m, gint mstr)

/*
 * fill in a bunch of monomials up to order N
 */
  
{
  gdouble pw[64] ;
  gint n, i, j, k ;
  
  m[0] = 1.0 ;

  pw[0] = pw[1] = pw[2] = 1.0 ;

  for ( n = 1 ; n <= N ; n ++ ) {

    pw[3*n+0] = pw[3*(n-1)+0]*x[0] ;
    pw[3*n+1] = pw[3*(n-1)+1]*x[1] ;
    pw[3*n+2] = pw[3*(n-1)+2]*x[2] ;
    
    for ( i = 0 ; i <= n ; i ++ ) {
      for ( j = 0 ; j <= n-i ; j ++ ) {
	k = n - i - j ;
	m[bsi_poly_coefficient_3d_index(i,j,k)*mstr] =
	  pw[3*i+0]*pw[3*j+1]*pw[3*k+2] ;
      }
    }
  }
  
  return 0 ;
}

gint bsi_polynomial_fit(GtsVertex *x, gint i, GtvVolume *v,
			gdouble *f, gint nf, gint fstr, 
			bsi_interpolator_t *interp, gdouble *work)

{
  gint n, info, j, k, lwork, lda ;
  gdouble *A, xi[3], iwork[4096] ;

  lwork = 4096 ;
  
  g_assert(fstr >= nf) ;
  g_assert(nf <= interp->nf) ;
  /* g_assert(interp->m < 3) ; */

  interp->rbf = NULL ;
  
  bsi_node_neighbours(x, i, v, interp->idx, &(interp->nn), interp->nnmax) ;

  A = work ;

  n = interp->nn ;

  interp->np = bsi_poly_coefficient_3d_offset(interp->m+1) ;
  
  /* g_assert(n >= interp->np) ; */

  lda = MIN(n, interp->np) ;
  
  memset(interp->c, 0, nf*(interp->np)*sizeof(gdouble)) ;
  interp->xc = &(x[i]) ;

  for ( j = 0 ; j < interp->nn ; j ++ ) {
    gts_vector_init(xi,
		    GTS_POINT(interp->xc),
		    GTS_POINT(&(x[interp->idx[j]]))) ;
    /*FORTRAN indexing in setting up polynomial fit matrix*/
    bsi_monomials_3d(interp->m, xi, &(A[j]), lda) ;

    for ( k = 0 ; k < nf ; k ++ ) {
      interp->c[k*n+j] = f[interp->idx[j]*fstr+k] ;
    }
  }

  dgels_("N", &n, &lda, &nf, A, &n, interp->c, &(interp->np), iwork,
	 &lwork, &info) ;
  
  return 0 ;
}

gint bsi_node_neighbours(GtsVertex *x, gint i,
			 GtvVolume *v,
			 gint *nbrs, gint *nnbrs, gint nnmax)

/*
 * find neighbours of node x[i]
 */
  
{
  GSList *cells, *vertices, *k ;
  GtsVertex *vt[4] ;
  gint j ;
  
  *nnbrs = 0 ;
  
  cells = gtv_vertex_cells(&(x[i]), v) ;

  bsi_cell_neighbours(GTV_CELL(cells->data), x, 4, nbrs, nnbrs, nnmax) ;

  return 0 ;
  
  for ( k = cells ; k != NULL ; k = k->next ) {
    cells = g_slist_concat(gtv_cell_neighbours(GTV_CELL(k->data), v),
			   cells) ;
  }
  
  vertices = NULL ;
  for ( ; cells != NULL ; cells = cells->next ) {
    gtv_tetrahedron_vertices(GTV_TETRAHEDRON(cells->data),
			     &(vt[0]), &(vt[1]), &(vt[2]), &(vt[3])) ;
    for ( j = 0 ; j < 4 ; j ++ ) {
      if ( /* (vt[j] != &(x[i])) && */
	   (g_slist_find(vertices, vt[j]) == NULL) ) {
	vertices = g_slist_prepend(vertices, vt[j]) ;
      }
    }
  }

  for ( ; (vertices != NULL) && ((*nnbrs) < nnmax) ;
	vertices = vertices->next ) {
    nbrs[(*nnbrs)] = GTS_VERTEX(vertices->data) - &(x[0]) ;
    (*nnbrs) ++ ;
  }
  
  return 0 ;
}

static gboolean in_list(gint *list, gint n, gint j)

{
  gint i ;

  for ( i = 0 ; i < n ; i ++ ) {
    if ( list[i] == j ) return TRUE ;
  }
  
  return FALSE ;
}

gint bsi_cell_neighbours(GtvCell *c, GtsVertex *x, gint d,
			 gint *nbrs, gint *nnbrs, gint nnmax)

/*
 * find neighbour vertices of cell c
 *
 * first four entries are the indices of the cell vertices and
 * subsequent entries are the indices of vertices in successive
 * "layers" around the cell up to a depth d
 */
  
{
  GSList *i ;
  GtsVertex *v1, *v2, *v3, *v4 ;
  gint j, k, k0, k1 ;

  g_assert(nnmax > 4) ;
  g_assert(d >= 0) ;
  
  gtv_tetrahedron_vertices(GTV_TETRAHEDRON(c), &v1, &v2, &v3, &v4) ;

  nbrs[0] = v1 - &(x[0]) ; nbrs[1] = v2 - &(x[0]) ;
  nbrs[2] = v3 - &(x[0]) ; nbrs[3] = v4 - &(x[0]) ;

  *nnbrs = 4 ;
  k0 = 0 ; k1 = 4 ;

  for ( ; d > 0 ; d -- ) {
    for ( k = k0 ; k < k1 ; k ++ ) {
      v1 = &(x[nbrs[k]]) ;
      for ( i = v1->segments ; i != NULL ; i = i->next ) {
	j = GTS_SEGMENT(i->data)->v1 - &(x[0]) ;
	if ( in_list(nbrs, *nnbrs, j) == FALSE ) {
	  nbrs[(*nnbrs)] = j ; (*nnbrs) ++ ;
	}
	if ( (*nnbrs) == nnmax ) return 0 ;

	j = GTS_SEGMENT(i->data)->v2 - &(x[0]) ;
	if ( in_list(nbrs, *nnbrs, j) == FALSE ) {
	  nbrs[(*nnbrs)] = j ; (*nnbrs) ++ ;
	}
	if ( (*nnbrs) == nnmax ) return 0 ;
      }
    }
    k0 = k1 ; k1 = (*nnbrs) ;
  }
  
  return 0 ;
}

bsi_polynomial_t *bsi_polynomial_new(gint order, gint nf)

{
  bsi_polynomial_t *p ;

  p = (bsi_polynomial_t *)g_malloc0(sizeof(bsi_polynomial_t)) ;

  p->omax  = order ; p->order = -1 ;
  p->nfmax = nf ; p->nf = -1 ;
  
  p->c = (gdouble *)g_malloc0
    (nf*bsi_poly_coefficient_3d_offset(order+1)*sizeof(gdouble)) ;
  
  return p ;
}

bsi_polynomial_workspace_t *bsi_polynomial_workspace_new(gint order,
							 gint nf,
							 gint nnbrs)

{
  bsi_polynomial_workspace_t *w ;
  gint nc ;
  
  w = (bsi_polynomial_workspace_t *)
    g_malloc0(sizeof(bsi_polynomial_workspace_t)) ;

  /* nc = (order+1)*(order+2)*(order+3)/6 ; */
  nc = bsi_poly_coefficient_3d_offset(order+1) ;
  
  w->omax = order ;  w->order = -1 ;
  w->nnmax = nnbrs ; w->nnbrs = -1 ;
  w->nfmax = nf ;    w->nf    = -1 ;
  w->nbrs = (gint *)g_malloc0(nnbrs*sizeof(gint)) ;
  w->work = (gdouble *)g_malloc0(((3+nf)*nnbrs+nc*nf)*sizeof(gdouble)) ;
  
  w->mp = mop_polynomial_alloc(nnbrs, nf, order) ;
  w->wp = mop_polynomial_workspace_alloc(nnbrs, nf, order) ;
  
  return w ;
}

static gint tetrahedron_centroid(GtvTetrahedron *t, gdouble *c)

{
  GtsVertex *v1, *v2, *v3, *v4 ;

  gtv_tetrahedron_vertices(t, &v1, &v2, &v3, &v4) ;

  c[0] = 0.25*(GTS_POINT(v1)->x + GTS_POINT(v2)->x +
	       GTS_POINT(v3)->x + GTS_POINT(v4)->x) ;
  c[1] = 0.25*(GTS_POINT(v1)->y + GTS_POINT(v2)->y +
	       GTS_POINT(v3)->y + GTS_POINT(v4)->y) ;
  c[2] = 0.25*(GTS_POINT(v1)->z + GTS_POINT(v2)->z +
	       GTS_POINT(v3)->z + GTS_POINT(v4)->z) ;
  
  return 0 ;
}

static gint merge_interpolant(mop_polynomial_t *p, gdouble *c, gint nc,
			      gdouble *pm)

{
  gint n, i, j, k, idx ;
  gdouble pc ;
  
  n = mop_polynomial_order(p) ;

  memset(pm, 0, nc*(n+1)*(n+2)*(n+3)/6*sizeof(gdouble)) ;

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

gint bsi_polynomial_cell_interpolant(GtvCell *c, GtsVertex *x,
				     gdouble *f, gint nf, gint fstr,
				     gint depth,
				     bsi_polynomial_t *p, gint order,
				     bsi_polynomial_workspace_t *w)

{
  gpointer udata[4] ;
  gdouble *xi, *fi, *ci, tol ;
  gint i, j, k ;
  
  g_assert(nf <= p->nfmax) ;
  g_assert(nf <= w->nfmax) ;
  g_assert(order <= p->omax) ;
  g_assert(order <= w->omax) ;

  tol = 1e-6 ;
  
  bsi_cell_neighbours(c, x, depth, w->nbrs, &(w->nnbrs), w->nnmax) ;

  tetrahedron_centroid(GTV_TETRAHEDRON(c), p->xc) ;
  
  udata[0] = x ; udata[1] = p->xc ;
  g_qsort_with_data(&(w->nbrs[4]), w->nnbrs-4, sizeof(gint),
		    compare_distance, udata) ;

  xi = w->work ; fi = &(xi[3*w->nnbrs]) ; ci = &(fi[nf*w->nnbrs]) ;
  
  for ( i = 0 ; i < w->nnbrs ; i ++ ) {
    j = w->nbrs[i] ;
    xi[3*i+0] = GTS_POINT(&(x[j]))->x - p->xc[0] ;
    xi[3*i+1] = GTS_POINT(&(x[j]))->y - p->xc[1] ;
    xi[3*i+2] = GTS_POINT(&(x[j]))->z - p->xc[2] ;
    for ( k = 0 ; k < nf ; k ++ ) fi[nf*i+k] = f[fstr*j+k] ;
  }
  
  mop_polynomial_set_points(w->mp, xi, NULL, w->nnbrs) ;
  /* mop_polynomial_basis_points(w->mp, order, tol, w->wp) ; */
  mop_polynomial_basis_power(w->mp, order, tol, w->wp) ;

  mop_polynomial_make(w->mp, w->wp) ;

  mop_polynomial_normalize(w->mp, w->wp) ;

  mop_polynomial_transform(w->mp, fi, nf, ci, w->wp) ;

  merge_interpolant(w->mp, ci, nf, p->c) ;

  p->nf = nf ; p->order = order ;
  
  return 0 ;
}

gint bsi_polynomial_evaluate(bsi_polynomial_t *p, gdouble *x, gdouble *f)

{
  gdouble xc[3] ;

  xc[0] = x[0]-p->xc[0] ; xc[1] = x[1]-p->xc[1] ; xc[2] = x[2]-p->xc[2] ;

  bsi_poly_evaluate_3d(p->c, p->nf, p->order, xc, f) ;
  
  return 0 ;
}

gint bsi_polynomial_shift_origin(bsi_polynomial_t *p, gdouble *x,
				 bsi_polynomial_workspace_t *w)

{
  gint n ;
  gdouble dx[3] ;

  dx[0] = x[0]-p->xc[0] ; dx[1] = x[1]-p->xc[1] ; dx[2] = x[2]-p->xc[2] ;
  
  n = bsi_poly_coefficient_3d_offset(p->order+1) ;
  memcpy(w->work, p->c, n*p->nf*sizeof(gdouble)) ;

  bsi_poly_shift_3d(w->work, p->nf, p->order, dx, p->c) ;
  
  p->xc[0] = x[0] ; p->xc[1] = x[1] ; p->xc[2] = x[2] ;

  return 0 ;
}

gint bsi_mesh_interp_linear(GtsVertex *x, GtvVolume *v,
			    gdouble *f, gint nf, gint fstr,
			    GtsPoint *xi, gdouble *fi)

{
  GtvCell *c ;
  gdouble w[4] ;
  gint i1, i2, i3, i4, j ;
  
  c = gtv_point_locate(GTS_POINT(xi), v, NULL) ;

  if ( c == NULL ) return -1 ;

  gtv_tetrahedron_point_barycentric(GTV_TETRAHEDRON(c), xi, w) ;

  bsi_cell_indices(c, x, &i1, &i2, &i3, &i4) ;

  for ( j = 0 ; j < nf ; j ++ ) {
    fi[j] = w[0]*f[i1*fstr+j] + w[1]*f[i2*fstr+j] +
      w[2]*f[i3*fstr+j] + w[3]*f[i4*fstr+j] ;
  }
  
  return 0 ;
}
