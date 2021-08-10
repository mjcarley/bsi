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
#include <math.h>
#include <string.h>

#include <gtv.h>

#include <bsi.h>

static void refine_cell(GtvCell *c, gpointer *rdata)

{
  GtsVertex *nodes = rdata[0] ;
  gdouble *data = rdata[1] ;
  gint nd = *(gint *)rdata[2] ;
  gint *nn = rdata[3] ;
  gint nnmax = *(gint *)rdata[4] ;
  gdouble tol = *(gdouble *)rdata[5] ;
  gpointer idata = rdata[6] ;
  bsi_source_func_t sfunc = rdata[7] ;
  GtsVertex *v1, *v2, *v3, *v4, *vt ;
  gdouble fi[8], xt, yt, zt ;
  gint i1, i2, i3, i4, i ;
  
  if ( *nn >= nnmax ) return ;
  
  gtv_tetrahedron_vertices(GTV_TETRAHEDRON(c), &v1, &v2, &v3, &v4) ;

  /* i1 = v1 - &(nodes[0]) ; */
  /* i2 = v2 - &(nodes[0]) ; */
  /* i3 = v3 - &(nodes[0]) ; */
  /* i4 = v4 - &(nodes[0]) ; */

  bsi_cell_indices(c, nodes, &i1, &i2, &i3, &i4) ;
  
  /* g_assert(&(nodes[i1]) == v1) ; */
  for ( i = 0 ; i < nd ; i ++ ) {
    fi[i] = 0.25*(data[i1*nd+i] + data[i2*nd+i] +
		  data[i3*nd+i] + data[i4*nd+i]) ; 
  }

  xt = 0.25*(GTS_POINT(v1)->x + GTS_POINT(v2)->x +
	     GTS_POINT(v3)->x + GTS_POINT(v4)->x) ;
  yt = 0.25*(GTS_POINT(v1)->y + GTS_POINT(v2)->y +
	     GTS_POINT(v3)->y + GTS_POINT(v4)->y) ;
  zt = 0.25*(GTS_POINT(v1)->z + GTS_POINT(v2)->z +
	     GTS_POINT(v3)->z + GTS_POINT(v4)->z) ;

  vt = gts_vertex_new(gts_vertex_class(), xt, yt, zt) ;

  sfunc(&(GTS_POINT(vt)->x), &(data[(*nn)*nd]), nd, idata) ;

  for ( i = 0 ; i < nd ; i ++ ) {
    if ( fabs(data[(*nn)*nd+i] - fi[i]) > tol ) {
      memcpy(&(nodes[(*nn)]), vt, sizeof(GtsVertex)) ;
      (*nn) ++ ;
      /* gts_object_destroy(GTS_OBJECT(vt)) ; */
      return ;
    }
  }
  
  /* gts_object_destroy(GTS_OBJECT(vt)) ; */

  return ;
}

gint bsi_distribution_refine(GtvVolume *v, GtsVertex *x, gdouble *f,
			     gint nf, gint nx, gint nxmax,
			     bsi_source_func_t sfunc, gpointer sdata,
			     gdouble tol)

{
  gint i ;
  gpointer rdata[8] ;
  gboolean stop ;
  
  for ( i = 0 ; i < nx ; i ++ ) {
    sfunc(&(GTS_POINT(&(x[i]))->x), &(f[i*nf]), nf, sdata) ;
  }

  rdata[0] = x ;
  rdata[1] = f ;
  rdata[2] = &nf ;
  rdata[3] = &nx ;
  rdata[4] = &nxmax ;
  rdata[5] = &tol ;
  rdata[6] = sdata ;
  rdata[7] = sfunc ;
  
  stop = FALSE ;
  while ( stop == FALSE ) {
    i = nx ;

    gtv_volume_foreach_cell(v, (GtsFunc)refine_cell, rdata) ;
    
    if ( i == nx ) stop = TRUE ;
    if ( nx == nxmax ) stop = TRUE ;

    for ( ; i < nx ; i ++ ) {
      gtv_delaunay_add_vertex(v, &(x[i]), NULL) ;
    }
  }
  
  return nx ;
}

GtvVolume *bsi_make_delaunay_volume(GtsVertex *x, gint nx)

{
  GtvVolume *v ;
  gdouble len ;
  gint i ;
  GSList *j, *cells ;
  GtsVertex *v1, *v2, *v3, *v4, *p ;
  GtvCell *c ;
  
  len = 1e6 ;
  
  v = gtv_volume_new(gtv_volume_class(),
		     gtv_cell_class(),
		     gtv_facet_class(),
		     gts_edge_class(),
		     gts_vertex_class()) ;

  c = GTV_CELL(gtv_tetrahedron_large((GtvTetrahedronClass *)
				     gtv_cell_class(),
				     gtv_facet_class(),
				     gts_edge_class(),
				     gts_vertex_class(),
				     len)) ;
  gtv_volume_add_cell(v, c) ;

  gtv_tetrahedron_vertices(GTV_TETRAHEDRON(c), &v1, &v2, &v3, &v4) ;

  for ( i = 0 ; i < nx ; i ++ ) {
    p = &(x[i]) ;
    if ( gtv_delaunay_add_vertex(v, p, NULL) == GTV_VERTEX_NOT_IN_VOLUME )
      fprintf(stderr,
	      "vertex (%lg,%lg,%lg) not inside convex hull\n",
	      GTS_POINT(p)->x, GTS_POINT(p)->y, GTS_POINT(p)->z) ;
  }

  cells = gtv_vertex_cells(v1, v) ;
  for ( j = cells ; j != NULL ; j = j->next )
    gtv_volume_remove_cell(v, GTV_CELL(j->data)) ;

  cells = gtv_vertex_cells(v2, v) ;
  for ( j = cells ; j != NULL ; j = j->next )
    gtv_volume_remove_cell(v, GTV_CELL(j->data)) ;

  cells = gtv_vertex_cells(v3, v) ;
  for ( j = cells ; j != NULL ; j = j->next )
    gtv_volume_remove_cell(v, GTV_CELL(j->data)) ;
  
  cells = gtv_vertex_cells(v4, v) ;
  for ( j = cells ; j != NULL ; j = j->next )
    gtv_volume_remove_cell(v, GTV_CELL(j->data)) ;

  return v ;
}

gint bsi_cell_indices(GtvCell *c, GtsVertex *x,
		      gint *i1, gint *i2, gint *i3, gint *i4)

{
  GtsVertex *v1, *v2, *v3, *v4 ;
  
  gtv_tetrahedron_vertices(GTV_TETRAHEDRON(c), &v1, &v2, &v3, &v4) ;

  *i1 = v1 - &(x[0]) ;
  *i2 = v2 - &(x[0]) ;
  *i3 = v3 - &(x[0]) ;
  *i4 = v4 - &(x[0]) ;

  g_assert(&(x[(*i1)]) == v1) ;
  g_assert(&(x[(*i2)]) == v2) ;
  g_assert(&(x[(*i3)]) == v3) ;
  g_assert(&(x[(*i4)]) == v4) ;
  
  return 0 ;
}

gint bsi_points_origin_width(GtsVertex *x, gint nx,
			     gdouble *xmin, gdouble *xmax,
			     gdouble *D, gboolean init)

{
  gint i ;

  if ( init ) {
    xmin[0] = xmin[1] = xmin[2] =  G_MAXDOUBLE ;
    xmax[0] = xmax[1] = xmax[2] = -G_MAXDOUBLE ;
  }

  for ( i = 0 ; i < nx ; i ++ ) {
    xmin[0] = MIN(xmin[0], GTS_POINT(&(x[i]))->x) ;
    xmin[1] = MIN(xmin[1], GTS_POINT(&(x[i]))->y) ;
    xmin[2] = MIN(xmin[2], GTS_POINT(&(x[i]))->z) ;
    xmax[0] = MAX(xmax[0], GTS_POINT(&(x[i]))->x) ;
    xmax[1] = MAX(xmax[1], GTS_POINT(&(x[i]))->y) ;
    xmax[2] = MAX(xmax[2], GTS_POINT(&(x[i]))->z) ;
  }    

  *D = xmax[0] - xmin[0] ;
  *D = MAX(*D, xmax[1] - xmin[1]) ;
  *D = MAX(*D, xmax[2] - xmin[2]) ;
  
  return 0 ;
}
