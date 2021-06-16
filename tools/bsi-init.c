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

gchar *progname ;

static GtvVolume *make_delaunay_volume(GtsVertex *nodes, gint nv)

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

  for ( i = 0 ; i < nv ; i ++ ) {
    p = &(nodes[i]) ;
    if ( gtv_delaunay_add_vertex(v, p, NULL) == GTV_VERTEX_NOT_IN_VOLUME )
      fprintf(stderr,
	      "%s: vertex (%lg,%lg,%lg) not inside convex hull\n",
	      progname,
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

gint init_box(GtsVertex *nodes,
	      gdouble xmin, gdouble xmax, gint nx,
	      gdouble ymin, gdouble ymax, gint ny,
	      gdouble zmin, gdouble zmax, gint nz,
	      gint nnmax)
{
  gint nn, ix, iy, iz ;
  gdouble x, y, z ;
  GtsVertex *v ;
  
  nn = 0 ;

  for ( ix = 0 ; ix <= nx ; ix ++ ) {
    x = xmin + (xmax - xmin)*ix/nx ;
    for ( iy = 0 ; iy <= ny ; iy ++ ) {
      y = ymin + (ymax - ymin)*iy/ny ;
      for ( iz = 0 ; iz <= nz ; iz ++ ) {
	z = zmin + (zmax - zmin)*iz/nz ;
	v = gts_vertex_new(gts_vertex_class(), x, y, z) ;
	memcpy(&(nodes[nn]), v, sizeof(GtsVertex)) ;
	gts_object_destroy(GTS_OBJECT(v)) ;
	nn ++ ;
      }
    }
  }
  
  return nn ;
}

gint init_func(GtsVertex *v, gdouble *f, gint nf, gpointer data)

{
  gdouble s, r, r2, r0, th, om ;

  s = 0.25 ; r0 = 0.9 ;
  
  r = sqrt(GTS_POINT(v)->x*GTS_POINT(v)->x +
	   GTS_POINT(v)->y*GTS_POINT(v)->y) ;
  th = atan2(GTS_POINT(v)->y, GTS_POINT(v)->x) ;

  r2 = (r-r0)*(r-r0) + GTS_POINT(v)->z*GTS_POINT(v)->z ;

  om = exp(-r2/s/s) ;

  f[0] = -om*sin(th) ; 
  f[1] =  om*cos(th) ; 
  f[2] =  0.0 ;
  
  return 0 ;
}

static void refine_cell(GtvCell *c, gpointer *rdata)

{
  GtsVertex *nodes = rdata[0] ;
  gdouble *data = rdata[1] ;
  gint nd = *(gint *)rdata[2] ;
  gint *nn = rdata[3] ;
  gint nnmax = *(gint *)rdata[4] ;
  gdouble tol = *(gdouble *)rdata[5] ;
  gpointer idata = rdata[6] ;
  GtsVertex *v1, *v2, *v3, *v4, *vt ;
  gdouble fi[8], xt, yt, zt ;
  gint i1, i2, i3, i4, i ;
  
  if ( *nn >= nnmax ) return ;
  
  gtv_tetrahedron_vertices(GTV_TETRAHEDRON(c), &v1, &v2, &v3, &v4) ;

  i1 = v1 - &(nodes[0]) ;
  i2 = v2 - &(nodes[0]) ;
  i3 = v3 - &(nodes[0]) ;
  i4 = v4 - &(nodes[0]) ;

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

  init_func(vt, &(data[(*nn)*nd]), nd, idata) ;

  for ( i = 0 ; i < nd ; i ++ ) {
    if ( fabs(data[(*nn)*nd+i] - fi[i]) > tol ) {
      memcpy(&(nodes[(*nn)]), vt, sizeof(GtsVertex)) ;
      (*nn) ++ ;
      gts_object_destroy(GTS_OBJECT(vt)) ;
      return ;
    }
  }
  
  gts_object_destroy(GTS_OBJECT(vt)) ;

  return ;
}

gint refine_nodes(GtvVolume *v, GtsVertex *nodes, gdouble *data, gint nd,
		  gint nn, gint nnmax, gdouble tol, gpointer idata)

{
  gint i ;
  gpointer rdata[8] ;
  gboolean stop ;
  
  for ( i = 0 ; i < nn ; i ++ ) {
    init_func(&(nodes[i]), &(data[i*nd]), nd, idata) ;
  }

  rdata[0] = nodes ;
  rdata[1] = data ;
  rdata[2] = &nd ;
  rdata[3] = &nn ;
  rdata[4] = &nnmax ;
  rdata[5] = &tol ;
  rdata[6] = idata ;
  
  stop = FALSE ;
  while ( stop == FALSE ) {
    i = nn ;

    gtv_volume_foreach_cell(v, (GtsFunc)refine_cell, rdata) ;
    
    if ( i == nn ) stop = TRUE ;
    if ( nn == nnmax ) stop = TRUE ;

    for ( ; i < nn ; i ++ ) {
      gtv_delaunay_add_vertex(v, &(nodes[i]), NULL) ;
    }
  }
  
  return nn ;
}

gint main(gint argc, gchar **argv)

{
  gdouble xmin, xmax, ymin, ymax, zmin, zmax, *data, tol ;
  gint nnmax, nn, nd, nx, ny, nz ;
  GtsVertex *nodes ;
  GtvVolume *v ;
  FILE *output ;

  progname = g_strdup(g_path_get_basename(argv[0])) ;
  
  nd = 3 ; nnmax = 32768*8 ;
  nx = 9 ; ny = 9 ; nz = 5 ;
  tol = 1e-2 ;

  output = stdout ;
  
  nodes = (GtsVertex *)g_malloc0(nnmax*sizeof(GtsVertex)) ;
  data  = (gdouble *)g_malloc0(nnmax*nd*sizeof(gdouble)) ;

  xmin = -2.0 ; xmax = 2.0 ; 
  ymin = -2.0 ; ymax = 2.0 ; 
  zmin = -1.0 ; zmax = 1.0 ; 
  
  nn = init_box(nodes, xmin, xmax, nx, ymin, ymax, ny, zmin, zmax, nz,
		nnmax) ;

  v = make_delaunay_volume(nodes, nn) ;

  fprintf(stderr, "%s: refining node distribution\n", progname) ;
  nn = refine_nodes(v, nodes, data, nd, nn, nnmax, tol, NULL) ;
  fprintf(stderr, "%s: distribution refined to %d nodes\n",
	  progname, nn) ;

  gtv_volume_write(v, output) ;
  
  return 0 ;
}
