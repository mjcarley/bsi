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

#include <glib.h>

#include <bsi.h>

gchar *progname ;

static gint read_nodes(FILE *f, GtsVertex **nodes, gdouble **data,
		       gint *nnodes, gint *nd)

{
  gint i, j ;
  gdouble x, y, z ;
  GtsVertex *v ;

  fscanf(f, "%d", nnodes) ;
  fscanf(f, "%d", nd) ;

  *nodes = (GtsVertex *)g_malloc0((*nnodes)*sizeof(GtsVertex)) ;
  *data  = (gdouble   *)g_malloc0((*nnodes)*(*nd)*sizeof(gdouble)) ;

  for ( i = 0 ; i < *nnodes ; i ++ ) {
    fscanf(f, "%lg %lg %lg", &x, &y, &z) ;
    v = gts_vertex_new(gts_vertex_class(), x, y, z) ;
    memcpy(&((*nodes)[i]), v, sizeof(GtsVertex)) ;
    gts_object_destroy(v) ;
    for ( j = 0 ; j < (*nd) ; j ++ )
      fscanf(f, "%lg", &((*data)[i*(*nd)+j])) ;
  }
  
  return 0 ;
}

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

gint main(gint argc, gchar **argv)

{
  FILE *input, *output ;
  GtvVolume *v ;
  GtsVertex *nodes ;
  gdouble *data ;
  gint nnodes, nd ;
  
  progname = g_strdup(g_path_get_basename(argv[0])) ;

  input = stdin ;
  output = stdout ;

  fprintf(stderr, "reading nodes\n") ;
  read_nodes(input, &nodes, &data, &nnodes, &nd) ;

  fprintf(stderr, "tetrahedralizing nodes\n") ;
  v = make_delaunay_volume(nodes, nnodes) ;
  
  /*data lookup test*/
  /* gint i, idx ; */
  /* for ( i = 0 ; i < nnodes ; i ++ ) { */
  /*   idx = &(nodes[i]) - &(nodes[0]) ; */
  /*   fprintf(stderr, "%d %d\n", i, idx) ; */
  /* } */
       
  fprintf(stderr, "writing volume\n") ;
  gtv_volume_write(v, output) ;

  return 0 ;
}
