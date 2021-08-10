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
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

#include <glib.h>

#include <bsi.h>

gchar *progname ;

static gint init_box(GtsVertex *nodes,
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

gint main(gint argc, gchar **argv)

{
  gdouble xmin, xmax, ymin, ymax, zmin, zmax, *data, tol ;
  gint nnmax, nn, nd, nx, ny, nz, i, j ;
  GtsVertex *nodes ;
  GtvVolume *v ;
  FILE *output ;
  gboolean write_tetgen ;
  gchar ch ;
  
  progname = g_strdup(g_path_get_basename(argv[0])) ;

  nd = 3 ; nnmax = 32768*16 ;
  nx = 33 ; ny = 33 ; nz = 9 ;
  tol = 1e-2 ;

  write_tetgen = FALSE ;
  
  output = stdout ;
  
  while ( (ch = getopt(argc, argv, "Tt:")) != EOF ) {
    switch ( ch ) {
    default: g_assert_not_reached() ; break ;
    case 't': tol = atof(optarg) ; break ;
    case 'T': write_tetgen = TRUE ; break ;
    }
  }
  
  nodes = (GtsVertex *)g_malloc0(nnmax*sizeof(GtsVertex)) ;
  data  = (gdouble *)g_malloc0(nnmax*nd*sizeof(gdouble)) ;

  xmin = -2.0 ; xmax = 2.0 ;
  ymin = -2.0 ; ymax = 2.0 ;
  zmin = -1.0 ; zmax = 1.0 ;

  /* xmin = -1.5 ; xmax = 1.5 ;  */
  /* ymin = -1.5 ; ymax = 1.5 ;  */
  /* zmin = -0.5 ; zmax = 0.5 ;  */
  
  nn = init_box(nodes, xmin, xmax, nx, ymin, ymax, ny, zmin, zmax, nz,
		nnmax) ;

  v = bsi_make_delaunay_volume(nodes, nn) ;

  fprintf(stderr, "%s: refining node distribution\n", progname) ;
  nn = bsi_distribution_refine(v, nodes, data, nd, nn, nnmax,
			       bsi_source_func_ring_gaussian, NULL, tol) ;
  fprintf(stderr, "%s: distribution refined to %d nodes\n",
	  progname, nn) ;

  if ( !write_tetgen ) {
    fprintf(output, "%d %d\n", nn, nd) ;

    for ( i = 0 ; i < nn ; i ++ ) {
      fprintf(output, "%1.16e %1.16e %1.16e",
	      GTS_POINT(&(nodes[i]))->x,
	      GTS_POINT(&(nodes[i]))->y,
	      GTS_POINT(&(nodes[i]))->z) ;
      for ( j = 0 ; j < nd ; j ++ ) {
	fprintf(output, " %1.16e", data[i*nd+j]) ;
      }
      fprintf(output, "\n") ;
    }

    return 0 ;
  }

  fprintf(output, "%d %d 0 0\n", nn, nd) ;

  for ( i = 0 ; i < nn ; i ++ ) {
    fprintf(output, "%d %1.16e %1.16e %1.16e",
	    i,
	    GTS_POINT(&(nodes[i]))->x,
	    GTS_POINT(&(nodes[i]))->y,
	    GTS_POINT(&(nodes[i]))->z) ;
    for ( j = 0 ; j < nd ; j ++ ) {
      fprintf(output, " %1.16e", data[i*nd+j]) ;
    }
    fprintf(output, "\n") ;
  }

    
  
  /* gtv_volume_write(v, output) ; */
  
  return 0 ;
}
