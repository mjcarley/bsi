/* This file is part of Tetquad, a set of functions for evaluating
 * integrals on tetrahedra.
 *
 * Copyright (C) 2021 Michael Carley
 *
 * Tetquad is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version. Tetquad is distributed in the
 * hope that it will be useful, but WITHOUT ANY WARRANTY; without even
 * the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 * PURPOSE.  See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Tetquad.  If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef _TETQUAD_H_INCLUDED_
#define _TETQUAD_H_INCLUDED_

#include <glib.h>

#define BETWEEN(_a,_b,_c)					\
  (((_a)<=(_b) && (_a)>=(_c)) || ((_a)>=(_b) && (_a)<=(_c)))

#define tq_vector_init(_c,_a,_b)			       \
  ((_c)[0] = (_b)[0] - (_a)[0],				       \
   (_c)[1] = (_b)[1] - (_a)[1],				       \
   (_c)[2] = (_b)[2] - (_a)[2])

#define tq_vector_scalar(_a,_b)				\
  ((_a)[0]*(_b)[0] + (_a)[1]*(_b)[1] + (_a)[2]*(_b)[2])

#define tq_vector_cross(_c,_a,_b)			\
  do {							\
 (_c)[0] = (_a)[1]*(_b)[2] - (_a)[2]*(_b)[1]  ;		\
 (_c)[1] = (_a)[2]*(_b)[0] - (_a)[0]*(_b)[2]  ;		\
 (_c)[2] = (_a)[0]*(_b)[1] - (_a)[1]*(_b)[0]  ;		\
  } while (0) 

#define tq_vector_length2(_a) (tq_vector_scalar((_a),(_a)))

#define tq_vector_distance2(_a,_b)		\
  (((_a)[0] - (_b)[0])*((_a)[0] - (_b)[0]) +	\
   ((_a)[1] - (_b)[1])*((_a)[1] - (_b)[1]) +	\
   ((_a)[2] - (_b)[2])*((_a)[2] - (_b)[2]))

/*a.(c-b)*/
#define tq_vector_diff_scalar(_a,_b,_c)				\
  (((_a)[0])*((_c)[0]-(_b)[0]) +				\
   ((_a)[1])*((_c)[1]-(_b)[1]) +				\
   ((_a)[2])*((_c)[2]-(_b)[2]))


typedef gint (*tq_tetquad_func_t)(gdouble *y, gdouble *s,
				  gdouble R, gdouble wt,
				  gdouble *q, gint nq, gpointer data) ;

gint tq_tet_quad(gdouble *x1, gdouble *x2, gdouble *x3, gdouble *x4,
		 gdouble *qph, gint qpstr, gdouble *wph, gint wpstr,
		 gint nph,
		 gdouble *qth, gint qtstr, gdouble *wth, gint wtstr,
		 gint nth,
		 gdouble *qr, gint qrstr, gdouble *wr, gint wrstr,
		 gint nr,
		 tq_tetquad_func_t qfunc, gpointer qdata,		     
		 gdouble *q, gint nq) ;
gint tq_transform_matrix(gdouble *x, gdouble *A) ;

gint tq_tet_dihedral_angles(gdouble *xa, gdouble *xb, gdouble *xc,
			    gdouble *oa, gdouble *ob, gdouble *oc,
			    gdouble *ab, gdouble *bc, gdouble *ca) ;

#endif /*_TETQUAD_H_INCLUDED_*/
