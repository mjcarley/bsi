#include <glib.h>

#include "havoc-private.h"

#define PRE_CHECK_TET

/* #define FULL_COORDINATE_CHECK */

#define EPSILON 1e-15

static gboolean intersect_triangle(gdouble orig[3], gdouble dir[3],
				   gdouble vert0[3], gdouble vert1[3], 
				   gdouble vert2[3], gdouble *t, 
				   gdouble *u, gdouble *v)
     /*
       Taken from Moller and Trumbore, Journal of Graphics Tools,
       1997, 2(1):21--28, http://jgt.akpeters.com/papers/MollerTrumbore97/
      */

{
  gdouble edge1[3], edge2[3], tvec[3], pvec[3], qvec[3]; 
  gdouble det, inv_det;

  /* find vectors for two edges sharing vert0 */ 
  vector_init(edge1, vert1, vert0); vector_init(edge2, vert2, vert0);

  /* begin calculating determinant - also used to calculate U parameter */ 
  vector_cross(pvec, dir, edge2);

  /* if determinant is near zero, ray lies in plane of triangle */ 
  det = vector_scalar(edge1, pvec);

  /*this is the in-plane check*/
  if (det > -EPSILON && det < EPSILON) return FALSE; 

  inv_det = 1.0 / det ;
  /* calculate distance from vert0 to ray origin */ 
  vector_init(tvec, orig, vert0) ;

  /* calculate U parameter and test bounds */ 
  *u = vector_scalar(tvec, pvec) * inv_det; 
  if (*u < 0.0 || *u > 1.0) return FALSE ;

  /* prepare to test V parameter */ 
  vector_cross(qvec, tvec, edge1);

  /* calculate V parameter and test bounds */ 
  *v = vector_scalar(dir, qvec) * inv_det; 
  if (*v < 0.0 || *u + *v > 1.0) return FALSE  ;

  /* calculate t, ray intersects triangle */ 
  *t = vector_scalar(edge2, qvec) * inv_det; 

  return TRUE ;
}


static inline gint face_interp(gdouble *w0, gdouble *w1, 
			       gdouble *w2, gdouble *w3,
			       gint face, gdouble u, gdouble v, gdouble *f)

{
  gdouble w ;
  switch (face) {
  default: g_assert_not_reached() ; break ;
  case 0: w = 1.0 - u - v ;
    f[0] = w*w0[0] + u*w1[0] + v*w2[0] ;
    f[1] = w*w0[1] + u*w1[1] + v*w2[1] ;
    f[2] = w*w0[2] + u*w1[2] + v*w2[2] ;
    break ;
  case 1: w = 1.0 - u - v ;
    f[0] = w*w1[0] + u*w2[0] + v*w3[0] ;
    f[1] = w*w1[1] + u*w2[1] + v*w3[1] ;
    f[2] = w*w1[2] + u*w2[2] + v*w3[2] ;
    break ;
  case 2: w = 1.0 - u - v ;
    f[0] = w*w2[0] + u*w3[0] + v*w0[0] ;
    f[1] = w*w2[1] + u*w3[1] + v*w0[1] ;
    f[2] = w*w2[2] + u*w3[2] + v*w0[2] ;
    break ;
  case 3: w = 1.0 - u - v ;
    f[0] = w*w3[0] + u*w0[0] + v*w1[0] ;
    f[1] = w*w3[1] + u*w0[1] + v*w1[1] ;
    f[2] = w*w3[2] + u*w0[2] + v*w1[2] ;
    break ;
  }

  return 0 ;
}

static gboolean tet_ray_cut(gdouble *x0, gdouble *x1, 
			    gdouble *x2, gdouble *x3,
			    gdouble *w0, gdouble *w1, 
			    gdouble *w2, gdouble *w3,
			    gdouble *x, gdouble *s,
			    gdouble h,
			    gdouble *R0, gdouble *R1,
			    gdouble *f0, gdouble *f1)

{
  gdouble u[4], v[4], t[4] ;
  gint face[4], i ;

#ifdef PRE_CHECK_TET
  gdouble r[3], R2, rds ;

/*   vector_init(r,x0,x) ; rds = vector_scalar(r,s) ; */
/*   R2 = r[0]*r[0] + r[1]*r[1] + r[2]*r[2] ; */
/*   if ( rds < 0 ) return FALSE ; */
/*   if ( rds*rds < (R2 - 2.0*h) && R2 > 0.1 ) return FALSE ; */
#endif /*PRE_CHECK_TET*/

  i = 0 ; 
#ifdef FULL_COORDINATE_CHECK
  /*if we want to check coordinate by coordinate*/
  if ( x[0] == x0[0] &&  x[1] == x0[1] &&  x[2] == x0[2] ) {
#else
  if ( x == x0 ) {
#endif /*FULL_COORDINATE_CHECK*/
    if ( intersect_triangle(x, s, x1, x2, x3, t, u, v) ) {
      if ( *t > 0 ) {
	*R0 = 0.0 ; *R1 = *t ; 
	f0[0] = w0[0] ; f0[1] = w0[1] ; f0[2] = w0[2] ;
	face_interp(w0, w1, w2, w3, 1, *u, *v, f1) ;
      } else {
	*R1 = 0.0 ; *R0 = *t ; 
	f1[0] = w0[0] ; f1[1] = w0[1] ; f1[2] = w0[2] ;
	face_interp(w0, w1, w2, w3, 1, *u, *v, f0) ;
      }      
      return TRUE ;
    }
    return FALSE ;
  }

#ifdef FULL_COORDINATE_CHECK
  /*if we want to check coordinate by coordinate*/
  if ( x[0] == x1[0] &&  x[1] == x1[1] &&  x[2] == x1[2] ) {
#else
  if ( x == x1 ) {
#endif /*FULL_COORDINATE_CHECK*/
    if ( intersect_triangle(x, s, x2, x3, x0, t, u, v) ) {
      if ( *t > 0 ) {
	*R0 = 0.0 ; *R1 = *t ; 
	f0[0] = w1[0] ; f0[1] = w1[1] ; f0[2] = w1[2] ;
	face_interp(w0, w1, w2, w3, 2, *u, *v, f1) ;
      } else {
	*R1 = 0.0 ; *R0 = *t ; 
	f1[0] = w1[0] ; f1[1] = w1[1] ; f1[2] = w1[2] ;
	face_interp(w0, w1, w2, w3, 2, *u, *v, f0) ;
      }      
      return TRUE ;
    }
    return FALSE ;
  }

#ifdef FULL_COORDINATE_CHECK
  /*if we want to check coordinate by coordinate*/
  if ( x[0] == x2[0] &&  x[1] == x2[1] &&  x[2] == x2[2] ) {
#else
  if ( x == x2 ) {
#endif /*FULL_COORDINATE_CHECK*/
    if ( intersect_triangle(x, s, x3, x0, x1, t, u, v) ) {
      if ( *t > 0 ) {
	*R0 = 0.0 ; *R1 = *t ; 
	f0[0] = w2[0] ; f0[1] = w2[1] ; f0[2] = w2[2] ;
	face_interp(w0, w1, w2, w3, 3, *u, *v, f1) ;
      } else {
	*R1 = 0.0 ; *R0 = *t ; 
	f1[0] = w2[0] ; f1[1] = w2[1] ; f1[2] = w2[2] ;
	face_interp(w0, w1, w2, w3, 3, *u, *v, f0) ;
      }      
      return TRUE ;
    }
    return FALSE ;
  }

#ifdef FULL_COORDINATE_CHECK
  /*if we want to check coordinate by coordinate*/
  if ( x[0] == x3[0] &&  x[1] == x3[1] &&  x[2] == x3[2] ) {
#else
    if ( x == x3 ) {
#endif /*FULL_COORDINATE_CHECK*/
    if ( intersect_triangle(x, s, x0, x1, x2, t, u, v) ) {
      if ( *t > 0 ) {
	*R0 = 0.0 ; *R1 = *t ; 
	f0[0] = w3[0] ; f0[1] = w3[1] ; f0[2] = w3[2] ;
	face_interp(w0, w1, w2, w3, 0, *u, *v, f1) ;
      } else {
	*R1 = 0.0 ; *R0 = *t ; 
	f1[0] = w3[0] ; f1[1] = w3[1] ; f1[2] = w3[2] ;
	face_interp(w0, w1, w2, w3, 0, *u, *v, f0) ;
      }      
      return TRUE ;
    }
    return FALSE ;
  }

  if ( intersect_triangle(x, s, x0, x1, x2, &t[i], &u[i], &v[i]) ) {
    face[i] = 0 ; i ++ ;
  }

  if ( intersect_triangle(x, s, x1, x2, x3, &t[i], &u[i], &v[i]) ) {
    face[i] = 1 ; i ++ ;
  }

  if ( intersect_triangle(x, s, x2, x3, x0, &t[i], &u[i], &v[i]) ) {
    face[i] = 2 ; i ++ ;
  }

  /*can't be a full intersection*/
  if ( i == 0 ) return FALSE ;

  if ( intersect_triangle(x, s, x3, x0, x1, &t[i], &u[i], &v[i]) ) {
    face[i] = 3 ; i ++ ;
  }

  /*only taking positive radii*/
/*   if ( t[0] < 0 || t[1] < 0 ) return FALSE ; */

  /*ray passes through one corner*/
  if ( i == 1 ) return FALSE ;
  
  if ( i == 3 ) {
    /*piercing a corner*/
    if ( t[0] == t[1] ) {
      t[0] = t[2] ; face[0] = face[2] ;
    }
    i = 2 ;
  }

  if ( t[0] == t[1] ) return FALSE ;
  g_assert(t[0] != t[1]) ;
  if ( i == 2 ) {
    if ( t[0] < t[1] ) {
      *R0 = t[0] ; *R1 = t[1] ;
      face_interp(w0, w1, w2, w3, face[0], u[0], v[0], f0) ;
      face_interp(w0, w1, w2, w3, face[1], u[1], v[1], f1) ;
    } else {
      *R0 = t[1] ; *R1 = t[0] ;
      face_interp(w0, w1, w2, w3, face[0], u[0], v[0], f1) ;
      face_interp(w0, w1, w2, w3, face[1], u[1], v[1], f0) ;
    }
    return TRUE ;
  }

  g_assert_not_reached() ;

  return TRUE ;
}


static gdouble tet_longest_side(gdouble *x0, gdouble *x1,
				gdouble *x2, gdouble *x3)

{
  gdouble lmax = 0.0 ;

  lmax = MAX(lmax,
	     (x0[0]-x1[0])*(x0[0]-x1[0]) +
	     (x0[1]-x1[1])*(x0[1]-x1[1]) +
	     (x0[2]-x1[2])*(x0[2]-x1[2])) ;

  lmax = MAX(lmax,
	     (x0[0]-x2[0])*(x0[0]-x2[0]) +
	     (x0[1]-x2[1])*(x0[1]-x2[1]) +
	     (x0[2]-x2[2])*(x0[2]-x2[2])) ;

  lmax = MAX(lmax,
	     (x0[0]-x3[0])*(x0[0]-x3[0]) +
	     (x0[1]-x3[1])*(x0[1]-x3[1]) +
	     (x0[2]-x3[2])*(x0[2]-x3[2])) ;

  lmax = MAX(lmax,
	     (x1[0]-x2[0])*(x1[0]-x2[0]) +
	     (x1[1]-x2[1])*(x1[1]-x2[1]) +
	     (x1[2]-x2[2])*(x1[2]-x2[2])) ;

  lmax = MAX(lmax,
	     (x1[0]-x3[0])*(x1[0]-x3[0]) +
	     (x1[1]-x3[1])*(x1[1]-x3[1]) +
	     (x1[2]-x3[2])*(x1[2]-x3[2])) ;

  lmax = MAX(lmax,
	     (x2[0]-x3[0])*(x2[0]-x3[0]) +
	     (x2[1]-x3[1])*(x2[1]-x3[1]) +
	     (x2[2]-x3[2])*(x2[2]-x3[2])) ;

  return lmax ;
}

gint havoc_tet_velocity(gdouble *x0, gdouble *x1,
			gdouble *x2, gdouble *x3,
			gdouble *w0, gdouble *w1,
			gdouble *w2, gdouble *w3,
			gint ns, gdouble *s, gdouble *wt,
			gint np, gdouble *x,
			gdouble *v, gint stride)

{
  gint i, j ;
  gdouble w[3], *sj, h ;
  gdouble R0, R1, f0[3], f1[3], df[3], dt, sgn ;

  h = tet_longest_side(x0, x1, x2, x3) ;

  for ( j = 0 ; j < ns ; j ++ ) {
    sj = &(s[3*j]) ;
    for ( i = 0 ; i < np ; i ++ ) {
      if ( tet_ray_cut(x0, x1, x2, x3, w0, w1, w2, w3,
		       &(x[3*i]), sj, h, &R0, &R1, f0, f1) ) {
	sgn = (R1 > 0.0 ? 1.0 : -1.0) ;
	w[0] = 0.5*(f0[0] + f1[0]) ;
	w[1] = 0.5*(f0[1] + f1[1]) ;
	w[2] = 0.5*(f0[2] + f1[2]) ;
	vector_cross(df,sj,w) ;
	dt = wt[j]*(R1-R0)*sgn ;
	v[i*stride+0] += df[0]*dt ; 
	v[i*stride+1] += df[1]*dt ; 
	v[i*stride+2] += df[2]*dt ;
      }
    }
  }

  return 0 ;
}
