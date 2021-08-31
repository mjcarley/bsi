#ifndef __DUFFY_H_INCLUDED__
#define __DUFFY_H_INCLUDED__

#define duffy_vector3(_c,_a,_b)			\
  do {						\
    (_c)[0] = (_b)[0] - (_a)[0] ;		\
    (_c)[1] = (_b)[1] - (_a)[1] ;		\
    (_c)[2] = (_b)[2] - (_a)[2] ;		\
  } while (0) 

#define duffy_length3sqr(_a)			\
  (((_a)[0])*((_a)[0]) + ((_a)[1])*((_a)[1]) + ((_a)[2])*((_a)[2]))

#define duffy_diff_scalar3(_a,_b,_c)				\
  (((_a)[0])*((_c)[0]-(_b)[0]) +				\
   ((_a)[1])*((_c)[1]-(_b)[1]) +				\
   ((_a)[2])*((_c)[2]-(_b)[2]))

#define duffy_scalar3(_a,_b)			\
  ((_a)[0]*(_b)[0] + (_a)[1]*(_b)[1] + (_a)[2]*(_b)[2])

#define duffy_cross3(_c,_a,_b)				\
  do {							\
    (_c)[0] = (_a)[1]*(_b)[2] - (_a)[2]*(_b)[1]  ;	\
    (_c)[1] = (_a)[2]*(_b)[0] - (_a)[0]*(_b)[2]  ;	\
    (_c)[2] = (_a)[0]*(_b)[1] - (_a)[1]*(_b)[0]  ;	\
  } while (0) 

#define duffy_distance3sqr(_a,_b)		\
  (((_a)[0] - (_b)[0])*((_a)[0] - (_b)[0]) +	\
   ((_a)[1] - (_b)[1])*((_a)[1] - (_b)[1]) +	\
   ((_a)[2] - (_b)[2])*((_a)[2] - (_b)[2]))

#define SIGN(_x) ((_x) < 0.0 ? -1 : ((_x) > 0.0 ? 1 : 0))

typedef gint (*duffy_func_t)(gdouble *y, gdouble *x, gdouble wt,
			     gdouble *q, gint nq, gpointer data) ;

gint duffy_ref_tet_quad_plain(gdouble *x1, gdouble *x2,
			      gdouble *x3, gdouble *x4,
			      gdouble V,
			      gdouble al, gdouble bt,
			      gdouble *qu, gint ustr,
			      gdouble *wu, gint wustr, gint nu,
			      gdouble *qv, gint vstr,
			      gdouble *wv, gint wvstr, gint nv,
			      gdouble *qw, gint wstr,
			      gdouble *ww, gint wwstr, gint nw,
			      duffy_func_t func, gpointer data,
			      gdouble *q, gint nq) ;
gint duffy_ref_tet_quad_ch(gdouble *x1, gdouble *x2,
			   gdouble *x3, gdouble *x4,
			   gdouble V,
			   gdouble al, gdouble bt,
			   gdouble *qu, gint ustr,
			   gdouble *wu, gint wustr, gint nu,
			   gdouble *qv, gint vstr,
			   gdouble *wv, gint wvstr, gint nv,
			   gdouble *qw, gint wstr,
			   gdouble *ww, gint wwstr, gint nw,
			   duffy_func_t func, gpointer data,
			   gdouble *q, gint nq) ;
gint duffy_tet_quad(gdouble *x1, gdouble *x2,
		    gdouble *x3, gdouble *x4,
		    gdouble al, gdouble bt,
		    gdouble *qu, gint ustr, gdouble *wu, gint wustr, gint nu,
		    gdouble *qv, gint vstr, gdouble *wv, gint wvstr, gint nv,
		    gdouble *qw, gint wstr, gdouble *ww, gint wwstr, gint nw,
		    duffy_func_t func, gpointer data,
		    gdouble *q, gint nq) ;

gint duffy_tet_quad_general(gdouble *x1, gdouble *x2,
			    gdouble *x3, gdouble *x4,
			    gdouble *x,
			    gdouble al,  gdouble bt,
			    gdouble *qu, gint ustr,
			    gdouble *wu, gint wustr, gint nu,
			    gdouble *qv, gint vstr,
			    gdouble *wv, gint wvstr, gint nv,
			    gdouble *qw, gint wstr,
			    gdouble *ww, gint wwstr, gint nw,
			    duffy_func_t func, gpointer data,
			    gdouble *q, gint nq) ;

#endif /*__DUFFY_H_INCLUDED__*/
