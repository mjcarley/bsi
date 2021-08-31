#ifndef __TETWRAP_H_INCLUDED__
#define __TETWRAP_H_INCLUDED__

#include <glib.h>

extern "C" gint tetwrap_tetrahedralize(gchar *switches,
				       gpointer in, gpointer out,
				       gpointer addin, gpointer bgmin) ;
extern "C" gpointer tetwrap_new(gdouble *nodes, gint nnodes) ;
extern "C" gint tetwrap_examine(gpointer ta) ;

#endif /*__TETWRAP_H_INCLUDED__*/
