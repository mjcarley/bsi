#include <memory>

#include <glib.h>

#include <tetgen.h>

#include "tetwrap.h"

gint tetwrap_tetrahedralize(gchar *switches, gpointer in, gpointer out,
			    gpointer addin, gpointer bgmin)

{
  tetrahedralize("",
		 (tetgenio *)in,
		 (tetgenio *)out,
		 (tetgenio *)addin,
		 (tetgenio *)bgmin) ;
  
  return 0 ;
}

gpointer tetwrap_new(gdouble *nodes, gint nnodes)

{
  tetgenio *a = (new tetgenio())  ;

  a->pointlist = nodes ;
  a->numberofpoints = nnodes ;
  
  return a ;
}

gint tetwrap_examine(gpointer ta)

{
  tetgenio *a = (tetgenio *)ta ;

  fprintf(stderr, "Anything\n") ;
  
  return 0 ;
}
