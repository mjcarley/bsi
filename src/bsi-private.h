#ifndef __BSI_PRIVATE_H_INCLUDED__
#define __BSI_PRIVATE_H_INCLUDED__

#define REAL gdouble
#define bool gboolean

// A "polygon" describes a simple polygon (no holes). It is not necessarily
//   convex. Each polygon contains a number of corners (points) and the same
//   number of sides (edges).  The points of the polygon must be given in
//   either counterclockwise or clockwise order and they form a ring, so 
//   every two consecutive points forms an edge of the polygon.
typedef struct {
  int *vertexlist;
  int numberofvertices;
} polygon;

//A "facet" describes a polygonal region possibly with holes, edges, and 
//  points floating in it.  Each facet consists of a list of polygons and
//  a list of hole points (which lie strictly inside holes).
typedef struct {
  polygon *polygonlist;
  int numberofpolygons;
  REAL *holelist;
  int numberofholes;
} facet;

// A "voroedge" is an edge of the Voronoi diagram. It corresponds to a
//   Delaunay face.  Each voroedge is either a line segment connecting
//   two Voronoi vertices or a ray starting from a Voronoi vertex to an
//   "infinite vertex".  'v1' and 'v2' are two indices pointing to the
//   list of Voronoi vertices. 'v1' must be non-negative, while 'v2' may
//   be -1 if it is a ray, in this case, the unit normal of this ray is
//   given in 'vnormal'.
typedef struct {
  int v1, v2;
  REAL vnormal[3];
} voroedge;

// A "vorofacet" is an facet of the Voronoi diagram. It corresponds to a
//   Delaunay edge.  Each Voronoi facet is a convex polygon formed by a
//   list of Voronoi edges, it may not be closed.  'c1' and 'c2' are two
//   indices pointing into the list of Voronoi cells, i.e., the two cells
//   share this facet.  'elist' is an array of indices pointing into the
//   list of Voronoi edges, 'elist[0]' saves the number of Voronoi edges
//   (including rays) of this facet.
typedef struct {
  int c1, c2;
  int *elist;
} vorofacet;


// Additional parameters associated with an input (or mesh) vertex.
//   These informations are provided by CAD libraries.
typedef struct {
  REAL uv[2];
  int tag;
  int type; // 0, 1, or 2.
} pointparam;

// Callback functions for meshing PSCs.
typedef REAL (* GetVertexParamOnEdge)(void*, int, int);
typedef void (* GetSteinerOnEdge)(void*, int, REAL, REAL*);
typedef void (* GetVertexParamOnFace)(void*, int, int, REAL*);
typedef void (* GetEdgeSteinerParamOnFace)(void*, int, REAL, int, REAL*);
typedef void (* GetSteinerOnFace)(void*, int, REAL*, REAL*);

// A callback function for mesh refinement.
typedef bool (* TetSizeFunc)(REAL*, REAL*, REAL*, REAL*, REAL*, REAL);


typedef struct {
  // Items are numbered starting from 'firstnumber' (0 or 1),
  // default is 0.
  int firstnumber;

  // Dimension of the mesh (2 or 3), default is 3.
  int mesh_dim;

  // Does the lines in .node file contain index or not, default is 1.
  int useindex;

  // 'pointlist':  An array of point coordinates.  The first point's x
  //   coordinate is at index [0] and its y coordinate at index [1], its
  //   z coordinate is at index [2], followed by the coordinates of the
  //   remaining points.  Each point occupies three REALs.
  // 'pointattributelist':  An array of point attributes.  Each point's
  //   attributes occupy 'numberofpointattributes' REALs.
  // 'pointmtrlist': An array of metric tensors at points. Each point's
  //   tensor occupies 'numberofpointmtr' REALs.
  // 'pointmarkerlist':  An array of point markers; one integer per point.
  // 'point2tetlist': An array of tetrahedra indices; one integer per point.
  REAL *pointlist;
  REAL *pointattributelist;
  REAL *pointmtrlist;
  int  *pointmarkerlist;
  int  *point2tetlist;
  pointparam *pointparamlist;
  int numberofpoints;
  int numberofpointattributes;
  int numberofpointmtrs;
 
  // 'tetrahedronlist':  An array of tetrahedron corners.  The first
  //   tetrahedron's first corner is at index [0], followed by its other
  //   corners, followed by six nodes on the edges of the tetrahedron if the
  //   second order option (-o2) is applied. Each tetrahedron occupies
  //   'numberofcorners' ints.  The second order nodes are ouput only.
  // 'tetrahedronattributelist':  An array of tetrahedron attributes.  Each
  //   tetrahedron's attributes occupy 'numberoftetrahedronattributes' REALs.
  // 'tetrahedronvolumelist':  An array of constraints, i.e. tetrahedron's
  //   volume; one REAL per element.  Input only.
  // 'neighborlist':  An array of tetrahedron neighbors; 4 ints per element.
  // 'tet2facelist':  An array of tetrahedron face indices; 4 ints per element.
  // 'tet2edgelist':  An array of tetrahedron edge indices; 6 ints per element.
  int  *tetrahedronlist;
  REAL *tetrahedronattributelist;
  REAL *tetrahedronvolumelist;
  int  *neighborlist;
  int  *tet2facelist;
  int  *tet2edgelist;
  int numberoftetrahedra;
  int numberofcorners;
  int numberoftetrahedronattributes;

  // 'facetlist':  An array of facets.  Each entry is a structure of facet.
  // 'facetmarkerlist':  An array of facet markers; one int per facet.
  facet *facetlist;
  int *facetmarkerlist;
  int numberoffacets;

  // 'holelist':  An array of holes (in volume).  Each hole is given by a
  //   seed (point) which lies strictly inside it. The first seed's x, y and z
  //   coordinates are at indices [0], [1] and [2], followed by the
  //   remaining seeds.  Three REALs per hole.
  REAL *holelist;
  int numberofholes;

  // 'regionlist': An array of regions (subdomains).  Each region is given by
  //   a seed (point) which lies strictly inside it. The first seed's x, y and
  //   z coordinates are at indices [0], [1] and [2], followed by the regional
  //   attribute at index [3], followed by the maximum volume at index [4].
  //   Five REALs per region.
  // Note that each regional attribute is used only if you select the 'A'
  //   switch, and each volume constraint is used only if you select the
  //   'a' switch (with no number following).
  REAL *regionlist;
  int numberofregions;

  // 'refine_elem_list': An array of tetrahedra to be refined.  The first
  //   tetrahedron's first corner is at index [0], followed by its other
  //   corners. Four integers per element.
  // 'refine_elem_vol_list':  An array of constraints, i.e. tetrahedron's
  //   volume; one REAL per element.
  int  *refine_elem_list;
  REAL *refine_elem_vol_list;
  int  numberofrefineelems;

  // 'facetconstraintlist':  An array of facet constraints.  Each constraint
  //   specifies a maximum area bound on the subfaces of that facet.  The
  //   first facet constraint is given by a facet marker at index [0] and its
  //   maximum area bound at index [1], followed by the remaining facet con-
  //   straints. Two REALs per facet constraint.  Note: the facet marker is
  //   actually an integer.
  REAL *facetconstraintlist;
  int numberoffacetconstraints;

  // 'segmentconstraintlist': An array of segment constraints. Each constraint
  //   specifies a maximum length bound on the subsegments of that segment.
  //   The first constraint is given by the two endpoints of the segment at
  //   index [0] and [1], and the maximum length bound at index [2], followed
  //   by the remaining segment constraints.  Three REALs per constraint.
  //   Note the segment endpoints are actually integers.
  REAL *segmentconstraintlist;
  int numberofsegmentconstraints;


  // 'trifacelist':  An array of face (triangle) corners.  The first face's
  //   three corners are at indices [0], [1] and [2], followed by the remaining
  //   faces.  Three ints per face.
  // 'trifacemarkerlist':  An array of face markers; one int per face.
  // 'o2facelist':  An array of second order nodes (on the edges) of the face.
  //   It is output only if the second order option (-o2) is applied. The
  //   first face's three second order nodes are at [0], [1], and [2],
  //   followed by the remaining faces.  Three ints per face.
  // 'face2tetlist':  An array of tetrahedra indices; 2 ints per face.
  // 'face2edgelist':  An array of edge indices; 3 ints per face.
  int *trifacelist;
  int *trifacemarkerlist;
  int *o2facelist;
  int *face2tetlist;
  int *face2edgelist;
  int numberoftrifaces;

  // 'edgelist':  An array of edge endpoints.  The first edge's endpoints
  //   are at indices [0] and [1], followed by the remaining edges.
  //   Two ints per edge.
  // 'edgemarkerlist':  An array of edge markers; one int per edge.
  // 'o2edgelist':  An array of midpoints of edges. It is output only if the
  //   second order option (-o2) is applied. One int per edge.
  // 'edge2tetlist':  An array of tetrahedra indices.  One int per edge.
  int *edgelist;
  int *edgemarkerlist;
  int *o2edgelist;
  int *edge2tetlist;
  int numberofedges;

  // 'vpointlist':  An array of Voronoi vertex coordinates (like pointlist).
  // 'vedgelist':  An array of Voronoi edges.  Each entry is a 'voroedge'.
  // 'vfacetlist':  An array of Voronoi facets. Each entry is a 'vorofacet'.
  // 'vcelllist':  An array of Voronoi cells.  Each entry is an array of
  //   indices pointing into 'vfacetlist'. The 0th entry is used to store
  //   the length of this array.
  REAL *vpointlist;
  voroedge *vedgelist;
  vorofacet *vfacetlist;
  int **vcelllist;
  int numberofvpoints;
  int numberofvedges;
  int numberofvfacets;
  int numberofvcells;


  // Variable (and callback functions) for meshing PSCs.
  void *geomhandle;
  GetVertexParamOnEdge getvertexparamonedge;
  GetSteinerOnEdge getsteineronedge;
  GetVertexParamOnFace getvertexparamonface;
  GetEdgeSteinerParamOnFace getedgesteinerparamonface;
  GetSteinerOnFace getsteineronface;

  // A callback function.
  TetSizeFunc tetunsuitable;

  /* // Input & output routines. */
  /* bool load_node_call(FILE* infile, int markers, int uvflag, char*); */
  /* bool load_node(char*); */
  /* bool load_edge(char*); */
  /* bool load_face(char*); */
  /* bool load_tet(char*); */
  /* bool load_vol(char*); */
  /* bool load_var(char*); */
  /* bool load_mtr(char*); */
  /* bool load_elem(char*); */
  /* bool load_poly(char*); */
  /* bool load_off(char*); */
  /* bool load_ply(char*); */
  /* bool load_stl(char*); */
  /* bool load_vtk(char*); */
  /* bool load_medit(char*, int); */
  /* bool load_neumesh(char*, int); */
  /* bool load_plc(char*, int); */
  /* bool load_tetmesh(char*, int); */
  /* void save_nodes(const char*); */
  /* void save_elements(const char*); */
  /* void save_faces(const char*); */
  /* void save_edges(char*); */
  /* void save_neighbors(char*); */
  /* void save_poly(const char*); */
  /* void save_faces2smesh(char*); */

  /* // Read line and parse string functions. */
  /* char *readline(char* string, FILE* infile, int *linenumber); */
  /* char *findnextfield(char* string); */
  /* char *readnumberline(char* string, FILE* infile, char* infilename); */
  /* char *findnextnumber(char* string); */
  
}   tetwrap_tetgenio_t ; // class tetgenio

gint tetwrap_tetrahedralize(gchar *switches,
			    tetwrap_tetgenio_t *in,
			    tetwrap_tetgenio_t *out,
			    gpointer addin, gpointer bgmin) ;
tetwrap_tetgenio_t *tetwrap_new(gdouble *nodes, gint nnodes) ;
gint tetwrap_examine(gpointer ta) ;

#endif /*__BSI_PRIVATE_H_INCLUDED__*/
