typedef struct particle{
  int iOrder, iAggOrder, iColor;
  double dMass;
  double dRadius;
  double dKn;
  double dKt;
  double dEpsN;
  double dEpsT;
  double dMuS;
  //Vec vNhat; //unit vector of relative velocity between particle and neighbor
  //Vec vOrient;
  //Mat mInertia;
  Vec r;
  Vec v;
  Vec w;
  Vec a;
  Vec vnOld;
  Vec vShear;
  
  //IFDEF POLY
  double dDensity; //Need to define density? or can we get this from the sphere equivalent
  Mat mInertia;
  tVertex vertices;
  tEdge edges;
  tFace faces;
  int iNumVertices;
  //struct couple *cpl; //Don't know if this is necessary _DEBUG_ -RB
} PARTICLE;


  /*SAT structures*/
typedef struct couple {
  /*This is a structure of particle p, so really only need to carry neighbor?*/
  PARTICLE *p;	/* First polyhedron of collision test. */
  PARTICLE *pn;	/* Second polyhedron of collision test. */
  bool plane_exists;	/* prospective separating plane flag */
  double pln_pnt1[3];	/* 1st point used to form separating plane. */
  double pln_pnt2[3];	/* 2nd point used to form separating plane. */
  int vert_indx[4][2]; /* cached points for distance algorithm. */
  int n;		/* number of cached points, if any. */
  } *COUPLE;











