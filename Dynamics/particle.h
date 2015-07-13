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
  
  tVertex vertices;
  tEdge edges;
  tFace faces;
} PARTICLE;













