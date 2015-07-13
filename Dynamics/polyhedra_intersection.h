/*Overlap Test*/
int bOverlapTestPoly(PARTICLE *p, PARTICLE *pn){
 
    int bOverlap = 0;
  
  
    Vec vRelDist; //relative distance
    double dRelDist;
  
    vecSub(p->r,pn->r,vRelDist);//Find relative Distance and Velocity
    dRelDist = vecMag(vRelDist);
  
    if ((p->dRadius+pn->dRadius)>dRelDist)
      bOverlap = 1;
    else
      bOverlap = 0;

    return bOverlap;
}


int WhichSide(PointSet S, PointD, Point P){
  //S vertices are projected to the for of P+t*D. Return value is +1 if all t > 0
  // -1 if all t < 0, 0 otherwise, in which case the line splits the polygon
  positive = 0; negative = 0;
  for (i=0; i< C.N; i++){ //C.N? number of vertices?
    t=Dot(D, S.V(i)-P)
    if (t > 0)
      positive++;
    else if (t < 0)
      negative++;
    if (positive && negative) return 0;
    }
  return (positive ? +1 : -1)
}

bool TestIntersection3d(ConvexPolyhedron C0, ConvexPolyhedron C1){

  /*Test faces of C0 for seperation. because of counterclockwise ordering,
    the projection interval for C0 is [m,0] where m<=0. Only try to determine
    if C1 is on the positive side of the line
  */
  for (i=0; i < C0.L; i++){ //CO.L is number of faces?
    D=C0.F(i).normal; //outward pointing normal to face -- need to define these?
                      //Is this a unit vector?
    if (WhichSide(C1.V, D, C0.F(i).vertex) > 0)
      //C1 is entirely on 'positive' side of line C0.F(i).vertex + t*D
      return FALSE;
    }
  
  for (i=0; i < C1.L; i++){
    D=C1.F(i).normal; //outward pointing
    if (WhichSide(C0.V, D, C1.F(i).vertex) > 0)
      //C0 is entirely on positive side of line C1.F(i).vertex + t*D
      return FALSE;
    }
  
  //Test cross product of pairs of edges one from each polyhedron
  for (i=0; i < C0.M; i++){
  
    for (j=0; j < C1.M; j++){
      
      D = Cross(C0.E(i), C1.E(j));
      
      int side0 = WhichSide(C0.V, D, C0.E(i).vertex);
      if (side0==0)
        continue;
      
      int side1 = WhichSide(C1.V, D, C0.E(i).vertex);
      if (side1==0)
        continue;

      if (side0*side1<0)
        return FALSE;
      }
    }
  
    return TRUE;
  }





/*Intersection Construction*/






/*Moment of Inertia Calculation*/