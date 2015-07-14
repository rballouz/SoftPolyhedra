tVertex PolyDualUnion(tVertex p_vd, tVertex pn_vd);
tVertex PolyConstructDual(PARTICLE *p);
void PolyMOIcompute(PARTICLE *p);
void RotatePoly(PARTICLE *p, double dDelta);
void SpringForce (PARTICLE *p, PARTICLE *pn);
int bOverlapTest(PARTICLE *p, PARTICLE *pn);
void matPrint(Mat a); //print out a 3x3 matrix
void FileInput(PARTICLE *p, PARTICLE *pn);
void FileOutput(PARTICLE *p, PARTICLE *pn, char *achOutFile);

/*
void PolyConstructIntersection(PARTICLE *p, PARTICLE *pn)
{
  PARTICLE *p_d,*pn_d,*p_u;
  p_d=malloc(sizeof(PARTICLE));
  pn_d=malloc(sizeof(PARTICLE));
  p_u=malloc(sizeof(PARTICLE)); //need to allocate memory to pointer

  //Construct Duals of p and pn: p_d and pn_d
  
  //create union of p_d and pn_d: p_u
  
  //find convex hull of p_u (only a list of vertices)

  //Construct Duals of of convex hull of p_u
*/

tVertex PolyDualUnion(tVertex p_vd, tVertex pn_vd)
{
  tVertex v_u;
  v_u=p_vd; //link to head of first list
  v_u->prev->next=pn_vd; //link end of first list to head of second list
  v_u->prev->next->prev=v_u->prev; //link head of second to end of first
  v_u->prev=pn_vd->prev; //link head of first to end of second
  v_u->prev->next=v_u; //link end of second to head of first

  return v_u;
  }

tVertex PolyConstructDual(PARTICLE *p) //Create vertex data structure for dual? in particle structure? or return vertex structure?
{
  tFace f;
  tVertex vertices, v;
  Vec v0,v1,v2,n,temp0,temp1;
  int vnum;

  vertices=NULL;
  f=p->faces;
  
  do{
    vecSet(v0, f->vertex[0]->v[X], f->vertex[0]->v[Y], f->vertex[0]->v[Z]);
    vecSet(v1, f->vertex[1]->v[X], f->vertex[1]->v[Y], f->vertex[1]->v[Z]);
    vecSet(v2, f->vertex[2]->v[X], f->vertex[2]->v[Y], f->vertex[2]->v[Z]);

    //Solve Plane equation for normal vector of 3 non-collinear points
    vecSub(v0,v1,temp0);
    vecSub(v0,v2,temp1);
    vecCross(temp0,temp1,n);
    
    //CreateVertexList
    NEW(v, tsVertex);
    ADD(vertices, v);
    v->v[X]=n[X];
    v->v[Y]=n[Y];
    v->v[Z]=n[Z];
    v->vnum = vnum++;
    f=f->next;
    } while (f != p->faces);
  
  return vertices;
}
  
#define PolyMOIsubexp(w0, w1, w2, f1, f2, f3, g0, g1, g2)\
{\
  double temp0, temp1, temp2;\
  temp0=w0+w1;\
  f1=temp0+w2;\
  temp1=w0*w0;\
  temp2=temp1+w1*temp0;\
  f2=temp2+w2*f1;\
  f3=w0*temp1+w1*temp2+w2*f2;\
  g0=f2+w0*(f1+w0);\
  g1=f2+w1*(f1+w1);\
  g2=f2+w2*(f1+w2);\
  }

void PolyMOIcompute(PARTICLE *p)
{
  double one,x,y,z,x_2,y_2,z_2,xy,yz,zx;
  one=0;x=0;y=0;z=0;x_2=0;y_2=0;z_2=0;xy=0;yz=0;zx=0;
  double x0, y0, z0, x1, y1, z1, x2, y2, z2;
  double a1, a2, b1, b2, c1, c2, d0 ,d1, d2;
  double f1x, f2x, f3x, g0x, g1x, g2x;
  double f1y, f2y, f3y, g0y, g1y, g2y;
  double f1z, f2z, f3z, g0z, g1z, g2z;
  double dDensity = p->dDensity;
  double dMass;
  /*loop through faces*/
  tFace f;
  f=p->faces;
  do{
  
    x0=f->vertex[0]->v[X];
    y0=f->vertex[0]->v[Y];
    z0=f->vertex[0]->v[Z];
    x1=f->vertex[1]->v[X];
    y1=f->vertex[1]->v[Y];
    z1=f->vertex[1]->v[Z];
    x2=f->vertex[2]->v[X];
    y2=f->vertex[2]->v[Y];
    z2=f->vertex[2]->v[Z];
  
    a1=x1-x0;
    a2=x2-x0;
    b1=y1-y0;
    b2=y2-y0;
    c1=z1-z0;
    c2=z2-z0;
    d0=(b1*c2)-(b2*c1);
    d1=(a2*c1)-(a1*c2);
    d2=(a1*b2)-(a2*b1);
    
    //Compute Integral Terms
    PolyMOIsubexp(x0, x1, x2, f1x, f2x, f3x, g0x, g1x, g2x);
    PolyMOIsubexp(y0, y1, y2, f1y, f2y, f3y, g0y, g1y, g2y);
    PolyMOIsubexp(z0, z1, z2, f1z, f2z, f3z, g0z, g1z, g2z);
    
    //Update Integrals
    one += d0*f1x;
    x   += d0*f2x;
    y   += d1*f2y;
    z   += d2*f2z;
    x_2  += d0*f3x;
    y_2  += d1*f3y;
    z_2  += d2*f3z;
    xy  += d0*(y0*g0x+y1*g1x+y2*g2x);
    yz  += d1*(z0*g0y+z1*g1y+z2*g2y);
    zx  += d2*(x0*g0z+x1*g1z+x2*g2z);
    f=f->next;
    } while (f != p->faces);
    
    //Multiply by coefficients
    one*=1./6.;x*=1./24.;y*=1./24.;z*=1./24.;x_2*=1./60.;y_2*=1./60.;z_2*=1./60.;xy*=1./120.;yz*=1./120.;zx*=1./120.;
  
    //Define Mass
    p->dMass = one*dDensity;
  
    //Center of Mass
    p->r[0]=x/dMass;
    p->r[1]=y/dMass;
    p->r[2]=z/dMass;
  
    //Inertia Tensor relative to center of Mass
    p->mInertia[0][0]=y_2+z_2-p->dMass*(p->r[1]*p->r[1]+p->r[2]*p->r[2]);
    p->mInertia[1][1]=x_2+z_2-p->dMass*(p->r[2]*p->r[2]+p->r[0]*p->r[0]);
    p->mInertia[2][2]=x_2+y_2-p->dMass*(p->r[0]*p->r[0]+p->r[1]*p->r[1]);
    p->mInertia[0][1]=-(xy-p->dMass*p->r[0]*p->r[1]);
    p->mInertia[1][0]=p->mInertia[0][1];
    p->mInertia[1][2]=-(yz-p->dMass*p->r[1]*p->r[2]);
    p->mInertia[2][1]=p->mInertia[1][2];
    p->mInertia[0][2]=-(zx-p->dMass*p->r[2]*p->r[0]);
    p->mInertia[2][0]=p->mInertia[0][2];
  }

/*Print vertex pairs of edges*/
void RotatePoly(PARTICLE *p, double dDelta)
{
  Vec xAxis, yAxis, zAxis;
  xAxis[0]=1;
  xAxis[1]=0;
  xAxis[2]=0;
  
  yAxis[0]=0;
  yAxis[1]=1;
  yAxis[2]=0;
  
  zAxis[0]=0;
  zAxis[1]=0;
  zAxis[2]=1;
  
  tVertex  v;
  v = p->vertices;
  do{
    /*Rotation Transformation Goes Here*/
    vecRotate(v->v, xAxis, p->w[0]*dDelta);
    vecRotate(v->v, yAxis, p->w[1]*dDelta);
    vecRotate(v->v, zAxis, p->w[2]*dDelta);
    
    v=v->next;
    } while(v != p->vertices);
}

/*Spring Force*/
void SpringForce (PARTICLE *p, PARTICLE *pn){

  Vec vRelDist, vNormal, vTmp; //relative distance
  double dRelDist, x; //Magnitude of Rel Distance, and Amount of Overlap
  double dForce;
  vecSub(pn->r,p->r,vRelDist);//Find relative Distance and Velocity
  dRelDist = vecMag(vRelDist);
  vecScale(vRelDist, 1.0/dRelDist,vNormal); //vNormal points from particle to neighbor
  
  
  /*Calculate the degree of Overlap, then apply accelerations to the particles*/
  x=(p->dRadius+pn->dRadius) - dRelDist;
  /* Normal Deformation */
  dForce = (-p->dKn * x); //Acceleration from normal force: Mag
  vecScale(vNormal, dForce, vTmp); //Scaling with direction
 
  /*Updating acceleration*/
  vecScale(vTmp, 1.0/p->dMass, p->a);
  /*Newton's 3rd law*/
  vecScale(vTmp, -1, vTmp);
  vecScale(vTmp, 1.0/pn->dMass, pn->a); //Applying to neighbor
}


/*Overlap Test*/
int bOverlapTest(PARTICLE *p, PARTICLE *pn){
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


/*File Input Function*/
void FileInput(PARTICLE *p, PARTICLE *pn){
  FILE *fp;
  fp=fopen("initial.bt","r");

  fscanf(fp,"%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d",&p->iOrder,&p->iAggOrder,&p->dMass,&p->dRadius,&p->r[0],&p->r[1],&p->r[2],&p->v[0],&p->v[1],&p->v[2],&p->w[0],&p->w[1],&p->w[2],&p->iColor);
  fscanf(fp,"%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d",&pn->iOrder,&pn->iAggOrder,&pn->dMass,&pn->dRadius,&pn->r[0],&pn->r[1],&pn->r[2],&pn->v[0],&pn->v[1],&pn->v[2],&pn->w[0],&pn->w[1],&pn->w[2],&pn->iColor);

  (void) fclose(fp);
}

/*File Output Function*/
void FileOutput(PARTICLE *p, PARTICLE *pn, char *achOutFile){
  FILE *fp;
  fp=fopen(achOutFile,"w");

  fprintf(fp,"%d %d %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %d\n",p->iOrder,p->iAggOrder,p->dMass,p->dRadius,p->r[0],p->r[1],p->r[2],p->v[0],p->v[1],p->v[2],p->w[0],p->w[1],p->w[2],p->iColor);
  fprintf(fp,"%d %d %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %d",pn->iOrder,pn->iAggOrder,pn->dMass,pn->dRadius,pn->r[0],pn->r[1],pn->r[2],pn->v[0],pn->v[1],pn->v[2],pn->w[0],pn->w[1],pn->w[2],pn->iColor);

  (void) fclose(fp);
}


/*Print Matrix*/
void matPrint(Mat a){

  int i;
  
  for (i=0;i<3;i++)
    printf("%f  %f  %f\n",a[i][0],a[i][1],a[i][2]);
}