void RotatePoly(PARTICLE *p, double dDelta);
void SpringForce (PARTICLE *p, PARTICLE *pn);
int bOverlapTest(PARTICLE *p, PARTICLE *pn);
void matPrint(Mat a); //print out a 3x3 matrix
void FileInput(PARTICLE *p, PARTICLE *pn);
void FileOutput(PARTICLE *p, PARTICLE *pn, char *achOutFile);


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