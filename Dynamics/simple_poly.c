#include<stdio.h>
#include<stdlib.h>
#include<assert.h>
#include<math.h>
#include"linear_algebra.h"
#include"macros.h"
#include"my_chull.h"
#include"particle.h"
#include"my_chull_func.h"
#include"mylib.h"
#include"smoothfcn_eq.h"

int main(){

  PARTICLE *p,*pn;
  p=malloc(sizeof(PARTICLE)); //need to allocate memory to pointer
  pn=malloc(sizeof(PARTICLE));
  FILE *fp;

	double dDelta,t,tmax;
  double CnPreFac, CtPreFac, dLnEpsN, dLnEpsT, dLnEpsNsq, dLnEpsTsq, a;;
  int i;
  int nSteps,iStep=0;
  int iOutFreq;
  int bOverlap=0;
  char achOutFile[20];

  //Need to initialize dKn for each particle
  /*TODO: move this to an initialization function to clean up*/
  p->dKn = 1.06e-15;
  p->dKt = 2.*p->dKn/7.;
  p->dEpsN = 0.8;
  p->dEpsT = 0.8;
  p->dMuS = 0.4;
  p->vertices=NULL;
  p->edges=NULL;
  p->faces=NULL;
  
  
  pn->dKn = p->dKn;
  pn->dKt = p->dKt;
  pn->dEpsN = p->dEpsN;
  pn->dEpsT = p->dEpsT;
  pn->dMuS = p->dMuS;
  pn->vertices=NULL;
  pn->edges=NULL;
  pn->faces=NULL;
  
  dDelta=1.47e-10;
  nSteps = 40000;
  iOutFreq = nSteps / 50;
  t=0;
  tmax=dDelta * nSteps;
 
 
  const double pi_sq = M_PI*M_PI;
  static int bCnCtPreFacCalculated = 0; /* since dEpsN and dEpsT are identical for all particles in sim */

  if (!bCnCtPreFacCalculated) {
    /* damping term: normal */
    dLnEpsN = log(p->dEpsN);
    dLnEpsNsq = dLnEpsN*dLnEpsN;
    a = pi_sq + dLnEpsNsq;
    CnPreFac = -sign(dLnEpsN)*sqrt(dLnEpsNsq*p->dKn/a);
    CnPreFac += CnPreFac;

    /* damping term: tangential */
    dLnEpsT = log(p->dEpsT);
    dLnEpsTsq = dLnEpsT*dLnEpsT;
    a = pi_sq + dLnEpsTsq;
    CtPreFac = -sign(dLnEpsT)*sqrt(dLnEpsTsq*p->dKt/a);
    CtPreFac += CtPreFac;
    bCnCtPreFacCalculated = 1;
    }

  /*Take input from FILE here - Initialization*/
  FileInput(p, pn);

  /*Construct Hull for each particle - Edges, Faces based on shape information (i.e. Vertex location about center of particle - store information in memory)*/
  /*Put these 3 / 4 calls into a single function that takes in the particle struct*/
  ReadVertices(p);
  DoubleTriangle(p);
  ConstructHull(p);
  PrintFaces(p);
  PrintVertices(p);
  //p->w[2]=0.5;
  //dDelta=0.5;
  //RotatePoly(p, dDelta);
  
  
  return 0;
}

//K: loop over i (particles) and k (components)
//D: loop over i (particles) and k (components)
///--->compute acceleration
//K: loop over i (particles) and k (components)

