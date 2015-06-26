#include<stdio.h>
#include<stdlib.h>
#include<assert.h>
#include<math.h>
#include"linear_algebra.h"
#include"particle.h"
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
  p->dKn = 1.06e-15;
  p->dKt = 2.*p->dKn/7.;
  p->dEpsN = 0.8;
  p->dEpsT = 0.8;
  p->dMuS = 0.4;
  
  pn->dKn = p->dKn;
  pn->dKt = p->dKt;
  pn->dEpsN = p->dEpsN;
  pn->dEpsT = p->dEpsT;
  pn->dMuS = p->dMuS;
  
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

  
  /*Simulation loop here*/
  while (t<tmax){
  
    /* Output to file Here*/
    if ((t == 0) || (iStep % iOutFreq == 0)){
      /*Define File Name based on time step*/
      sprintf(achOutFile,"ss.%09d.bt",iStep);
      /*Output function*/
		  FileOutput(p, pn, achOutFile);
    }
    /* Output-Function Ends Here*/
    

    /* LeapFrog */
    for (i=0;i<3;i++){
      p->v[i]=p->v[i]+(0.5*dDelta*p->a[i]); //kick
      pn->v[i]=pn->v[i]+(0.5*dDelta*pn->a[i]);
    }

    for (i=0;i<3;i++){
  	  p->r[i]=p->r[i]+(dDelta*p->v[i]); //drift
  	  pn->r[i]=pn->r[i]+(dDelta*pn->v[i]); //
    }

    //Overlap Test
        bOverlap = bOverlapTest(p,pn);
    if (bOverlap){
      DoDEM(p, pn, dDelta, CnPreFac, CtPreFac);
      printf("%d\n",bOverlap);
    }
    else{
      vecZero(p->a); //Is this correct?
      vecZero(pn->a);
    }
    
    for (i=0;i<3;i++){
      p->v[i]=p->v[i]+(0.5*dDelta*p->a[i]); //kick
      pn->v[i]=pn->v[i]+(0.5*dDelta*pn->a[i]);
    }
    /* LeapFrog ENDS*/


    t+=dDelta;
    iStep++;
  }
  
  return 0;
}

//K: loop over i (particles) and k (components)
//D: loop over i (particles) and k (components)
///--->compute acceleration
//K: loop over i (particles) and k (components)

