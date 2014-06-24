#include<stdio.h>
#include<stdlib.h>
#include<assert.h>
#include<math.h>
#include"linear_algebra.h"
#include"mylib.h"

typedef struct {
  double mass;
  double rad;
  double k;
  double dOverlap;
  Vec vNhat; //unit vector of relative velocity between particle and neighbor
  Vec vOrient;
  Mat mInertia;
  Vec r;
  Vec v;
  Vec a;
} PARTICLE;

double OverlapTest(PARTICLE *p, PARTICLE *pn);
double OverlapTest(PARTICLE *p, PARTICLE *pn){
    Vec vRelDist, vRelVel; //relative distance and velocity between two particles
    vecSub(p->r,pn->r,vRelDist);//Find relative Distance and Velocity
    vecSub(p->v,pn->v,vRelVel);
    vecScale(vRelDist,1/vecMag(vRelDist),p->vNhat);
    vecScale(p->vNhat,-1,pn->vNhat);
    p->dOverlap=vecMag(vRelDist)-(p->rad+pn->rad);//check overlap distance
    pn->dOverlap=p->dOverlap;
}

int main(){

  //PARTICLE pd,pnd;
  PARTICLE *p,*pn;
  p=malloc(sizeof(PARTICLE)); //need to allocate memory to pointer
  pn=malloc(sizeof(PARTICLE));


	double dDelta,t,tmax, theta, phi;
  int i,j=10;
  Vec vRelDist, vRelVel, vNhat, vnNhat; //relative distance and velocity between two particles
  dDelta=0.01;
  t=0;
  tmax=15;
 
  vecZero(p->r);
  vecSet(p->v,0.5,0.,0);
  vecZero(p->a);
  p->mass=1.;
  p->rad=0.5;

  vecSet(pn->r,5.,-5,0.);
  vecSet(pn->v,0,0.5,0);
  vecZero(pn->a);
  pn->mass=1.;
  pn->rad=0.5;
  vecZero(vRelDist);
  vecZero(vRelVel);

  while (t<tmax){
    //print to stdout
    if (j==10){
		  printf("%f  %f  %f  %f  %f  %f  %f\n",t,p->r[0],pn->r[0],p->r[1],pn->r[1],p->r[2],pn->r[2]);
      j=0;
    }
      j++;
  
    for (i=0;i<3;i++){
      p->v[i]=p->v[i]+(0.5*dDelta*p->a[i]); //kick
      pn->v[i]=pn->v[i]+(0.5*dDelta*pn->a[i]);
    }

    for (i=0;i<3;i++){
  	  p->r[i]=p->r[i]+(dDelta*p->v[i]); //drift
  	  pn->r[i]=pn->r[i]+(dDelta*pn->v[i]); //
    }
    //Overlap Test
    OverlapTest(p,pn);
    springforce(p->mass,p->dOverlap, p->vNhat, p->a);
    springforce(pn->mass,pn->dOverlap, pn->vNhat, pn->a);
    
    for (i=0;i<3;i++){
      p->v[i]=p->v[i]+(0.5*dDelta*p->a[i]); //kick
      pn->v[i]=pn->v[i]+(0.5*dDelta*pn->a[i]);
    }
    t+=dDelta;
  }

  return 0;
}

//K: loop over i (particles) and k (components)
//D: loop over i (particles) and k (components)
///--->comput acceleration
//K: loop over i (particles) and k (components)

