#include<stdio.h>
#include<stdlib.h>
#include<assert.h>
#include<math.h>
#include"linear_algebra.h"
#include"macros.h"
#include"dynamics_defs.h"
#include"force.h"
#include"polyhedra.h"



int main(){

  //PARTICLE pd,pnd;
  PARTICLE *p,*pn;
  p=malloc(sizeof(PARTICLE)); //need to allocate memory to pointer

  tVertex v;

  ReadVertices(p);
  ReadEdges(p);

 
    float     x0, y0, z0, x1, y1, z1;
    tEdge   e;
    
    e = p->edges;
    do{
        x0=e->endpts[0]->v[X];
        y0=e->endpts[0]->v[Y];
        z0=e->endpts[0]->v[Z];
        x1=e->endpts[1]->v[X];
        y1=e->endpts[1]->v[Y];
        z1=e->endpts[1]->v[Z];
        printf("%.2f  %.2f  %.2f  %.2f  %.2f  %.2f\n",x0,y0,z0,x1,y1,z1);
        e=e->next;
    } while(e != p->edges);

/*
	double dDelta,t,tmax;
  int i;
  dDelta=0.01;
  t=0;
  tmax=15;

  vecZero(p->r);
  vecSet(p->v,0.5,0.,0);
  vecZero(p->a);
  p->mass=1.;
  p->rad=0.5;

  while (t<tmax){

		  //printf("%f  %f  %f  %f\n",t,p->r[0],p->r[1],p->r[2]);
  
    for (i=0;i<3;i++){
      p->v[i]=p->v[i]+(0.5*dDelta*p->a[i]); //kick
    }

    for (i=0;i<3;i++){
  	  p->r[i]=p->r[i]+(dDelta*p->v[i]); //drift
    }
    
    for (i=0;i<3;i++){
      p->v[i]=p->v[i]+(0.5*dDelta*p->a[i]); //kick
    }
    t+=dDelta;
  }
*/
  return 0;
}
