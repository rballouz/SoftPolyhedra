#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"macros.h"
#include"chull_3d_defs.h"
#include"chull_3d_func.h"

int main()
{
    ReadVertices();
    DoubleTriangle();
    ConstructHull();
    PrintEdges();
    
    /*
    tVertex v;
    v=vertices;
    do  {
    
        printf("%d  %d  %d\n",v->v[X],v->v[Y],v->v[Z]);
        v=v->next;
    } while (v!=vertices);
    */
    
    return 0;
}