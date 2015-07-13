/*Polyhedra Related functions*/
void ReadVertices(PARTICLE *p);
void ReadEdges(PARTICLE *p);

/*Reading in Polyhedra Shape*/
void ReadVertices(PARTICLE *p)
{
    FILE *fp;
    fp=fopen("hull_vertices.in","r");
    if (fp == NULL) {
        fprintf(stderr, "Can't open input file!\n");
        exit(1);
    }

    tVertex     v;
    float       x, y, z;
    int         vnum=0;
    
    while (fscanf(fp,"%f %f %f",&x, &y,&z) != EOF) {
      NEW( v, tsVertex );   
      ADD( p->vertices, v );
      v->v[X] = x;
      v->v[Y] = y;
      v->v[Z] = z;
      v->vnum = vnum++;
    }
}

/*Print vertex pairs of edges*/
void ReadEdges(PARTICLE *p)
{
    FILE *fp;
    fp=fopen("hull_edges.in","r");
    if (fp == NULL) {
        fprintf(stderr, "Can't open input file!\n");
        exit(1);
    }

    float   x0, y0, z0, x1, y1, z1;
    tEdge   e;    

    while (fscanf(fp,"%f %f %f %f %f %f",&x0, &y0, &z0, &x1, &y1, &z1) != EOF) {

        NEW( e, tsEdge );
        ADD( p->edges, e );

        NEW( e->endpts[0], tsVertex);
        NEW( e->endpts[1], tsVertex);
//        e->adjface[0] = e->adjface[1] = NULL;

        e->endpts[0]->v[X]=x0;
        e->endpts[0]->v[Y]=y0;
        e->endpts[0]->v[Z]=z0;
        e->endpts[1]->v[X]=x1;
        e->endpts[1]->v[Y]=y1;
        e->endpts[1]->v[Z]=z1;
    }

}

