/*Function Prototypes*/
tVertex MakeNullVertex(PARTICLE *p);
tEdge MakeNullEdge(PARTICLE *p);
tFace MakeNullFace(PARTICLE *p);
void DoubleTriangle(PARTICLE *p);
bool Collinear(tVertex a, tVertex b, tVertex c);
tFace MakeFace( tVertex v0, tVertex v1, tVertex v2, tFace fold, PARTICLE *p);
void ConstructHull(PARTICLE *p);
bool AddOne(tVertex pnt, PARTICLE *p);
int VolumeSign( tFace f, tVertex pnt);
double Volume( tFace f, tVertex pnt );
tFace MakeConeFace( tEdge e, tVertex pnt, PARTICLE *p);
void MakeCcw( tFace f, tEdge e, tVertex pnt);
void CleanUp( tVertex *pvnext, PARTICLE *p );
void CleanVertices( tVertex *pvnext, PARTICLE *p  );
void CleanEdges( PARTICLE *p  );
void CleanFaces( PARTICLE *p );
void ReadVertices(PARTICLE *p, char *achShapeFile);
void PrintEdges(PARTICLE *p);
void PrintFaces(PARTICLE *p);
void PrintVertices(PARTICLE *p);
void PrintVertList(tVertex vertices);

/*Functions*/
tVertex	MakeNullVertex( PARTICLE *p )
{
   tVertex  v;
   
   NEW( v, tsVertex );
   v->duplicate = NULL;
   v->onhull = !ONHULL;
   v->mark = !PROCESSED;
   ADD( p->vertices, v );

   return v;
}
/*---------------------------------------------------------------------
 DoubleTriangle builds the initial double triangle.  It first finds 3 
 noncollinear points and makes two faces out of them, in opposite order.
 It then finds a fourth point that is not coplanar with that face.  The  
 vertices are stored in the face structure in counterclockwise order so 
 that the volume between the face and the point is negative. Lastly, the
 3 newfaces to the fourth point are constructed and the data structures
 are cleaned up. 
---------------------------------------------------------------------*/
void    DoubleTriangle( PARTICLE *p )
{
   tVertex  v0, v1, v2, v3, t;
   tFace    f0, f1 = NULL;
   tEdge    e0, e1, e2, s;
   double    vol;
	
   /* Find 3 noncollinear points. */
   v0 = p->vertices;
   while ( Collinear( v0, v0->next, v0->next->next ) )
      if ( ( v0 = v0->next ) == p->vertices )
         printf("DoubleTriangle:  All points are Collinear!\n"), exit(0);
   v1 = v0->next;
   v2 = v1->next;
	
   /* Mark the vertices as processed. */
   v0->mark = PROCESSED;
   v1->mark = PROCESSED;
   v2->mark = PROCESSED;
   
   /* Create the two "twin" faces. */
   f0 = MakeFace( v0, v1, v2, f1, p );
   f1 = MakeFace( v2, v1, v0, f0, p );

   /* Link adjacent face fields. */
   f0->edge[0]->adjface[1] = f1;
   f0->edge[1]->adjface[1] = f1;
   f0->edge[2]->adjface[1] = f1;
   f1->edge[0]->adjface[1] = f0;
   f1->edge[1]->adjface[1] = f0;
   f1->edge[2]->adjface[1] = f0;
	
   /* Find a fourth, noncoplanar point to form tetrahedron. */
   v3 = v2->next;
   vol = VolumeSign( f0, v3 );
   while ( !vol )   {
      if ( ( v3 = v3->next ) == v0 ) 
         printf("DoubleTriangle:  All points are coplanar!\n"), exit(0);
      vol = VolumeSign( f0, v3 );
   }
	
   /* Insure that v3 will be the first added. */
   p->vertices = v3;

	
}

  
/*---------------------------------------------------------------------
ConstructHull adds the vertices to the hull one at a time.  The hull
vertices are those in the list marked as onhull.
---------------------------------------------------------------------*/
void	ConstructHull( PARTICLE *p )
{
   tVertex  v, vnext;
   double 	  vol;
   bool	    changed;	/* T if addition changes hull; not used. */

   v = p->vertices;
   do {
      vnext = v->next;
      if ( !v->mark ) {
         v->mark = PROCESSED;
        changed = AddOne( v, p );
        CleanUp( &vnext, p ); /* Pass down vnext in case it gets deleted. */
        }
      v = vnext;
   } while ( v != p->vertices );
}
/*---------------------------------------------------------------------
AddOne is passed a vertex.  It first determines all faces visible from 
that point.  If none are visible then the point is marked as not 
onhull.  Next is a loop over edges.  If both faces adjacent to an edge
are visible, then the edge is marked for deletion.  If just one of the
adjacent faces is visible then a new face is constructed.
---------------------------------------------------------------------*/
bool 	AddOne( tVertex pnt, PARTICLE *p )
{
   tFace  f; 
   tEdge  e, temp;
   double vol;
   bool	  vis = FALSE;

   /* Mark faces visible from p. */
   f = p->faces;
   do {
      vol = VolumeSign( f, pnt );
      if ( vol < 0 ) {
	 f->visible = VISIBLE;  
	 vis = TRUE;                      
      }
      f = f->next;
   } while ( f != p->faces );

   /* If no faces are visible from p, then p is inside the hull. */
   if ( !vis ) {
      pnt->onhull = !ONHULL;
      return FALSE; 
   }

   /* Mark edges in interior of visible region for deletion.
      Erect a newface based on each border edge. */
   e = p->edges;
   do {
      temp = e->next;
      if ( e->adjface[0]->visible && e->adjface[1]->visible )
        /* e interior: mark for deletion. */
        e->delete = REMOVED;
      else if ( e->adjface[0]->visible || e->adjface[1]->visible ) 
        /* e border: make a new face. */
        e->newface = MakeConeFace( e, pnt, p );
      e = temp;
   } while ( e != p->edges );
   return TRUE;
}

/*---------------------------------------------------------------------
VolumeSign returns the sign of the volume of the tetrahedron determined by f
and p.  VolumeSign is +1 iff p is on the negative side of f,
where the positive side is determined by the rh-rule.  So the volume 
is positive if the ccw normal to f points outside the tetrahedron.
The final fewer-multiplications form is due to Bob Williamson.
---------------------------------------------------------------------*/
int  VolumeSign( tFace f, tVertex pnt)
{
   double  vol;
   double  ax, ay, az, bx, by, bz, cx, cy, cz;

   ax = f->vertex[0]->v[X] - pnt->v[X];
   ay = f->vertex[0]->v[Y] - pnt->v[Y];
   az = f->vertex[0]->v[Z] - pnt->v[Z];
   bx = f->vertex[1]->v[X] - pnt->v[X];
   by = f->vertex[1]->v[Y] - pnt->v[Y];
   bz = f->vertex[1]->v[Z] - pnt->v[Z];
   cx = f->vertex[2]->v[X] - pnt->v[X];
   cy = f->vertex[2]->v[Y] - pnt->v[Y];
   cz = f->vertex[2]->v[Z] - pnt->v[Z];

   vol =   ax * (by*cz - bz*cy)
         + ay * (bz*cx - bx*cz)
         + az * (bx*cy - by*cx);


   /* The volume sign should be an integer.*/
   /* DEBUG: should probably use some epsilon instead of 0.5*/
   if      ( vol >  0.0 )  return  1;
   else if ( vol <  0.0 )  return -1;
   else                    return  0;
}
/*---------------------------------------------------------------------
Same computation, but computes using doubles, and returns the actual volume.
DEBUG: Double check that this really is the volume of a polyhedron defined 
by face and non-collinear vertex
---------------------------------------------------------------------*/
double  Volume( tFace f, tVertex pnt )
{
   double  vol;
   double  ax, ay, az, bx, by, bz, cx, cy, cz;

   ax = f->vertex[0]->v[X] - pnt->v[X];
   ay = f->vertex[0]->v[Y] - pnt->v[Y];
   az = f->vertex[0]->v[Z] - pnt->v[Z];
   bx = f->vertex[1]->v[X] - pnt->v[X];
   by = f->vertex[1]->v[Y] - pnt->v[Y];
   bz = f->vertex[1]->v[Z] - pnt->v[Z];
   cx = f->vertex[2]->v[X] - pnt->v[X];
   cy = f->vertex[2]->v[Y] - pnt->v[Y];
   cz = f->vertex[2]->v[Z] - pnt->v[Z];

   vol =  (ax * (by*cz - bz*cy)
         + ay * (bz*cx - bx*cz)
         + az * (bx*cy - by*cx));

   return vol;
}



/*---------------------------------------------------------------------
MakeConeFace makes a new face and two new edges between the 
edge and the point that are passed to it. It returns a pointer to
the new face.
---------------------------------------------------------------------*/
tFace	MakeConeFace( tEdge e, tVertex pnt, PARTICLE *p )
{
   tEdge  new_edge[2];
   tFace  new_face;
   int 	  i, j;

   /* Make two new edges (if don't already exist). */
   for ( i=0; i < 2; ++i ) 
      /* If the edge exists, copy it into new_edge. */
      if ( !( new_edge[i] = e->endpts[i]->duplicate) ) {
	 /* Otherwise (duplicate is NULL), MakeNullEdge. */
	 new_edge[i] = MakeNullEdge(p);
	 new_edge[i]->endpts[0] = e->endpts[i];
	 new_edge[i]->endpts[1] = pnt;
	 e->endpts[i]->duplicate = new_edge[i];
      }

   /* Make the new face. */
   new_face = MakeNullFace(p);
   new_face->edge[0] = e;
   new_face->edge[1] = new_edge[0];
   new_face->edge[2] = new_edge[1];
   MakeCcw( new_face, e, pnt );
        
   /* Set the adjacent face pointers. */
   for ( i=0; i < 2; ++i )
      for ( j=0; j < 2; ++j )  
	 /* Only one NULL link should be set to new_face. */
	 if ( !new_edge[i]->adjface[j] ) {
	    new_edge[i]->adjface[j] = new_face;
	    break;
	 }
        
   return new_face;
}

/*---------------------------------------------------------------------
MakeCcw puts the vertices in the face structure in counterclock wise 
order.  We want to store the vertices in the same 
order as in the visible face.  The third vertex is always pnt.

Although no specific ordering of the edges of a face are used
by the code, the following condition is maintained for each face f:
one of the two endpoints of f->edge[i] matches f->vertex[i]. 
But note that this does not imply that f->edge[i] is between
f->vertex[i] and f->vertex[(i+1)%3].  (Thanks to Bob Williamson.)
---------------------------------------------------------------------*/
void	MakeCcw( tFace f, tEdge e, tVertex pnt )
{
   tFace  fv;   /* The visible face adjacent to e */
   int    i;    /* Index of e->endpoint[0] in fv. */
   tEdge  s;	/* Temporary, for swapping */
      
   if  ( e->adjface[0]->visible )      
        fv = e->adjface[0];
   else fv = e->adjface[1];
       
   /* Set vertex[0] & [1] of f to have the same orientation
      as do the corresponding vertices of fv. */ 
   for ( i=0; fv->vertex[i] != e->endpts[0]; ++i )
      ;
   /* Orient f the same as fv. */
   if ( fv->vertex[ (i+1) % 3 ] != e->endpts[1] ) {
      f->vertex[0] = e->endpts[1];  
      f->vertex[1] = e->endpts[0];    
   }
   else {                               
      f->vertex[0] = e->endpts[0];   
      f->vertex[1] = e->endpts[1];      
      SWAP( s, f->edge[1], f->edge[2] );
   }
   /* This swap is tricky. e is edge[0]. edge[1] is based on endpt[0],
      edge[2] on endpt[1].  So if e is oriented "forwards," we
      need to move edge[1] to follow [0], because it precedes. */
   
   f->vertex[2] = pnt;
}
 
/*---------------------------------------------------------------------
MakeNullEdge creates a new cell and initializes all pointers to NULL
and sets all flags to off.  It returns a pointer to the empty cell.
---------------------------------------------------------------------*/
tEdge 	MakeNullEdge( PARTICLE *p )
{
   tEdge  e;

   NEW( e, tsEdge );
   e->adjface[0] = e->adjface[1] = e->newface = NULL;
   e->endpts[0] = e->endpts[1] = NULL;
   e->delete = !REMOVED;
   ADD( p->edges, e );
   return e;
}

/*--------------------------------------------------------------------
MakeNullFace creates a new face structure and initializes all of its
flags to NULL and sets all the flags to off.  It returns a pointer
to the empty cell.
---------------------------------------------------------------------*/
tFace 	MakeNullFace( PARTICLE *p )
{
   tFace  f;
   int    i;

   NEW( f, tsFace);
   for ( i=0; i < 3; ++i ) {
      f->edge[i] = NULL;
      f->vertex[i] = NULL;
   }
   f->visible = !VISIBLE;
   ADD( p->faces, f );
   return f;
}

/*---------------------------------------------------------------------
MakeFace creates a new face structure from three vertices (in ccw
order).  It returns a pointer to the face.
---------------------------------------------------------------------*/
tFace   MakeFace( tVertex v0, tVertex v1, tVertex v2, tFace fold, PARTICLE *p )
{
   tFace  f;
   tEdge  e0, e1, e2;

   /* Create edges of the initial triangle. */
   if( !fold ) {
     e0 = MakeNullEdge(p);
     e1 = MakeNullEdge(p);
     e2 = MakeNullEdge(p);
   }
   else { /* Copy from fold, in reverse order. */
     e0 = fold->edge[2];
     e1 = fold->edge[1];
     e2 = fold->edge[0];
   }
   e0->endpts[0] = v0;              e0->endpts[1] = v1;
   e1->endpts[0] = v1;              e1->endpts[1] = v2;
   e2->endpts[0] = v2;              e2->endpts[1] = v0;
	
   /* Create face for triangle. */
   f = MakeNullFace(p);
   f->edge[0]   = e0;  f->edge[1]   = e1; f->edge[2]   = e2;
   f->vertex[0] = v0;  f->vertex[1] = v1; f->vertex[2] = v2;
	
   /* Link edges to face. */
   e0->adjface[0] = e1->adjface[0] = e2->adjface[0] = f;
	
   return f;
}

/*---------------------------------------------------------------------
CleanUp goes through each data structure list and clears all
flags and NULLs out some pointers.  The order of processing
(edges, faces, vertices) is important.
---------------------------------------------------------------------*/
void	CleanUp( tVertex *pvnext, PARTICLE *p )
{
   CleanEdges(p);
   CleanFaces(p);
   CleanVertices( pvnext, p );
}

/*---------------------------------------------------------------------
CleanEdges runs through the edge list and cleans up the structure.
If there is a newface then it will put that face in place of the 
visible face and NULL out newface. It also deletes so marked edges.
---------------------------------------------------------------------*/
void	CleanEdges( PARTICLE *p )
{
   tEdge  e;	/* Primary index into edge list. */
   tEdge  t;	/* Temporary edge pointer. */
		
   /* Integrate the newface's into the data structure. */
   /* Check every edge. */
   e = p->edges;
   do {
      if ( e->newface ) { 
	 if ( e->adjface[0]->visible )
	    e->adjface[0] = e->newface; 
	 else	e->adjface[1] = e->newface;
	 e->newface = NULL;
      }
      e = e->next;
   } while ( e != p->edges );

   /* Delete any edges marked for deletion. */
   while ( p->edges && p->edges->delete ) {
      e = p->edges;
      DELETE( p->edges, e );
   }
   e = p->edges->next;
   do {
      if ( e->delete ) {
	 t = e;
	 e = e->next;
	 DELETE( p->edges, t );
      }
      else e = e->next;
   } while ( e != p->edges );
}

/*---------------------------------------------------------------------
CleanFaces runs through the face list and deletes any face marked visible.
---------------------------------------------------------------------*/
void	CleanFaces( PARTICLE *p  )
{
   tFace  f;	/* Primary pointer into face list. */
   tFace  t;	/* Temporary pointer, for deleting. */
	

   while ( p->faces && p->faces->visible ) {
      f = p->faces;
      DELETE( p->faces, f );
   }
   f = p->faces->next;
   do {
      if ( f->visible ) {
	 t = f;
	 f = f->next;
	 DELETE( p->faces, t );
      }
      else f = f->next;
   } while ( f != p->faces );
}

/*---------------------------------------------------------------------
CleanVertices runs through the vertex list and deletes the 
vertices that are marked as processed but are not incident to any 
undeleted edges. 
The pointer to vnext, pvnext, is used to alter vnext in
ConstructHull() if we are about to delete vnext.
---------------------------------------------------------------------*/
void	CleanVertices( tVertex *pvnext, PARTICLE *p )
{
   tEdge    e;
   tVertex  v, t;
	
   /* Mark all vertices incident to some undeleted edge as on the hull. */
   e = p->edges;
   do {
      e->endpts[0]->onhull = e->endpts[1]->onhull = ONHULL;
      e = e->next;
   } while (e != p->edges);
	
   /* Delete all vertices that have been processed but
      are not on the hull. */
   while ( p->vertices && p->vertices->mark && !p->vertices->onhull ) {
      /* If about to delete vnext, advance it first. */
      v = p->vertices;
      if( v == *pvnext )
         *pvnext = v->next;
      DELETE( p->vertices, v );
   }
   v = p->vertices->next;
   do {
      if ( v->mark && !v->onhull ) {    
	 t = v; 
	 v = v->next;
	 if( t == *pvnext )
         *pvnext = t->next;
	 DELETE( p->vertices, t );
      }
      else v = v->next;
   } while ( v != p->vertices );
	
   /* Reset flags. */
   v = p->vertices;
   do {
      v->duplicate = NULL; 
      v->onhull = !ONHULL; 
      v = v->next;
   } while ( v != p->vertices );
}

/*---------------------------------------------------------------------
Collinear checks to see if the three points given are collinear,
by checking to see if each element of the cross product is zero.
---------------------------------------------------------------------*/
bool	Collinear( tVertex a, tVertex b, tVertex c )
{
   return 
         ( c->v[Z] - a->v[Z] ) * ( b->v[Y] - a->v[Y] ) -
         ( b->v[Z] - a->v[Z] ) * ( c->v[Y] - a->v[Y] ) == 0
      && ( b->v[Z] - a->v[Z] ) * ( c->v[X] - a->v[X] ) -
         ( b->v[X] - a->v[X] ) * ( c->v[Z] - a->v[Z] ) == 0
      && ( b->v[X] - a->v[X] ) * ( c->v[Y] - a->v[Y] ) -
         ( b->v[Y] - a->v[Y] ) * ( c->v[X] - a->v[X] ) == 0  ;
}

/*Reads in x,y,z from command line*/
void ReadVertices(PARTICLE *p, char *achShapeFile)
{
    FILE *fp;
    fp=fopen(achShapeFile,"r");
    if (fp == NULL) {
        fprintf(stderr, "Can't open input file!\n");
        exit(1);
    }

    tVertex     v;
    double       x, y, z;
    int         vnum=0;
    
    while (fscanf(fp,"%lf %lf %lf",&x, &y,&z) != EOF) {
        v = MakeNullVertex(p);
        v->v[X] = x;
        v->v[Y] = y;
        v->v[Z] = z;
        //vecScale(v->v,p->dRadius,v->v); /*Scale by Particle Radius*/
        v->vnum = vnum++;
    }
}

/*Print vertex pairs of edges*/
void PrintEdges(PARTICLE *p)
{
  double x0, y0, z0, x1, y1, z1;
  tEdge e;

  e = p->edges;
  do{
    x0=e->endpts[0]->v[X];
    y0=e->endpts[0]->v[Y];
    z0=e->endpts[0]->v[Z];
    x1=e->endpts[1]->v[X];
    y1=e->endpts[1]->v[Y];
    z1=e->endpts[1]->v[Z];
    printf("%.2g  %.2g  %.2g  %.2g  %.2g  %.2g\n",x0,y0,z0,x1,y1,z1);
    e=e->next;
    } while(e != p->edges);
}

void PrintFaces(PARTICLE *p)
{
  double x0, y0, z0, x1, y1, z1, x2, y2, z2;
  tFace f;

  f = p->faces;
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

    printf("%.2g  %.2g  %.2g  %.2g  %.2g  %.2g  %.2g  %.2g  %.2g\n",x0,y0,z0,x1,y1,z1,x2,y2,z2);
    f=f->next;
    } while(f != p->faces);
    printf("------\n");
}

void PrintVertices(PARTICLE *p)
{
  double x0, y0, z0;
  tVertex v;

  v = p->vertices;
  do{
    x0=v->v[X];
    y0=v->v[Y];
    z0=v->v[Z];
    
    printf("%.2g  %.2g  %.2g\n",x0,y0,z0);
    v=v->next;
    } while(v != p->vertices);
    printf("------\n");
}

void PrintVertList(tVertex vertices)
{
  double x0, y0, z0;
  tVertex v;

  v = vertices;
  do{
    x0=v->v[X];
    y0=v->v[Y];
    z0=v->v[Z];
    
    printf("%.2g  %.2g  %.2g\n",x0,y0,z0);
    v=v->next;
    } while(v != vertices);
    printf("------\n");
}
