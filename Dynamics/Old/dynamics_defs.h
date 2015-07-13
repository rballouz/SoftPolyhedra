/*Xor Function*/
bool Xor( bool x, bool y)
{
    return !x ^ !y;
}

/*Vertex, Edge, and Face data structures*/
typedef struct tVertexStructure tsVertex;
typedef tsVertex *tVertex;

typedef struct tEdgeStructure tsEdge;
typedef tsEdge *tEdge;

typedef struct tFaceStructure tsFace;
typedef tsFace *tFace;

struct tVertexStructure{
    float     v[3];
    int     vnum;
    tVertex next, prev;
};

struct tEdgeStructure{
//    tFace   adjface[2];
    tVertex endpts[2];
    tEdge   next, prev;
};

struct tFaceStructure{
    tEdge   edge[3];
    tVertex vertex[3];
    tFace   next, prev;
};


/* Particle Definition*/

/*typedef struct tPolyhedron {
  tVertex vertices;
  tEdge edges;
  tFace faces;
} POLYHEDRON;
*/

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


  tVertex vertices;
  tEdge edges;
  tFace faces;

} PARTICLE;


/*Global Head Pointers
tVertex vertices    =   NULL;
tEdge edges       =   NULL;
tFace faces         =   NULL;
*/
