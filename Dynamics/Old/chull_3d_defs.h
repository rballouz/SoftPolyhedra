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
    tEdge   duplicate; /*pointer to incident cone edge (or NULL)*/
    bool    onhull; /*T iff point on hull*/
    bool    mark; /*T iff point already processed*/
    tVertex next, prev;
};

struct tEdgeStructure{
    tFace   adjface[2];
    tVertex endpts[2];
    tFace   newface; /*pointer to incident cone face*/
    bool    delete; /*T iff edge should be deleted*/
    tEdge   next, prev;
};

struct tFaceStructure{
    tEdge   edge[3];
    tVertex vertex[3];
    bool    visible; /*T iff face visible from new point*/
    tFace   next, prev;
};

/*Xor Function*/
bool Xor( bool x, bool y)
{
    return !x ^ !y;
}

/*Global Head Pointers*/
tVertex vertices    =   NULL;
tEdge edges       =   NULL;
tFace faces         =   NULL;

