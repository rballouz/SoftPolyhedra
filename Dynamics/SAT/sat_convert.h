#include <math.h>

typedef long	       Boolean;
/*Boolean type - defined in macros.h*/
typedef enum {FALSE,TRUE} bool;

#define FUZZ	       (0.00001)
#define TRUE	       (1)
#define FALSE	       (0)
#define MAX_VERTS      (1026)

#define ABS(x)	       ( (x) > 0 ? (x) : -(x) )
#define EQZ(x)	       ( ABS((x)) < FUZZ ? TRUE : FALSE )

#define DOT3(u,v)      ( u[0]*v[0] + u[1]*v[1] + u[2]*v[2])
#define VECADD3(r,u,v) { r[0]=u[0]+v[0]; r[1]=u[1]+v[1]; r[2]=u[2]+v[2]; }
#define VECADDS3(r,a,u,v){r[0]=a*u[0]+v[0]; r[1]=a*u[1]+v[1]; r[2]=a*u[2]+v[2];}
#define VECSMULT3(a,u) { u[0]= a * u[0]; u[1]= a * u[1]; u[2]= a * u[2]; }
#define VECSUB3(r,u,v) { r[0]=u[0]-v[0]; r[1]=u[1]-v[1]; r[2]=u[2]-v[2]; }
#define CPVECTOR3(u,v) { u[0]=v[0];	 u[1]=v[1];	 u[2]=v[2]; }
#define VECNEGATE3(u)  { u[0]=(-u[0]);	 u[1]=(-u[1]);	 u[2]=(-u[2]); }

#define GET(u,i,j,s) (*(u+i*s+j))
#define GET3(u,i,j,k,s) (*(u+i*(s*2)+(j*2)+k))

/*This information is defined in the particle structure?-not really needed keeping this here for reference*/
typedef struct polyhedron {
   double   verts[MAX_VERTS][3]; /* 3-D vertices of polyhedron. */
   int	    m;			 /* number of 3-D vertices.  */
   double   trn[3];		 /* translational position in world coords. */
   double   itrn[3];		 /* inverse of translational position. */
} *Polyhedron;

/* Couples two polyhedron, to tests for intersection (repeatedly and using spatial coherence[?]). The backbone of the intersection test, maybe make this an attribute of a particle?*/

typedef struct couple {
   /*This is a structure of particle p, so really only need to carry neighbor?*/
   //PARTICLE  *p;	/* First polyhedron of collision test. */
   PARTICLE  *pn;	/* Second polyhedron of collision test. */
   bool     plane_exists;	/* prospective separating plane flag */
   double      pln_pnt1[3];	/* 1st point used to form separating plane. */
   double      pln_pnt2[3];	/* 2nd point used to form separating plane. */
   int	       vert_indx[4][2]; /* cached points for distance algorithm. */
   int	       n;		/* number of cached points, if any. */
} COUPLE;

/* Evaluate the support (H) and contact (C) function at A for a given polytope.*/
/* This function should take in the polyhedra and vector A, and 
return: max_value[P_V_i dot A], and the index/vector P_V_i*/

double PolySupportContact(tVertex PolyVertices, Vec A, Vec vContactP, int *dContactIndex)
{
  int		 i=0;
  double val, max_val=0;
  tVertex v;
  
  /*looping through each vertex of polytope P finding the vector that gives largest value when dotted with A*/
  v=PolyVertices;
  do{
    val = vecDot(v->v, A);
    if (val > max_val) {
      *dContactIndex = i;
      max_val = val;
      }
    i++;
    v=v->next;
    
    } while (v != PolyVertices);
  vecCopy(PolyVertices[*dContactIndex].v,vContactP);

  return max_val;
}
/* Evaluate the support (H) and contact (C) function at A for a set difference of two polytopes.*/
/* This function should take in two polyhedra and a vector A, and
return: Cs (contact function set difference), Polyedra indices for solution, and Hs (support function difference)*/

double PolySupportContactDiff(tVertex Poly1, tVertex Poly2, Vec A, Vec vContactS, int *dPolyIndex1, int *dPolyIndex2)
{
  Vec vContact1, vContact2, negA;
  double vH1, vH2;
  vH1 = PolySupportContact( Poly1, A, vContact1, dPolyIndex1);
  vecScale(A, -1, negA);
  vH2 = PolySupportContact( Poly2, negA, vContact2, dPolyIndex2);
  
  vecSub(vContact1, vContact2, vContactS);

  return (vH1+vH2);
}

/*Alternate function to compute the point in a polytope closest to the
 *   origin in 3-space. The polytope size m is restricted to 1 < m <= 4.
 *   This function is called only when comp_sub_dist fails.*/

int PolySubDistBack(tVertex Poly, int stop_index, double D_P, double Di_P, int Is, int c2)
{
  bool	  first, pass;
  int		  i, k, s, is, best_s;
  float	  sum, v_aff, best_v_aff;

  first = TRUE;  best_s = -1;
  for (s = 0; s < stop_index; s++) {
    pass = TRUE;
    if (D_P[s] > 0.0) {
      for (i = 1; i <= GET(Is,s,0,c2); i++) {
        is = GET(Is,s,i,c2);
        if (Di_P[s][is] <= 0.0)
        pass = FALSE;
        }
      }
    else
      pass = FALSE;

    if (pass) {
      /*** Compute equation (33) in Gilbert ***/
      k = GET(Is,s,1,c2);
      sum = 0;
      for (i = 1; i <= GET(Is, s, 0, c2); i++) {
        is = GET(Is,s,i,c2);
        sum += Di_P[s][is] * DOT3(P[is],P[k]);
        }
      v_aff = sqrt(sum / D_P[s]);
      if (first) {
        best_s = s;
        best_v_aff = v_aff;
        first = FALSE;
        }
      else {
        if (v_aff < best_v_aff) {
         best_s = s;
         best_v_aff = v_aff;
          }
        }
      }
    }
  if (best_s == -1) {
    printf("backup failed\n");
    exit(0);
    }
  return best_s;
}























/*Function to compute the minimum distance between two convex polytopes */
/* returns Vector different(minVec), points nearest one another (near_indx), Lambda (needed to compute nearest points), and m3 update number of indices for polyhedra in near_indx*/
void PolyMinDist(tVertex Poly1, tVertex Poly2, Vec minVec, int near_indx[4][2], double lambda[], int *dNumInd)
{
   bool pass;
   int  set_size, I[4], i, j, i_tab[4], j_tab[4], P1_i, P2_i, k;
   double ContactDiff, Pk[4][3], Pk_subset[4][3], Vk[3], neg_Vk[3], Cp[3], Gp;

   if ((*dNumInd) == 0) {	     /** if *m3 == 0 use single point initialization **/
      set_size = 1;
      VECSUB3(P1[0], P2[0], Pk[0]);	 /** first elementary polytope **/
      i_tab[0] = j_tab[0] = 0;
   }
   else {				 /** else use indices from near_indx **/
      for (k = 0; k < (*m3); k++) {
	 i = i_tab[k] = near_indx[k][0];
	 j = j_tab[k] = near_indx[k][1];
	 VECSUB3(Pk[k], P1[i], P2[j]);	 /** first elementary polytope **/
      }
      set_size = *m3;
   }

   pass = FALSE;
   while (!pass) {

      /** compute Vk **/

      if (set_size == 1) {
	 CPVECTOR3(Vk, Pk[0]);
	 I[0] = 0;
      }
      else
	 set_size = sub_dist(Pk, set_size, Vk, I, lambda);

      /** eq. (13) **/

      CPVECTOR3(neg_Vk, Vk);	  VECNEGATE3(neg_Vk);
      Gp = DOT3(Vk, Vk) + Hs(P1, m1, P2, m2, neg_Vk, Cp, &P1_i, &P2_i);

      /** keep track of indices for P1 and P2 **/

      for (i = 0; i < set_size; i++) {
	 j = I[i];
	 i_tab[i] = i_tab[j];	  /** j is value from member of some Is **/
	 j_tab[i] = j_tab[j];	  /** j is value from member of some Is **/
      }

      if (EQZ(Gp))		  /** Do we have a solution **/
	 pass = TRUE;
      else {
	 for (i = 0; i < set_size; i++) {
	    j = I[i];
	    CPVECTOR3(Pk_subset[i], Pk[j]);  /** extract affine subset of Pk **/
	 }
	 for (i = 0; i < set_size; i++)
	    CPVECTOR3(Pk[i], Pk_subset[i]);  /** load into Pk+1 **/

	 CPVECTOR3(Pk[i], Cp);		     /** Union of Pk+1 with Cp **/
	 i_tab[i] = P1_i;  j_tab[i] = P2_i;
	 set_size++;
      }
   }

   CPVECTOR3(VP, Vk);			  /** load VP **/
   *m3 = set_size;
   for(i = 0; i < set_size; i++) {
      near_indx[i][0] = i_tab[i];	  /** set indices of near pnt. in P1 **/
      near_indx[i][1] = j_tab[i];	  /** set indices of near pnt. in P2 **/
   }
}







