typedef long	       Boolean;

#define FUZZ	       (0.00001)
#define MAX_VERTS      (1026)

#define ABS(x)	       ( (x) > 0 ? (x) : -(x) )
#define EQZ(x)	       ( ABS((x)) < FUZZ ? TRUE : FALSE )

#define GET(u,i,j,s) (*(u+i*s+j))
#define GET3(u,i,j,k,s) (*(u+i*(s*2)+(j*2)+k))
#define GET_VINDEX(verts,i) {if (verts->vnum == i)\
      ;\
    else{\
      do {\
        verts=verts->next;\
      } while (verts->vnum != i);}\
    }
/* Couples two polyhedron, to tests for intersection (repeatedly and using spatial coherence[?]). The backbone of the intersection test, maybe make this an attribute of a particle?*/
/* Evaluate the support (H) and contact (C) function at A for a given polytope.*/
/* This function should take in the polyhedra and vector A, and 
return: max_value[P_V_i dot A], and the index/vector P_V_i*/

double PolySupportContact(double PolyVerts[][3], int r, Vec A, Vec vContactP, int *dContactIndex)
{
  int		 i=0;
  double val, max_val;
  *dContactIndex = 0;
  
  max_val = vecDot(PolyVerts[0], A);
  *dContactIndex = 0;

  for (i = 1; i < r; i++) {
    val = vecDot(PolyVerts[i], A);
    if (val > max_val) {
      *dContactIndex = i;
      max_val = val;
      }
    }
  vecCopy(PolyVerts[*dContactIndex], vContactP);
  
  return max_val;
}
/* Evaluate the support (H) and contact (C) function at A for a set difference of two polytopes.*/
/* This function should take in two polyhedra and a vector A, and
return: Cs (contact function set difference), Polyedra indices for solution, and Hs (support function difference)*/

double PolySupportContactDiff(double P1[][3], int m1, double P2[][3], int m2, Vec A, Vec vContactS, int *dPolyIndex1, int *dPolyIndex2)
{
  Vec vContact1, vContact2, negA;
  double vH1, vH2;

  vH1 = PolySupportContact( P1, m1, A, vContact1, dPolyIndex1);
  vecScale(A, -1, negA);
  vH2 = PolySupportContact( P2, m2, negA, vContact2, dPolyIndex2);
  vecSub(vContact1, vContact2, vContactS);
  return (vH1+vH2);
}

/*Alternate function to compute the point in a polytope closest to the
 *   origin in 3-space. The polytope size m is restricted to 1 < m <= 4.
 *   This function is called only when comp_sub_dist fails.*/

int PolySubDistBack(double Poly[][3], int stop_index, double D_P[], double Di_P[][4], int *Is, int c2)
{
  Boolean	  first, pass;
  int		  i, k, s, is, best_s;
  float	  sum, v_aff, best_v_aff;

  first = TRUE;
  best_s = -1;
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
        sum += Di_P[s][is] * vecDot(Poly[is],Poly[k]);
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


/*   Function to compute the point in a polytope closest to the origin in
 *   3-space. The polytope size m is restricted to 1 < m <= 4.*/
int PolyCompSubDist(double Poly[][3], int m, int *jo, int *Is, int *IsC, double near_pnt[], int near_indx[], double lambda[])
{
  Boolean pass, fail;
  int i, j, k, isp, is, s, row, col, stop_index, c1, c2;
  double D_P[15], x[3], Dj_P, Di_P[15][4];
  static int combinations[5] = {0,0,3,7,15};
  Vec vTmp;
  
  c1 = m;
  c2 = m + 1;		    /** row offsets for IsC and Is **/

  /** Initialize Di_P for singletons **/

  Di_P[0][0] = Di_P[1][1] = Di_P[2][2] = Di_P[3][3] = 1.0;
  s = 0;
  pass = FALSE;

  while ((s < stop_index) && (!pass)) {	/* loop through each subset */
    D_P[s] = 0.0;  fail = FALSE;
    for (i = 1; i <= GET(Is,s,0,c2); i++) {	/** loop through all Is **/
      is = GET(Is,s,i,c2);
      if (Di_P[s][is] > 0.0)			/** Condition 2 Theorem 2 **/
        D_P[s] += Di_P[s][is];		/** sum from eq. (16)	  **/
      else
        fail = TRUE;
      }

    for (j = 1; j <= GET(IsC,s,0,c1); j++) {	/** loop through all IsC **/
      Dj_P = 0;  k = GET(Is,s,1,c2);
      isp = GET(IsC,s,j,c1);
      for (i = 1; i <= GET(Is,s,0,c2); i++) {
        is = GET(Is,s,i,c2);
        vecSub(Poly[k], Poly[isp],x);		  /** Wk - Wj  eq. (18) **/
        Dj_P += Di_P[s][is] * vecDot(Poly[is], x); /** sum from eq. (18) **/
        }
      row = GET3(jo,s,isp,0,c1);
      col = GET3(jo,s,isp,1,c1);
      Di_P[row][col] = Dj_P;			 /** add new cofactors	**/

      if (Dj_P > 0.00001)		/*RB debug this; scale by radius?*/	 /** Condition 3 Theorem 2 **/
        fail = TRUE;
      }

    if ((!fail) && (D_P[s] > 0.0))  /** Conditions 2 && 3 && 1 Theorem 2  **/
      pass = TRUE;
    else
      s++;
    }

  if (!pass) {
    printf("*** using backup procedure in sub_dist\n");
    s = PolySubDistBack(Poly, stop_index, D_P, Di_P, Is, c2);
    }

  near_pnt[0] = near_pnt[1] = near_pnt[2] = 0.0;
  j = 0;
  for (i = 1; i <= GET(Is,s,0,c2); i++) {	 /** loop through all Is **/
    is = GET(Is,s,i,c2);
    near_indx[j] = is;
    lambda[j] = Di_P[s][is] / D_P[s];		       /** eq. (17)  **/
    
    vecScale(Poly[is], lambda[j], vTmp);
    vecAdd(vTmp, near_pnt, near_pnt);         /** eq. (17)  **/
    j++;
    }

  return (i-1);
}

/*   Function to compute the point in a polytope closest to the origin in
 *   3-space. The polytope size m is restricted to 1 < m <= 4.*/

int PolySubDist(double Poly[][3], int m, double near_pnt[], int near_indx[], double lambda[])
{
  int size;

/*
*
*  Tables to index the Di_P cofactor table in comp_sub_dist.  The s,i
*  entry indicates where to store the cofactors computed with Is_C.
*
*/

  static int jo_2[2][2][2]  =
            { {{0,0}, {2,1}},
              {{2,0}, {0,0}}};

  static int jo_3[6][3][2]  =
            { {{0,0}, {3,1}, {4,2}},
              {{3,0}, {0,0}, {5,2}},
              {{4,0}, {5,1}, {0,0}},
              {{0,0}, {0,0}, {6,2}},
              {{0,0}, {6,1}, {0,0}},
              {{6,0}, {0,0}, {0,0}}};

  static int jo_4[14][4][2] =
            { { {0,0}, {4,1}, {5,2}, {6,3}},
              { {4,0}, {0,0}, {7,2}, {8,3}},
              { {5,0}, {7,1}, {0,0}, {9,3}},
              { {6,0}, {8,1}, {9,2}, {0,0}},
              { {0,0}, {0,0},{10,2},{11,3}},
              { {0,0},{10,1}, {0,0},{12,3}},
              { {0,0},{11,1},{12,2}, {0,0}},
              {{10,0}, {0,0}, {0,0},{13,3}},
              {{11,0}, {0,0},{13,2}, {0,0}},
              {{12,0},{13,1}, {0,0}, {0,0}},
              { {0,0}, {0,0}, {0,0},{14,3}},
              { {0,0}, {0,0},{14,2}, {0,0}},
              { {0,0},{14,1}, {0,0}, {0,0}},
              {{14,0}, {0,0}, {0,0}, {0,0}}};


/*
*  These tables represent each Is.  The first column of each row indicates
*  the size of the set.
*
*/
  static int Is_2[3][3] =
            { {1,0,0},
              {1,1,0},
              {2,0,1}};

  static int Is_3[7][4] =
            { {1,0,0,0},
              {1,1,0,0},
              {1,2,0,0},
              {2,0,1,0},
              {2,0,2,0},
              {2,1,2,0},
              {3,0,1,2}};

  static int Is_4[15][5] =
            { {1,0,0,0,0},
              {1,1,0,0,0},
              {1,2,0,0,0},
              {1,3,0,0,0},
              {2,0,1,0,0},
              {2,0,2,0,0},
              {2,0,3,0,0},
              {2,1,2,0,0},
              {2,1,3,0,0},
              {2,2,3,0,0},
              {3,0,1,2,0},
              {3,0,1,3,0},
              {3,0,2,3,0},
              {3,1,2,3,0},
              {4,0,1,2,3}};

/*
*  These tables represent each Is complement. The first column of each row
*  indicates the size of the set.
*
*/
  static int IsC_2[3][2] = { {1,1}, {1,0}, {0,0}};

  static int IsC_3[7][3] =
            { {2,1,2},
              {2,0,2},
              {2,0,1},
              {1,2,0},
              {1,1,0},
              {1,0,0},
              {0,0,0}};

  static int IsC_4[15][4] =
            { {3,1,2,3},
              {3,0,2,3},
              {3,0,1,3},
              {3,0,1,2},
              {2,2,3,0},
              {2,1,3,0},
              {2,1,2,0},
              {2,0,3,0},
              {2,0,2,0},
              {2,0,1,0},
              {1,3,0,0},
              {1,2,0,0},
              {1,1,0,0},
              {1,0,0,0},
              {0,0,0,0}};

/** Call comp_sub_dist with appropriate tables according to size of P **/

  switch (m) {
    case 2:
      size = PolyCompSubDist(Poly, m, jo_2[0][0], Is_2[0], IsC_2[0], near_pnt, near_indx, lambda);
      break;
    case 3:
      size = PolyCompSubDist(Poly, m, jo_3[0][0], Is_3[0], IsC_3[0], near_pnt, near_indx, lambda);
      break;
    case 4:
      size = PolyCompSubDist(Poly, m, jo_4[0][0], Is_4[0], IsC_4[0], near_pnt, near_indx, lambda);
      break;
    }

  return size;
}

/*   Function to compute the minimum distance between two convex polytopes*/
void PolyMinDist(double P1[][3], int m1, double P2[][3], int m2, double VP[], int near_indx[][2], double lambda[], int *m3)
{
  Boolean pass;
  int set_size, I[4], i, j, i_tab[4], j_tab[4], P1_i, P2_i, k;
  double Pk[4][3], Pk_subset[4][3], Vk[3], neg_Vk[3], Cp[3], Gp;
  
  if ((*m3) == 0) {	     /** if *m3 == 0 use single point initialization **/
    set_size = 1;
    vecSub(P1[0],P2[0],Pk[0]);
    i_tab[0] = j_tab[0] = 0;
    }
  else {				 /** else use indices from near_indx **/
    for (k = 0; k < (*m3); k++) {
      i = i_tab[k] = near_indx[k][0];
      j = j_tab[k] = near_indx[k][1];
      vecSub(P1[i],P2[j],Pk[k]);
      }
    set_size = *m3;
    }
  
  pass = FALSE;
  while (!pass) {
  /** compute Vk **/
    if (set_size == 1) {
      vecCopy(Pk[0],Vk);
      I[0] = 0;
      }
    else
      set_size = PolySubDist(Pk, set_size, Vk, I, lambda);

  /** eq. (13) **/
    vecScale(Vk, -1, neg_Vk);
    Gp = vecDot(Vk, Vk) + PolySupportContactDiff(P1, m1, P2, m2, neg_Vk, Cp, &P1_i, &P2_i);
   
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
        vecCopy(Pk[j],Pk_subset[i]); /** extract affine subset of Pk **/
        }
      for (i = 0; i < set_size; i++)
        vecCopy(Pk_subset[i],Pk[i]);
      vecCopy(Cp, Pk[i]);
      i_tab[i] = P1_i;
      j_tab[i] = P2_i;
      set_size++;
      }
    }


  vecCopy(Vk, VP);
  *m3 = set_size;
  for(i = 0; i < set_size; i++) {
    near_indx[i][0] = i_tab[i];	  /** set indices of near pnt. in P1 **/
    near_indx[i][1] = j_tab[i];	  /** set indices of near pnt. in P2 **/
    }
}

/****************************************************************************
*   Function to compute a proper separating plane between a pair of
 *   polytopes.	 The plane will be a support plane for polytope 1.
 *
 *   On Entry:
 *	couple - couple structure for a pair of polytopes.
 *
 *   On Exit:
 *	couple - containing new proper separating plane, if one was
 *		 found.
 *
 *   Function Return :
 *	result of whether a separating plane exists, or not.
 *
 ****************************************************************************/

Boolean PolyGetNewPlane(COUPLE couple)
{
  tVertex poly1, poly2;
  Boolean plane_exists;
  double pnts1[MAX_VERTS][3], pnts2[MAX_VERTS][3], dist, u[3], v[3], lambda[4], VP[3];
  int i, k, m1, m2;
  Vec vTmp;
  plane_exists = FALSE;
  poly1 = couple->p->vertices;
  poly2 = couple->pn->vertices;
  m1 = couple->p->iNumVertices;
  for (i = 0; i < m1; i++) {
    GET_VINDEX(poly1,i);
    vecCopy(poly1->v,pnts1[i]);
    vecAdd(pnts1[i],couple->p->r,pnts1[i]);
    }

  m2 = couple->pn->iNumVertices;
  for (i = 0; i < m2; i++) {
    GET_VINDEX(poly2,i);
    vecCopy(poly2->v,pnts2[i]);
    vecAdd(pnts2[i],couple->pn->r,pnts2[i]);
    }
  
/** solve eq. (1) for two polytopes **/
  PolyMinDist(pnts1, m1, pnts2, m2, VP, couple->vert_indx, lambda, &couple->n);
  dist = sqrt(vecDot(VP,VP));   /** distance between polytopes **/

  if (!EQZ(dist)) {	       /** Does a separating plane exist **/
    plane_exists = TRUE;
    u[0] = u[1] = u[2] = v[0] = v[1] = v[2] = 0.0;
    for (i = 0; i < couple->n; i++) {
      k = couple->vert_indx[i][0];
      vecScale(pnts1[k], lambda[i], vTmp);
      vecAdd(vTmp, u, u);
      
      k = couple->vert_indx[i][1];
      vecScale(pnts2[k], lambda[i], vTmp);
      vecAdd(vTmp, v, v);
      }
/** Store separating plane in P1's local coordinates **/
    vecScale(couple->p->r,-1,vTmp); /*itrn*/
    vecAdd(u, vTmp, u);
    vecAdd(v, vTmp, v);
/** Place separating plane in couple data structure **/
    vecCopy(u, couple->pln_pnt1);
    vecCopy(u, couple->pln_pnt2);
    }
  
  return plane_exists;
}

/*** RJR 05/26/93 ***********************************************************
 *
 *   Function to detect if two polyhedra are intersecting.
 *
 *   On Entry:
 *	couple - couple structure for a pair of polytopes.
 *
 *   On Exit:
 *
 *   Function Return :
 *	result of whether polyhedra are intersecting or not.
 *
 ****************************************************************************/

Boolean PolyCollisionTest(COUPLE couple)
{
  tVertex poly1, poly2;
  Vec vTmp;
//Polyhedron	  polyhedron1, polyhedron2;
  Boolean collide, loop;
  double u[3], v[3], norm[3], d;
  int i, m;
  
  collide = FALSE;
  poly1 = couple->p->vertices;
  poly2 = couple->pn->vertices;
  
//polyhedron1 = couple->polyhdrn1;	polyhedron2 = couple->polyhdrn2;
  
  if (couple->plane_exists) {

/** Transform proper separating plane to P2 local coordinates.   **/
/** This avoids the computational cost of applying the	       **/
/** transformation matrix to all the vertices of P2.	       **/
    vecCopy(couple->pln_pnt1, u);
    vecCopy(couple->pln_pnt2, v);
    vecAdd(u, couple->p->r, u);
    vecAdd(v, couple->p->r, v);
    vecScale(couple->pn->r, -1, vTmp);
    vecAdd(u, vTmp, u);
    vecAdd(v, vTmp, v);
    vecSub(v,u, norm);

    m = couple->pn->iNumVertices;
    i = 0;
    loop = TRUE;
    while ((i < m) && (loop)) {

      /** Evaluate plane equation **/
      GET_VINDEX(poly2,i)
      vecSub(poly2->v, u, v);
      d = vecDot(v, norm);
      if (d <= 0.0) {		 /** is P2 in opposite half-space **/
        loop = FALSE;
        if (!PolyGetNewPlane(couple)) {
          collide = TRUE;		 /** Collision **/
          couple->plane_exists = FALSE;
          }
        }
      i++;
      }
    }
  else {
    if (PolyGetNewPlane(couple)){
      couple->plane_exists = TRUE;
    }    /** No Collision **/
    else {
      collide = TRUE;
    } /** Collision **/
  }
  return collide;
}
