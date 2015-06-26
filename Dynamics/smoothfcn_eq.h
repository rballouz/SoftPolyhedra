void DoDEM(PARTICLE *p, PARTICLE *pn, double dDelta, double CnPreFac, double CtPreFac);

void DoDEM(PARTICLE *p, PARTICLE *pn, double dDelta, double CnPreFac, double CtPreFac){


/* Put Variable Declaration here */
  double m1, m2, r1, r2, r1sq, r2sq, dReducedMass, dSqrtReducedMass, x, d;
  double l1, l2, b, Un, Ut, Cn, Ct, a1, i1, b1, a2, i2, b2, N, T, F, F2, a;
  int i = 0;
  double dTangentialSpringDrag = 0;
  Vec vRelDist, n, t, s1, s2, v, s, u, un, ut;
  Vec vDum, vDum2, Fn, Ft, Ftn;

  /* set mass and radius and reduced mass*/
  m1 = p->dMass;
  m2 = pn->dMass;
  r1 = p->dRadius;
  r2 = pn->dRadius;
  
  r1sq=r1*r1;
  if (r1 == r2)
    r2sq = r1sq;
  else
    r2sq = r2*r2;
  if (m1 == m2) /* save time if masses are equal */
    dReducedMass = 0.5*m1;
  else
    dReducedMass = m1*m2/(m1 + m2);
  
  dSqrtReducedMass = sqrt(dReducedMass);


  /*Seperation == sqrt(fDist2); lever arm calculated ; overlap calculated*/
  //d = sqrt(nnList[j].fDist2);
  vecSub(pn->r,p->r,vRelDist);
  d = vecMag(vRelDist);
  
  if (r1 == r2) /* save time if radii are equal */
    l1 = l2 = 0.5*d;
  else {
    l1 = (r1sq - r2sq + (d*d))/(d + d);
    l2 = d - l1;
    //mdlassert(smf->pkd->mdl,l1 >= 0.0 && l2 >= 0.0);
    }

  x = r1 + r2 - d; /* probably want to be checking x/d for avg and max, warning if max is large... */
  //mdlassert(smf->pkd->mdl,x >= 0.0);

  /* a = 1 / seperation */
  a = 1./d;
  /*Set n = n-hat : seperation/mag(seperation)*/
  //vectorSet(n,-nnList[j].dx*a,-nnList[j].dy*a,-nnList[j].dz*a); /* sign: n[i] = unit vector of [rn - r], points from particle's center toward neighbor's center */
  
  
  vecScale(vRelDist, a, n);
  /* use predicted w since we're in post-drift pre-kick stage (already half-stepped v's and w's)*/
  vecCross(p->w,n,s1);
  vecCross(pn->w,n,s2);
  /*calculating component for u, l(n-hat x w)*/
  vecScale(s1,l1,s1);
  vecScale(s2,-l2,s2);

  /*magnitudes of Un and Ut*/
  Un = Ut = 0.;

  /* Calculate u, un, ut */
  /*why not use linalg here?*/
  for (i=0;i<3;i++)
    v[i] = pn->v[i] - p->v[i];
  
  for (i=0;i<3;i++) {
    s[i] = s2[i] - s1[i];
    u[i] = v[i] + s[i];
    Un += u[i]*n[i];
    }

  for (i=0;i<3;i++) {
    un[i] = Un*n[i];
    ut[i] = u[i] - un[i];
    Ut += ut[i]*ut[i];
    }
  Ut = sqrt(Ut);
  
   /* unit tangential velocity vector */
  if (Ut != 0.0)
    vecScale(ut,1./Ut,t);
  else
    vecZero(t); /* avoid dividing by zero */


  /* PreFac's calculated earlier, before calculating reduced masses ... */
  /* damping term: normal */
  Cn = CnPreFac*dSqrtReducedMass;
  /* damping term: tangential */
  Ct = CtPreFac*dSqrtReducedMass; /* ...relating dEpsN,dEpsT to Cn,Ct comes from Cleary etal. '98 */

  /*a1= inverse mass; i1 = mom. of inertia for sphere; b1 = inverse of MoI*/
  a1 = 1./m1;
  i1 = 0.4*m1*r1sq;
  b1 = 1./i1;

  if (m1 == m2) { /* save time if masses are equal */
    a2 = a1;
    if (r1 == r2) { /* save time if both masses and radii are equal */
      i2 = i1;
      b2 = b1;
      }
    else {
      i2 = 0.4*m2*r2sq;
      b2 = 1./i2;
      }
    }
  else {
    a2 = 1./m2;
    i2 = 0.4*m2*r2sq;
    b2 = 1./i2;
    }


  /* Tracking motion of contact point, vnOld is previous steps n-hat?*/
  if (vecMagSq(p->vnOld) != 0.0) {
    //mdlassert(smf->pkd->mdl,pe->liOverlapCounter > 1); /* should not be first step in overlap */

    /*Vector addition of spin's into vDum*/
    vecAdd(p->w,pn->w,vDum);
    /*Vector addition of n-hat and old n-hat to vDum2*/
    vecAdd(p->vnOld,n,vDum2);
    /*Unit Vector of addition */
    vecNorm(vDum2);
    /*vShear essentially S vector? : rotating the S vector about n-hat*/
    vecRotate(p->vShear,vDum2,0.5*dDelta*vecDot(vDum,vDum2));
    
    /*vector orthogonal to n-hat and old n-hat*/
    vecCross(p->vnOld,n,vDum); // axis of rotation (not normalized)
    /*if n-hat changes over the last time-step*/
    if ((a = vecMagSq(vDum)) != 0.0) {
      a = sqrt(a);
      vecScale(vDum,1./a,vDum);
      /*rotate about the orthogonal vector*/
      vecRotate2(p->vShear,vDum,a,vecDot(n,p->vnOld)); // converting to an angle and using vectorRotate() failed when recomputing the sin and cos
      }
    }

  /* set vnOld to current n for use in next step in the case that the overlap persists to the next step */
  vecCopy(n,p->vnOld);

  N = T = 0.0; /* normal, tangential forces (scalar) */

  for (i=0;i<3;i++) { /*DEBUG: optimize this function for cases when !dMuS)*/
    /*This is the integral portion of S's calculation eq. 7 in Schwartz et al (2012)*/
    p->vShear[i] += ut[i]*dDelta; /* this step's contribution to the tangential spring */
    Fn[i] = -p->dKn*x*n[i] + Cn*un[i]; /* check to ensure that all other forces are computed first (e.g. gravity) */
    Ft[i] = p->dKt*p->vShear[i] + Ct*ut[i]; /*      ""      */
    N += Fn[i]*Fn[i]; /* square of normal force */
    T += Ft[i]*Ft[i]; /* square of tangential force */
    }

  F2 = N*p->dMuS*p->dMuS; /* square of max allowed static friction */
  N = sqrt(N); /* N is now defined as the normal force (not squared) */

  /*minimizing the tangential Force*/
  if (F2 < T) {

    F = N*p->dMuS; /* max allowed static friction (not squared) */
    /*When is this not equal to zero?*/
    if (dTangentialSpringDrag != 0.0) {

    /* compute the new magnitude of vShear, but if tangential kintetic friction alone exceeds static friction, vShear will be zeroed */
    a = F - Ct*Ut;
    b = a <= 0. ? 0. : vecMag(p->vShear);
    a = b <= 0. ? 0. : dTangentialSpringDrag*a/(p->dKt*b);

    /* set vShear */
    vecScale(p->vShear,a,p->vShear); /* vShear for this step */
    vecScale(t,min(F,Ct*Ut),vDum); /* the force due to tangential kinetic friction */
    vecScale(p->vShear,p->dKt,vDum2);
    vecAdd(vDum,vDum2,vDum); /* total tangential force due to both static and kinetic friction */

    //mdlassert(smf->pkd->mdl,T > 0.);
    if (a)
      vecScale(Ft,a*b/sqrt(T),p->vShear);
    else
      vecZero(p->vShear);
    
    vecCopy(vDum,Ft);
    }
    else {
      vecZero(p->vShear);
      vecScale(t,min(F,Ct*Ut),Ft);
      }
    }
  
  /*Torque*/
  vecCross(Ft,n,Ftn); /* tangential DEM induced torque per unit length */
  
  /*Updating acceleration and spins based on forces and torques*/
  for (i=0;i<3;i++) {
    p->a[i] = (Fn[i] + Ft[i])*a1;
    p->w[i] -= b1*(l1*(Ftn[i]))*0.5*dDelta;
    pn->a[i] = -(Fn[i] + Ft[i])*a2;
    pn->w[i] -= b2*(l2*(Ftn[i]))*0.5*dDelta; /*shouldn't this be plus?*/

    /*p->dDeltaAccel[i] += (Fn[i] + Ft[i])*a1;
    /p->wDot[i] -= b1*(l1*(Ftn[i] - vRft[i]) - vTfn[i]);
    /pn->dDeltaAccel[i] -= (Fn[i] + Ft[i])*a2;
    /pn->wDot[i] -= b2*(l2*(Ftn[i] + vRft[i]) + vTfn[i]);*/
    }
}
  /* end loop over nSmooth */


