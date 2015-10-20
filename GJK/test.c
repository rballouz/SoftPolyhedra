/*Main in gjkdemo.c*/
int main( int argc, char ** argv)
{
  /* vertex arrays: both objects, and a centred form of the second object */
  REAL initial_points1[MAX_POINTS][3], initial_points2[MAX_POINTS][3],
  transformed_points1[MAX_POINTS][3], transformed_points2[MAX_POINTS][3],
  (*ppoints1)[3], (*ppoints2)[3];
  /* edge ring arrays for both objects */
  /*VertexID is an int*/
  /*MAX_RING_SIZE_MULTIPLIER=8, MAX_POINTS=1000*/
  VertexID fixed_rings1[MAX_RING_SIZE_MULTIPLIER * MAX_POINTS];
  VertexID fixed_rings2[MAX_RING_SIZE_MULTIPLIER * MAX_POINTS];
  VertexID * rings1, * rings2;
  REAL (*tr1)[DIM+1] = 0, (*tr2)[DIM+1] = 0;
  unsigned long this_seed, original_seed; /* random number seeds */
  int num_pts, repeat, num_repeats, instance, num_instances;
  long run, num_runs, num_errors, total_g_test, total_simplices,
  total_vertices, total_dp, total_sdp, total_mult_divide_only,
  total_ops, num_zero_runs, tick_frequency;
  double aver_verts, aver_dp, aver_sdp, aver_ops, aver_time;
  double aver_time_g, aver_time_simp, aver_time_dp, aver_time_op;
  REAL dist, sqdist, dp_error_val, total_dist, aver_dist, aver_zeros;
  REAL wit1[3], wit2[3], one_to_two[3], two_to_one[3], hull2_shift[3];
  struct simplex_point simplex_record, initial_simplex_record;
  double start_time, total_time, initial_time, elapsed_time;
  char outbuff[256]; /* line output buffer */
  int i, d, nbad1, nbad2;
  REAL err, err1, err2;
  /* various flags */
  int quiet = 0, rotating = 0, use_polyballs = 1, allow_hill_climbing = 1;
  int pass_transforms = 1;
  /* various transformations */
  /* definition of local type for a rigid-body transformation
  struct transf { double matrix[DIM][DIM+1]; };*/

  struct transf initial_trf1, initial_trf2, trf2,
  * ptrf1, * ptrf2, delta_rot,current_rot;
  struct Object_structure obj1, obj2;
  /* Object structure: holds basic information about each object
  Definition transferred from gjk.h ...
  struct Object_structure {
  int numpoints;
  REAL (* vertices)[DIM];
  int * rings;
  };
  typedef struct Object_structure * Object;*/

  num_instances = 10;
  num_repeats = 100;
  tick_frequency = 0;
  prog_name = argv[0];

  /* now deal with setting the intial value of the seed. */
  original_seed = (unsigned long) time( NULL);

  /* check whether we are doing hill-climbing */
  if ( allow_hill_climbing ) {
    rings1 = fixed_rings1; /*rings are integer pointers, fixed_rings is a an array of ints. */
    rings2 = fixed_rings2; /* hill climbing implies some kind of iterative method. Maybe we're allocate some static memory here?*/
    }
  else {
    rings1 = rings2 = 0;
    }

  /* passing some initial transformation of the object?*/
  /*same deal as before, memory allocation*/
  if ( pass_transforms ) {
    ptrf1 = &initial_trf1;
    ptrf2 = &        trf2;
    tr1   = ptrf1->matrix;
    tr2   = ptrf2->matrix;
    ppoints1 = initial_points1;
    ppoints2 = initial_points2;
    }
  else {
    tr1 = tr2 = 0;
    ptrf1 = ptrf2 = 0;
    ppoints1 = transformed_points1; /*These are arrays (type REAL = double) with dimensions [MAX_POINTS][3]*/
    ppoints2 = transformed_points2;
    }

  /* Now the main loop to run GJK itself */

  repeat = 0;  instance = 0;  run = 0;   seed = original_seed;

  while ( run<num_runs ) {
    
    /**********************************************************
    STEP 1: Initialize Vertex Structures for Shapes
            & Apply Transformations 
    **********************************************************/
    if ( instance == 0 ) { /* then create new shapes */
      this_seed = seed;
      /*Create the two objects, each with num_vertices=num_pts, rings (edge structure) required for hill_climbing algorithm*/
      #ifdef QHULL
      if ( !use_polyballs ) {
        create_hull( num_pts, initial_points1, rings1);
        create_hull( num_pts, initial_points2, rings2);
        }
       else
      #endif
      {
        create_polyball( num_pts, initial_points1, rings1);
        create_polyball( num_pts, initial_points2, rings2);
        }

      /* define initial rotations and translations*/
      /*
      generate_small_random_rotation( &initial_trf1);
      generate_small_random_rotation( &initial_trf2);
      mkid( &current_rot);
      generate_initial_random_shift( hull2_shift);
      */
      /*mkid = function to return Identity Transformation: 3x4 Matrix*/
      mkid( &initial_trf1);
      mkid( &initial_trf2);
      mkid( &current_rot);
      generate_initial_random_shift( hull2_shift);
      hull2_shift[0] = hull2_shift[1] = 0; /*shifting randomly in the z-direction*/

      /* If we are not going to pass transformations, we
      can transform all the points in the first hull now.
      */
      if ( !pass_transforms )
        transform_hull( num_pts, initial_points1, ppoints1, &initial_trf1);
      } /****End of Instance==0******/
    
    else  /* If not the first instance: shift the second hull, and rotate if necessary */
    {
      if ( rotating ) {
        generate_small_random_rotation( &delta_rot);
        compose_transformations( &current_rot, &delta_rot, &current_rot);
        }
      add_delta_random_shift( hull2_shift);
      hull2_shift[0] = hull2_shift[1] = 0;
      }

    /* Now construct appropriate transformation for second hull */
    compose_transformations( &trf2, &current_rot, &initial_trf2);
    for ( i=0 ; i<3 ; i++ )
      trf2.matrix[i][3] += hull2_shift[i];
	 
    /* transform the point in the second hull
    if we are not passing transformations */
    if ( !pass_transforms )
      transform_hull( num_pts, initial_points2, ppoints2, &trf2);
    /**********Done With Step 1: Applying Transformations******************************************/
    
    /**********************************************************
    STEP 2: (Re-)Initialize Objects (after applying transformations)
            , Simplex Structure, then Compute Distance (Call gjk_disance routine)
    **********************************************************/
    /* save the current value of the simplex record */
    initial_simplex_record = simplex_record; /*simplex_point structure defined in gjk.h , holds a lot of information, mainly the "error" i.e. error in the shortest distance vector - should be close to zero*/
    
    obj1.numpoints = num_pts; /*object structure defined in gjk.h (has 3 traits: num_pts, vertex pointer, ring pointer)*/
    obj1.vertices = ppoints1;
    obj1.rings = rings1;

    obj2.numpoints = num_pts;
    obj2.vertices = ppoints2;
    obj2.rings = rings2;

    /* now time num_repeats calls to GJK */
    start_time = clocktime();

    for ( repeat = 0 ; repeat < num_repeats ; repeat++ ) {
      simplex_record = initial_simplex_record;
      /* the call to GJK itself */
      sqdist = gjk_distance( &obj1, tr1, &obj2, tr2, wit1, wit2, &simplex_record, (instance>0));
      }
    total_time += clocktime() - start_time;
    dist = sqrt( sqdist);
    /* Compute the expected width of the error bound */
    err = ( dist>0 ) ? ( simplex_record.error / dist ) : 0.0 ;

    total_dist += dist;
    /**********Done With Step 2: computing gjk_distance******************************************/
  
    /**********************************************************
    STEP 3: Testing the Answer
    **********************************************************/
    /* Now to test the answers.  Compute the displacement vectors
    Don't assume that the witness points were returned. */
    gjk_extract_point( &simplex_record, 1, wit1); /*Extracts the witness points from the simplex_record structure*/
    gjk_extract_point( &simplex_record, 2, wit2);
    for ( d=0 ; d<3 ; d++ ) {
      two_to_one[d] = wit1[d] - wit2[d];
      one_to_two[d] = - two_to_one[d];
      }

    /* dp_error_val should be zero, to the accuracy of our arithmetic */
    dp_error_val = sqdist - dot_product( two_to_one, two_to_one);

    /* now check, for each hull, how many of its points lie to the
    wrong side of a plane that lies through the witness point
    and perpendicular to the displacement vector
    */
    if ( sqdist>err*err ) {
      nbad1 = num_further( num_pts, err, &err1, ppoints1, ptrf1,
      one_to_two, wit1);
      nbad2 = num_further( num_pts, err, &err2, ppoints2, ptrf2,
      two_to_one, wit2);

      }
    else {
      /* zero was returned, so check that the witness point
      lies within both hulls. Requires QHULL
      */
      #ifdef QHULL
      nbad1 = num_outside( num_pts, err, &err1, ppoints1, ptrf1, wit1);
      nbad2 = num_outside( num_pts, err, &err2, ppoints2, ptrf2, wit2);
      #else
      nbad1 = nbad2 = 0; /* can't test, assume OK */
      #endif /* QHULL */
      }

    instance = ( instance+1 ) % num_instances;
    run++;
    }
      /* end of main loop that calls GJK */

   elapsed_time = clocktime() - initial_time;
   /**********Done With Step 3: Testing the Answer******************************************/

   
   exit( 0);
   }


/******************************************************************************************************************************
******  The actual gjk_distance algorithm in gjk.c
******************************************************************************************************************************/
REAL gjk_distance( Object obj1, Transform tr1, Object obj2, Transform tr2, REAL wpt1[DIM], REAL wpt2[DIM], struct simplex_point * simplex, int use_seed) {

   VertexID v, p, maxp, minp;
   REAL minus_minv, maxv, sqrd, g_val;
   REAL displacementv[DIM], reverse_displacementv[DIM];
   REAL local_witness1[DIM], local_witness2[DIM];
   REAL local_fdisp[DIM], local_rdisp[DIM], trv[DIM];
   REAL * fdisp, * rdisp;
   struct simplex_point local_simplex;
   int d, compute_both_witnesses, use_default, first_iteration, max_iterations;
   double oldsqrd;

   assert( NumVertices( obj1)>0 && NumVertices( obj2)>0 );

   use_default = first_iteration = 1;
#ifdef CONSTRUCT_TABLES
   initialise_simplex_distance();
		/* will return immediately if already initialised */
#else
   assert( PRE_DEFINED_TABLE_DIM >= DIM );
#endif /* CONSTRUCT_TABLES */

   compute_both_witnesses = ( wpt1!=0 ) || ( wpt2!=0 ) ||
                            (  tr1!=0 ) || (  tr2!=0 );

   if ( wpt1==0 )
       wpt1 = local_witness1;

   if ( wpt2==0 )
       wpt2 = local_witness2;

   fdisp = IdentityTransform( tr1) ?         displacementv : local_fdisp;
   rdisp = IdentityTransform( tr2) ? reverse_displacementv : local_rdisp;

   if ( simplex==0 ) {
      use_seed = 0;
      simplex = & local_simplex;
   }

   if ( use_seed==0 ) {
      simplex->simplex1[0] = 0;    simplex->simplex2[0] = 0;
      simplex->npts = 1;           simplex->lambdas[0] = ONE;
      simplex->last_best1 = 0;     simplex->last_best2 = 0;
      add_simplex_vertex( simplex, 0,
			  obj1, FirstVertex( obj1), tr1,
			  obj2, FirstVertex( obj2), tr2);
   }
   else {
      /* If we are being told to use this seed point, there
         is a good chance that the near point will be on
         the current simplex.  Besides, if we don't confirm
         that the seed point given satisfies the invariant
         (that the witness points given are the closest points
         on the current simplex) things can and will fall down.
         */
      for ( v=0 ; v<simplex->npts ; v++ )
	add_simplex_vertex( simplex, v,
          obj1, simplex->simplex1[v], tr1,
	  obj2, simplex->simplex2[v], tr2);
   }

   /* Now the main loop.  We first compute the distance between the
      current simplicies, the check whether this gives the globally
      correct answer, and if not construct new simplices and try again.
      */

   max_iterations = NumVertices( obj1)*NumVertices( obj2) ;
      /* in practice we never see more than about 6 iterations. */

   /* Counting the iterations in this way should not be necessary;
      a while( 1) should do just as well. */
   while ( max_iterations-- > 0 ) {

     if ( simplex->npts==1 ) { /* simple case */
       simplex->lambdas[0] = ONE;
     }
     else { /* normal case */
       compute_subterms( simplex);
       if ( use_default ) { 
	 use_default = default_distance( simplex);
       }
       if ( !use_default ) {
	 backup_distance( simplex);
       }
     }

     /* compute at least the displacement vectors given by the
	simplex_point structure.  If we are to provide both witness
	points, it's slightly faster to compute those first.
     */
     if ( compute_both_witnesses ) {
       compute_point( wpt1, simplex->npts, simplex->coords1,
		      simplex->lambdas);
       compute_point( wpt2, simplex->npts, simplex->coords2,
		      simplex->lambdas);
      
       overd( d) {
	 displacementv[ d]         = wpt2[d] - wpt1[d];
	 reverse_displacementv[ d] = - displacementv[d];
       }
     }
     else {
       overd( d) {
	 displacementv[d] = 0;
	 for ( p=0 ; p<simplex->npts ; p++ )
	   displacementv[d] +=
	     DO_MULTIPLY( simplex->lambdas[p],
			  simplex->coords2[p][d] - simplex->coords1[p][d]);
	 reverse_displacementv[ d] = - displacementv[d];
       }
     }
	 
     sqrd = OTHER_DOT_PRODUCT( displacementv, displacementv);

     /* if we are using a c-space simplex with DIM_PLUS_ONE
	points, this is interior to the simplex, and indicates
	that the original hulls overlap, as does the distance 
	between them being too small. */
     if ( sqrd<EPSILON ) {
       simplex->error = EPSILON;
       return sqrd;                 
     }

     if ( ! IdentityTransform( tr1) )
       ApplyInverseRotation( tr1,         displacementv, fdisp);

     if ( ! IdentityTransform( tr2) )
       ApplyInverseRotation( tr2, reverse_displacementv, rdisp);

     /* find the point in obj1 that is maximal in the
	direction displacement, and the point in obj2 that
	is minimal in direction displacement;
     */
     maxp = support_function(
		  obj1,
		  ( use_seed ? simplex->last_best1 : InvalidVertexID),
		  &maxv, fdisp
		  );

     minp = support_function(
		  obj2,
		  ( use_seed ? simplex->last_best2 : InvalidVertexID),
		  &minus_minv, rdisp
		  );

     /* Now apply the G-test on this pair of points */

     INCREMENT_G_TEST_COUNTER;

     g_val = sqrd + maxv + minus_minv;

     if ( ! IdentityTransform( tr1) ) {
       ExtractTranslation( tr1, trv);
       g_val += OTHER_DOT_PRODUCT(         displacementv, trv);
     }

     if ( ! IdentityTransform( tr2) ) {
       ExtractTranslation( tr2, trv);
       g_val += OTHER_DOT_PRODUCT( reverse_displacementv, trv);
     }

     if ( g_val < 0.0 )  /* not sure how, but it happens! */
       g_val = 0;

     if ( g_val < EPSILON ) {
       /* then no better points - finish */
       simplex->error = g_val;
       return sqrd;
     }

     /* check for good calculation above */
     if ( (first_iteration || (sqrd < oldsqrd))
	  && (simplex->npts <= DIM ) ) {
       /* Normal case: add the new c-space points to the current
	  simplex, and call simplex_distance() */
       simplex->simplex1[ simplex->npts] = simplex->last_best1 = maxp;
       simplex->simplex2[ simplex->npts] = simplex->last_best2 = minp;
       simplex->lambdas[ simplex->npts] = ZERO;
       add_simplex_vertex( simplex, simplex->npts,
			   obj1, maxp, tr1,
			   obj2, minp, tr2);
       simplex->npts++;
       oldsqrd = sqrd;
       first_iteration = 0;
       use_default = 1;
       continue;
     }

     /* Abnormal cases! */ 
     if ( use_default ) {
       use_default = 0;
     }
     else { /* give up trying! */
       simplex->error = g_val;
       return sqrd;
     }
   } /* end of `while ( 1 )' */

   return 0.0; /* we never actually get here, but it keeps some
                  fussy compilers happy.  */
}

