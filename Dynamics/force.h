void springforce (float mass, double x, double *n, double *a);
void matPrint(Mat a); //print out a 3x3 matrix

void springforce (float mass, double x, double *n, double *a){

  double a_mag;
  if (x<0) //Determine magnitude of Force
	  a_mag=-x/mass;
  else
    a_mag=0;
  vecScale(n,a_mag,a);//Assign Force magnitudes depending on Angles
}

void matPrint(Mat a){

  int i;
  for (i=0;i<3;i++)
    printf("%f  %f  %f\n",a[i][0],a[i][1],a[i][2]);
}

