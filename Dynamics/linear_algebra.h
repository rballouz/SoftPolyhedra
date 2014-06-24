#include<assert.h>
#include<math.h>
typedef double Vec[3];
typedef Vec Mat[3]; 

int sign(double A);
#define sign(A) ((A)<0.?-1:(A)>0.?1:0)
double min(double A,double B);
#define min(A,B) ((A) > (B) ? (B) : (A))
double max(double A,double B);
#define max(A,B) ((A) > (B) ? (A) : (B))

static void vecSet(Vec v,double x,double y,double z);
void vecZero(Vec v); 
void vecCopy(const Vec u,Vec v);
void vecScale(const Vec u,double scalar,Vec v);
void vecAdd(const Vec v1,const Vec v2,Vec v);
void vecSub(const Vec v1,const Vec v2,Vec v);
double vecDot(const Vec v1,const Vec v2);
double vecMagSq(const Vec v);
double vecMag(const Vec v);
void vecNorm(Vec v);
void vecCross(const Vec v1,const Vec v2,Vec v);
void vecGetBasis(Vec a,Vec b,Vec c);
void vecRotate(Vec vPoint,const Vec vAxis,double dAngle);
void vecRotate2(Vec vPoint,const Vec vAxis,double sa,double ca);
void vecGetClosestPointOnLine(const Vec A,const Vec B,const Vec P,Vec Q);
void matCopy(Mat a,Mat b);
void matTransform(Mat m,const Vec u,Vec v);
void matTranspose(Mat a,Mat b);
void matSwapRows(Mat m,int row1,int row2);
void matIdentity(Mat m);
void matDiagonal(const Vec v,Mat m);
double matSumOffDiagElem(Mat m);
double matSumAbsOffDiagElem(Mat m);
void matScale(Mat a,double s,Mat b);
void matMultiply(Mat a,Mat b,Mat c);
void matInverse(Mat mat_in,Mat mat_out);
void jacobi(Mat m,Vec eig_vals,Mat eig_vecs);

static void vecSet(Vec v,double x,double y,double z)
{
	v[0] = x;
	v[1] = y;
	v[2] = z;
	}

void vecZero(Vec v)
{
	vecSet(v,0.0,0.0,0.0);
	}

void vecCopy(const Vec u,Vec v)
{
	v[0] = u[0];
	v[1] = u[1];
	v[2] = u[2];
	}

void vecScale(const Vec u,double scalar,Vec v)
{
	v[0] = u[0]*scalar;
	v[1] = u[1]*scalar;
	v[2] = u[2]*scalar;
	}

void vecAdd(const Vec v1,const Vec v2,Vec v)
{
	v[0] = v1[0] + v2[0];
	v[1] = v1[1] + v2[1];
	v[2] = v1[2] + v2[2];
	}

void vecSub(const Vec v1,const Vec v2,Vec v)
{
	v[0] = v1[0] - v2[0];
	v[1] = v1[1] - v2[1];
	v[2] = v1[2] - v2[2];
	}

double vecDot(const Vec v1,const Vec v2)
{
	return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
	}

double vecMagSq(const Vec v)
{
	return vecDot(v,v);
	}

double vecMag(const Vec v)
{
	return sqrt(vecMagSq(v));
	}

void vecNorm(Vec v)
{
	double mag = vecMag(v);
	assert(mag > 0.0);
	vecScale(v,1.0/mag,v);
	}

void vecCross(const Vec v1,const Vec v2,Vec v)
{
	v[0] = v1[1]*v2[2] - v1[2]*v2[1];
	v[1] = v1[2]*v2[0] - v1[0]*v2[2];
	v[2] = v1[0]*v2[1] - v1[1]*v2[0];
	}

void vecGetBasis(Vec a,Vec b,Vec c)
{
	/* Given vec "a", this routine returns orthonormal basis (a,b,c) */

	Vec v,ctmp;
	Mat I;
	double proj;

	matIdentity(I); /* unit mat */

	/* Get spanning set...first guess: choose I[1] (y-hat) and I[2] (z_hat) as 2nd & 3rd vecs */

	vecCopy(I[1],b);
	vecCopy(I[2],c);

	/* If "a" is actually null, set it to I[0] (x-hat) and return */

	if (a[0] == 0.0 && a[1] == 0.0 && a[2] == 0.0) {
		vecCopy(I[0],a);
		return;
	}

	/*
	** If "a" does not have an x component, make 2nd vec I[0] (x-hat). If in
	** addition "a" does not have a y component, make 3rd vec I[1] (y-hat).
	** Now a, b, and c span 3-space.
	*/

	if (a[0] == 0.0) {
		vecCopy(I[0],b);
		if (a[1] == 0.0)
			vecCopy(I[1],c);
	}

	/* Construct orthonormal basis using the Gram-Schmidt orthonormalization process */

	vecNorm(a); /* first basis vec */

	/* Construct second basis vec */

	proj = vecDot(a,b);
	vecScale(a,proj,v);
	vecSub(b,v,b);
	vecNorm(b);

	/* Construct third basis vec */

	proj = vecDot(a,c);
	vecScale(a,proj,v);
	vecSub(c,v,ctmp);
	proj = vecDot(b,c);
	vecScale(b,proj,v);
	vecSub(ctmp,v,c);
	vecNorm(c);
	}

void vecRotate(Vec vPoint,const Vec vAxis,double dAngle)
{
	/*
	** Rotates vPoint around vAxis (a unit vec) by dAngle radians.
	** vPoint is overwritten.  Based on algorithm by Glenn Murray at
	** http://inside.mines.edu/~gmurray/ArbitraryAxisRotation/.
	*/

	double sa,ca;

	sa = sin(dAngle);
	ca = cos(dAngle);
	vecRotate2(vPoint,vAxis,sa,ca);
	}

void vecRotate2(Vec vPoint,const Vec vAxis,double sa,double ca)
{
	/*
	** Rotates vPoint around vAxis (a unit vec) by the angle that
	** is given by arccosine(ca) [and arcsin(sa)].
	** vPoint is overwritten.  Based on algorithm by Glenn Murray at
	** http://inside.mines.edu/~gmurray/ArbitraryAxisRotation/.
	*/

	double dot,x,y,z,u,v,w;

	dot = vecDot(vPoint,vAxis);
	x = vPoint[0];
	y = vPoint[1];
	z = vPoint[2];
	u = vAxis[0];
	v = vAxis[1];
	w = vAxis[2];
	vecSet(vPoint,
			  u*dot + (x*(v*v + w*w) - u*(v*y + w*z))*ca + (-w*y + v*z)*sa,
			  v*dot + (y*(u*u + w*w) - v*(u*x + w*z))*ca + ( w*x - u*z)*sa,
			  w*dot + (z*(u*u + v*v) - w*(u*x + v*y))*ca + (-v*x + u*y)*sa);
	}

void vecGetClosestPointOnLine(const Vec A,const Vec B,const Vec P,Vec Q)
{
	/*
	** Given points A & B on a line, and point P, this routine
	** returns the point Q on the line closest to P.
	*/

	Vec v1,v2;
	double T;

	vecSub(A,B,v1);
	vecSub(P,B,v2);
	T = vecDot(v1,v2)/vecMagSq(v1);
	vecScale(v1,T,v1);
	vecAdd(B,v1,Q);
}

void vecGetClosestPointOnSegment(const Vec A,const Vec B,const Vec P,Vec Q)
{
	/*
	** Given segment (finite line) AB and the point P, this
	** routine returns the point Q on the segment closest to P.
	*/

	Vec v1,v2;
	double T;

	vecSub(A,B,v1);
	vecSub(P,B,v2);
	T = vecDot(v1,v2)/vecMagSq(v1);
	T = min(T,1.0);
	T = max(T,0.0);
	vecScale(v1,T,v1);
	vecAdd(B,v1,Q);
}

void matCopy(Mat a,Mat b)
{
	vecCopy(a[0],b[0]);
	vecCopy(a[1],b[1]);
	vecCopy(a[2],b[2]);
	}

void matTransform(Mat m,const Vec u,Vec v)
{
	v[0] = vecDot(m[0],u);
	v[1] = vecDot(m[1],u);
	v[2] = vecDot(m[2],u);
	}

void matTranspose(Mat a,Mat b)
{
	b[0][0] = a[0][0];
	b[0][1] = a[1][0];
	b[0][2] = a[2][0];
	b[1][0] = a[0][1];
	b[1][1] = a[1][1];
	b[1][2] = a[2][1];
	b[2][0] = a[0][2];
	b[2][1] = a[1][2];
	b[2][2] = a[2][2];
	}

void matSwapRows(Mat m,int row1,int row2)
{
	Vec tmp;

	vecCopy(m[row2],tmp);
	vecCopy(m[row1],m[row2]);
	vecCopy(tmp,m[row1]);
	}

void matIdentity(Mat m)
{
	m[0][0] = 1.0;
	m[0][1] = 0.0;
	m[0][2] = 0.0;
	m[1][0] = 0.0;
	m[1][1] = 1.0;
	m[1][2] = 0.0;
	m[2][0] = 0.0;
	m[2][1] = 0.0;
	m[2][2] = 1.0;
	}

void matDiagonal(const Vec v,Mat m)
{
	m[0][0] = v[0];
	m[0][1] = 0.0;
	m[0][2] = 0.0;
	m[1][0] = 0.0;
	m[1][1] = v[1];
	m[1][2] = 0.0;
	m[2][0] = 0.0;
	m[2][1] = 0.0;
	m[2][2] = v[2];
	}

double matSumOffDiagElem(Mat m)
{
	return m[0][1] + m[0][2] + m[1][0] +
		m[1][2] + m[2][0] + m[2][1];
	}

double matSumAbsOffDiagElem(Mat m)
{
	return fabs(m[0][1]) + fabs(m[0][2]) + fabs(m[1][0]) +
		fabs(m[1][2]) + fabs(m[2][0]) + fabs(m[2][1]);
	}

void matScale(Mat a,double s,Mat b)
{
	vecScale(a[0],s,b[0]);
	vecScale(a[1],s,b[1]);
	vecScale(a[2],s,b[2]);
	}

void matMultiply(Mat a,Mat b,Mat c)
{
	c[0][0] = a[0][0]*b[0][0] + a[0][1]*b[1][0] + a[0][2]*b[2][0];
	c[0][1] = a[0][0]*b[0][1] + a[0][1]*b[1][1] + a[0][2]*b[2][1];
	c[0][2] = a[0][0]*b[0][2] + a[0][1]*b[1][2] + a[0][2]*b[2][2];
	c[1][0] = a[1][0]*b[0][0] + a[1][1]*b[1][0] + a[1][2]*b[2][0];
	c[1][1] = a[1][0]*b[0][1] + a[1][1]*b[1][1] + a[1][2]*b[2][1];
	c[1][2] = a[1][0]*b[0][2] + a[1][1]*b[1][2] + a[1][2]*b[2][2];
	c[2][0] = a[2][0]*b[0][0] + a[2][1]*b[1][0] + a[2][2]*b[2][0];
	c[2][1] = a[2][0]*b[0][1] + a[2][1]*b[1][1] + a[2][2]*b[2][1];
	c[2][2] = a[2][0]*b[0][2] + a[2][1]*b[1][2] + a[2][2]*b[2][2];
	}

void matInverse(Mat mat_in,Mat mat_out)
{
	int row_to_pivot; /* Actual row that will be pivoted upon. */
	int pivot_row;    /* Current row of the mat.  We will make row_to_pivot
						 equal to this by moving that row to this position. */
	double max;        /* The maximum value in the correct place.  The row we
						 pivot is based on the row with the maximum value. */
	int x,y;
	double scale,custom_scale;
	Mat m;

	matCopy(mat_in,m);
	matIdentity(mat_out);

	/* Which row of the mat are we on? */
	for (pivot_row = 0; pivot_row < 3; ++pivot_row) {
		/* First, we identify the largest element in the column in question and
		 *  then move it to be in the current position. */

		/* Start looking at the current row. */
		row_to_pivot = pivot_row;
		max = fabs(m[pivot_row][pivot_row]);
		for (y = pivot_row + 1; y < 3; ++y) {
			/* If element is large, mark it as maximum. */
			if (fabs(m[y][pivot_row]) > max) {
				row_to_pivot = y;
				max = fabs(m[y][pivot_row]);
				}
			}

		/* Okay, we now know what row to pivot.  Move it to the right place. */
		if (row_to_pivot != pivot_row) {
			matSwapRows(m,row_to_pivot,pivot_row);
			matSwapRows(mat_out,row_to_pivot,pivot_row);
			}

		/* Row is in place, now we move on to the actual pivot. */

		/* First we compute how much we need to scale the pivot row by, which is
		 *  1 / pivot_element. */
		scale = 1.0 / m[pivot_row][pivot_row];

		/* Next, we compute the base pivot_row.  We will just leave this in the
		 *  mat at the appropriate spot. */
		for (x = 0; x < 3; ++x) {
			m[pivot_row][x] *= scale;
			mat_out[pivot_row][x] *= scale;
			}

		/* Finally, for all rows except the pivot row, we subtract the pivot row 
		 *  of the mat.  This is the actual pivot.  We will need a custom scale
		 *  factor for each row to eliminate the element. */
		for (y = 0; y < 3; ++y) {
			if (y != pivot_row) {
				/* Get the custom scale for this row. */
				custom_scale = m[y][pivot_row];

				/* Now pivot. */
				for (x = 0; x < 3; ++x) {
					m[y][x] -= custom_scale * m[pivot_row][x];
					mat_out[y][x] -= custom_scale * mat_out[pivot_row][x];
					}
				}
			}
		}
	}
#define MAX_JACOBI_SWEEPS 50
#define JACOBI_N 3 /* hard-wired for 3D */

void jacobi(Mat m,Vec eig_vals,Mat eig_vecs)
{
	int y,x;                  /* Row, column of mat for element to be 
							   *  eliminated. */
	int j;                    /* Column of mat. */
	int sweep_count;          /* Current number of sweeps made */
	double off_diag_sum;       /* Sum of off diagonal elements. */
	double threshold;          /* Threshold for performing the rotation.  For the
							   *  first three sweeps, this is set equal to 1/5
							   *  the off diagonal sum divided by n^2.  After
							   *  four sweeps, this is 0. */
	Mat a;                 /* Copy of input mat. */
	double s;                  /* Rotation angle sine */
	double c;                  /* Rotation angle cosine */
	double t;
	double tau;
	double theta; 
	double temp;

	/* Make copy of input mat. */
	matCopy(m,a);

	/* Initialize result. */
	matIdentity(eig_vecs);
 
	/* Do Jacobi rotations */
	for (sweep_count = 0; sweep_count < MAX_JACOBI_SWEEPS; ++sweep_count) {
		/* First, check to see if we are done. */
		off_diag_sum = matSumAbsOffDiagElem(a);
		if (off_diag_sum == 0.0) { /* Make that underflow... */
			/* Store results */
			for (x = 0; x < JACOBI_N; ++x)
				eig_vals[x] = a[x][x];
			/* Eig vecs is ready */
			return;
			}

		/* Not done yet, so determine the threshold value for performing the
		 *  rotation. */
		if (sweep_count < 3)
			threshold = 0.2 * off_diag_sum / (JACOBI_N * JACOBI_N);
		else
			threshold = 0.0;
  
		/* Do a sweep; this means do a rotation for each off diag element in the 
		 *  mat.  Here, (x,y) are the coordinates for the element being 
		 *  "eliminated."  We only need eliminate the top triangle of the symmetric
		 *  mat. */
		for (y = 0; y < JACOBI_N - 1; ++y)
			for (x = y + 1; x < JACOBI_N; ++x) {
    
				/* Now, we check to see if we should bother with the rotation.  If the
				 *  element being "eliminated" is too small, then we don't bother.  Too
				 *  small in this case is when the element * 100 is small enough that
				 *  adding it to either the diagonal element in its row or the diagonal
				 *  element in its column is not enough to change the diagonal element
				 *  within machine precision.  We only do this check if we are on the
				 *  fourth or higher sweep.  If we "don't bother" we set the element to
				 *  0 (eliminating it) and then move on. */
				if (sweep_count > 3)
					if (((100.0 * a[y][x] + fabs(a[x][x])) == fabs(a[x][x])) &&
						((100.0 * a[y][x] + fabs(a[y][y])) == fabs(a[y][y])))
						{
						a[y][x] = a[x][y] = 0.0;
						continue;
						}

				/* Now we check the threshold value. */
				if (fabs(a[y][x]) <= threshold)
					continue;
     
				/* Okay, we must perform the rotation.
				 * theta -> t -> c, s, tau -> rotation */

				/* theta, t */

				/* Start theta calculation.  Check size; if too big, get t without
				 *  the need to square theta, otherwise, get t the normal way. */
				theta = a[x][x] - a[y][y];
				/* If theta too big... */
				if ((fabs(theta) + 100.0 * a[y][x]) == fabs(theta)) {
					/* Do t = 1/2/theta (note that two formulas are combined here) */
					t = a[y][x] / theta;
					}
				else { /* theta just right... */
					theta = 0.5 * theta / a[y][x];
					/* t = sgn(theta)/( |theta| + sqrt(theta^2 + 1) ). */
					t = 1.0/(fabs(theta) + sqrt(1.0+theta*theta));
					if (theta < 0.0)
						t = -t;
					}
				/* c, s, tau */
				c = 1.0 / sqrt(1 + t*t);
				s = t * c;
				tau = s / (1.0 + c);

				/* Really do rotation, finally! */
				a[y][y] -= t * a[y][x];
				a[x][x] += t * a[y][x];
				a[y][x] = a[x][y] = 0.0;

				for (j = 0; j < JACOBI_N; ++j)
					if (j != x && j != y) {
						temp = a[j][x];

						/* Do column w/ index x */
						a[j][x] += s*(a[j][y] - tau * a[j][x]);
						a[x][j] = a[j][x];

						/* Do column w/ index y */
						a[j][y] -= s*(temp + tau * a[j][y]);
						a[y][j] = a[j][y];
						}

				/* Update result. */
				for (j = 0; j < JACOBI_N; ++j) {
					temp = eig_vecs[j][x];

					/* Do column w/ index x */
					eig_vecs[j][x] += s*(eig_vecs[j][y] - tau * eig_vecs[j][x]);

					/* Do column w/ index y */
					eig_vecs[j][y] -= s*(temp + tau * eig_vecs[j][y]);  
					}
				}
		}

	/* Should never get here.  Means we had too many sweeps. */
	assert(0);
	}


