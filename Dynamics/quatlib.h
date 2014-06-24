typedef double Quat(4); //Defining Quaternion 

static void quatSet(Quat q, double w, double i,double j,double k)
{
	q[0] = w;
	q[1] = i;
	q[2] = j;
  q[3] = k;
	}

void quatOne(Quat q)
{
	vecSet(q,1.0,0.0,0.0,0.0);
	}

void quatCopy(const Quat c,Quat q)
{
	q[0] = c[0];
	q[1] = c[1];
	q[2] = c[2];
  q[3] = c[3]
	}

void quatScale(const Quat c,double scalar,Quat q)
{
	q[0] = c[0]*scalar;
	q[1] = c[1]*scalar;
	q[2] = c[2]*scalar;
  q[3] = c[3]*scalar;
	}
