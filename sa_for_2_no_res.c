#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

double ff(int x,int y,int kk); //The energy function.
double temp(int m); //The cooling schedule (temperature).



int main()
{
//Parameters:
        int kk = 2; //Right side of the diophantine equation.
        int lbx = -pow(10,4); //Lower bound for the x coordinate of a state.
        int ubx = 0; //Upper bound for the x coordinate.
        int lby = 0; //Lower bound for the y coordinate of a state.
        int uby = pow(10,4); //Upper bound for the y coordinate.
        double thr = pow(10,-5); //Solution threshold. Used to decide when to check whether a solution has been found.

//Variables:
	time_t t0 = time(NULL); //Time used to seed the RNG.
	int x,y,xn,yn; //Coordinates of the current and new states.
	double z,zn; //Energy values for the above.
	int m; //Time.
	int a,b,c,d; //Variables used in the generation of a new state.
	double prob; //Probability for accepting a transition to the already generated state.
	int zt; //Auxiliary integer varialbe.
	double r; //Auxiliary double variable.

//Seeding the RNG:
	srand48(t0);

//Random initial state and its energy:
	x = lround(drand48()*(ubx-lbx)+lbx); y = lround(drand48()*(uby-lby)+lby);
	z = ff(x,y,kk);

	if(z <= thr) //Solution check.
	{
		zt = lround(cbrt(kk-pow(x,3)-pow(y,3)));
		if(pow(x,3)+pow(y,3)+pow(zt,3) == kk)
		{
			printf("{%d,%d,%d}\n",x,y,zt);
			goto end;
		}
	}

//Iterations:
	m = 1; while(1)
	{
//Confinement:
		a = -10; b = 10; c = -10; d = 10;
		if(x+a <= lbx) { a = lbx-x; }
		if(x+b >= ubx) { b = ubx-x; }
		if(y+c <= lby) { c = lby-y; }
		if(y+d >= uby) { d = uby-y; }

//State generation:
		roll:
		xn = x+lround(drand48()*(b-a)+a);
		yn = y+lround(drand48()*(d-c)+c);
		if((xn == x) && (yn == y)) { goto roll; }
		zn = ff(xn,yn,kk);

//State transition:
		if(zn <= z)
		{
			x = xn;
			y = yn;
			z = zn;
		}
		else
		{
			prob = exp((z-zn)/temp(m));
			r = drand48();
			if(r <= prob)
			{
				x = xn;
				y = yn;
				z = zn;
			}
		}

		if(z <= thr) //Solution check.
		{
			zt = lround(cbrt(kk-pow(x,3)-pow(y,3)));
			if(pow(x,3)+pow(y,3)+pow(zt,3) == kk)
			{
				printf("{%d,%d,%d}\n",x,y,zt);
				goto end;
			}
		}
		m = m+1;
	}

	end:

	return 0;
}



double ff(int x,int y,int kk)
{
	double z;

	z = fabs(cbrt(kk-pow(x,3)-pow(y,3))-lround(cbrt(kk-pow(x,3)-pow(y,3))));

	return z;
}



double temp(int m)
{
	double z;

	z = 1.0/(log(m)+0.01);

	return z;
}
