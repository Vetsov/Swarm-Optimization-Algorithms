#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

double ff(int x, int y, int kk); //The fitness function.
int sort_and_break(double q[], int n); //A function which chooses the number of a best particle breaking ties randomly.



int main()
{
//Parameters:
        int kk = 2; //Right side of the diophantine equation.
        int lb[2] = {-pow(10,4),0}; //Lower bounds for the search space.
        int ub[2] = {0,pow(10,4)}; //Upper bounds for the search space.
        double thr = pow(10,-5); //Solution threshold. Used to decide when to check whether a solution has been found.
        int s = 50; //Number of particles.
        int dp = (s/5); //Number of particles that need to be at the best swarm best position in order to disperse all particles.

//Variables:
	time_t t0 = time(NULL); //Time used to seed the RNG.
	int pos[s][2]; //Particle positions.
	double fitf[s]; //Fitness functions for the respective particles.
	int nbrs[s][3]; //Neighbours of each particle.
	int sb[2],bsb[2],nb[s][2]; //Position of the swarm best particle, best swarm best position and position of the neighbours best particles.
	double sbfitf,bsbfitf,nbfitf[s]; //Fitness functions for the above.
        int dc; //Dispersion counter to decide when to disperse the particles.
	double nbrsfitf[3],r; //Auxiliary double variables.
        int i,j,zt; //Auxiliary integer variables.

//Seeding the RNG:
	srand48(t0);

//Calculating the neighbours of each particle:
	i = 0; while(i < s)
	{
		nbrs[i][0] = (s+i-1)%s;
		nbrs[i][1] = (s+i)%s;
		nbrs[i][2] = (s+i+1)%s;
		i = i+1;
	}

//Initialisation:
	i = 0; while(i < s)
	{
		pos[i][0] = lround(drand48()*(ub[0]-lb[0])+lb[0]);
		pos[i][1] = lround(drand48()*(ub[1]-lb[1])+lb[1]);
		fitf[i] = ff(pos[i][0],pos[i][1],kk);
		i = i+1;
	}

//Determining the swarm best and best swarm best particles:
	i = sort_and_break(fitf, s);
	sb[0] = pos[i][0]; sb[1] = pos[i][1]; sbfitf = fitf[i];
	bsb[0] = sb[0]; bsb[1] = sb[1]; bsbfitf = sbfitf;

	if(bsbfitf <= thr) //Solution check.
	{
		zt = lround(cbrt(kk-pow(bsb[0],3)-pow(bsb[1],3)));
		if(pow(bsb[0],3)+pow(bsb[1],3)+pow(zt,3) == kk)
		{
			printf("{%d,%d,%d}\n",bsb[0],bsb[1],zt);
			goto end;
		}
	}

//Calculating the best particle among the neighbours:
	i = 0; while(i < s)
	{
		nbrsfitf[0] = fitf[nbrs[i][0]];
		nbrsfitf[1] = fitf[nbrs[i][1]];
		nbrsfitf[2] = fitf[nbrs[i][2]];
		j = sort_and_break(nbrsfitf, 3);
		nb[i][0] = pos[(s+i+j-1)%s][0];
		nb[i][1] = pos[(s+i+j-1)%s][1];
		nbfitf[i] = fitf[(s+i+j-1)%s];
		i = i+1;
	}

//Iterations:
	while(1)
	{
//Position update:
		if(dc >= dp)
		{
			i = 0; while(i < s)
			{
				j = 0; while(j < 2)
				{
					r = drand48();
					pos[i][j] = lround(r*(ub[j]-lb[j])+lb[j]);
					j = j+1;
				}
				i = i+1;
			}
		}
		else
		{
			i = 0; while(i < s)
			{
				j = 0; while(j < 2)
				{
					r = drand48();
					pos[i][j] = lround(nb[i][j]+(r/2.0)*(sb[j]+bsb[j]-2.0*pos[i][j]));
					j = j+1;
				}
				i = i+1;
			}
		}

//Confinement and dispersion condition check:
		i = 0; dc = 0; while(i < s)
		{
			if((pos[i][0] < lb[0]) || (pos[i][0] > ub[0]) || (pos[i][1] < lb[1]) || (pos[i][1] > ub[1]))
			{
				pos[i][0] = lround(drand48()*(ub[0]-lb[0])+lb[0]);
				pos[i][1] = lround(drand48()*(ub[1]-lb[1])+lb[1]);
			}
			fitf[i] = ff(pos[i][0],pos[i][1],kk);
			if((pos[i][0] == bsb[0]) && (pos[i][1] == bsb[1])) { dc = dc+1; }
			i = i+1;
		}

//Swarm best, best swarm best and neighbours best updates:
		i = sort_and_break(fitf, s);
		sb[0] = pos[i][0]; sb[1] = pos[i][1]; sbfitf = fitf[i];
		r = drand48();
		if(sbfitf < bsbfitf)
		{
			bsb[0] = sb[0]; bsb[1] = sb[1]; bsbfitf = sbfitf;

			if(bsbfitf <= thr) //Solution check.
			{
				zt = lround(cbrt(kk-pow(bsb[0],3)-pow(bsb[1],3)));
				if(pow(bsb[0],3)+pow(bsb[1],3)+pow(zt,3) == kk)
				{
					printf("{%d,%d,%d}\n",bsb[0],bsb[1],zt);
					goto end;
				}
			}

		}

		if((sbfitf > bsbfitf) && (r <= (1-((sbfitf-bsbfitf)/0.5))))
		{
			bsb[0] = sb[0]; bsb[1] = sb[1]; bsbfitf = sbfitf;

			if(bsbfitf <= thr) //Solution check.
			{
				zt = lround(cbrt(kk-pow(bsb[0],3)-pow(bsb[1],3)));
				if(pow(bsb[0],3)+pow(bsb[1],3)+pow(zt,3) == kk)
				{
					printf("{%d,%d,%d}\n",bsb[0],bsb[1],zt);
					goto end;
				}
			}

		}

//Calculating the best particle among the neighbours:
		i = 0; while(i < s)
		{
			nbrsfitf[0] = fitf[nbrs[i][0]];
			nbrsfitf[1] = fitf[nbrs[i][1]];
			nbrsfitf[2] = fitf[nbrs[i][2]];
			j = sort_and_break(nbrsfitf, 3);
			nb[i][0] = pos[(s+i+j-1)%s][0];
			nb[i][1] = pos[(s+i+j-1)%s][1];
			nbfitf[i] = fitf[(s+i+j-1)%s];
			i = i+1;
		}
	}

	end:

	return 0;
}



double ff(int x, int y, int kk)
{
	double z;

	z = fabs(cbrt(kk-pow(x,3)-pow(y,3))-lround(cbrt(kk-pow(x,3)-pow(y,3))));

	return z;
}

int sort_and_break(double q[], int n)
{
        int k,i;
        int c = 0;
        int b[n];
        double gv = q[0];

	i = 0; while(i < n) { b[i] = 0; i = i+1; }

        i = 1; while(i < n)
        {
                if(q[i] == gv)
                {
                        c = c+1;
                        b[c] = i;
                }
                if(q[i] < gv)
                {
                        gv = q[i];
                        c = 0;
                        b[c] = i;
                }
                i = i+1;
        }

        k = lround(drand48()*c);

        return b[k];
}
