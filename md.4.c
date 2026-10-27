#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define kB 0.08617718  /*  kB= 1./11.604=0.086 meV/K  */
//#define CUTOFF 3.0
/* UNITS: meV, ps, Ang -- mass is in units of ~1.602e-26 kg ~= 9.648 amu  */

double CUTOFF = 3.0;
int N;             //number of atoms
char *atomtype;    //atomic species (1 char)
double dt;         //timestep
double mass;       //mass
double kspring;    //spring constant
double d0;         //equilibrium spring distance
double c_rep;      //repulsive interaction of scatterers
double kbox;       //spring constant of box wall
double side,halfside;//side of the box
long int nsteps;   //number of simulation steps
long int nstep;    //actual step, changing during simulation
long int print;    //print every these steps
double *x,*y,*z,*Vx,*Vy,*Vz,*Fx,*Fy,*Fz;
double *xprev,*yprev,*zprev,*xtmp,*ytmp,*ztmp; // to store the old coordinates for Verlet
double epot,ekin,etot,Tm,Tw;
char inputxyz[30],trajfile[30];
FILE *ftraj;

//===================
// CLOSEST NEIGHBOURS
//===================
struct myNeighbours{
	int num; // total neighbours
	int *closest; // array of idxs of neighbours
};
struct myNeighbours *nlist = NULL; // array of all neighbours

double square(double val){ return val*val; }
double cube(double val){ return val*val*val; }

void read_parameters(char *filename){
	char dummy[100];
	FILE *fin;
	void error() { fprintf(stderr,"ERROR in input parameters file %s\n",filename); exit(1); }

	if((fin=fopen(filename,"r"))==NULL) error();
	if( fscanf(fin,"%s %ld ", dummy, &nsteps  ) != 2 ) error();
	if( fscanf(fin,"%s %ld ", dummy, &print   ) != 2 ) error();
	if( fscanf(fin,"%s %lf ", dummy, &mass    ) != 2 ) error();
	if( fscanf(fin,"%s %lf ", dummy, &dt      ) != 2 ) error();
	if( fscanf(fin,"%s %lf ", dummy, &kspring ) != 2 ) error();
	if( fscanf(fin,"%s %lf ", dummy, &d0      ) != 2 ) error();
	if( fscanf(fin,"%s %s  ", dummy, inputxyz ) != 2 ) error();
	if( fscanf(fin,"%s %s  ", dummy, trajfile ) != 2 ) error();
	if( fscanf(fin,"%s %lf ", dummy, &c_rep   ) != 2 ) error();
	if( fscanf(fin,"%s %lf ", dummy, &side    ) != 2 ) error();
	if( fscanf(fin,"%s %lf ", dummy, &kbox    ) != 2 ) error();

	//printout parameters to check
	fprintf(stderr,"nsteps %ld\n" ,nsteps);
	fprintf(stderr,"print %ld\n"  ,print);
	fprintf(stderr,"mass %g\n"    ,mass);
	fprintf(stderr,"timestep %g\n",dt);
	fprintf(stderr,"kspring %g\n" ,kspring);
	fprintf(stderr,"d0 %g\n"      ,d0);
	fprintf(stderr,"inputxyz %s\n",inputxyz);
	fprintf(stderr,"trajfile %s\n",trajfile);
	fprintf(stderr,"c_rep %g\n"   ,c_rep);
	fprintf(stderr,"side %g\n"    ,side);
	fprintf(stderr,"kbox %g\n"    ,kbox);

	fclose(fin);
	halfside=side/2.;

	fprintf(stderr,"\nAll parameters read from file\n");
}

void readxyz(){
	int i,n=0,readCharCount;
	char dummy[300];
	FILE *fin;
	void error() { fprintf(stderr,"ERROR in input structure file %s, %c %g %g %g\n",inputxyz,atomtype[i],x[i],y[i],z[i]); exit(1); }

	if((fin=fopen(inputxyz,"r"))==NULL) error();

	//the first line should contain the number of atoms
	fgets(dummy,sizeof(dummy),fin);
	if( sscanf(dummy,"%d",&N) != 1 ) error();

	if( (x=(double *)  calloc(N, sizeof(double))) == NULL ) n++;
	if( (y=(double *)  calloc(N, sizeof(double))) == NULL ) n++;
	if( (z=(double *)  calloc(N, sizeof(double))) == NULL ) n++;
	if( (xprev=(double *) calloc(N, sizeof(double))) == NULL ) n++;
	if( (yprev=(double *) calloc(N, sizeof(double))) == NULL ) n++;
	if( (zprev=(double *) calloc(N, sizeof(double))) == NULL ) n++;
	if( (xtmp=(double *) calloc(N, sizeof(double))) == NULL ) n++;
	if( (ytmp=(double *) calloc(N, sizeof(double))) == NULL ) n++;
	if( (ztmp=(double *) calloc(N, sizeof(double))) == NULL ) n++;
	if( (Vx=(double *) calloc(N, sizeof(double))) == NULL ) n++;
	if( (Vy=(double *) calloc(N, sizeof(double))) == NULL ) n++;
	if( (Vz=(double *) calloc(N, sizeof(double))) == NULL ) n++;
	if( (Fx=(double *) calloc(N, sizeof(double))) == NULL ) n++;
	if( (Fy=(double *) calloc(N, sizeof(double))) == NULL ) n++;
	if( (Fz=(double *)  calloc(N, sizeof(double))) == NULL ) n++;
	if( (atomtype=(char *) calloc(N,sizeof(char))) == NULL ) n++;
	if(n>0){fprintf(stderr,"Allocation error %d\n",n); exit(1);}

	// the second line is a comment, we skip it
	fgets(dummy,sizeof(dummy),fin);

	for(i=0;i<N;i++){
		//read positions
		fgets(dummy,sizeof(dummy),fin);
		if( sscanf(dummy,"%s %lf %lf %lf%n", &atomtype[i],&x[i],&y[i],&z[i],&readCharCount) !=4 ) error();
		//read velocities if present, otherwise set them to zero
		if( sscanf(dummy+readCharCount,"%lf %lf %lf", &Vx[i],&Vy[i],&Vz[i]) !=3 ) Vx[i]=Vy[i]=Vz[i]=0.;
	}
	fclose(fin);
	fprintf(stderr,"Simulation with %d atoms\n",N);
}

double distx(int i, int j){
	double dX= x[i]-x[j];
	if(kbox>0.) return dX;              //no PBC
	else return dX-round(dX/side)*side; //yes PBC
}

double disty(int i, int j){
	double dY= y[i]-y[j];
	if(kbox>0.) return dY;              //no PBC
	else return dY-round(dY/side)*side; //yes PBC
}

double distz(int i, int j){
	double dZ= z[i]-z[j];
	if(kbox>0.) return dZ;              //no PBC
	else return dZ-round(dZ/side)*side; //yes PBC
}

double dist(int i, int j){
	return sqrt( square(distx(i,j)) + square(disty(i,j)) + square(distz(i,j)) );
}

// neighbouring list


void build_closest(){
	double D0 = CUTOFF + 0.5;
	// allocate empty
	if(nlist == NULL){
		calloc(N,sizeof(struct myNeighbours));
	}

	// loop over atoms
	for(int i=0;i<N;i++){
		// count neighbours
		int counter = 0;
		for(int j=0;j<N;j++){
			if(j==i) continue;
			
		}
		// allocate right size array
	}
}



void calc_forces(){
	int i,j; double effe,d,d4,d12;
	//reset the forces
	for (i=0;i<N;i++) { Fx[i]=Fy[i]=Fz[i]=0; }

	//chain
	for (i=0;i<N-1;i++)
	    if(atomtype[i]=='B' && atomtype[i+1]=='B') {
		    effe = -kspring*(dist(i,i+1) - d0)/dist(i,i+1);
		    Fx[i] += effe*distx(i,i+1); Fx[i+1] += -effe*distx(i,i+1);
		    Fy[i] += effe*disty(i,i+1); Fy[i+1] += -effe*disty(i,i+1);
		    Fz[i] += effe*distz(i,i+1); Fz[i+1] += -effe*distz(i,i+1);
	    }

	//scatterers
	for (i=0;i<N;i++)
	    for(j=0;j<i;j++){
		if(atomtype[i]=='B' && atomtype[j]=='B') continue;

		d   = dist(i,j);
		// cutoff 
		if(d > CUTOFF) continue;

		d4  = d*d*d*d;
		d12 = d4*d4*d4;
		effe = 10.*c_rep/d12;

		Fx[i] +=  effe*distx(i,j);  Fx[j] += -effe*distx(i,j);
		Fy[i] +=  effe*disty(i,j);  Fy[j] += -effe*disty(i,j);
		Fz[i] +=  effe*distz(i,j);  Fz[j] += -effe*distz(i,j);
		}

	//box (if no PBC)
	if(kbox>0.) for (i=0;i<N;i++){
	  if ( x[i] < -halfside ) Fx[i] -= kbox*(x[i] + halfside);
	  if ( x[i] >  halfside ) Fx[i] -= kbox*(x[i] - halfside);
	  if ( y[i] < -halfside ) Fy[i] -= kbox*(y[i] + halfside);
	  if ( y[i] >  halfside ) Fy[i] -= kbox*(y[i] - halfside);
	  if ( z[i] < -halfside ) Fz[i] -= kbox*(z[i] + halfside);
	  if ( z[i] >  halfside ) Fz[i] -= kbox*(z[i] - halfside);
       }

}

void calc_energies(){
	int i,j,Nm=0,Nw=0; double d,d5,d10;

	ekin=Tm=Tw=0;
	for(i=0;i<N;i++)
	    if(atomtype[i]=='B'){ Tm+=( Vx[i]*Vx[i] + Vy[i]*Vy[i] + Vz[i]*Vz[i] ); Nm++; }
	    else                { Tw+=( Vx[i]*Vx[i] + Vy[i]*Vy[i] + Vz[i]*Vz[i] ); Nw++; }

	ekin = 0.5 * mass * ( Tm + Tw );
	Tm= mass*Tm/(3.*Nm*kB);
	Tw= mass*Tw/(3.*Nw*kB);

	epot=0;

	//chain
	for (i=0;i<N-1;i++)
	    if(atomtype[i]=='B' && atomtype[i+1]=='B')
		epot += 0.5*kspring*square( dist(i,i+1) - d0 );

        //scatterers
	for (i=0;i<N;i++)
	    for(j=0;j<i;j++){
		if(atomtype[i]=='B' && atomtype[j]=='B') continue;
		d   = dist(i,j);
		// cutoff 
		if(d > CUTOFF) continue;
		d5  = d*d*d*d*d;
		d10 = d5*d5;
		epot += c_rep/d10;
	    }

	//box (if no PBC)
	if(kbox>0.) for (i=0;i<N;i++){
	  if ( x[i] < -halfside ) epot += 0.5*kbox*square(x[i] + halfside);
	  if ( x[i] >  halfside ) epot += 0.5*kbox*square(x[i] - halfside);
	  if ( y[i] < -halfside ) epot += 0.5*kbox*square(y[i] + halfside);
	  if ( y[i] >  halfside ) epot += 0.5*kbox*square(y[i] - halfside);
	  if ( z[i] < -halfside ) epot += 0.5*kbox*square(z[i] + halfside);
	  if ( z[i] >  halfside ) epot += 0.5*kbox*square(z[i] - halfside);
       }

	etot=ekin+epot;
}

void velocity_Verlet_integrator(){
	int i; static int once=1;
	if(once){ once=0; calc_forces(); }

	// update velocities(1)
	for (i=0;i<N;i++){
		Vx[i] += dt*Fx[i]/(2*mass);
		Vy[i] += dt*Fy[i]/(2*mass);
		Vz[i] += dt*Fz[i]/(2*mass);
	}

	//update positions
	for (i=0;i<N;i++){
		x[i] += Vx[i]*dt;
		y[i] += Vy[i]*dt;
		z[i] += Vz[i]*dt;
	}

	calc_forces(); // update the forces

	// update velocities(2)
	for (i=0;i<N;i++){
		Vx[i] += dt*Fx[i]/(2*mass);
		Vy[i] += dt*Fy[i]/(2*mass);
		Vz[i] += dt*Fz[i]/(2*mass);
	}
}

void write_traj(){
	static int once=1;
	int i;
	if(once){
		once=0;
		if((ftraj=fopen(trajfile,"w"))==NULL) { fprintf(stderr,"ERROR while saving trajectory file\n"); exit(1); }
	}

	fprintf(ftraj,"%d\nLattice=\"%g 0.0 0.0  0.0 %g 0.0  0.0 0.0 %g\"\tTIME: %g\n",N,side,side,side,dt*nstep);
	for(i=0;i<N;i++){
		fprintf(ftraj,"%c\t%.6f\t%.6f\t%.6f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%6e\n" ,atomtype[i],x[i],y[i],z[i],Vx[i],Vy[i],Vz[i],Fx[i],Fy[i],Fz[i]);
	}
	fflush(ftraj);
}

void printout(){
	static int once=1;
	if(once){
		once=0;
		printf("#time\tekin\t\tepot\t\tetot\t\tchain_length\tTmolecule\tTwater\n");
	}
	calc_energies();
	printf("%g\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",dt*nstep,ekin,epot,etot,dist(0,N-1),Tm,Tw );
	write_traj();
	fflush(stdout);
}


int main(int argc, char *argv[])
{
   // read parameter file
   if (argc < 2) {fprintf(stderr,"Provide parameter file\n"); exit(1);}
   read_parameters(argv[1]);
	 printf("argc = %d",argc);
	 if(argc >=2 && argv[2]){
		CUTOFF = atof(argv[2]);
		fprintf(stderr,"Using CUTOFF radius: %.6f\n",CUTOFF);
	 }


   // read initial conditions
   readxyz();

   // write trajectory and print some information
   printout();

   // dynamics loop
   for (nstep=1;nstep<=nsteps;nstep++)
   {

      velocity_Verlet_integrator(); // one time of dynamics

      // write trajectory and print some information
      if ( nstep % print == 0 ) printout();
   }

   fclose(ftraj);

   return 0;
}


