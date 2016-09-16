#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<malloc.h>
#include</usr/include/complex.h>
#include<time.h>
#include"blas.h"
#include<fftw3.h>


double dipole_integral;
double n, Energy, l;
double pi = 3.14159265358979, m = 1, hbar = 1, ec=1;
double h=6.626070040e-34;
double Lau = 5.2917721092e-11;
double L=2e-9/Lau;
// 1 Atomic unit of time in SI (s)
double time_au = 2.418884326505e-17;
// 1 atomic unit electric field in SI (V/m)
double field_au = 5.14220652e11;
double E_attraction;
//
//Fourier transform arrays
fftw_complex *dm_corr_func;

// // Misc
double dt;
double kOhb;
double omega, tau, start, tim, total_time;
int lasertime, itermax, ftlength, extra_pts;
int dm_fftw_iter;
int nocc, nvirt, next, nact;

// Arrays for storing electric field from FDTD output
double *fdtd_time;
double **fdtd_field;

// Scale electric field by this quantity
double E0 = sqrt(1e15*377)/field_au;

void Psin(double *XA, double *PSIA, int n, int numpts);

// Calculates dipole integral along X, Y and Z planes
double TransitionDipole_X(int nx, int nxp, int ny, int nyp, int nz, int nzp);
double TransitionDipole_Y(int ny, int nyp, int nx, int nxp, int nz, int nzp);
double TransitionDipole_Z(int nz, int nzp, int nx, int nxp, int ny, int nyp);

double En(int n);
double DeltaE(int fin, int init);
double GroundtstateEnergyCalc ( int numorbs, int *E);
double ExcitedstateEnergyCalc ( int orbital_i, int orbital_a, double energy_total, int *E);
void OrderPsis(int norbs, int *E, int **MO);
void HdotC(int end, int iter, double m, double *H0, double *mu, double *C, double *HdC);
void Propagate(int sup, int ncis, double *H0, double *mu, int **bas, int *E, int *init);
void PrintFrontier(int iter, int ncis, FILE *fp, double complex *c, int **bas, int *E);
double EFieldTime(int numlam, double tim, double *Pl, double *w_au);
void IncandescentLight(int numlam, double lam1, double lam2, double *Pl, double *w_au);
void buildBasis(int no, int nv, int **abaa,int **ibaa);
void OrderPsis3D(int norbs, int *E, int **MO);
void AbsorptionSpectrum(int iter, int extra_pts);
int ReadElectricField(int component, char *filename, double *fdtd_dt, double *fdtd_total_time);

int num, *E, **MO;
int **SQMAT_2_INT(int dim);
int **MAT_INT(int dim1, int dim2); // 2d array
int *VEC_INT(int dim);
double *VEC_DOUBLE(int dim);
char *VEC_CHAR(int dim);
double complex *VEC_CDOUBLE(int dim);
double **MAT_DOUBLE(int dim1, int dim2);

double *mx, *my, *mz;
// We will also keep a complex copy of the mu vector to use with the ZDOTC function 
// to compute the dipole moment from Trace(mu*D)
complex double *cmx, *cmy, *cmz;
int main() {

        int nmax, Nels;
        int *E, **MO, **Frontier;

        // it is possible to initialize the wavefunction as a superposition... these two variables define the initial superposition
        // e.g. if sup=1 and init[0] = 1., then the wavefunction is initialized as the ground state configuration... all other init[i] = 0.
        // on the otherhand, if sup=2, then the first two values in the init[] array can be nonzero... e.g. init[0] = sqrt(1/2) and init[1] = sqrt(1/2)
        // we will assume that only the first sup configurations have nonzero weight
        int sup, *init;

        int i, j, a, b, ncis, nstates, idx, jidx;
        int nxi, nxf, nyi, nyf, nzi, nzf;
        double mux, muy, muz;
        double *H0, *mu;

        // Ground-state energy and excited-state energies
        double Ei, Eg;

        // nx ny nz equal to three (3)
        nmax = 10;
        // number of electrons
        Nels = 1104;
        // number of occupied orbitals
        nocc =  Nels/2;
        // number of total unoccupied orbitals
        nvirt = pow(nmax,3) - nocc;
        // nuclear attraction energy
        E_attraction = 0.;
        // Number of active occupied orbitals to excite out of
        nact = 70;
        // Number of external orbitals to excite into
        next = 80;
        // Number of single excited configurations
        ncis = nact*next;
        // ground state plus single excited configurations
        nstates = ncis+1;
        printf("  nstates is %i\n",nstates);

        printf("  nvirt is %i\n",nvirt);
        printf("  nocc is %i\n",nocc);
        fdtd_time = VEC_DOUBLE(10000000);
        fdtd_field = MAT_DOUBLE(10000000,3);

        double fdtd_dt, fdtd_total_time;
        //  Read the electric fields from Lumerical output files!  Get the total time & the time-step in atomic units
        //  and also the total number of field snapshots (itermax)
        itermax = ReadElectricField(0, "FIELDS/TFSF_TiO2_121.1nm_0_331_Z_2nm_Pt.txt", &fdtd_dt, &fdtd_total_time);
        itermax = ReadElectricField(1, "FIELDS/TFSF_TiO2_121.1nm_0_331_Y_2nm_Pt.txt", &fdtd_dt, &fdtd_total_time);
        itermax = ReadElectricField(2, "FIELDS/TFSF_TiO2_121.1nm_0_331_X_2nm_Pt.txt", &fdtd_dt, &fdtd_total_time);

        total_time = fdtd_total_time;
        dt = fdtd_dt;
        kOhb = dt/hbar;        
        //itermax = (int)(total_time/dt);
        printf("  itermax is %i, total time is %f, dt is %f\n",itermax, total_time, dt);
        fflush(stdout);

        extra_pts = 1000;
        ftlength=itermax+extra_pts;

        dm_corr_func = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(ftlength));



        E = VEC_INT(3*nmax*nmax*nmax);
        MO = MAT_INT(3*nmax*nmax*nmax,3);

        printf("  ncis is %i\n",ncis); // what is ncis?

        // Allocate more space than we need for init vec!
        init = VEC_INT(nstates+1);

        // Allocate the other arrays based on system size
        H0 = VEC_DOUBLE(nstates*nstates);
        //mu = VEC_DOUBLE(nstates*nstates);

        // New array based on system size for X, Y & Z
        mx = VEC_DOUBLE(nstates*nstates);
        my = VEC_DOUBLE(nstates*nstates);
        mz = VEC_DOUBLE(nstates*nstates);

        // complex copy of mu vectors
        cmx = VEC_CDOUBLE(nstates*nstates);
        cmy = VEC_CDOUBLE(nstates*nstates);
        cmz = VEC_CDOUBLE(nstates*nstates);


        //  For each electronic state, this vector stores the index of the frontier orbital indices nx ny and nz
        Frontier = MAT_INT(nstates,3); // 1st, 2nd and 3rd rows

        //  Initialize wavefunction as ground state
        // THis defines the initial C vector

        sup = nstates;
        init[0] = 1.;
        init[1] = 0.;
        init[2] = 0.;
        init[3] = 0.;
        init[4] = 0.;
        init[5] = 0.;


        OrderPsis3D(nmax, E, MO);


        //  < Phi_0 | H | Phi_0 >
        Eg = GroundtstateEnergyCalc (  nocc, E);
        H0[0] = Eg-Eg;
        Frontier[0][0] = MO[nocc][0];
        Frontier[0][1] = MO[nocc][1];
        idx=1;
        // Dipole elements between reference and singly-excited states
        // < Phi_0 | H | Phi_i^a > = mu_ia,  < Phi_i^a | H | Phi_0 > = mu_ia

        for (i=nocc; i>(nocc-nact); i--) {
          for (a=nocc+1; a<=(next+nocc); a++) {
            // Get the nx ny and nz orbital indices from the composite orbital indices
            // i is the initial composite index and a is the final composite index
            nxi = MO[i-1][0];
            nyi = MO[i-1][1];
            nxf = MO[a-1][0];
            nyf = MO[a-1][1];
            nzi = MO[i-1][2];
            nzf = MO[a-1][2];

            // Get the transition dipole moments along x, y and z
            mux = TransitionDipole_X(nxi, nxf, nyi, nyf, nzi, nzf);
            muy = TransitionDipole_Y(nyi, nyf, nxi, nxf, nzi, nzf);
            muz = TransitionDipole_Z(nzi, nzf, nxi, nxf, nyi, nyf);

  /*          printf("  <P_0 | mu | P(%i,%i)>\n",i,a);
            printf("  Psi(%i,%i,%i) -> Psi(%i,%i,%i)\n",nxi,nyi,nzi,nxf,nyf,nzf);
            printf("  mux:  %f\n",mux);
            printf("  muy:  %f\n",muy);
            printf("  muz:  %f\n",muz); 
            printf("  Filling elements %i and %i\n",idx,idx*nstates);
            fflush(stdout);
*/
          //  printf("mux %f muy %f muz %f \n",mux,muy,muz);

            // Store them in the array
            // Real first
            mx[idx*nstates] = mux;
            mx[idx] = mux;

            my[idx*nstates] = muy;
            my[idx] = muy;

            mz[idx*nstates] = muz;
            mz[idx] = muz;


            cmx[idx] = mux + 0.*I;

            cmy[idx*nstates] = muy + 0.*I;
            cmy[idx] = muy + 0.*I;

            cmz[idx*nstates] = muz + 0.*I;
            cmz[idx] = muz + 0.*I;



            idx++;

	  }
        }

        // Now <Phi_i^a | H | \Phi_j^b >
        idx = 1.;
        for (i=nocc; i>(nocc-nact); i--) {

          for (a=nocc+1; a<=(nocc+next); a++) {



            jidx = 1.;
            for (j=nocc; j>(nocc-nact); j--) {

              for (b=nocc+1; b<=(nocc+next); b++) {

                // So only these values matter for energy.
                if (i==j && a==b) {

                  // These are the < Phi_i^a | T | Phi_i^a > terms

                  Ei = ExcitedstateEnergyCalc(i,a, Eg, E);
                  H0[idx*nstates+jidx] = Ei-Eg;
                  Frontier[idx][0] = MO[a-1][0];
                  Frontier[idx][1] = MO[a-1][1];
                  //printf(  "  Phi(%i,%i) has energy %f, idx is %i, jidx is %i, and idx*ncis+idx is %i \n",i,a, Ei,idx,jidx,idx*nstates+jidx);

                }

                // Slaters rules... <Phi_i^a | mu | Phi_i^b > = mu_ab
                if (i==j && a!=b) {

                  nxi = MO[a-1][0];
                  nxf = MO[b-1][0];
                  nyi = MO[a-1][1];
                  nyf = MO[b-1][1];
                  nzi = MO[a-1][2];
                  nzf = MO[b-1][2];

                  mux = TransitionDipole_X(nxi, nxf, nyi, nyf, nzi, nzf);
                  muy = TransitionDipole_Y(nyi, nyf, nxi, nxf, nzi, nzf);
                  muz = TransitionDipole_Z(nzi, nzf, nxi, nxf, nyi, nyf);

                  /*printf("  <P(%i,%i) | mu | P(%i,%i)>\n",i,a,j,b);
                  printf("  Psi(%i,%i,%i) -> Psi(%i,%i,%i)\n",nxi,nyi,nzi,nxf,nyf,nzf);
                  printf("  mux:  %f\n",mux);
                  printf("  muy:  %f\n",muy);
                  printf("  muz:  %f\n",muz);
                  printf("  filling element %i\n",idx*nstates+jidx);
                  fflush(stdout);
*/
                  // Real first
                  mx[idx*nstates + jidx] += mux;
                  my[idx*nstates + jidx] += muy;
                  mz[idx*nstates + jidx] += muz;

                  // Complex copy
                  cmx[idx*nstates + jidx] += (mux + 0.*I);
                  cmy[idx*nstates + jidx] += (muy + 0.*I);
                  cmz[idx*nstates + jidx] += (muz + 0.*I);
 


              }  

                // Slaters rules... <Phi_i^a | mu | Phi_j^a > = -mu_ij (note the negative sign!)
                if (i!=j && a==b) {

                  nxi = MO[i-1][0];
                  nxf = MO[j-1][0];
                  nyi = MO[i-1][1];
                  nyf = MO[j-1][1];
                  nzi = MO[i-1][2];
                  nzf = MO[j-1][2];

                  mux = TransitionDipole_X(nxi, nxf, nyi, nyf, nzi, nzf);
                  muy = TransitionDipole_Y(nyi, nyf, nxi, nxf, nzi, nzf);
                  muz = TransitionDipole_Z(nzi, nzf, nxi, nxf, nyi, nyf);

/*                  printf("  <P(%i,%i) | mu | P(%i,%i)>\n",i,a,j,b);
                  printf("  Psi(%i,%i,%i) -> Psi(%i,%i,%i)\n",nxi,nyi,nzi,nxf,nyf,nzf);
                  printf("  mux:  %f\n",mux);
                  printf("  muy:  %f\n",muy);
                  printf("  muz:  %f\n",muz);
                  printf("  filling element %i\n",idx*nstates+jidx);
*/
                  // Real first
                  mx[idx*nstates + jidx] -= mux;
                  my[idx*nstates + jidx] -= muy;
                  mz[idx*nstates + jidx] -= muz;

                  // Complex copy
                  cmx[idx*nstates + jidx] -= (mux + 0.*I);
                  cmy[idx*nstates + jidx] -= (muy + 0.*I);
                  cmz[idx*nstates + jidx] -= (muz + 0.*I);

                }   

                jidx++;

              }
            }
            idx++;
          }
        }

        //  Now the Hamiltonian matrix is formed...uncomment if you want to print it for some perverse reason

         
/*        for (i=0; i<nstates; i++) {
          for (j=0; j<nstates; j++) {
             //printf(" %6.4e ",H0[i*nstates+j]+mu[i*nstates+j]);
             printf(" %6.4e ",H0[i*nstates+j]+mx[i*nstates+j]+my[i*nstates+j]+mz[i*nstates+j]);
          }
          printf("\n");
        }
        
*/
        FILE *cisabfp;
 
        cisabfp = fopen("CIS_Absorption.txt","w");

        double DeltaE;
        idx=1;
        for (i=nocc; i>(nocc-nact); i--) {

          for (a=nocc+1; a<=(nocc+next); a++) {

            DeltaE = ExcitedstateEnergyCalc(i,a, Eg, E) - Eg;
            fprintf(cisabfp,"  %f  %f\n",45.6/DeltaE,fabs(mx[idx]+my[idx]+mz[idx]));
            idx++;

          }
        }

        fclose(cisabfp);



        // Now solve TDSE!
        Propagate(sup, nstates, H0, mu, Frontier, E, init);
        printf("  Computing Absorption Spectrum!\n");
        //AbsorptionSpectrum(itermax, extra_pts, mur, mui);

        free(H0);
        free(E);

        free(init);

        free(mx);
        free(my);
        free(mz);
   
        free(cmx);
        free(cmy);
        free(cmz);

        for (i=0; i<nmax*nmax+nmax; i++) {
           free(MO[i]);
        }
        free(MO);

  return 0;
  exit(0);
}

// Begin Functions
double En(int n) {


 return n*n*h*h*pi*pi/(2*m*L*L);
}

double GroundtstateEnergyCalc ( int numorbs, int *E) {
  double energy_sum, fac;
  int i;

  fac = hbar*hbar*pi*pi/(2.*m*L*L);
  energy_sum = 0.;

	for (i=0; i<numorbs; i++) {
	  energy_sum += (E[i] + E_attraction) ;
        }
	double energy_total = 2*energy_sum*fac;

return energy_total;
}

double ExcitedstateEnergyCalc ( int orbital_i, int orbital_a, double energy_total, int *E) {
        double fac = hbar*hbar*pi*pi/(2.*m*L*L);
	double excited_energy = energy_total- (fac*E[orbital_i-1]+E_attraction) + (fac*E[orbital_a-1]+E_attraction);
	return excited_energy;
}

double TransitionDipole_X(int nx, int nxp, int ny, int nyp, int nz, int nzp) {
  double dipole_integral;

if (ny==nyp && nz==nzp) {

  double integ1 = pow(L,2)*(pi*(nx-nxp)*sin(pi*(nx-nxp))+cos(pi*(nx-nxp))-1)/(pi*pi*pow((nx-nxp),2));
  //printf("  integ1 is %f\n",integ1);
  double integ2 = pow(-1*L,2)*(pi*(nx+nxp)*sin(pi*(nx+nxp))+cos(pi*(nx+nxp))-1)/(pi*pi*pow((nx+nxp),2));
  dipole_integral =(2/L)*(0.5)*(integ1-integ2);
}
else {
  dipole_integral = 0;
}
return dipole_integral;
}

// If X init and X final are equal, meaning from the orbitals.. 1,2.. etc, then calc. dipole integral.
double TransitionDipole_Y(int ny, int nyp, int nx, int nxp, int nz, int nzp) {
double dipole_integral;

if (nx==nxp && nz==nzp) {
  double integ1 = pow(L,2)*(pi*(ny-nyp)*sin(pi*(ny-nyp))+cos(pi*(ny-nyp))-1)/(pi*pi*pow((ny-nyp),2));
  double integ2 = pow(-1*L,2)*(pi*(ny+nyp)*sin(pi*(ny+nyp))+cos(pi*(ny+nyp))-1)/(pi*pi*pow((ny+nyp),2));
  dipole_integral =(2/L)*0.5*(integ1-integ2);
}
else {
  dipole_integral = 0.;
}
return dipole_integral;
}

double TransitionDipole_Z(int nz, int nzp, int nx, int nxp, int ny, int nyp) {
  double dipole_integral;

if (nx==nxp && ny==nyp) {

  double integ1 = pow(L,2)*(pi*(nz-nzp)*sin(pi*(nz-nzp))+cos(pi*(nz-nzp))-1)/(pi*pi*pow((nz-nzp),2));
  //printf("  integ1 is %f\n",integ1);
  double integ2 = pow(-1*L,2)*(pi*(nz+nzp)*sin(pi*(nz+nzp))+cos(pi*(nz+nzp))-1)/(pi*pi*pow((nz+nzp),2));
  dipole_integral =(2/L)*(0.5)*(integ1-integ2);
}
else {
  dipole_integral = 0.;
}
return dipole_integral;
}



void Psin(double *XA, double *PSIA, int n, int numpts) {

 int i;

 for(i=0; i<numpts; i++){

   printf("  %f\n", XA[i]);
}

}



void OrderPsis3D(int norbs, int *E, int **MO) {

  int i, j, k, l, c, d, swap, idx, tx, ty, tz, ix, iy, iz;
  int **N, tidx[3];
  int cond, Ecur;
  N = MAT_INT(2*(norbs+1)*(norbs+1)*(norbs+1),3);

  // Going to start with i and j=0 because that's the beginning
  // of the array... nx=i+1, ny=j+1

  for (i=0; i<norbs; i++) {
    for (j=0; j<norbs; j++) {
      for (k=0; k<norbs; k++) {

        idx = i*norbs*norbs+j*norbs+k;
        // l is related to nx^2 + ny^2 + nz^2, aka, (i+1)^2 + (j+1)^2 (k+1)^2
        l = (i+1)*(i+1) + (j+1)*(j+1) + (k+1)*(k+1);
        E[idx] = l;
        // index is and energy is ...
        printf("  Index is %i and Energy[%i,%i,%i] is %i\n",idx,i+1,j+1,k+1,l);
        // element N[k][0] is nx = i+1
        N[l][0] = i+1;
        // element N[k][1] is ny = j+1
        N[l][1] = j+1;
        // element N[k][2] is nz = k+1
        N[l][2] = k+1;

      }
    }
  }

  for (c = 0 ; c < ( norbs*norbs*norbs-1 ); c++)
  {
    for (d = 0 ; d < norbs*norbs*norbs - c - 1; d++)
    {
      if (E[d] > E[d+1]) /* For decreasing order use < */
      {
        swap       = E[d];
        E[d]   = E[d+1];
        E[d+1] = swap;
      }
    }
  }

// print all energy values
for (i=0; i<(norbs*norbs*norbs); i++) {
  printf(" E[%i] is %i \n",i,E[i]);
}


  c=0;
  do {

    Ecur = E[c];

    i=0;
    do {
      i++;
      j=0;

      do {
        j++;
        k=0;

        do {

          k++;

          cond=Ecur-(i*i+j*j+k*k);

          if (cond==0) {

            MO[c][0] = i;
            MO[c][1] = j;
            MO[c][2] = k;
            c++;

          }

        }while( Ecur==E[c] && k<norbs);

      }while( Ecur==E[c] && j<norbs);

    }while (Ecur==E[c] && i<norbs);

  }while(c<norbs*norbs*norbs);

  printf(" exit successful \n");


  for (i=0; i<(norbs*norbs*norbs); i++) {
    printf("  Psi( %i , %i, %i ) %i\n",MO[i][0],MO[i][1],MO[i][2],E[i]);
  }

}

double DeltaE(int fin, int init) {

                     return ((h*h*pi*pi)/(2*m*L*L))*(fin*fin-init*init);

}



int **MAT_INT(int dim1, int dim2){
  int i,j,k;
  double sum=0.0;
  int **M;
  M = (int **)malloc(dim1*sizeof(int *));
  if (M==NULL) {
     printf("\n\nMAT_INT: Memory allocation error\n\n");
     exit(0);
  }
  for (i=0; i<dim1; i++){
      M[i] = (int *)malloc(dim2*sizeof(int));
      if (M[i]==NULL) {
         printf("\n\nMAT_INT: Memory allocation error\n\n");
         exit(0);
      }
      for (j=0; j<dim2; j++){
          M[i][j] = 0;
      }
  }
  return M;
}


double **MAT_DOUBLE(int dim1, int dim2){

  int i,j,k;
  double sum=0.0;
  double **M;
  M = (double **)malloc(dim1*sizeof(double *));
  if (M==NULL) {
     printf("\n\nMAT_INT: Memory allocation error\n\n");
     exit(0);
  }
  for (i=0; i<dim1; i++){
      M[i] = (double *)malloc(dim2*sizeof(double));
      if (M[i]==NULL) {
         printf("\n\nMAT_INT: Memory allocation error\n\n");
         exit(0);
      }
      for (j=0; j<dim2; j++){
          M[i][j] = 0;
      }
  }
  return M;


}

char *VEC_CHAR(int dim){
  char *v;
  v = (char *)malloc(dim*sizeof(char));
  if (v==NULL) {
     printf("\n\nVEC_CHAR: Memory allocation error\n\n");
     exit(0);
  }
  return v;
}


int **SQMAT_2_INT(int dim){
  int i,j;
  int **M;
  M = (int **)malloc(dim*sizeof(int *));
  if (M==NULL) {
     printf("\n\nSQMAT_2_INT: Memory allocation error\n\n");
     exit(0);
  }
  for (i=0; i<dim; i++){
      M[i] = (int *)malloc(dim*sizeof(int));
      if (M[i]==NULL) {
         printf("\n\nSQMAT_2_INT: Memory allocation error\n\n");
         exit(0);
      }
      for (j=0; j<dim; j++)
          M[i][j] = 0;
  }
  return M;
}

double *VEC_DOUBLE(int dim){
  int i;
  double *v;
  v = (double *)malloc(dim*sizeof(double));
  if (v==NULL) {
     printf("\n\nVEC_DOUBLE: Memory allocation error\n\n");
     exit(0);
  }
  for (i=0; i<dim; i++) v[i] = 0.0;
  return v;
}


int *VEC_INT(int dim){
  int *v,i;
  v = (int *)malloc(dim*sizeof(int));
  if (v==NULL) {
     printf("\n\nVEC_INT: Memory allocation error\n\n");
     exit(0);
  }
  for (i=0; i<dim; i++) v[i] = 0;
  return v;
}


void Propagate(int sup, int ncis, double *H0, double *mu, int **bas, int *E, int *init) {
  int i, j, a, iter, piter, rdmiter, itermax;
  double fac, cfac, csum, c, C, pfac, rho, rdp, tntf=1./3924.;
  double *factorB, *factorb, *q, *Q, *p, *P, *HdP, *HdQ;
  double ei;
  double complex dmx, dmy, dmz, dmmag, alpha, beta;
  double complex *D, *cn, *ct;


  int II, AA, JJ, IDX, ni, li, mi;
  double *POP, *DPOP, C2, FAC;
  int *DIDX, *DPOP_I, DNUM;
  double *DEN;
  FAC = hbar*hbar*pi*pi/(2.*m*L*L);
    



  POP = VEC_DOUBLE(nocc+next+1);
  DPOP = VEC_DOUBLE(nocc+next+1);
  DIDX = VEC_INT(nocc+next+1);
  DEN  = VEC_DOUBLE(nocc+next+1);
  DPOP_I = VEC_INT(nocc+next+1);

  II=nocc-nact;
  int IIEN = E[II];
  DIDX[0] = II;
  DEN[0] = E[II]* hbar*hbar*pi*pi/(2*m*L*L) + E_attraction;
  
  for (i=nocc-nact; i<=nocc+next; i++) {
  //for (i=1; i<(nocc+nvirt); i++) {

    // did we find a non-degenerate orbital?
    if (E[i]>IIEN) {

      II++;
      IIEN=E[i];
    }

    DIDX[i] = II;
    DEN[II] = IIEN*hbar*hbar*pi*pi/(2*m*L*L)+E_attraction;
    printf("  DIDX[%i] is %i\n",i,II);
  }

  DNUM=II;

  for (i=nocc-nact; i<DNUM; i++) {

    printf("  DEN[%i] is %i\n",i,DEN[i]);

  }



  FILE *dmfp, *orbfp, *popfp;
  FILE *vfp, *poptrace;

  dm_fftw_iter = 0;

  alpha = 1.0 + 0.*I;
  beta = 0. + 0.*I;

  //Pl = VEC_DOUBLE(2*numlam);
  //w_au = VEC_DOUBLE(2*numlam);

  //IncandescentLight(numlam, 100e-9, 1100e-9, Pl, w_au);
  // Allocate memory
  factorB = VEC_DOUBLE(6);
  factorb = VEC_DOUBLE(6);

  D  = VEC_CDOUBLE(ncis*ncis);

  cn  = VEC_CDOUBLE(ncis*ncis);
  ct  = VEC_CDOUBLE(ncis*ncis);
  Q   = VEC_DOUBLE(ncis);
  P   = VEC_DOUBLE(ncis);
  q   = VEC_DOUBLE(ncis);
  p   = VEC_DOUBLE(ncis);
  HdP = VEC_DOUBLE(ncis);
  HdQ = VEC_DOUBLE(ncis);

  factorB[0] = (642.+sqrt(471.))*tntf;
  factorB[1] = 121.*(12.-sqrt(471.))*tntf;
  factorB[2] = 1. - 2.*(factorB[0]+factorB[1]);
  factorB[3] = factorB[1];
  factorB[4] = factorB[0];

  factorb[0] = 6./11.;
  factorb[1] = .5 - factorb[0];
  factorb[2] = factorb[1];
  factorb[3] = factorb[0];
  factorb[4] = 0.;

  cfac = 1./sqrt(sup);
  csum=0;
  for (i=0; i<sup; i++) {
    a = init[i];
    cn[i] = cfac*a + 0.*I;
    csum += cfac*cfac*a*a;
  }

  for (i=0; i<sup; i++) {
    cn[i] /= sqrt(csum);
  }

  pfac = sqrt(2);
  printf("  Printing initial wavefunction vector!\n");
  for (i=0; i<ncis; i++) {
    q[i] = pfac*creal(cn[i]);
    Q[i] = pfac*creal(cn[i]);
    printf("  %f  %f\n",creal(cn[i]),cimag(cn[i]));
  }



    C2 = 1.;
    for (II=nocc-nact; II<nocc; II++) {

     DPOP_I[DIDX[II]] += C2;

    }



  c = 0.;
  C = 0.;

  fac = 1./(sqrt(2));

  itermax = (int)(total_time/dt);

  piter=0;

  rdmiter=1;
  dmfp  = fopen("DATA/DipoleMoment_121.1nm_0_331_Pt.txt","w");
  orbfp = fopen("DATA/Orbital_121.1nm_0_331_Pt.txt","w");
  vfp   = fopen("DATA/Dipole_Vector_121.1nm_0_331_Pt.txt","w");
  popfp = fopen("DATA/Populations_121.1nm_0_331_Pt.txt","w");
  poptrace = fopen("DATA/Population_Trace_121.1nm_0_331_Pt.txt","w");

  for (iter=0; iter<itermax; iter++) {

    c = 0.;
    C = 0.;
    for (a=0; a<5; a++) {

      HdotC(ncis, iter, c, H0, mu, P, HdP);

      for (i=0; i<ncis; i++) {

        Q[i] += kOhb*factorB[a]*HdP[i];

      }

      c += factorb[a];
      C += factorB[a];

      HdotC(ncis, iter, C, H0, mu, Q, HdQ);

      for (i=0; i<ncis; i++) {

        P[i] -= kOhb*factorb[a]*HdQ[i];

      }
    }
    printf("  Printing updated wavefunction\n");
    for (i=0; i<ncis; i++) {

      q[i] = Q[i];
      p[i] = P[i];

      cn[i] = fac*q[i] + I*fac*p[i];
      ct[i] = fac*q[i] - I*fac*p[i];
      //printf("  %12.10e  %12.10e\n",creal(cn[i]),cimag(cn[i]));

    }

    rho = 0.;
    dmx  = 0. + 0.*I;
    dmy  = 0. + 0.*I;
    dmz  = 0. + 0.*I;
    // Form density matrix!
    F_ZGEMM( 'N', 'T', ncis, ncis, ncis, alpha, cn, ncis, ct, ncis, beta, D, ncis);

    // Print trace of DM
    rho=0.;
    for (i=0; i<ncis; i++) {
      rho+= creal(D[i*ncis+i]);
    }
    printf("  rho is %16.12f\n",rho);

    // Get orbital population densities
    C2 = creal(D[0]);

    // Initialize
    for (II=nocc-nact; II<(nocc+next); II++) {
      POP[II] = 0.;
      DPOP[II] = 0.;
    }
    // Fill up HF reference orbitals
    for (II=nocc-nact; II<nocc; II++) {

      POP[II] += C2;
      DPOP[DIDX[II]] += C2;
    }
    IDX=1;
    for (II=nocc; II>(nocc-nact); II--) {
    for (AA=nocc+1; AA<=(nocc+next); AA++) {

      C2 = creal(D[IDX*ncis+IDX]);
      // Go through all occupied orbitals in reference again
      for (JJ=nocc-nact; JJ<nocc; JJ++) {

        POP[JJ] += C2;
        DPOP[DIDX[JJ]] += C2;
      }
      // Now correct only for orbitals that were excited in this configuration
      POP[II-1] -= C2;
      POP[AA-1] += C2;
      DPOP[DIDX[II-1]] -= C2;
      DPOP[DIDX[AA-1]] += C2;
      IDX++;

    }
  }

  fprintf(popfp,"\n\n#%i\n",iter);
  for (II=nocc-nact; II<DNUM; II++) {
     fprintf(popfp,"  %f  %f\n",DEN[II],DPOP[II]-DPOP_I[II]);
  }
  //fprintf(poptrace,"\n  %12.10f ",iter*dt*time_au/1e-15);
  fprintf(poptrace,"\n %i ",iter);
  for (II=nocc-nact; II<DNUM; II++) {
    fprintf(poptrace," %12.10e ",DPOP[II]-DPOP_I[II]);
  }


  
    dmx = F_ZDOTC(ncis*ncis, D, 1, cmx, 1);
    dmy = F_ZDOTC(ncis*ncis, D, 1, cmy, 1);
    dmz = F_ZDOTC(ncis*ncis, D, 1, cmz, 1);
 
    fprintf(vfp,"\n\n#%i\n",iter);
    fprintf(vfp,"%12.10f  %12.10f  %12.10f  %12.10f  %12.10f  %12.10f \n",0.,0.,0.,dmx,dmy,dmz);
    fprintf(dmfp,"%12.10e  %12.10e %12.10e\n",iter*dt,creal(dmx+dmy+dmz),cimag(dmx+dmy+dmz));
    //mutr[iter] = creal(dmx + dmy + dmz);
    //muti[iter] = cimag(dmx + dmy + dmz);
    dm_corr_func[dm_fftw_iter] = dmx + dmy + dmz;
    dm_fftw_iter++;


    //if (iter%5==0) {

    //PrintFrontier(piter, ncis, orbfp, cn, bas, E);
    //piter++;
   // }
  }
  fclose(orbfp);
  fclose(dmfp);

  free(factorB);

  free(D);
  free(cn);
  free(ct);
  free(Q);
  free(P);
  free(q);
  free(p);
  free(HdP);
  free(HdQ);
}


void HdotC(int end, int iter, double m, double *H0, double *mu, double *C, double *HdC) {
  int i, j;
  double sum, field_x, field_y, field_z;
  double alpha, beta;
  double current_time, stencil_time;
  double current_fieldx, stencil_fieldx, current_fieldy, stencil_fieldy, current_fieldz, stencil_fieldz;
  double fieldx_slope, fieldy_slope, fieldz_slope;
  int n, lda, incx, incy;

  n = end; /* Size of Row ( the number of columns ) */
  lda = end; /* Leading dimension of 5 * 4 matrix is 5 */
  incx = 1;
  incy = 1;
  alpha = 1;
  beta = 0;

  //  Hdc = alpha*H*C + beta*HdC: 
  F_DGEMV( 'T', n, n, alpha, H0, lda, C, 1, beta, HdC, 1);

  // Adapted for when field is read from FDTD output and time is fixed by 
  // FDTD output... need to interpolate fields for the stencil in 
  // time demanded by Sanz-Serna algorithm

  // remember that we offset the fields and time arrays by 1 so that we can always take centered finite differences
  if (1) {
    // This is the time at the current iteration
    current_time = fdtd_time[iter+1];
    // This is the time at the current stencil point
    stencil_time = fdtd_time[iter+1] + m*dt;

    // These are the field values at the time of the current iteration
    current_fieldx = fdtd_field[iter+1][0];
    current_fieldy = fdtd_field[iter+1][1];
    current_fieldz = fdtd_field[iter+1][2];

    // These are the slopes of the fields at the current iteration
    fieldx_slope = (fdtd_field[iter+2][0] - fdtd_field[iter][0])/(fdtd_time[iter+2]-fdtd_time[iter]);
    fieldy_slope = (fdtd_field[iter+2][1] - fdtd_field[iter][1])/(fdtd_time[iter+2]-fdtd_time[iter]);
    fieldz_slope = (fdtd_field[iter+2][2] - fdtd_field[iter][2])/(fdtd_time[iter+2]-fdtd_time[iter]);


    // These are the fields values at the time at the current stencil point
    stencil_fieldx = current_fieldx + fieldx_slope*m*dt; 
    stencil_fieldy = current_fieldy + fieldy_slope*m*dt;
    stencil_fieldz = current_fieldz + fieldz_slope*m*dt;
   
    // Scale the field values 
    stencil_fieldx*=E0;
    stencil_fieldy*=E0;
    stencil_fieldz*=E0;

    //  HdC -= lpulse_x*mux*C
    F_DGEMV( 'T', n, n, -1*stencil_fieldx, mx, lda, C, 1, 1, HdC, 1);
    
    // HdC -= lpulse_y*muy*C
    F_DGEMV( 'T', n, n, -1*stencil_fieldy, my, lda, C, 1, 1, HdC, 1);

    // HdC -= lpulse_z*muz*C
    F_DGEMV( 'T', n, n, -1*stencil_fieldz, mz, lda, C, 1, 1, HdC, 1);

  }

}



void PrintFrontier(int iter, int ncis, FILE *fp, double complex *c, int **bas, int *E) {
  int i, j, k, a, numpts;
  int nx, ny;

  numpts = 100;
  double dx = L/numpts;
  double dy = L/numpts;
  double fac;

  fac = hbar*hbar*pi*pi/(2*m*L*L);

  double wfn_real, wfn_imag, wfn_mag, en, xval, yval;

  fprintf(fp,"\n\n#%i  \n",iter+1);
  for (i=0; i<numpts; i++) {

    wfn_real = 0.;
    wfn_imag = 0.;
    xval = i*dx;

    for (j=0; j<numpts; j++) {

      yval = j*dy;
      wfn_real = 0;
      wfn_imag = 0;
      for (k=0; k<ncis; k++) {

        en =  E[k]*fac;
        nx = bas[k][0];
        ny = bas[k][1];
        wfn_real += creal(c[k]*(2./L)*sin( (nx*pi/L)*xval)*sin( (ny*pi/L)*yval)*cexp(I*en*iter*dt));
        wfn_imag += cimag(c[k]*(2./L)*sin( (nx*pi/L)*xval)*sin( (ny*pi/L)*yval)*cexp(I*en*iter*dt));

      }

      fprintf(fp,"  %12.10e  %12.10e  %12.10e  %12.10e  %12.10e\n",xval,yval,wfn_real, wfn_imag,(wfn_real*wfn_real+wfn_imag*wfn_imag));

    }
    fprintf(fp,"\n");
  }
}

double complex *VEC_CDOUBLE(int dim) {
  int i;
  double complex *v;
  v = (double complex *)malloc(dim*sizeof(double complex));
  if (v==NULL) {
     printf("\n\nVEC_CDOUBLE:  Memory allocation error\n\n");
     exit(0);
  }
  for (i=0; i<dim; i++) v[i] = 0.0 + I*0.;
  return v;
}

double EFieldTime(int numlam, double tim, double *Pl, double *w_au) {
  int i;
  double Eval;

  Eval = 0.0;

  for (i=0; i<=numlam; i++) {

    Eval += Pl[i]*sin(w_au[i]*tim);

  }
  printf("  Field is %f\n",Eval);
  return Eval;
}



void IncandescentLight(int numlam, double lam1, double lam2, double *Pl, double *w_au) {
     int i, j;
     double c, kb, T, lam, nu, omega, nu_au, omega_au;
     double dlam, pmax;
     FILE *fp;

     fp = fopen("Incandescent.txt","w");
     c=299792458;
     kb=1.380642e-23;
     T=10000;

     pmax=-1000;
     omega_au = 4.134147e+16;
     nu_au = omega_au/(2*pi);

    dlam = (lam2-lam1)/numlam;

     for (i=0; i<=numlam; i++) {
       lam=lam1+dlam*i;
       nu=c/lam;
       omega=2*pi*nu;
       w_au[i] = omega/omega_au;
       Pl[i] = (2*h*pow(nu,3)/(c*c))*(1/(exp(h*nu/(kb*T))-1));
       if (Pl[i]>pmax) pmax=Pl[i];
     }

     for (i=0; i<=numlam; i++) {
       Pl[i] = Pl[i]/pmax;
       fprintf(fp,"  %12.10e  %12.10e\n",w_au[i],Pl[i]);
     }

     fclose(fp);
}

void buildBasis(int no, int nv, int **abaa,int **ibaa) {
   int n, i, a;

   a = no-1;
   i = no-1;
   n = 0;
   abaa[n][0] = i;
   abaa[n][1] = a;
   ibaa[i][a] = n;
   ibaa[a][i] = n;

   n = 1;
   for (a=no; a<(no+nv); a++) {
     for (i=no-1; i>=0; i--) {
          abaa[n][0] = i;
          abaa[n][1] = a;
          ibaa[a][i] = n;
          ibaa[i][a] = n;
          n++;
      }
   }
}

void AbsorptionSpectrum(int iter, int extra_pts, double *mutr, double *muti) {

  fftw_plan p;
  FILE *fp;

  long int i;
  double val, valr, vali, twopi;
  double w, w_si, nu_si, lambda;
  double *spec_real, *spec_imag;

  double c=299792458;
  double omega_au = 4.134147e+16;
  double nu_au = omega_au/(2*pi);
  twopi = 2.*pi;

  fp = fopen("AbsorptionSpectrum.txt","w");
  for (i=0; i<extra_pts; i++) {

    dm_corr_func[iter+i] = 0. + I*0.;

  }
  p = fftw_plan_dft_1d((extra_pts+iter),dm_corr_func,dm_corr_func,FFTW_BACKWARD,FFTW_ESTIMATE);
  fftw_execute(p);


  for (i=1; i<(iter+extra_pts); i++) {
    w = twopi*i/((iter+extra_pts)*dt);
    w_si = w*omega_au;
    nu_si = w_si/twopi;
    lambda=c/nu_si;
    valr = creal(dm_corr_func[i]/(iter));
    val = 2./3. * w * sqrt(valr*valr);
    fprintf(fp,"  %12.10lf  %12.12lf \n", lambda*1e9, val );

  }

  fclose(fp);
}

int ReadElectricField(int component, char *filename, double *fdtd_dt, double *fdtd_total_time) {
  FILE *fp;
  int i;
  double tval,tval_SI, tval_au, eval_SI, eval_au;
  fp = fopen(filename, "r");

  // We are going to offset these arrays by 1 so that we can take centered-finite differences 
  // of the field starting at the first timepoint... we need to do that 
  // because we need the field values at arbitrary times for the propagation step, so 
  // we will need the derivative of the field at a given time to estimate the field at some nearby time
  i=1;
  do {

    // read time (in fs)
    fscanf(fp,"%lf",&tval);
    // convert to atomic units of time 
    tval_SI = 1e-15*tval;
    tval_au = tval_SI/time_au;
    // store value in atomic units
    fdtd_time[i] = tval_au;

    // read field value in SI units
    fscanf(fp,"%lf",&eval_SI);

    // convert to atomic units of electric field
    //eval_au = eval_SI/1.0;
    fdtd_field[i][component] = eval_SI;
    printf("  component %i is %e\n",component, eval_SI); 
    i++;

  }while(!feof(fp));

  i--;
  *fdtd_dt = fdtd_time[2]-fdtd_time[1];
  *fdtd_total_time = fdtd_time[i];
 
  return i;

  fclose(fp);
}




/*
void AbsorptionSpectrum(int iter, int extra_pts, double *mutr, double *muti) {
  FILE *fp;
  long int i;
  double val, valr, vali, twopi;
  double w, w_si, nu_si, lambda;
  double *spec_real, *spec_imag;

  double c=299792458;
  double omega_au = 4.134147e+16;
  double nu_au = omega_au/(2*pi);



  fp = fopen("Absorption.txt", "w");
  twopi = 2.*pi;

  spec_real = VEC_DOUBLE(iter+extra_pts);
  spec_imag = VEC_DOUBLE(iter+extra_pts);

  for (i=0; i<extra_pts; i++) {
    mutr[iter+i] = 0.;
    muti[iter+i] = 0.;
  }



  compute_dft(mutr, muti, spec_real, spec_imag, iter+extra_pts);
  for (i=1; i<(iter+extra_pts); i++) {
    w = twopi*i/((iter+extra_pts)*dt);
    w_si = w*omega_au;
    nu_si = w_si/twopi;
    lambda=c/nu_si;
    fprintf(fp,"  %20.12lf  %20.12lf \n",lambda*1e9,sqrt(spec_real[i]*spec_real[i]+spec_imag[i]*spec_imag[i]));

  }

  fclose(fp);
  free(spec_real);
  free(spec_imag);
}

*/


