#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<malloc.h>

double pi = 3.14159265358979, m = 1, hbar = 1, ec=1;
double h=6.626070040e-34;
double Lau = 5.2917721092e-11;
// 1 Atomic unit of time in SI (s)
double time_au = 2.418884326505e-17;
// 1 atomic unit electric field in SI (V/m) 
double field_au = 5.14220652e11;
double eVtoAU = 0.03674930882;
void OrderPsis3D(int norbs, int *E, int **MO);
double *VEC_DOUBLE(int dim);
int **MAT_INT(int dim1, int dim2);
double L = 2e-9/Lau;;

int main() {

  int i, *E, **MO, nmax, norbs;

  // Fermi energy of gold is 5.53 eV
  // Fermi energy of silver is 5.49 eV
  // Fermi energy of iron is 11.1 eV
// Fermi energy of platinum is 9.75 eV
  double Ef = 9.75*eVtoAU;
  // 20 nm in atomic units
  // Energy prefactor in atomic units
  double fac = pi*pi/(2*L*L);
  nmax = 8;
  norbs = nmax*nmax*nmax;

  E = (int *)malloc(norbs*sizeof(int));
  MO = MAT_INT(norbs,3);

  OrderPsis3D(nmax, E, MO);

  double Ecurr=0.;
  int die=1;
  i=0;
  do {

    Ecurr = fac*E[i];
    
    if ( Ecurr>Ef) die=0;

    i++;
  }while( die==1 && i<norbs);

  printf("  Ecurr is %f\n",Ecurr);
  printf("  Orbital Number is %i\n",i);

return 0;
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
        printf("  Index is %i and Energy[%i,%i,%i] is %f\n",idx,i+1,j+1,k+1,pi*pi/(2*L*L)*l);
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
    printf("  Psi( %i , %i, %i ) %i %f\n",MO[i][0],MO[i][1],MO[i][2],i,pi*pi/(2*L*L)*E[i]);
  }

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

