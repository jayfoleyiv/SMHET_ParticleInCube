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
double L = 2e-9/Lau;;

  // Fermi energy of gold is 5.53 eV
  // Fermi energy of silver is 5.49 eV
  // Fermi energy of iron is 11.1 eV
  // Fermi energy of platinum is 9.75 eV


int main() {

double FermiEnergy=5.53*eVtoAU;
double term1, term2;
double fac = pi*pi/(2*L*L);
int electrons;
int nf;


term1 = (L*L*L)/(3*pi*pi);
term2 = pow(2*FermiEnergy, 3./2.);

electrons=(int)term1*term2;

nf=(int)sqrt(FermiEnergy*(2*L*L)/(pi*pi*1));

double eF= 0.5*pow( (3*pi*pi*(electrons)/(L*L*L)), 2./3.);

printf ("%f  from hbar^2/2m * (3pi^2 * Ne/Vnc)^2/3 \n", 27.211*eF);

printf("%i   from, Vnc/3pi^2 * (2mEf/hbar^2)^3/2  \n", electrons);

printf("%i   from sqrt(Ef*2*m*L^2/pi^2 hbar^2)  \n", nf);
//nf=5;
printf("%f   from hbar^2 pi^2/2mL^2( nx^2 + ny^2 + nz^2)\n",(27.211*fac*(nf*nf)));

return 0;
}


