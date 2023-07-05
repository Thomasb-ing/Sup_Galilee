#include <mpi.h>
# include <cstdlib>
# include <cmath>
# include <algorithm>
# include <iostream>
# include <cassert>
# include <stdint.h>
# include "Chronometer.hpp"

////
//// Init
////
void init( int* Ndim_tab, int nfic, double* temp0, double* kdiff, double* x , double* y, double* dx, int rank, int nbp )
{
  const double lx = 5./nbp;
  const double ly = 1;

  const double mu0 =0.0000004;   // coeff diffusion beton
  const double mu1 =0.0001;      // coeff diffusion cuivre

  dx[0] = lx/(Ndim_tab[0]-2*nfic);
  dx[1] = ly/(Ndim_tab[1]-2*nfic);

  const double x0 = 0 + lx*rank;
  const double y0 = 0;

  const double xinit = 2.;
  const double yinit = 0.95;

  for (int64_t i = 0; i < Ndim_tab[0] ; ++i ){
     x[i] = (i-1)*dx[0] + x0;
  for (int64_t j = 0; j < Ndim_tab[1] ; ++j ){
     y[j] = (j-1)*dx[1] + y0;

     int l = j*Ndim_tab[0]+ i;

     temp0[l]= 280;

     if( (x[i] <= xinit or x[i] >=xinit+1) and y[j] > yinit )  {kdiff[l]=mu0;}
     else                                                      {kdiff[l]=mu1;}
  }
  }
}

////
//// mise a jour
////
void mise_a_jour( int* ndim_rank, double* temp0, double* temp1, double* bilan, const double dt, int step )
{
 double cte_rk =0.5;
 if(step==0) { cte_rk = 1;}

 for (int64_t j = 1; j < ndim_rank[1]-1 ; ++j ){ 
  for (int64_t i = 1; i < ndim_rank[0]-1 ; ++i ){ 
    
     int    l = j*ndim_rank[0]+ i;
 
     temp1[l]    = temp0[l] - dt*cte_rk*bilan[l]; 
   }
  }
}

void diffusion( int* ndim,   double* temp, double* bilan, double* dx, double* kdiff, int step  )
{
  double cte_rk =1;
  if(step==0) { cte_rk = 0;}

    for (int64_t j = 1; j < ndim[1]-1 ; ++j ) {
      for (int64_t i = 1; i < ndim[0]-1 ; ++i ) { 

         int    l = j*ndim[0]+ i;// (i  , j  )
         int    l1= l+1;              // (i+1, j  )
         int    l2= l-1;              // (i-1, j  )
         int    l3= l+ndim[0];   // (i  , j+1)
         int    l4= l-ndim[0];   // (i  , j-1)

         double mup = 0.5*( kdiff[l] + kdiff[l1]);
         double mum = 0.5*( kdiff[l] + kdiff[l2]);

         bilan[l] = bilan[l]*cte_rk- ( mup*(temp[l1]-temp[l]) - mum*(temp[l]-temp[l2]) )/(dx[0]*dx[0]);

         mup = 0.5*( kdiff[l] + kdiff[l3]);
         mum = 0.5*( kdiff[l] + kdiff[l4]);

         bilan[l] += -( mup*(temp[l3]-temp[l]) - mum*(temp[l]-temp[l4]) )/(dx[1]*dx[1]);
      }
    }
}

int main( int nargc, char* argv[])
{

  int rank, nbp;
  MPI_Init(&nargc, &argv);
  MPI_Comm comm;
  MPI_Comm_dup(MPI_COMM_WORLD, &comm);
  MPI_Comm_size(comm, &nbp);
  MPI_Comm_rank(comm, &rank);
  
  
  char fileName[255];
  FILE* out;


  Chronometer chrono;
  
  sprintf(fileName, "Sortie%05d.txt", rank);
  out = fopen(fileName, "w");
 
  
  // //
  // //Determination la taille des grilles dans les direction X ( Ndim_tab[1]) et Y ( Ndim_tab[1]) 
  // //
  // int dim[2]; dim[0] =500; dim[1]=500;
  // int nfic     =  1;
  // int Ndim_tab[2];
  // Ndim_tab[0] = dim[0]+2*nfic;  // on suppose que tous les bloc on la meme taille
  // Ndim_tab[1] = dim[1]+2*nfic;  // on suppose que tous les bloc on la meme taille

  int dim[2]; dim[0] =500/nbp; dim[1]=500;  // on divise en 2 sur les x pour faire sur 2 process
  int nfic     =  1;
  int Ndim_tab[2];
  Ndim_tab[0] = dim[0]+2*nfic;  // on suppose que tous les bloc on la meme taille
  Ndim_tab[1] = dim[1]+2*nfic; 

  double *x,*y, *temp1, *temp0,  *bilan, *kdiff, *buffer, *buffer_s;
  x       = new double[Ndim_tab[0]];
  y       = new double[Ndim_tab[1]];
  bilan   = new double[Ndim_tab[0]*Ndim_tab[1]]; 
  temp1   = new double[Ndim_tab[0]*Ndim_tab[1]]; 
  temp0   = new double[Ndim_tab[0]*Ndim_tab[1]]; 
  kdiff   = new double[Ndim_tab[0]*Ndim_tab[1]]; 
  buffer = new double[nfic*2*(Ndim_tab[0])];
  buffer_s = new double[nfic*2*(Ndim_tab[0])];

  for (int64_t j = 0 ; 1 <nfic*2*Ndim_tab[0] ;j++) {buffer[j]=-40000 ; buffer_s[j]=40000;}

  chrono.click();
  double dx[2];

  init(  Ndim_tab, nfic,  temp0, kdiff, x, y, dx, rank , nbp); // changer pour intitier selon le processus
  fprintf(out, "dim blocX =  %d, dim blocY =  %d, dx= %f, dy= %f \n",Ndim_tab[0], Ndim_tab[1],  dx[0], dx[1] );

  for (int64_t j = 0; j < Ndim_tab[1] ; ++j ){
    for (int64_t i = 0; i < Ndim_tab[0] ; ++i ){ 

    int    l = j*Ndim_tab[0]+ i;
    fprintf(out, " Init: %f %f %f \n", x[i],y[j], temp0[l]); 
   }
    fprintf(out, " Init: \n"); 
  }



  const double dt       =0.01;  // pas de temps
  const double Tambiant =300.;  // Temperature ambiante en Kelvin
 
  int Nitmax      = 5000;
  int Stepmax     = 2;

  //Boucle en temps
  for (int64_t nit = 0; nit < Nitmax ; ++nit )
  { 
    //Boucle Runge-Kutta
    double *Tin;
    double *Tout;
    double *Tdiff;
    for (int64_t step = 0; step < Stepmax ; ++step )
    { 
       if(Stepmax==2){
         if(step==0) { Tin = temp0; Tout= temp1; Tdiff=temp0; }
         else        { Tin = temp0; Tout= temp0; Tdiff=temp1;}
       }
       else{ Tin = temp0; Tout= temp0; Tdiff=temp0; }

       diffusion(Ndim_tab, Tdiff, bilan,  dx, kdiff, step);
       //mise a jour point courant
       mise_a_jour(Ndim_tab, Tin, Tout, bilan,  dt, step);


       //conditions aux limites
       for (int64_t ific = 0; ific < nfic ; ++ific )
       {  
          //adiabatique en Jmin
          for (int64_t i = 0; i < Ndim_tab[0]  ; ++i )
          {  
           //Jmin
           int l0   = ific +i; 
           int l1   = Ndim_tab[0] +ific +i;

           Tout[l0] = Tout[l1];

           //paroi isotherme en Jmax 
           l0   = Ndim_tab[0]*(Ndim_tab[1]-nfic) +ific +i;
           Tout[l0] = Tambiant;
          }
       }
       for (int64_t ific = 0; ific < nfic ; ++ific )
       { 
          //paroi adiabatique en Imin
          for (int64_t j = 0; j < Ndim_tab[1]  ; ++j )
          {  
           //Imin
           int l0   = ific +j*Ndim_tab[0]; 
           int l1   = l0 +1;
           Tout[l0] = Tout[l1];

           //Imax
            l0   = ific + (j+1)*Ndim_tab[0] - 1;
            l1   = l0 - 1;
           Tout[l0] = Tout[l1];
          }
       }

     }  // Nstepmax
    }  // Nitmax

 for (int64_t i = 1; i < Ndim_tab[0]-1 ; ++i ){ 
  for (int64_t j = 1; j < Ndim_tab[1]-1 ; ++j ){ 

    int    l = j*Ndim_tab[0]+ i;
    fprintf(out, " Final %f %f %.8f \n", 0.5*(x[i]+x[i+1]),0.5*(y[j]+y[j+1]), temp0[l]); 
    //fprintf(out, " Final %f %f %.8f \n", x[i],y[j], kdiff[l]); 
   }
    fprintf(out, " Final \n"); 
  }

  double tensDt = chrono.click();
  fprintf(out, "Temps initialisation   %f  \n", tensDt);

  fclose(out);

  delete [] temp0;  delete [] temp1; delete [] bilan; delete [] x;

  // fin MPI
  MPI_Finalize();
  return EXIT_SUCCESS;
}

