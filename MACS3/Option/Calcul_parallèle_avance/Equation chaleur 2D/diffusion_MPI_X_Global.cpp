# include <mpi.h>
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
void init( int rank, int max_rankX, int max_rankY, int* ndim_rank, int* ndim, int nfic, double* temp0, double* kdiff, double* x , double* y, double* dxy,double* xy0  )
{
  const double mu0 =0.0000004;   // coeff diffusion beton
  const double mu1 =0.0001;      // coeff diffusion cuivre

  const double x0 = xy0[0];
  const double y0 = xy0[1];

  const double xinit = 2.;
  const double yinit = 0.95;

  for (int64_t i = 0; i < ndim_rank[0] ; ++i ){
     x[i] = (i-1)*dxy[0] + x0;
  for (int64_t j = 0; j < ndim_rank[1] ; ++j ){
     y[j] = (j-1)*dxy[1] + y0;

     int l = j*ndim_rank[0]+ i;

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

void diffusion( int* ndim_rank,   double* temp, double* bilan, double* dx, double* kdiff, int step  )
{
  double cte_rk =1;
  if(step==0) { cte_rk = 0;}

    for (int64_t j = 1; j < ndim_rank[1]-1 ; ++j ) {
      for (int64_t i = 1; i < ndim_rank[0]-1 ; ++i ) { 

         int    l = j*ndim_rank[0]+ i;// (i  , j  )
         int    l1= l+1;              // (i+1, j  )
         int    l2= l-1;              // (i-1, j  )
         int    l3= l+ndim_rank[0];   // (i  , j+1)
         int    l4= l-ndim_rank[0];   // (i  , j-1)

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
  MPI_Comm globalComm;
  MPI_Request request[2];
  MPI_Status  status;

  int rank, Max_rank;
  char fileName[255];
  FILE* out;

  MPI_Init(&nargc, &argv);
  MPI_Comm_dup(MPI_COMM_WORLD, &globalComm);
  MPI_Comm_rank(globalComm, &rank);
  MPI_Comm_size(globalComm, &Max_rank);

  int dim[2]; dim[0] =500; dim[1]=500;
  //int dim[2]; dim[0] =1000; dim[1]=1000;
  int nfic     =  1;

  //if ( nargc > 1 ) dim = atoi(argv[1]);

  Chronometer chrono;
  
  sprintf(fileName, "Sortie%05d.txt", rank);
  out = fopen(fileName, "w");
 
  //
  //
  //Determination partage de la grille entre les processus MPI
  //
  //
    //Si nombre processus pair, on decoupe la grille en 2 dans la direction Y
    int Max_rankY=1;

    //on deduit le nombre de processus répartis dans la directionX
    int Max_rankX=Max_rank/Max_rankY;


    // on determine la position du rank dans cette grille (I,J) de processus    
    //
    // Exemple
    //
    //      X
    // 4   5   6   7   
    // 0   1   2   3   Y
    //
    //  le rank 6 se situe en 2eme ligne et 3eme colonne
    //
    int Pos_rankY = rank/Max_rankX;
    int Pos_rankX = (rank-Pos_rankY*Max_rankX);

  fprintf(out, "Hello world from %d inside %d processus!\n", rank, Max_rank);
  fprintf(out, "rank= %d POsX= %d PosY= %d \n", rank, Pos_rankX, Pos_rankY );
  //
  //
  //Determination la taille des grilles dans les direction X ( Ndim_rank[0]) et Y ( Ndim_rank[1]) pour chacun des ranks
  //
  //
  int Ndim_rank[2];
  int NcellX = dim[0]/Max_rankX;
  int reste  = dim[0]-Max_rankX*NcellX;

  double xy0[2];
  xy0[1]=0;
  double lx =5;
  double ly =1;
  double dxy[2];
  dxy[0] = lx/ dim[0];
  dxy[1] = ly/ dim[1];

  if (rank < reste){ Ndim_rank[0] = NcellX+1+2*nfic; xy0[0]=( NcellX+1)*dxy[0]*rank; }
  else{ Ndim_rank[0] = NcellX+2*nfic; xy0[0]= (reste+(NcellX)*rank)*dxy[0];}

  Ndim_rank[1] = dim[1]/Max_rankY+2*nfic; 


  double *x,*y, *temp1, *temp0,  *bilan, *buffer, *buffer_s, *kdiff;
  x       = new double[Ndim_rank[0]];
  y       = new double[Ndim_rank[1]];
  bilan   = new double[Ndim_rank[0]*Ndim_rank[1]]; 
  temp1   = new double[Ndim_rank[0]*Ndim_rank[1]]; 
  temp0   = new double[Ndim_rank[0]*Ndim_rank[1]]; 
  kdiff   = new double[Ndim_rank[0]*Ndim_rank[1]]; 
  buffer  = new double[nfic*2*Ndim_rank[1]]; 
  buffer_s= new double[nfic*2*Ndim_rank[1]]; 
  
  double moy_loc[2];
  double moy[2];
  
  for (int64_t j = 0; j < nfic*2*(Ndim_rank[0]+Ndim_rank[1])  ; ++j ){ buffer[j]=-40000; buffer_s[j]=40000;}
  
  chrono.click();

  init(  rank, Max_rankX, Max_rankY, Ndim_rank, dim, nfic,  temp0, kdiff, x, y, dxy, xy0);
  fprintf(out, "dim blocX =  %d, dim blocY =  %d, dx= %f, dy= %f, x0= %f \n",Ndim_rank[0], Ndim_rank[1],  dxy[0], dxy[1], xy0[0] );

/*
  for (int64_t j = 0; j < Ndim_rank[1] ; ++j ){
    for (int64_t i = 0; i < Ndim_rank[0] ; ++i ){ 

    int    l = j*Ndim_rank[0]+ i;
    fprintf(out, " Init: %f %f %f \n", x[i],y[j], temp0[l]); 
   }
    fprintf(out, " Init: \n"); 
  }

*/

  //const double dt       =0.018;  // pas de temps
  const double dt       =0.01;  // pas de temps
  const double Tambiant =300.;  // Temperature ambiante en Kelvin
 
  //int Nitmax      = 150000;
  int Nitmax      =   1000;

  int Stepmax     = 2; //RK2
  //int Stepmax     = 1; //RK1

  //Boucle en temps
  for (int64_t nit = 0; nit < Nitmax ; ++nit )
  { 
    //Boucle Runge-Kutta
    double *Tin;
    double *Tout;
    double *Udiff;
    for (int64_t step = 0; step < Stepmax ; ++step )
    { 
       if(Stepmax==2){
         if(step==0) { Tin = temp0; Tout= temp1; Udiff=temp0; }
         else        { Tin = temp0; Tout= temp0; Udiff=temp1;}
       }
       else{ Tin = temp0; Tout= temp0; Udiff=temp0; }

       diffusion(Ndim_rank, Udiff, bilan,  dxy, kdiff, step);
       //mise a jour point courant
       mise_a_jour(Ndim_rank, Tin, Tout, bilan,  dt, step);

      //Application Condition limite
      //adiabatique en Jmin
      for (int64_t i = 0; i < Ndim_rank[0]  ; ++i )
      {  
        //Jmin
        int l0   = i; 
        int l1   = Ndim_rank[0] +i;

        Tout[l0] = Tout[l1];
      }
      //periodicité en Jmax 
      for (int64_t i = 0; i < Ndim_rank[0]  ; ++i )
      {  
        int l0   = Ndim_rank[0]*(Ndim_rank[1]-nfic)  +i;
        Tout[l0] = Tambiant;
      }
      if (Pos_rankX==0) 
      { 
          //paroi adiabatique en Imin
          for (int64_t j = 0; j < Ndim_rank[1]  ; ++j )
          {  
           //Imin
           int l0   = j*Ndim_rank[0]; 
           int l1   = l0 +1;

           Tout[l0] = Tout[l1];
          }
      }
      if (Pos_rankX==Max_rankX-1) 
      { 
          //paroi adiabatique en Imax 
          for (int64_t j = 0; j < Ndim_rank[1]  ; ++j )
          {  
           int l0   = (j+1)*Ndim_rank[0] - 1;
           int l1   = l0 - 1;
           Tout[l0] = Tout[l1];
          }
      }
      if (Max_rank!=1)

      { int shift_Imin=0;
        int shift_Imax=Ndim_rank[1];

        //nbre de raccord
        int nracX=0;
        if (Max_rankX!=1)
          {  
             if (Pos_rankX==0 )
             {
               nracX=1;
               //Remplissage buffer envoi imax
               for (int64_t j = 0; j < Ndim_rank[1] ; ++j )
               {  
                int l1                    = (j+1)*Ndim_rank[0] - 2;
                buffer_s[ j + shift_Imax ]= Tout[l1];
               }
             }
             else if (Pos_rankX==Max_rankX-1)
             {
               nracX=1;
               //Remplissage buffer envoi imin
               for (int64_t j = 0; j < Ndim_rank[1] ; ++j )
               {  
                int l1                    = j*Ndim_rank[0]+1;
                buffer_s[ j + shift_Imin ]= Tout[l1]; 
               }
             }
             else
             { nracX=2;
               //Remplissage buffer envoi
               for (int64_t j = 0; j < Ndim_rank[1] ; ++j )
               {  
                int l1   = j*Ndim_rank[0]+1;
                buffer_s[ j +shift_Imin ]= Tout[l1]; //imin

                l1   = (j+1)*Ndim_rank[0] - 2;
                buffer_s[ j+ shift_Imax  ]= Tout[l1];//imax
               }
             }
          }

        int rank_donor, rank_dest, shift;
        //
        //Reception non bloquante direction X
        //
        for (int64_t rac = 0; rac < nracX ; ++rac )
        {
           if(rac        == 0 )
             { rank_donor = rank -1;
               shift =shift_Imin;
               if(Pos_rankX == 0) {rank_donor = rank + 1; shift = shift_Imax;}
             }
           else
             { rank_donor = rank +1;
               if(Pos_rankX == Max_rankX-1) printf("Cas IMPOSSIBLE X \n");
               shift =shift_Imax;
             }
           int etiquette= rank + 10000*rank_donor;
           int size     = Ndim_rank[1];

           MPI_Irecv( buffer+shift, size, MPI_DOUBLE,  rank_donor,
                      etiquette,
                      globalComm,
                   &request[rac]);
        }
        //
        //Envoi  bloquant direction X
        //
        for (int64_t rac = 0; rac < nracX ; ++rac )
        {
           if(rac        == 1 )
             { rank_dest  = rank +1;
               shift = shift_Imax;
             }
           else
             { if(Pos_rankX == 0) { rank_dest  = rank +1; shift = shift_Imax; }
               else               { rank_dest  = rank -1; shift = shift_Imin; }
             }

           int etiquette= 10000*rank + rank_dest;
           int size     = Ndim_rank[1];
           MPI_Send(buffer_s+shift, size, MPI_DOUBLE,  rank_dest, etiquette, globalComm);

        }

        //barrier MPI
        for (int64_t rac = 0; rac < nracX ; ++rac )
        {
         int flag=0;
         MPI_Test(&request[rac], &flag, &status);
         while (!flag) { MPI_Test(&request[rac], &flag, &status); } 
        }

        //mise a jour a partir des buffer
        //
        //
        if (Max_rankX!=1) 
        { 
          if(Pos_rankX == 0)
          {
            for (int64_t j = 0; j < Ndim_rank[1]  ; ++j )
            {  
             int l1   = (j+1)*Ndim_rank[0] - 1;   //imax
             Tout[l1] = buffer[j + shift_Imax ];
            }
          }
          else if(Pos_rankX == Max_rankX-1)
          {
            for (int64_t j = 0; j < Ndim_rank[1]  ; ++j )
            {  
             int l0   =     j*Ndim_rank[0];       //imin
             Tout[l0] = buffer[j  + shift_Imin    ];
            }
          }
          else
          {
          //periodicité en Imax et Imin
            for (int64_t j = 0; j < Ndim_rank[1]  ; ++j )
            {  
             //Imin
             int l0   =     j*Ndim_rank[0];       //imin
             int l1   = (j+1)*Ndim_rank[0] - 1;   //imax

             //fprintf(out, "bufferRRR %f %f  \n", buffer[j], buffer[j+Ndim_rank[1]]  );
             Tout[l0] = buffer[j + shift_Imin ];
             Tout[l1] = buffer[j + shift_Imax ];
            }
          }
        }



      } // sequentiell/Parallel(Mpi=1 ou N)
     }  // Nstepmax

      // calcul temperature moyenne et conductivitÃ© 
       moy_loc[0] =0;
       moy_loc[1] =0;
          for (int64_t j = 0; j < Ndim_rank[1] ; ++j ){ 
            for (int64_t i = 0; i < Ndim_rank[0] ; ++i ){ 

             int l = j*Ndim_rank[0]+ i;
             moy_loc[0] +=Tout[l]; 
             moy_loc[1] +=kdiff[l]; 
            }
          }
       moy_loc[0] = moy_loc[0]/float(Ndim_rank[1]*Ndim_rank[0]);
       moy_loc[1] = moy_loc[1]/float(Ndim_rank[1]*Ndim_rank[0]);
       MPI_Reduce( moy_loc, moy, 2 , MPI_DOUBLE, MPI_SUM, 0, globalComm);
    
       if (rank==0)  fprintf(out, "Tmoy  %f %f \n", nit*dt, moy[0]/float(Max_rank) );
    }  // Nitmax

  for (int64_t i = 1; i < Ndim_rank[0]-1 ; ++i ){ 
    for (int64_t j = 1; j < Ndim_rank[1]-1 ; ++j ){ 

    int    l = j*Ndim_rank[0]+ i;
    fprintf(out, " Final %f %f %.8f \n", 0.5*(x[i]+x[i+1]),0.5*(y[j]+y[j+1]), temp0[l]); 
   }
    fprintf(out, " Final \n"); 
  }
/*
*/

  double tensDt = chrono.click();
  fprintf(out, "Temps initialisation   %f  \n", tensDt);

  fclose(out);

  //delete [] temp0;  delete [] temp1; delete [] bilan; delete [] x; delete[] request;
  delete [] temp0;  delete [] temp1; delete [] bilan; delete [] x;

  MPI_Finalize();

  return EXIT_SUCCESS;
}
