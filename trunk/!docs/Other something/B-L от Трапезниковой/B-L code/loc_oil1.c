/*    plane-filtration      */

#include <stdio.h>
#include <math.h>

typedef char * LPSTR;

/*     . . . 0 < omega < 1         */
#define OMEGAI(XX) (2./( 1. + sqrt(1.+(XX))))

/*     . . . 1 < omega < 2         */
#define OMEGAR(XX) (2./( 1. + sqrt(1.-(XX))))

/*     . . . limits                */
#define OMEGA0(XX) 1.e-4
#define OMEGA2(XX) (2.-1.e-4)

#define clrscr() system("cls")

#define K 0.5e-12     //абсолютная проницаемость
#define M 0.2          //пористость
#define MUW 1.2e-13    //вязкость воды
#define MUO 6.e-13     //вязкость нефти
#define PRESS 100.     //начальное давление
#define N1 2           //степени в функциях проницаемости
#define N2 2           //степени в функциях проницаемости
#define N3 0.5		//степени в функциях проницаемости
#define H 10.		//толщина пласта
#define SZ 0.8		//Критическая водонасыщенность
#define S_wr 0.1		//Связанная водонасыщенность (начальная?)
#define S1 0.70324	//Вспомогательное значение водонасыщенности. Используется в вычислениях проницаемости
#define QN 100.     //Мощность источников
#define ALEN 2.		//Связано с геометрией расположения скважин

#define NIH 5			//Связано с геометрией расположения скважин
#define NHALF 6        //Связано с геометрией расположения скважин
#define NJH 9			//Связано с геометрией расположения скважин	
#define NSIDE 4		//Связано с геометрией расположения скважин
#define Nx NIH*NHALF+3	//Размеры сетки
#define Ny NJH*NSIDE+3	//Размеры сетки
#define Nx1 Nx+1	//Размеры сетки
#define Ny1 Ny+1	//Размеры сетки
#define NT 1			//Количество подобластей (для параллельного варианта программы)
#define NJSIZE Ny1/NT	//Размеры сетки по одному из направлений в подобласти

#define nprint 200  //счетчик для записи результатов в файлы 

int lrpar(double ccs,double ccw,double cce,double ccn,double cosnx2,
          double cosny2,double cosnxy,int i,int j);

void relaxrb(double cc, double yy[Nx1][Ny1]);

double FKW(double S);
double DFKW(double S);
double FKO(double S);
double FW(double S);
double FK(double S);
void OIL1();

/*    . . . 0 < omega < 2         */
double OMEGAC(double xxr, double xxi)
{
   return (2./( 1. + sqrt(1. - xxr + xxi/(1.-pow(xxr,0.333333333333333)))));
}


FILE * xgrid;
FILE * ygrid;
FILE * rec;
FILE * out;     /* file-description for user */

static double CW[Nx1][Ny1];           //Коэффициенты разностной схемы для расчёта давления
static double CS[Nx1][Ny1];           // 
static double CN[Nx1][Ny1];       //   
static double CE[Nx1][Ny1];       //
static double F[Nx1][Ny1];        //Правая часть для этой разностной схемы
static double Sw_old[Nx1][Ny1];	//Водонасыщенность на предыдущем слое
static double Sw_new[Nx1][Ny1];	//Водонас. на новом временном слое
static double P[Nx1][Ny1];		//Давление (вычисляется через водонасыщенность??)
static double Q[Nx1][Ny1];		//Накачиваемая в скважины вода или извлекаемая жидкость 
static double SWRES[Nx1][Ny1];	//Результирующая водонасыщенность (параллельный алгоритм??)
static double PRES[Nx1][Ny1];		//Результирующее давление
static double F_bl[Nx1][Ny1];		//Функция Баклея-Леверетта ??
static double omg[Nx1][Ny1];		//Параметры локальной релаксации
static double SS[Nx1][NJSIZE];	//Дополнитеьные массивы для решения уравнения Баклея-Леверетта
static double SSS[Nx1][NJSIZE];  //
static double FN[Nx1][NJSIZE];   //
static double FN1[Nx1][NJSIZE];  // 
//static double P0[NI1][NJ1];       
static double QTR[Nx1][NJSIZE];	//Дебиты на скважинах
static double FW0TR[Nx1][NJSIZE];//Значени функции Баклея-Леверетта на скважинах
//static double omg[Nx1][Ny1];
static double X[Nx1];         /* x-coordinate */
static double Y[Ny1];         /* y-coordinate */

double cosnx, cosny, cosnxy, cosnx2, cosny2;      

double EPSGAM;
double EPS;
double HX,HY,HT,AK,A1,A2;

int N,NIS1,NJS1,NIS2,NJS2,MAXDEL,MAXGAM;
int I,J,III,INDIC;

double QHALF,QDOUBLE,TIME,DPDX1,DPDX2,DPDY1,DPDY2,SUM;
double A,B,C,ADASH,BDASH,flag;
double AMAXK,AMAXPX,PX,AMAXPY,PY,AMAXP;

char ch,namefl[13];


double FW(double S)
{
      return(FKW(S)/MUW)/((FKW(S)/MUW)+(FKO(S)/MUO));
//		  return FW(S);
}

double FKW(double S)
{
      if((S_wr <= S) && (S < S1)) return(pow(((S-S_wr)/(SZ-S_wr)), N2));
      else
		  if((S1 <= S) && (S <= 1.)) return(pow(0.8*(S-S_wr)/(SZ-S_wr) ,N3));
      else
		  if((0. <= S) && (S < S_wr)) return(0.);
      else
		  if(S < 0.) return(0.);
      else
		  if(S > 1.) return(pow( 8.*(S-S_wr)/(SZ-S_wr) ,N3));
//		  return FKW(S);
      }

double DFKW(double S)
{
      if((S_wr <= S) && (S < S1)) return(pow( N2*(S-S_wr),(N2-1))/pow((SZ-S_wr),N2));
      else
		  if((S1 <= S) && (S <= 1.)) return(pow(8.*N3*(S-S_wr),(N3-1))/pow((SZ-S_wr),N3));
      else
		  if((0. <= S) && (S < S_wr)) return(0.);
      else
		  if(S < 0.) return(0.);
      else
	      if(S > 1.) return(pow(8.*N3*(S-S_wr),(N3-1))/pow((SZ-S_wr),N3));
//		  return DFKW(S);
}

double FKO(double S)
{
      if((S_wr <= S) && (S <= SZ)) return(pow((SZ-S)/(SZ-S_wr),N1));
      else
		  if((0. <= S) && (S < S_wr)) return(1.);
      else 
		  if((SZ < S) && (S <= 1.)) return(0.);
      else 
		  if(S < 0.) return(1.);
      else 
		  if(S > 1.) return(0.);
//		  return FKO(S);
}

double FO(double S)
{
      return((FKO(S)/MUO)/((FKW(S)/MUW)+(FKO(S)/MUO)));
//	  		  return FO(S);
}

double FK(double S)
{
      return(-K*(FKW(S)/MUW+FKO(S)/MUO));
//		  return FK(S);
}   
    
/******************************************************/
/*                                                    */
/*            MAIN - PROGRAM                          */
/******************************************************/

void main()
{

 NIS1=Nx-1;
 NIS2=Nx-2;
 NJS1=Ny-1;
 NJS2=Ny-2;
 EPS=1.e-4;

 HX=5.;
 HY=5.;

  X[1]=0.;                 /* initial coordinates */
  Y[1]=0.;

  for(I=2; I<=Nx; I++)     /* grid in x-direction */
     {
     X[I]=X[I-1]+HX;
     }

  for(J=2; J<=Ny; J++)     /* grid in y-direction */
     {
     Y[J]=Y[J-1]+HY;
     }

//-----------initial zero values----------

  for(J=0; J<=Ny; J++)
        for(I=0; I<=Nx; I++)
        {
         CW[I][J]=0.;
         CS[I][J]=0.;
         CN[I][J]=0.;
         CE[I][J]=0.;
         F[I][J]=0.;
		 omg[I][J]=0.;
		}  

      cosnx = cos( 3.1415926/(Nx-3) );     /* const for over-relaxation  */
      cosnx2 = cosnx*cosnx;
      cosny = cos( 3.1415926/(Ny-3) );
      cosny2 = cosny*cosny;
      cosnxy =cosnx*cosny;

//------------------------------------

       OIL1();

//------------------------------------


}

/*********************************************
*                                            *
*           OIL - PREDICTION                 *
*                                            *
*********************************************/

void OIL1()

{


  for(J=0; J<=Ny; J++)
        for(I=0; I<=Nx; I++)
        {
         Q[I][J]=0.;
         F_bl[I][J]=0.;
		}

      QHALF=0.5*QN;
      QDOUBLE=2.*QN;
      Q[NIH+2][2]=-QHALF;
      Q[NIH*3+2][2]=QN;
      Q[NIH*5+2][2]=-QHALF;
      Q[2][NJH+2]=QN;
      Q[NIH*2+2][NJH+2]=-QN;
      Q[NIH*4+2][NJH+2]=-QN;
      Q[NIH*6+2][NJH+2]=QN;
      Q[NIH+2][NJH*2+2]=-QN;
      Q[NIH*3+2][NJH*2+2]=QDOUBLE;
      Q[NIH*5+2][NJH*2+2]=-QN;
      Q[2][NJH*3+2]=QN;
      Q[NIH*2+2][NJH*3+2]=-QN;
      Q[NIH*4+2][NJH*3+2]=-QN;
      Q[NIH*6+2][NJH*3+2]=QN;
      Q[NIH+2][NJH*4+2]=-QHALF;
      Q[NIH*3+2][NJH*4+2]=QN;
      Q[NIH*5+2][NJH*4+2]=-QHALF;

  for(J=0; J<=Ny; J++)
        for(I=0; I<=Nx; I++)
        {
		if ( Q[I][J] > 0.)
         {
          F_bl[I][J]=FW(SZ); 
         }         
		}
  for(J=0; J<=Ny; J++)
        for(I=0; I<=Nx; I++)
        {
         QTR[I][J]=Q[I][J];
         FW0TR[I][J]=F_bl[I][J];
		}


/*---------------- initial values ------------------------*/

  for(J=0; J<=Ny; J++)
        for(I=0; I<=Nx; I++)
        {
         Sw_old[I][J]=S_wr;
		 P[I][J]=PRESS;
		}

     AK=K/(M*MUW);
  

label25: N++;
		 TIME+=HT;

		 printf(" \n");
		 printf("Step N = %d, Time Step HT = %6.4e, TIME = %6.4e \n",N,HT,TIME);

  III=1;
  for(J=0; J<=Ny; J++)
        for(I=0; I<=Nx; I++)
        {
         SS[I][J]=Sw_old[I][J];
		}

label702:  for(J=0; J<=Ny; J++)
        for(I=0 ; I<=Nx; I++)
        {
      
       if ( QTR[I][J] < 0. ) FW0TR[I][J]=FW(SS[I][J]);       
		}

  for(J=2; J<=NJS1; J++)
  {
        for(I=2; I<=NIS1; I++)
        {
        DPDX1=(FK(SS[I+1][J])+FK(SS[I][J]))/2.*(P[I+1][J]-P[I][J])/HX;
        DPDX2=(FK(SS[I][J])+FK(SS[I-1][J]))/2.*(P[I][J]-P[I-1][J])/HX;
        DPDY1=(FK(SS[I][J+1])+FK(SS[I][J]))/2.*(P[I][J+1]-P[I][J])/HY;
        DPDY2=(FK(SS[I][J])+FK(SS[I][J-1]))/2.*(P[I][J]-P[I][J-1])/HY;

        A1=(((DPDX1+fabs(DPDX1))/2.-(DPDX2-fabs(DPDX2))/2.)*FW(SS[I][J])-
			 (DPDX2+fabs(DPDX2))/2.*FW(SS[I-1][J])+(DPDX1-fabs(DPDX1))/
               2.*FW(SS[I+1][J]))/HX;
        A2=(((DPDY1+fabs(DPDY1))/2.-(DPDY2-fabs(DPDY2))/2.)*FW(SS[I][J])-
			 (DPDY2+fabs(DPDY2))/2.*FW(SS[I][J-1])+(DPDY1-fabs(DPDY1))/
              2.*FW(SS[I][J+1]))/HX;
           if(III == 1) 
		   {
           Sw_new[I][J]=SS[I][J]+0.01;
           FN[I][J]=M/HT*(SS[I][J]-Sw_old[I][J])-
			       QTR[I][J]/(H*HX*HY)*FW0TR[I][J]+A1+A2;
		   }
	       else
			   if  (fabs(SS[I][J]-SSS[I][J]) > 1.e-3)
			   { 
				FN[I][J]=M/HT*(SS[I][J]-Sw_old[I][J])-QTR[I][J]/(H*HX*HY)*
                   FW0TR[I][J]+A1+A2;
                Sw_new[I][J]=SS[I][J]-(SS[I][J]-SSS[I][J])*FN[I][J]/
                     (FN[I][J]-FN1[I][J]);
		   
			   }
		}
        Sw_new[1][J]=Sw_new[2][J];
        Sw_new[Nx][J]=Sw_new[NIS1][J];
        }

        for(I=1; I<=Nx; I++)
		{
         Sw_new[I][1]=Sw_new[I][2];
         Sw_new[I][Ny]=Sw_new[I][NJS1];
        }
      SUM=0.;
      for(J=1; J<=Ny; J++)
        for(I=1; I<=Nx; I++)
		{
         SUM+=pow((Sw_new[I][J]-SS[I][J]),2)*HX*HY;
        }

      SUM=sqrt(fabs(SUM));

      if(SUM < 1.e-4) goto label705;
	  
      INDIC=0;

      III=III+1;

      for(J=1; J<=Ny; J++)
        for(I=1; I<=Nx; I++)
		{
         SSS[I][J]=SS[I][J];
         SS[I][J]=Sw_new[I][J];
         FN1[I][J]=FN[I][J];
		}
      goto label702;
	  
label705:   INDIC=1;

	printf("END of SWNEW:     SW_center = %6.4e \n",Sw_new[NIH*3+2][NJH*2+2]);

      for(J=2; J<=NJS1; J++)
        for(I=2; I<=NIS1; I++)
		{
            A = (FK(Sw_new[I][J])+FK(Sw_new[I-1][J]))/(2.*HX*HX);
            B = (FK(Sw_new[I+1][J])+FK(Sw_new[I][J]))/(2.*HX*HX);
            ADASH = (FK(Sw_new[I][J])+FK(Sw_new[I][J-1]))/(2.*HY*HY);
            BDASH = (FK(Sw_new[I][J+1])+FK(Sw_new[I][J]))/(2.*HY*HY);

            C = A+B+ADASH+BDASH;

            F[I][J]=-QTR[I][J]/(H*HX*HY)/C;

            CW[I][J] = A/C;
            CE[I][J] = B/C;
            CS[I][J] = ADASH/C;
            CN[I][J] = BDASH/C;

                 flag=lrpar(CS[I][J], CW[I][J], CE[I][J], CN[I][J],
                        cosnx2, cosny2, cosnxy, I, J);

                 if(flag!=0)
                    {
                       printf("Omega evaluation error \n");
                       exit(1);
                    }
			}

//----------------------------------------------------

         relaxrb(C,P);	  

//----------------------------------------------------

      for(J=1; J<=Ny; J++)
        for(I=1; I<=Nx; I++)
         Sw_old[I][J]=Sw_new[I][J];


      printf("			TIME STEP IS OVER \n");

/*-----------records in files------------------------------*/

	  if((N % nprint)==0)
		{
		 printf("Attention! I am writting!\n");

		 /*-----------two-dimensional grafics --------------*/

		 sprintf(namefl, "sw_%d.dat",N);
         out=fopen(namefl,"w");
		 fprintf(out,"TITLE =  \"Water Saturation\" \n");
		 fprintf(out,"VARIABLES = \"X\",\"Y\",\"SWNEW\" \n");
         fprintf(out,"ZONE T = \"BIG ZONE\", I=%d  ,J=%d, F = POINT \n", Nx, Ny);

			for(J=1; J<=Ny; J++)
		       {			
			for(I=1; I<=Nx; I++)
			   {
			     fprintf(out,"%6.4e %6.4e %6.4e \n",X[I],Y[J],Sw_new[I][J]);
				     }
			   }
		 fclose(out);

		 sprintf(namefl, "p_%d.dat",N);
         out=fopen(namefl,"w");
		 fprintf(out,"TITLE =  \"Pressure\" \n");
		 fprintf(out,"VARIABLES = \"X\",\"Y\",\"P\" \n");
         fprintf(out,"ZONE T = \"BIG ZONE\", I=%d  ,J=%d, F = POINT \n", Nx, Ny);

			for(J=1; J<=Ny; J++)
		       {			
			for(I=1; I<=Nx; I++)
			   {
			     fprintf(out,"%6.4e %6.4e %6.4e \n",X[I],Y[J],P[I][J]);
				     }
			   }
		  fclose(out);
		}


  if(N <= 40000) goto label25;

//      goto label25;


}


int lrpar (double ccs,double ccw,double cce,double ccn,
           double cosnx2,double cosny2,double cosnxy,int i,int j )

/*=====================================================================C
C     Local relaxation parameter evaluation based on:                  C
C                                                                      C
C     1. E.F.F. Botta, A.E.P. Veldman.                                 C
C        On Local Relaxation Methods and Their Application to          C
C        Convection-Diffusion Equations.                               C
C        J.Comput. Phys., 1982, v.48., N 1, pp.127-149.                C
C     2. L.W. Ehrlich.                                                 C
C        An Ad Hoc SOR Method.                                         C
C        J.Comput. Phys., 1981, v.44., N 1, pp.31-45.                  C
C     3. L.W. Ehrlich.                                                 C
C        The Ad-Hoc SOR Method: a Local Relaxation Scheme.             C
C        Elliptic Problem Solvers II                                   C
C        (Eds. G. Birkhoff and A. Schoenstadt), pp.257-269.            C
C        Academic Press, Orlando, 1984.                                C
C                                                                      C
C     for the solution of the following eqn:                           C
C                                                                      C
C         - CCS*Y(i,j-1) - CCW*Y(i-1,j) + Y(i,j)                       C
C         - CCE*Y(i+1,j) - CCN*Y(i,j+1) = FC(i,j)                      C
C                                                                      C
C----------------------------------------------------------------------*/
{
      double cwce,cscn,c;
      double rmur2, rmui2;
      int flag;

      flag = 0;
      cwce = ccw*cce;
      cscn = ccs*ccn;
      c = cwce*cscn;


      if ( c < 0. )
         {
         /*  . . . complex value . . .  */
/*         printf("c= %12.7e \n",c);*/

         if ( cwce > 0. )
            {
             rmur2 =  4.0*cwce*cosnx2;
             rmui2 = -4.0*cscn*cosny2;
            }
         else
            {
             rmur2 =  4.0*cscn*cosny2;
             rmui2 = -4.0*cwce*cosnx2;
            }
         if ( rmur2 >= 1. )
            omg[i][j] = OMEGA0( 0. );
         else
            omg[i][j] = OMEGAC(rmur2,rmui2);
         return flag;
         }
      if ( c == 0. )
         {
         /*  . . . zero coefficients . . . */
/*         printf("c= %12.7e \n",c);  */

         if ( cwce < 0. )
            {
             rmui2 =  -4.0*cwce*cosnx2;
             omg[i][j] = OMEGAI( rmui2 );
             return flag;
            }
         if ( cwce == 0. )
            {
             if ( cscn < 0. )
                {
                 rmui2 = -4.0*cscn*cosny2;
                 omg[i][j] = OMEGAI( rmui2 );
                 return flag;
                }
             if ( cscn == 0. )
                {
                 /*  . . . CWCE = CSCN = 0 !!!  */
                 flag = 1;
                 return flag;
                }
             if ( cscn > 0. )
                {
                 rmur2 =  4.0*cscn*cosny2;
                 if ( rmur2 >= 1.0 )
                     omg[i][j] = OMEGA2( 2. );
                 else
                     omg[i][j] = OMEGAR( rmur2 );
                 return flag;
                }
             }
         if ( cwce > 0. )
            {
             rmur2 =  4.0*cwce*cosnx2;
             if ( rmur2 >= 1.0 )
                 omg[i][j] = OMEGA2( 2. );
             else
                 omg[i][j] = OMEGAR( rmur2 );
             return flag;
            }
         }
      if ( c > 0. )
         {
          /*     . . . singular case . . .   */
/*         printf("c= %12.7e \n",c);
          printf("singular case! n\");  */
          if ( cwce > 0.0 )
             {
              rmur2 =  4.0*
                 ( cwce*cosnx2 + 2.0*sqrt(c)*cosnxy + cscn*cosny2 );
/*              printf("rmur2= %12.7e \n",rmur2);  */

              if ( rmur2 >= 1.0 )
                  omg[i][j] = OMEGA2( 2. );
              else
                  omg[i][j] = OMEGAR( rmur2 );
              }
          else
              {
               rmui2 = -4.0*
                  ( cwce*cosnx2 - 2.0*sqrt(c)*cosnxy + cscn*cosny2 );
               omg[i][j] = OMEGAI( rmui2 );
              }
          return flag;
          }
    return flag;
} /* end of lrpar  */


/**************************************************
 *                                                *
 *             OVER-RELAXATION                    *
 *      with red/black data ordering              *
 *                                                *
 **************************************************/

void relaxrb(double cc, double yy[Nx1][Ny1])
{
      int i,j,ind,itr;
      double d,work,omega;

      itr=0;
	  ind=0;


L10:  ++itr;
      for (j = 2; j <= NJS1; j++)
	  {
          if(j%2==0)
             {
             for (i = 2; i <= NIS1; i+=2)
			 {
               omega=omg[i][j];  
               yy[i][j]=(1.0-omega)*yy[i][j]+
                        omega*(CW[i][j]*yy[i-1][j]+CE[i][j]*yy[i+1][j]+
                        CS[i][j]*yy[i][j-1]+CN[i][j]*yy[i][j+1]+F[i][j]);
              }
        }
      else
        {
         for (i = 3; i <= NIS1; i+=2)
		 {
                  omega=omg[i][j];
                  yy[i][j]=(1.0-omega)*yy[i][j]+
                        omega*(CW[i][j]*yy[i-1][j]+CE[i][j]*yy[i+1][j]+
                        CS[i][j]*yy[i][j-1]+CN[i][j]*yy[i][j+1]+F[i][j]);
                  }
             }
          }

      for (j = 2; j <= NJS1; j++){
          if(j%2==0)
             {
             for (i = 3; i <= NIS1; i+=2){
                  omega=omg[i][j];
                  yy[i][j]=(1.0-omega)*yy[i][j]+
                        omega*(CW[i][j]*yy[i-1][j]+CE[i][j]*yy[i+1][j]+
                        CS[i][j]*yy[i][j-1]+CN[i][j]*yy[i][j+1]+F[i][j]);

                  }
             }
          else
             {
             for (i = 2; i <= NIS1; i+=2){
                  omega=omg[i][j];
                  yy[i][j]=(1.0-omega)*yy[i][j]+
                        omega*(CW[i][j]*yy[i-1][j]+CE[i][j]*yy[i+1][j]+
                         CS[i][j]*yy[i][j-1]+CN[i][j]*yy[i][j+1]+F[i][j]);

                  }
             }
          }
      /*-----------boundary conditions --------------------------------*/

      for (i = 1; i <= Nx; i++)
      {
           yy[i][Ny]=yy[i][NJS1];
           yy[i][1]=yy[i][2];
      }

      for (j = 1; j <= Ny; j++)
      {
           yy[Nx][j]=yy[NIS1][j];
           yy[1][j]=yy[2][j];
      }


      /*------------------ residual -----------------------------------*/

      d=0.;
      for (j = 2; j <= NJS1; j++){
           for (i = 2; i <= NIS1; i++){
                work=CW[i][j]*yy[i-1][j]-yy[i][j]+CE[i][j]*yy[i+1][j]
                +CS[i][j]*yy[i][j-1]+CN[i][j]*yy[i][j+1]+F[i][j];
                d+=work*work*HX*HY;
           }
      }

      d = sqrt(fabs(d*cc*cc));

      if(d > EPS)
        {
         ind=0;
//         printf("Itr= %d:  Indic= %d,   del=%12.10e \n",itr,ind,d);  
         goto L10;
        }
         ind=1;
/*         printf("Indic= %d  \n",ind);*/
	 /*         printf("END:  Itr= %d  \n",itr); */


   printf("END of LOCREL_RB: P_center  = %6.4e, niter = %d,  del = %6.4e \n",
	                                                  yy[NIH*3+2][NJH*2+2],itr,d);   

   return; 

} /* end of relaxrb */

