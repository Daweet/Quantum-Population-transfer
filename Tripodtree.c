/***************************************************************************************/
/*--FourLevel TRIPOD system for Tree     --*/
/*--Integrating 4x4 FSTIRAP using predictor corrector algorithm;    --*/
/*--The system is taken to be off  resonace               --*/
/*--Same intensity for both P and S pulse taken            --*/
/*--P is for pump pulse, and S is for stokes pulse in EIT's   --*/
/**********************************************/
/***************************************************************************************/
/*--Initial population at |0> , state |3> is the highest level    --*/
/*Relaxation term included in the detuining that is the environment*/
/*--The system is taken to be ON resonace NON-DISSIPATIVE  ..2 photon resonace             --*/
/*--State |3> couples with |0> (Pump0), with |1> (stokes0 and pump1 at d/t times), with |2> (stokes1)     --*/
/*--Enerergy leveles arranged as E3 > E2 > E1 > E0            --*/
/*--Const detuning ....... Coherences Long time to    --*/
/***************************************************************************************/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
//#include "gasdev.c"
/*  #include <fftw3.h> */

#define PI 3.14159265
int const  N=4;

/***************************** parameters  ****************************************************************/
double const sigma=0.413e12  ;
double const E0=(50.0/(0.413e12));   /*E0=0.1*(20/sigma)*/
double const gamma1=-((1.0)/((16.4)*0.413e12))  ;/*(1.0/(67.768e11))*/
double const gamma2=-((1.0)/(0.413e19))/*(1.0)/(41.32e17)19*/  ;
//double const gamma1=-1.0/(67.77e11)  ;
//double const gamma2=-1.0/(41.32e17)  ;
double const dt=10.0e8;/*works for very fast modulation 10.0e3*/
double const tmax=94.5e12 ; /*64.5e13*/

/**********************************************************************************************************/

void predic ( double dt, int N, double A[] ,double A_1[] ,double A_2[] ,double A_3[] ,double A_4[],
             double B[] ,double B_1[] ,double B_2[] ,double B_3[] ,double B_4[]);
void corre ( double dt, int N, double A[] ,double A_1[] ,double A_2[] ,double A_3[] ,double A_4[],
            double B[] ,double B_1[] ,double B_2[] ,double B_3[] ,double B_4[], double H[N][N]);

void initialize ( int N, double A[] ,double A_1[] ,double A_2[] ,double A_3[] ,double A_4[],
                 double B[] ,double B_1[] ,double B_2[] ,double B_3[] ,double B_4[], double z[], double phi[]);

double pump,stokes;
double E_field0(double t),pump0,E_fild0(double t),stokes0;
double E_field1(double t),pump1,E_fild1(double t),stokes1;
double E_field2(double t),pump2,E_fild2(double t),stokes2;
double Fluctn (double t),Fluc;
double C_1,C_2,C_3,totp,dr;
double delta,delta1,delta2;double n,val,val1,val2;
int tr;

/* void FFT(int steps); */
int main(void)
{
    
    //long idum; /*idum is seed*/
    //idum=-123456789;
    //printf("Value: %ld",&idum);
    
    double  norm, Ex, Ey, Ez, end_pulse=12.5*(sigma) ;
    int     i=0,j=0, steps=0 , steps0=0, steps1=0, flag=0, flag2=0,counter;
    char    st[130];
    double  tauc ;
    double  A[N] , A_1[N] , A_2[N] , A_3[N] , A_4[N], B[N] ,B_1[N] , B_2[N] , B_3[N] , B_4[N], z[N], phi[N],  ran_phi[N];
    double  H[N][N], H0[N][N],Aerr[N],Berr[N];
    
    
    
    
    FILE * Cpr31;
    Cpr31 = fopen("Cpr31.txt","w");
    
    
    FILE * ChrR31;
    ChrR31 = fopen("ChrR31.txt","w");
    
    FILE * ChrI31;
    ChrI31 = fopen("ChrI31.txt","w");
    
    FILE * radn1;
    radn1 = fopen("radn1.txt","w");
    
    
    FILE * Delta;
    Delta = fopen("Delta.txt","w");
    
    
    //for(tr=0;tr<1;tr++){
    double  t=-4.5e12 ;
    
    
    /*
     sprintf(st,"/home/scr/dawit/random3/result/pop_tr%d",tr);
     Cpr31 = fopen(st,"w");
     
     sprintf(st,"/home/scr/dawit/random3/result/Chr_tr%d",tr);
     Chr31 = fopen(st,"w");
     
     sprintf(st,"/home/scr/dawit/random3/result/fluc_tr%d",tr);
     radn1 = fopen(st,"w");
     
     sprintf(st,"/home/scr/dawit/random3/result/Detun_tr%d",tr);
     Delta = fopen(st,"w");
     */
    
    initialize ( N,  A , A_1 , A_2 , A_3 , A_4, B ,B_1 , B_2 , B_3 , B_4, z, phi );
    
    //idum += 123458 ;
    //printf("Value: %ld",&idum);
    //val=/*5.0+5.0**/gasdev(&idum);
    
    while(t<tmax){  /****** main loop *********/
        
        pump0=1.0*E_field0(t);
        stokes0=1.0*E_fild0(t);
        pump1=1.0*E_field1(t);/*0.6*/
        stokes1=1.0*E_fild1(t);
        pump2=1.0*E_field2(t);
        stokes2=1.0*E_fild2(t);
        pump=pump0+pump1+1.0*pump2;
        stokes=stokes0+stokes1+1.0*stokes2;
        /**************************************************/
        //val=Fluctn(t);
        
        double p,k;double h,ep;
        double  lambda,D; double  tcor=10.0e9;/*tcor=10.0e14,,D=10.0e13;*/
        lambda=(1.0)/(tcor);
        D=10.0e10;
        p= (double)rand()/(double)RAND_MAX;
        k= (double)rand()/(double)RAND_MAX;
        /*
         n=ran1(&idum);
         k=ran1(&idum);*/
        ep=exp(-dt/tcor);
        h=sqrt(-2.0*D*lambda*(1.0-ep*ep)*log(p))*cos(2.0*PI*k);
        val=(val*ep)+(h);
        /**************************************************/
        //val=0.0;
        //val=/*5.0+5.0**/gasdev(&idum);
        
        
        n=(0.2*(val-0.5));/*a=0.2*/
        //n=val;
        delta=0.0*(1.376e-12)/*1.376e-12 works too*/;
        
        delta1=(1.376e-12)*n;
        delta2=0.2*(1.376e-12)*n;/*check this value???*/
        /******** calculate H(t) ********/
        
        H[0][0]=H[0][1]=H[0][2]=0.0;
        H[0][3]=0.5*(pump0+stokes2);
        
        H[1][0]=H[1][2]=0.0;;
        H[1][1]=1.0*(delta);
        H[1][3]=0.5*(stokes0+pump1);
        
        
        H[2][0]=H[2][1]=0.0;
        H[2][2]=1.0*(delta);/*delta1/2*/
        H[2][3]=0.5*(stokes1+pump2);
        
        
        
        H[3][0]=0.5*(pump0+stokes2);
        H[3][1]=0.5*(stokes0+pump1);
        H[3][2]=0.5*(stokes1+pump2);H[3][3]=0.0*(delta)/*off two resonance*/;
        
        /*If H[3][3]=2.0*(delta) trasfer effiececy decreases*/
        
        
        /********************************************/
        
        //if (t>-3.0e12 && t<3.0e12){
        //if (fmod(t,10.0e9)==0.00000000000000000){/*printing every 10.0e8 */
        fprintf(Cpr31,"%e\t %17.17e\t %17.17e\t %17.17e\t %17.17e\n", t,A[0]*A[0]+B[0]*B[0],A[1]*A[1]+B[1]*B[1],A[2]*A[2]+B[2]*B[2],A[3]*A[3]+B[3]*B[3] );
        
        //fprintf(Chr31,"%e\t %17.17e\t %17.17e\t %17.17e\t %17.17e\t %17.17e\t %17.17e\n", t,0.5*(A[0]+B[3]),0.5*(A[1]+B[4]),0.5*(A[2]+B[5]),0.5*(B[0]-A[3]), 0.5*(B[1]-A[4]),0.5*(B[2]-A[5]) );
         /*************order of coherences 01,02,03,12,13,23*******************************/
        fprintf(ChrR31,"%e\t %17.17e\t %17.17e\t %17.17e\t %17.17e\t %17.17e\t %17.17e\n", t,A[0]*A[1]+B[0]*B[1],A[0]*A[2]+B[0]*B[2],A[0]*A[3]+B[0]*B[3],A[1]*A[2]+B[1]*B[2],A[1]*A[3]+B[1]*B[3],A[2]*A[3]+B[2]*B[3] );
        
        fprintf(ChrI31,"%e\t %17.17e\t %17.17e\t %17.17e\t %17.17e\t %17.17e\t %17.17e\n", t,A[0]*B[1]-B[0]*A[1],A[0]*B[2]-B[0]*A[2],A[0]*B[3]-B[0]*A[3],A[1]*B[2]-B[1]*A[2],A[1]*B[3]-B[1]*A[3],A[2]*B[3]-B[2]*A[3] );
        
        
        
        //}

        
        /********************************************/
        //if (t>-3.0e12 && t<3.0e12){
        //if (t>-3.0e8 && t<3.0e8){
        //if (fmod(t,10.0e9)==0.00000000000000000){/*printing every 10.0e8 */
        fprintf(radn1,"%e\t %17.17e\t %17.17e\t %17.17e\n", t,delta,pump,stokes );
        //}
        
        //if (t>-3.0e12 && t<3.0e12){
        //if (t>-3.0e8 && t<3.0e8){
        fprintf(Delta,"%e\t %17.17e\n",  t,(1.376e-12)*(0.1*(val)));
        //}
        
        
        predic (dt,  N,  A , A_1 , A_2 , A_3 , A_4, B ,B_1 , B_2 , B_3 , B_4 );
        corre  (dt,  N,  A , A_1 , A_2 , A_3 , A_4, B ,B_1 , B_2 , B_3 , B_4, H );
        //corre  (dt,  N,  A , A_1 , A_2 , A_3 , A_4, B ,B_1 , B_2 , B_3 , B_4, H );
        
        t=t+dt;
        
        
        
    }  /* for while t loop*/
    
    
    
    
    
    fclose(Cpr31);
    
    fclose(ChrR31);
    
    fclose(ChrI31);
    
    fclose(radn1);
    
    fclose(Delta);
    
    
    
    
    
    
}/*main*/
//} /*tr*/
/*************************************************************************/
void initialize ( int N, double A[] ,double A_1[] ,double A_2[] ,double A_3[] ,double A_4[],
                 double B[] ,double B_1[] ,double B_2[] ,double B_3[] ,double B_4[], double z[], double phi[]){
    
    int  i  ;
    
    
    for( i=0 ; i<N ; i++){
        
        
        A[i]=0.0;
        B[i]=0.0;
        
        
        A_1[i] = 0.0;
        A_2[i] = 0.0;
        A_3[i] = 0.0;
        A_4[i] = 0.0;
        
        
        B_1[i] = 0.0;
        B_2[i] = 0.0;
        B_3[i] = 0.0;
        B_4[i] = 0.0;
        
        A[0]=1.0;
        
    }
    
    
    
    
}/** initialize****/

/*****************************************************************************************************/
double E_field0 (double t){
    
    double E,  sigmaT, t0;
    
    
    sigmaT=sigma;
    
    t0=3.9*sigmaT  ;
    
    
    
    E= (E0)*sin(PI/2)*exp( -(t-t0)*(t-t0)/(sigmaT*sigmaT) );
    
    
    
    return(E);
    
    
}

/**************************************************************************************************/
double E_fild0 (double t){
    
    double E,  sigmaT, t0,t1;
    
    
    sigmaT=sigma;
    
    t0=3.9*sigma  ;
    
    t1=2.5*sigma  ;
    
    E= (E0*cos(PI/2)*exp( -(t-t0)*(t-t0)/(sigmaT*sigmaT) ))+(E0*exp( -(t-t1)*(t-t1)/(sigmaT*sigmaT) ));
    
    
    return(E);
    
    
}
/**************************************************************/
/*****************************************************************************************************/
double E_field1 (double t){
    
    double E,  sigmaT, t0;
    
    
    sigmaT=sigma;
    
    t0=18.9*sigmaT  ;
    
    
    E= (E0)*sin(PI/4)*exp( -(t-t0)*(t-t0)/(sigmaT*sigmaT) );
    
    
    
    return(E);
    
    
}

/**************************************************************************************************/
double E_fild1 (double t){
    
    double E,  sigmaT, t0,t1;
    
    
    sigmaT=sigma;
    
    t0=18.9*sigma  ;
    t1=17.5*sigma  ;
    
    E= (E0*cos(PI/4)*exp( -(t-t0)*(t-t0)/(sigmaT*sigmaT) ))+(E0*exp( -(t-t1)*(t-t1)/(sigmaT*sigmaT) ));
    
    
    return(E);
    
    
}

/**************************************************************/
/*****************************************************************************************************/
double E_field2 (double t){
    
    double E,  sigmaT, t0;
    
    
    sigmaT=sigma;
    
    t0=34.9*sigmaT  ;
    
    
    E= (E0)*sin(PI/2)*exp( -(t-t0)*(t-t0)/(sigmaT*sigmaT) );
    
    
    
    return(E);
    
    
}

/**************************************************************************************************/
double E_fild2 (double t){
    
    double E,  sigmaT, t0,t1;
    
    
    sigmaT=sigma;
    
    t0=34.9*sigma  ;
    
    t1=33.5*sigma  ;
    
    E= (E0*cos(PI/2)*exp( -(t-t0)*(t-t0)/(sigmaT*sigmaT) ))+(E0*exp( -(t-t1)*(t-t1)/(sigmaT*sigmaT) ));
    
    
    return(E);
    
    
}
/**************************************************************************************************/
/**************************************************************/
void predic ( double dt, int N, double A[] ,double A_1[] ,double A_2[] ,double A_3[] ,double A_4[],
             double B[] ,double B_1[] ,double B_2[] ,double B_3[] ,double B_4[]){
    
    int  i  ;
    double  c1, c2, c3, c4;
    
    c1 = dt;
    c2 = c1 * dt / 2.0;
    c3 = c2 * dt / 3.0;
    c4 = c3 * dt / 4.0;
    
    
    
    for( i=0 ; i<N ; i++){
        
        A[i]   = A[i]   + c1*A_1[i] + c2*A_2[i] + c3*A_3[i] + c4*A_4[i];
        A_1[i] = A_1[i] + c1*A_2[i] + c2*A_3[i] + c3*A_4[i];
        A_2[i] = A_2[i] + c1*A_3[i] + c2*A_4[i];
        A_3[i] = A_3[i] + c1*A_4[i];
        
        B[i]   = B[i]   + c1*B_1[i] + c2*B_2[i] + c3*B_3[i] + c4*B_4[i];
        B_1[i] = B_1[i] + c1*B_2[i] + c2*B_3[i] + c3*B_4[i];
        B_2[i] = B_2[i] + c1*B_3[i] + c2*B_4[i];
        B_3[i] = B_3[i] + c1*B_4[i];
        
    } /* for */
    
    
    
}  /* function predic  */

/********************************************************************************************/

void corre ( double dt, int N, double A[] ,double A_1[] ,double A_2[] ,double A_3[] ,double A_4[],
            double B[] ,double B_1[] ,double B_2[] ,double B_3[] ,double B_4[], double H[N][N]){
    
    
    int  i,j;
    double   A_1i,B_1i;
    double   CORA, CORB;
    double   C1, C2, C3, C4;
    double   CA, CA_2, CA_3, CA_4,CB, CB_2, CB_3, CB_4;
    double   GEAR0, GEAR2, GEAR3, GEAR4;
    
    
    double  SA[N] , SA_1[N] , SA_2[N] , SA_3[N] , SA_4[N], SB[N] ,SB_1[N] , SB_2[N] , SB_3[N] , SB_4[N];
    
    GEAR0 = 251.0 / 720.0  ;
    GEAR2 = 11.0 / 12.0    ;
    GEAR3 = 1.0 / 3.0     ;
    GEAR4 = 1.0 / 24.0    ;
    
    C1 = dt;
    C2 = C1 * dt / 2.0;
    C3 = C2 * dt / 3.0;
    C4 = C3 * dt / 4.0;
    
    CA   = GEAR0 * C1;
    CA_2 = GEAR2 * C1 / C2;
    CA_3 = GEAR3 * C1 / C3;
    CA_4 = GEAR4 * C1 / C4;
    
    CB   = GEAR0 * C1;
    CB_2 = GEAR2 * C1 / C2;
    CB_3 = GEAR3 * C1 / C3;
    CB_4 = GEAR4 * C1 / C4;
    
    for( i=0 ; i<N ; i++){
        SA[i]   = A[i];
        SA_1[i] = A_1[i];
        SA_2[i] = A_2[i];
        SA_3[i] = A_3[i];
        SA_4[i] = A_4[i];
        
        SB[i]   = B[i];
        SB_1[i] = B_1[i];
        SB_2[i] = B_2[i];
        SB_3[i] = B_3[i];
        SB_4[i] = B_4[i];
    }
    
    for( i=0 ; i<N ; i++){
        
        A_1i   = 0.0;
        for( j=0 ; j<N ; j++){
            
            A_1i= A_1i+ SB[j]*H[i][j] ;/*idA/dt=iHB*/
            
        }
        
        
        B_1i   = 0.0;
        for( j=0 ; j<N ; j++){
            
            B_1i= B_1i - SA[j]*H[i][j] ;/*dB/dt=-HA*/
            
        }
        
        CORA   = A_1i - SA_1[i];
        
        A[i]   = SA[i]   + CA   * CORA;
        A_1[i] = A_1i;
        A_2[i] = SA_2[i] + CA_2 * CORA;
        A_3[i] = SA_3[i] + CA_3 * CORA;
        A_4[i] = SA_4[i] + CA_4 * CORA;
        
        /******************************/
        
        
        
        CORB   = B_1i - SB_1[i];
        
        B[i]   = SB[i]   + CB   * CORB;
        B_1[i] = B_1i;
        B_2[i] = SB_2[i] + CB_2 * CORB;
        B_3[i] = SB_3[i] + CB_3 * CORB;
        B_4[i] = SB_4[i] + CB_4 * CORB;
        
    } /*for*/
    
    
}  /* correc function */

/**************************************************************************************************/
