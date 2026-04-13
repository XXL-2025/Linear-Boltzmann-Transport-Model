//#include "Pythia8/Pythia.h"

#include<iostream>
#include<fstream>
#include<iomanip>
#include<string>
#include <vector>

#include<cstdio>
#include<cstdlib>
#include<cmath>
#include<ctime>
#include<cstdlib>

//#include"TH1.h"
//#include"TH2.h"
//#include"TH3.h"
//#include"TFile.h"


using namespace std;
//using namespace Pythia8;

//#include "fastjet/ClusterSequence.hh"

//using namespace fastjet;

//root for interactive graphics.
//#include "TVirtualPad.h"
//#include "TApplication.h"
//#include "TSystem.h"

//...constant

const int NPartons = 10000;

// for heavy quark radiation table
const int HQener_gn=400;
const int t_gn=75;
const int temp_gn=100;

const double HQener_max=1000.0;
const double t_min=0.0;
const double t_max=15.0;
const double temp_max=0.65;
const double temp_min=0.15;

// for MC initialization of jet partons
const double ipTmin=0.0;
const double ipTmax=70.0;
const double eta_cut=0.5;
const int maxMC=2000000;


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////...HQupdate

const int N_p1=100;
const int N_T=50;
const int N_e2=75;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////...HQupdate


//xx
//25-12-14
const int N_p1XX = 500;
const int N_TXX = 60;
//xx

const int N_p1XXX= 1000; //??-1-12
const int N_TXXX = 60; //??-1-12




//  const double  pi=3.1415926;  //conflict with fastjet  ???
const double    CA=3.0; 
const double    sctr=0.197;	     // 1 GeV^-1 = 0.197 fm
const double    pi=3.1415926; 
//...constant
//const double    pi=3.1415926;
//const double    CA=3.0;
const double    CF=4.0/3.0;
//const double    sctr=0.1973;	     //fm to GeV^-1
const double    hydro_Tc=0.165;
const double    KPamp=2.0;
const double    KPsig=5.0;
const double    KTamp=2.0;
const double    KTsig=0.05*hydro_Tc;
const double    preKT=1.0;












class LBTclass{


public:




char bulk3D_route[1024];
char parton_route[1024];





//...input with current machine time
//...random number seed (any negative integer)
	  
   //long  NUM1 = -33;
                  
long  NUM1;  

//    scattering rate
double Rg[50][100];         //total gluon scattering rate as functions of initial energy and temperature 
double Rg1[50][100];        //gg-gg              CT1
double Rg2[50][100];        //gg-qqbar           CT2
double Rg3[50][100];        //gq-qg              CT3
double Rq[50][100];         //total gluon scattering rate as functions of initial energy and temperature
double Rq3[50][100];        //qg-qg              CT13
double Rq4[50][100];        //qiqj-qiqj          CT4
double Rq5[50][100];        //qiqi-qiqi          CT5
double Rq6[50][100];        //qiqibar-qjqjbar    CT6
double Rq7[50][100];        //qiqibar-qiqibar    CT7
double Rq8[50][100];        //qqbar-gg           CT8

double qhatLQ[50][100];
double qhatG[50][100];
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////...HQupdate



double RHQ[50][100];        //total scattering rate for heavy quark
double RHQ11[50][100];      //Qq->Qq
double RHQ12[50][100];      //Qg->Qg
double qhatHQ[50][100];     //qhat of heavy quark
double qhat_over_T3;        //qhat/T^3 for heavy quark as fnc of (T,p)
double D2piT;

// for heavy quark radiation table

double dNg_over_dt_c[75+2][100+1][400+1];
double dNg_over_dt_q[75+2][100+1][400+1];
double dNg_over_dt_g[75+2][100+1][400+1];
double max_dNgfnc_c[75+2][100+1][400+1];
double max_dNgfnc_q[75+2][100+1][400+1];
double max_dNgfnc_g[75+2][100+1][400+1];


double delta_tg;
double delta_temp;
double delta_HQener;

// for MC initialization of jet partons


double initMCX[2000000],initMCY[2000000];




///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////...HQupdate













//....Hydro profile
double temphydro[50][80][80][40];
double edhydro[50][80][80][40];
double fractionhydro[50][80][80][40];
double VXhydro[50][80][80][40];
double VYhydro[50][80][80][40];
double VZhydro[50][80][80][40];


double numjet[1000000]={0.0};
double Xnucleon[1000000]={0.0};
double Ynucleon[1000000]={0.0};
double numjettotal=0.0;

////////////////////////////////////////////////////////////


double dxh_r;
double dyh_r;
double detah_r;  
double dtauh_r;

double x0h_r;
double y0h_r;
double eta0h_r;  
double tau0h_r;

//....geometry profile
//double numjet_r[10001];
//double Xnucleon_r[10001];
//double Ynucleon_r[10001];

double numjettotal_r;
int numnucleon_r;

char readinhydro_route[1024];


////////////////////////////////////////////////////////////



//...input information



int switchsingle;
int switchmedium;

double temp0medium;
double VXmedium;
double VYmedium;
double VZmedium;

double mass;
double px0;
double py0;
double pz0;
double en0;

double Xinitial;
double Yinitial;
double Zinitial;






int    ncall;
int    nprint;		

int    nrank;                    //NOT NECESSARY

double dt;                       //dtime when tauswitch is turned off (0) / dtau when tauswitch is turned on (1)
double ti0;                      //initial time
double timend;                   //end time or tau RENAME
double al;                       //2 dimensional plots range NOT NECESSARY

double ener;                     //initial energy of the jet parton/heavy quark
double amss;                     //initial mass of the jet parton/heavy quark 
double temp;                    
double temp0;                  
double temp00; 
double Ecut;                     //energy cut of the recoiled partons
		
int    Kjet;                     //initial flavor of the jet parton
int	   Kqhat0;                   //Debye screening mass switch RENAME
int    Kalphas;                  //alphas switch		
int    KINT;                     //radiation switch
int    KINT0;                     //radiation switch
int    Kradiation;
int    Kprimary;
int    tauswitch;                //coordinate switch
int    corswitch;

int    LBTswitch;                //switch of the Linear Boltzmann Transport

int    switchCoLBT_Hydro;
		
//...time-tau 				
double dtau;		
double tau0;		
double tauend;		
double time0;

//...variables with switch		
double alphas;
double qhat0;                    //Debye mass RENAME
double qhat00;		

int np;                          //number of partons at this time step ti				
		
//...radiation block		
int    icl22;                   
int    icl23;                    //the numerical switch in colljet23
int    iclrad;                   //the numerical switch in radiation 	  
int    isp;                      //the splitting function switch
int    nj;                       //number of leading shower partons 

//...Hydroprofile
double tauhydroend;
double tauhydrofile;
double tauhydro0;
double dtauhydro;
int ntauhydro;
int nxhydro;
int nyhydro;
int netahydro;

//...global variable qhat	
int counth100;
  
//double qhat;	                 //transport parameter

double qhat;	                 //transport parameter

double dng0[101][101];	 //table of dn/dkperp2/dx 	  

double Vtemp[4];

//...time system
double tirad[NPartons];
double tiscatter[NPartons];
double tiform[NPartons];	 //pythia 


double Tint_lrf[NPartons];          //for heavy quark
double eGluon;
double nGluon;








	  
//....radiated gluon
double radng[NPartons];	  
double xwm[3];
double wkt2m;
	  
double vf[4];              //flow velocity in tau-eta coordinate
double vfcar[4];           //flow velocity in t-z coordinate

double vp[4];              //position of particle		
double vc0[4];             //flow velocity     

//...dimensions in subroutine colljet and twcoll
double vc[4];              
double pc0[4];
double pc2[4];
double pc3[4];
double pc4[4];		
double p0[4];
double p2[4];
double p3[4];		
double p4[4];

double pc00[4];
double pc30[4];	
		
double pc01[4];
double pb[4];
		
double Pj0[4];
		
double V[4][NPartons];           //parton position 
double P[4][NPartons];           //parton 4-momentum
double V0[4][NPartons];          //negative parton position
double P0[4][NPartons];          //negative parton 4-momentum

double Prad[4][NPartons];		
		
int NR[NPartons];                  //scattering rank 
int KATT1[NPartons];               //parton flavor  
int KATT10[NPartons];              //negative parton flavor

double PGm[4];
double tjp[NPartons];
double Vfrozen[4][NPartons];     //parton final 4 coordinate
double Vfrozen0[4][NPartons];    //negative parton final 4 coordinate 
double Ecmcut;                   //energy cut for free streaming   

//........................................................................................................NCUT

double VV[4][NPartons];
double VV0[4][NPartons];
double PP[4][NPartons];
double PP0[4][NPartons];
int CAT[NPartons];
int CAT0[NPartons];
		
int ncut;
int ncut0;
int Reachtauend;

//........................................................................................................NCUTEND

//////////////
//...test
        int n_coll22;
		
        //int n_coll22=0;		
        int n_coll23;
        int ng_coll23;
        int ng_nrad;
        int n_radiation;
        int ng_radiation;
        int n_gluon;
        int n_sp1;
        int n_sp2;
//      ofstream datEg("./Eg.dat");
//////////////	  


        int ntest22;
        int ntestrad;


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////...HQupdate

// Variables for HQ 2->2
double min_p1,min_T,min_e2,max_p1,max_T,max_e2,bin_p1,bin_T,bin_e2;

double distFncB[50][100][75],distFncF[50][100][75],distMaxB[50][100][75],distMaxF[50][100][75];
double distFncBM[50][100],distFncFM[50][100];

//int loopN=10000;
int loopN;


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////...HQupdate



/////////////////////////////////////////////////////////////////////////////////////////////////////////////
double WT[NPartons];
double WT0[NPartons];
/////////////////////////////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////////////////////////////////

double vcfrozen[4][NPartons];
double vcfrozen0[4][NPartons];


double Tfrozen[NPartons];
double Tfrozen0[NPartons];


/////////////////////////////////////////////////////////////////////////////////////////////////////////////

double tempHydro[NPartons];
double vxHydro[NPartons];
double vyHydro[NPartons];
double vzHydro[NPartons];
double fractionHydro[NPartons];








//xx-12-14
double striRTE;
double runKT;
double preKTXX;
double fixAlphas = 0.3;
double KPampXX = 2.0;
double KTampXX = 2.0;
void lamXX(int KATT0, double& RTE, double& striRTE, double E, double T, double& T1, double& T2, double& E1, double& E2, int& iT1, int& iT2, int& iE1, int& iE2);
double nHQgluonXX(int parID, double dtLRF, double& time_gluon, double& temp_med, double& HQenergy, double& max_Ng);
double striRHQ[60][20];
double striRHQ11[60][20];
double striRHQ12[60][20];
double striqhatHQ[60][20];
double striRTE1;
double striRTE2;
double striqhatTP;
double scaleSigma = 4.0;
void flavorXX(int& CT, int& KATT0, int& KATT2, int& KATT3, double RTE, double E, double T, double& T1, double& T2, double& E1, double& E2, int& iT1, int& iT2, int& iE1, int& iE2);
void linearXX(int KATT, double E, double T, double& T1, double& T2, double& E1, double& E2, int& iT1, int& iT2, int& iE1, int& iE2, double& RTEg, double& RTEg1, double& RTEg2, double& RTEg3, double& RTEq, double& RTEq3, double& RTEq4, double& RTEq5, double& RTEq6, double& RTEq7, double& RTEq8, double& RTEHQ, double& RTEHQ11, double& RTEHQ12, double& qhatTP);
void collHQ22Qq(int CT, double temp, double qhat, double v0[4], double p0[4], double p2[4], double p3[4], double p4[4], double& qt);
void collHQ22Qg(int CT, double temp, double qhat, double v0[4], double p0[4], double p2[4], double p3[4], double p4[4], double& qt);
const double a_md = 0.1;
const double b_md = 2.0;
double runalphas = 1.0;
double alphafactor = 1.0;
double stridistFncBM[N_TXX][N_p1XX], stridistFncFM[N_TXX][N_p1XX];
double stridistFncB[N_TXX][N_p1XX][N_e2], stridistFncF[N_TXX][N_p1XX][N_e2], stridistMaxB[N_TXX][N_p1XX][N_e2], stridistMaxF[N_TXX][N_p1XX][N_e2];
double Mqc2qc_prime(double s, double t, double M, double temp);
double Mgc2gc_prime(double s, double t, double M, double temp);
double msprime;
double mdprime;
const double as_ms = 0.0;
const double bs_ms = 0.1;
const double alphasprime = 0.29;
const double sigmaprime = 0.225;
const double    Nc = 3.0;
void collHQ23XX(int parID, double temp_med, double qhat0ud, double v0[4], double p0[4], double p2[4], double p3[4], double p4[4], double qt, int& ic, double Tdiff, double HQenergy, double max_Ng, double xLow, double xInt);
void radiationHQXX(int parID, double qhat0ud, double v0[4], double P2[4], double P3[4], double P4[4], double Pj0[4], int& ic, double Tdiff, double HQenergy, double max_Ng, double temp_med, double xLow, double xInt);
const double HQener_maxXX = 1000.0;
static const int HQener_gnXX = 500;
const double temp_maxXX = 0.65;
static const int temp_gnXX = 100;
double qhatGXX[60][20];
double qhatLQXX[60][20];
double RgXX[60][20];
double Rg1XX[60][20];
double Rg2XX[60][20];
double Rg3XX[60][20];
double RqXX[60][20];
double Rq3XX[60][20];
double Rq4XX[60][20];
double Rq5XX[60][20];
double Rq6XX[60][20];
double Rq7XX[60][20];
double Rq8XX[60][20];
double RHQXX[60][20];
double RHQ11XX[60][20];
double RHQ12XX[60][20];
double qhatHQXX[60][20];
double dNg_over_dt_cXX[t_gn + 2][temp_gnXX + 1][HQener_gnXX + 1] = { 0.0 };
double dNg_over_dt_qXX[t_gn + 2][temp_gnXX + 1][HQener_gnXX + 1] = { 0.0 };
double dNg_over_dt_gXX[t_gn + 2][temp_gnXX + 1][HQener_gnXX + 1] = { 0.0 };
double max_dNgfnc_cXX[t_gn + 2][temp_gnXX + 1][HQener_gnXX + 1] = { 0.0 };
double max_dNgfnc_qXX[t_gn + 2][temp_gnXX + 1][HQener_gnXX + 1] = { 0.0 };
double max_dNgfnc_gXX[t_gn + 2][temp_gnXX + 1][HQener_gnXX + 1] = { 0.0 };
double min_p1XX = 0.0;
double max_p1XX = 1000.0;
double bin_p1XX = (max_p1XX - min_p1XX) / N_p1XX;
double min_TXX = 0.1;
double max_TXX = 0.7;
double bin_TXX = (max_TXX - min_TXX) / N_TXX;
double distFncBMXX[N_T][N_p1], distFncFMXX[N_T][N_p1];
double distFncBXX[N_TXX][N_p1XX][N_e2], distFncFXX[N_TXX][N_p1XX][N_e2], distMaxBXX[N_TXX][N_p1XX][N_e2], distMaxFXX[N_TXX][N_p1XX][N_e2];
double delta_tempXX;
double delta_HQenerXX;
const double    epsilonXX = 1e-6;
double EjpXX;
double DebyeMass2XX(int& Kqhat0, double alphas, double temp0);
double ElabXX;
void bulklinearXX(double tau, double x, double y, double eta, double& temp, double& VX, double& VY, double& VZ, int& fraction);








//??-1-11
static const int dimParList = 2000000;
int color[dimParList] = { 0 };
int anticolor[dimParList] = { 0 };
int color0[dimParList] = { 0 };
int anticolor0[dimParList] = { 0 };
int fake[dimParList] = { 0 };
int fake0[dimParList] = { 0 };

//??-1-12
int maxColor ;
double RHQb[60][20];        
double RHQ11b[60][20]; 
double RHQ12b[60][20];   
double qhatHQb[60][20];
double dNg_over_dt_b[75 + 2][100 + 1][1000 + 1] = { 0.0 };
double max_dNgfnc_b[75 + 2][100 + 1][1000 + 1] = { 0.0 };
double distFncBMb[60][1000], distFncFMb[60][1000];
double distFncBb[60][1000][75],distFncFb[60][1000][75],distMaxBb[60][1000][75],distMaxFb[60][1000][75];
double KPampb = 5.0;
double KPsigb = 5.0;
double KTampb = 0.0;
double KTsigb = 0.05;
const double HQener_maxXXX = 2000.0;
static const int HQener_gnXXX = 1000;
double delta_HQenerXXX;
double min_p1XXX = 0.0;
double max_p1XXX = 2000.0;
double bin_p1XXX = (max_p1XXX - min_p1XXX) / N_p1XXX;
double min_TXXX = 0.1;
double max_TXXX = 0.7;
double bin_TXXX = (max_TXXX - min_TXXX) / N_TXXX;
void collHQ22XXX(int CT, double temp, double qhat, double v0[4], double p0[4], double p2[4], double p3[4], double p4[4], double& qt, int parID);






//...functions	  

void LinearBoltzmannTransport(int &n, double &ti);

void trans(double v[4],double p[4]);
void transback(double v[4],double p[4]);
  
float ran0(long *idum);
   
double alphas0(int &Kalphas,double temp0);
double DebyeMass2(int &Kqhat0,double alphas,double temp0);
	  
void lam(int KATT0,double &RTE,double E,double T,double &T1,double &T2,double &E1,double &E2,int &iT1,int &iT2,int &iE1,int &iE2);	  

void flavor(int &CT,int &KATT0,int &KATT2,int &KATT3,double RTE,double E,double T,double &T1,double &T2,double &E1,double &E2,int &iT1,int &iT2,int &iE1,int &iE2);
void linear(int KATT,double E,double T,double &T1,double &T2,double &E1,double &E2,int &iT1,int &iT2,int &iE1,int &iE2,double &RTEg,double &RTEg1,double &RTEg2,double &RTEg3,double &RTEq,double &RTEq3,double &RTEq4,double &RTEq5,double &RTEq6,double &RTEq7,double &RTEq8);
void twflavor(int &CT,int &KATT0,int &KATT2,double E,double T);
void twlinear(int KATT,double E,double T,double &RTEg1,double &RTEg2,double &RTEq6,double &RTEq7,double &RTEq8);
void colljet22(int CT,double temp,double qhat0ud,double v0[4],double p0[4],double p2[4],double p3[4],double p4[4],double &qt);
void twcoll(int CT,double qhat0ud,double v0[4],double p0[4],double p2[4]);
	  
void titau(double ti,double vf[4],double vp[4],double p0[4],double &Vx,double &Vy,double &Veta,double &Xtau);


//void colljet23( double temp, double qhat0ud, double v0[4], double p0[4],double p2[4], double p3[4], double p4[4], double qt, int &ic, double tint, double tisint, double tirint, double dtlastrad, double Elabint, double Ejpint);
//void radiation(double qhat0ud,double v0[4],double P2[4],double P3[4],double P4[4],double Pj0[4], int &ic, double tint, double tisint, double tirint, double dtlastrad, double Elabint, double Ejpint);

void colljet23( double temp, double qhat0ud, double v0[4], double p0[4],double p2[4], double p3[4], double p4[4], double qt, int &ic, double tint, double tisint, double tirint, double dtlastrad, double Elabint, double Ejpint);
void radiation(double qhat0ud,double v0[4],double P2[4],double P3[4],double P4[4],double Pj0[4], int &ic, double tint, double tisint, double tirint, double dtlastrad, double Elabint, double Ejpint);


void rotate(double px,double py,double pz,double p[4],int icc);
void dngintt(double dtint,double tint,double tisint,double tirint,double dtlastrad,double qhat0int,double Elabint,double Ejpint);
double angluon(double dtint,double tint,double tisint,double tirint,double qhat0int,double Elabint,double Ejpint,double dtlrf,double dtlastrad);
double dngrad(double xw0, double wkt20, double Elab0, double Ejp0, double tint0, double tis0, double tir0, double qhat0int0, double dtlastrad);
double dng(double x,double y);
	  
int KPoisson(double alambda);	  

void bulklinear(double tau, double x,double y,double eta,double &ed, double &temp,double &fraction,double &VX,double &VY,double &VZ);






//////////////////////////////////////////////////////////////////////////////////////////////////...NEW



void radiationHQ(int parID, double qhat0ud, double v0[4], double P2[4], double P3[4], double P4[4], double Pj0[4], int &ic, double Tdiff, double HQenergy, double max_Ng, double temp_med, double xLow, double xInt);  
void collHQ22(int CT,double temp,double qhat,double v0[4],double p0[4],double p2[4],double p3[4],double p4[4],double &qt);
double Mqc2qc(double s, double t, double M);
double Mgc2gc(double s, double t, double M);

void collHQ23(int parID, double temp_med, double qhat0ud, double v0[4], double p0[4], double p2[4], double p3[4], double p4[4], double qt, int &ic, double Tdiff, double HQenergy, double max_Ng, double xLow, double xInt);
double dNg_over_dxdydt(int parID, double x0g, double y0g, double HQenergy, double HQmass, double temp_med, double Tdiff);
double tau_f(double x0g, double y0g, double HQenergy, double HQmass);
double splittingP(int parID, double z0g);
double lambdas(double kTFnc);
double nflavor(double kTFnc);
double alphasHQ(double kTFnc, double tempFnc);
double nHQgluon(int parID,double dtLRF,double &time_gluon,double &temp_med,double &HQenergy,double &max_Ng);

void read_xyMC(int &numXY);
void jetInitialize(int numXY);

//////////////////////////////////////////////////////////////////////////////////////////////////...NEW



void LBTinitialize();

void LBTclear();


void rotatePsi(double Psi, double P1[3], double P2[3]);
double psi_std(double px, double py);



};











