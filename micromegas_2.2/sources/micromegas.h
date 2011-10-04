#ifndef  __MICROMEGAS__
#define  __MICROMEGAS__

#ifdef __cplusplus
extern "C" {
#endif 

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<unistd.h>
#include"../CalcHEP_src/include/num_out.h"
#include"../CalcHEP_src/c_source/model_aux/include/SLHAreader.h"
typedef struct numout
{
  void * handle;
  double ** link;
  double *Q,*SC,*GG;
  int init;
  CalcHEP_interface * interface; 
} numout;

extern int sortOddParticles(char * name);
extern double sWidth;

/*============================
     Particles
==============================*/

extern long pNum(char * name);  /* returns PDG code */
extern double pMass(char * name); /* returns particle mass */

/*======= Subprocesses ===========*/
  typedef struct txtListStr
  {  struct txtListStr * next;
     char  *txt;
  } txtListStr;

  typedef txtListStr * txtList;
  
extern txtList  makeDecayList(char * pname, int nx);
extern void massFilter(double M, txtList * List);
extern void gammaGluFilter(txtList* List);
extern void process2Lib(char * process,char * lib);
extern void cleanTxtList(txtList L);
extern double findBr(txtList L, char * pattern);
extern void printTxtList(txtList L, FILE *f);

/*=============================
     (1)2->2 matrix elements
=============================*/
extern double decay2Info(char * pname, FILE* f);
extern numout *  newProcess(char*Process,char* lib);
extern double cs22(numout * cc, int nsub,double P,  double cos1, double cos2 , int * err);
extern double pWidth2(numout * cc,int nsub);
extern int  procInfo1(numout*cc, int *nsub, int * nin, int *nout);
extern int procInfo2(numout*cc,int nsub,char**name,double*mass);
extern double Helicity[2];
extern double hCollider(double Pcm, int pp, char * name1,char *name2);
double pWidth(char *name, txtList *L,int *dim);


/*===================
      Variables 
=====================*/

extern double*varAddress(char * name); 
extern int    assignVal(char * name, double val);
extern int    findVal(char * name,double * val);
extern void   assignValW(char * name, double val);
extern double findValW(char * name);
extern int    readVar(char *fname);
  
/*===========================================
   Checking of parameters 
  ===========================================*/ 
extern void printVar(FILE *f);
extern void  printMasses(FILE * f, int sort);
extern double lopmass_(void);

/*=====================================
    Relic density evaluation 
  =====================================*/ 

extern double darkOmega(double *Xf,int Fast, double Beps);
extern double darkOmegaFO(double *Xf,int fast,double Beps);
extern double printChannels(double Xf,double cut,double Beps,int prcnt,FILE *f );   
extern void improveCrossSection(long n1,long n2,long n3,long n4,double Pcm, 
                                                            double * addr);
extern char * wimpAnnLib(void);
/*===============================================
    Annihilation spectra
=================================================*/
                                                                                
extern double calcSpectrum(double v ,int outP,double *tab,int*errcode);
/* input parameters:
       a)  v velocity in c=1 units, v is about 0.001.
       b)  outP   0-gamma; 1-positron; 2-antiproton; 3,4,5 -
         for neutrinos (electron, muon and tau correspondinly)
output parameters:
       a) returned value v*sigma/M_lop^2 in sm^3/(GeV^2*sec)
       b) array of 250 double elements to store spectra.
       c) errcode !=0 if  in process of decays of outgoing particles 
          we get a non-SM particle which has not 1->2 decays.  
*/
extern double zInterp(double x, double * tab);
/* input parameter:
        a) x=log(E/Nmass)
        b) double tab[250] table obtainaed by calcSpectrum;
   output:
         dN/dx where N is the number spectrum.
*/
extern void spectrInfo(double Xmin,double*tab, double * Ntot,double*Etot);
/* input parameters:
  a)    Xmin = Emin/mLSP - characterizes minimal energy under considerstion
  b)    double tab[250]  - table for spectrum distribution.
  output:
  a)    Ntot - total number of particles with energy E> Xmin*mLsp
  b)    Everage energy of particle  divided on mLsp
*/
                                                                                
extern int spectrTable(double*tab, char*fname,char*mess,double Xmin,int N);
/* writes into  the file fname a table of  dN/d(log_10(E)) distribution.
   input parameters:
   a) tab - a table for distrubution
   b)  fname - a file name
   c) mess - a title
   d) Emin/mLsp
   f) N<=300 -a number of points to write
                                                                                
   The file can be displayed  by
     ./sources/CalcHEP_src/bin/tab_view <fname
*/

extern int displaySpectrum(double*tab, char*mess,double Xmin);

extern double HaloFactor(double fi,double (*rhoQ)(double));
/*
  Intergates halo squred density along the ling of slight.
  fi - angle in radians
  rhoQ  presents squred density of Dark Matter in [GeV/cm^3]^2 units as 
  a function of cerner galactic distance in kpc units.
  
  Return value is done in  GeV^2/cm^5   
*/ 

extern double rhoQisotermic(double x);
/* 
   presents  an example of rhoQ for HaloFactor.
   It corresponds to modified isotermic density with absence of clumbs.
*/ 
                                                                                
extern char * outNames[6];

extern void mInterp(double Nmass,  int  CHin,int  CHout, double*tab);

/*  Direct Detection */
extern int QCDcorrections, Twist2On;
extern void setProtonFF(double *scalar, double *ps_vector, double * sigma);
extern void setNeutronFF(double *scalar, double *ps_vector, double * sigma);
extern void getScalarFF(double muDmd,double msDmd,double sigmaPiN,double sigma0,
                        double * protonFF, double * neutronFF);

extern double FeScLoop(double sgn, double mq,double msq,double mne);
extern int nucleonAmplitudes( double (*LF)(double,double,double,double), 
                  double*pA0,double*pA5,double*nA0,double*nA5); 

extern void  SetFermi(double C,double B, double a);
extern double FermiFF(int A, double Qfermi);

extern void SetfMaxwell(
double DV,    
double v1,    /* velocity of Earth respect to WINP  in km/s units     */
double vmax
);

extern double fDvMaxwell(double v);

extern void SetfDelta(double V0);
extern double fDvDelta(double v);

extern void setRecoilEnergyGrid(double step, int dim);

extern double nucleusRecoil(
double rho, /* WINP density near the Earth in GeV/cm^3 units */
double(*fDv)(double),   /*  f(v)/v where f(v) is velocity distribution, 
                            v in km/s  */
int A, int Z, double J, /* nuclear atomic number, charge and spin   */

      double(*S00)(double),double(*S01)(double),double(*S11)(double),
double(*LF)(double,double,double,double),     
double * dNdE   /* distribution of number of events respect to recoil 
                 energy  in 1/(Kg*Day*Kev) inits. Presinted as array 
                 with 200 elements with 1KeV step (from 0 to 199KeV) */
);
/* nucleusRecoil returns the total number of events for 1day*1Kg */
                  
extern double nucleusRecoil0(double rho, double (*vfv)(double),
 int A,int Z,double J,double Sp,double Sn,
 double(*LF)(double,double,double,double), double*dNdE);
 
extern double nucleusRecoilAux(
      double rho, double(*vfv)(double),
      int A, int Z, double J,
      double(*S00)(double),double(*S01)(double),double(*S11)(double),
      double Mwimp,
      double cs_SI_P,double cs_SI_N,  double cs_SD_P, double cs_SD_N,
      double * dNdE);

extern double nucleusRecoil0Aux(
      double rho, double(*vfv)(double),
      int A, int Z, double J,
      double Sp, double Sn,
      double Mwimp,
      double cs_SI_P,double cs_SI_N,  double cs_SD_P, double cs_SD_N,
      double * dNdE);

extern double MaxGapLim(double x, double mu);
/* S.Yellin, Phys.Rev. D66,032005(2002)                                                                                                                        
  Lut us  a theoretical model  predicts homogenious event distribution
  with everage number of events mu. Let experiment gets a gap bitween points
  where according to theory x point are expected. Then the theoretical model
  is non-confirmed with probability MaxGap
*/
                                                                                                                                    

extern int displayRecoilPlot(double * tab, char * text, double E1, double E2);

extern double cutRecoilResult(double *tab, double E1, double E2);

typedef double (double2double)(double);

extern  double2double     S00F19,S11F19,S01F19,       S00Na23,S11Na23,S01Na23,
S00Al27,S11Al27,S01Al27,  S00Si29,S11Si29,S01Si29,    S00K39,S11K39,S01K39,
S00Ge73,S11Ge73,S01Ge73,  S00Nb93,S11Nb93,S01Nb93,    S00Te125,S11Te125,S01Te125,
S00I127,S11I127,S01I127,  S00Xe129,S11Xe129,S01Xe129, S00Xe131,S11Xe131,S01Xe131,
S00Pb207,S11Pb207,S01Pb207;

extern  double2double     S00Na23A,S11Na23A,S01Na23A, S00Si29A,S11Si29A,S01Si29A, 
S00Te125A,S11Te125A,S01Te125A,   S00I127A,S11I127A,S01I127A,
S00Xe129A,S11Xe129A,S01Xe129A,   S00Xe131A,S11Xe131A,S01Xe131A,
S00Ge73A,S11Ge73A,S01Ge73A;
extern  double2double S00Xe131B,S11Xe131B,S01Xe131B;

/* for testing SD form factors */
extern  int PlotSS(double (*f)(double),  int A,char * title, double Emax);
extern  int Plot3SS(double xiP,double xiN,
         double (*S00)(double), double (*S01)(double),double (*S11)(double),
         int A, double J, char * title, double Emax);
extern  int Plot3SS0(double xiP,double xiN, double Sp, double Sn,
    int A, double J, char* title, double Emax);
     
      
#define Sp_H1     ( 0.5)
#define Sn_H1       0.
#define Sp_He3    (-0.081)
#define Sn_He3      0.552
#define Sp_F19    ( 0.4751)
#define Sn_F19    (-0.0087)
#define Sp_Na23   ( 0.2477)
#define Sn_Na23   ( 0.0198)
#define Sp_Te125  ( 0.001) 
#define Sn_Te125  ( 0.287)
#define Sp_I127     0.309 
#define Sn_I127     0.075
#define Sp_Xe129    0.028  
#define Sn_Xe129    0.359
#define Sp_Xe131  (-0.009)  
#define Sn_Xe131  (-0.227)
#define Sp_Al27   ( 0.343)
#define Sn_Al27     0.0296
#define Sp_Si29   (-0.0019)
#define Sn_Si29   ( 0.1334)
#define Sp_K39    (-0.184) 
#define Sn_K39    ( 0.054)
#define Sp_Ge73   ( 0.03)
#define Sn_Ge73   ( 0.378)
#define Sp_Nb93   ( 0.46)
#define Sn_Nb93   ( 0.08) 
#define Sp_Cs133  (-0.370)   /* Phis.Lett.B254,220,(1991)*/
#define Sn_Cs133  ( 0.003)  
#define Sp_Pb207  (-0.010)  
#define Sn_Pb207  (-0.149)


/*
http://www.nndc.bnl.gov/nudat2/indx_sigma.jsp
*/
#define J_H1    0.5
#define J_He3   0.5
#define J_F19   0.5
#define J_Na23  1.5
#define J_Al27  2.5
#define J_Si29  0.5
#define J_K39   1.5
#define J_Ge73  4.5
#define J_Nb93  4.5
#define J_Te125 0.5 
#define J_I127  2.5
#define J_Xe129 0.5
#define J_Xe131 1.5
#define J_Cs133 3.5
#define J_Pb207 0.5

#define Z_H    1
#define Z_He   2
#define Z_F    9
#define Z_Na  11
#define Z_Al  13
#define Z_Si  14
#define Z_K   19
#define Z_Ge  32
#define Z_Nb  41
#define Z_Te  52
#define Z_I   53
#define Z_Xe  54
#define Z_Pb  82
#define Z_Cs  55

#ifdef __cplusplus
}
#endif 

#endif
