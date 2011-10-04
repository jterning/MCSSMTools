#ifndef  __MICRO_FORT__
#define  __MICRO_FORT__


extern int sortoddparticles_(char * name,int len);
/* 
    INTEGER FUNCTION sortOddParticles(name)
    character *(*) name
*/


/*=============================
     MSSM 2->2 cross sections
=============================*/

extern void newprocess_(char*Process,char* lib, int* address,int len1,int len2);
/* 
   Subroutine newProcess(Process, lib, address)
   character *(*) Process
   character *(*) lib
   integer address(2)
*/

extern double cs22_(int * address, int * nsub,double * P,  double * cos1, 
   double * cos2 , int * err);
/* 
   REAL*8 FUNCTION cs22( address,nsub,P,cos1,cos2,err) 
   integer address(2),nsub,err
   real*8 P,cos1,cos2
*/   
      
extern double pwidth2_(int * address,int *nsub);
/* 
    REAL*8 FUNCTION pWidth2(address, nsub)
    integer address(2),nsub
*/

extern void  procinfo1_(int * address, int * ntot, int *nin, int *nout);
/*   
    subroutine procInfo1(address, ntot,nin,nout)
    integer address(2), ntot, nin,nout
*/      
 
extern void procinfo2_(int * address ,int *nsub,char*names,double*masses,int len);
/*  
    subroutine  procInfo2(address, names, masses)
    integer address(2)
    character *(*) names
    real*8 (*) masses
*/      

extern  double decay2info_(char * pname, int *N, int len);
/*   
    real*8 function decay2Info(pname, Nfile)
    character*(*) pname
    integer Nfile 
*/

extern void sethelicities(double *h1, double *h2);
/*
    subroutine setHelicities(h1,h2)
    real*8 h1,h2 
*/
/*===================
      Variables 
=====================*/


extern int  assignval_(char * name, double * val, int len);
/*  
      integer function assignVal(name, val)
      character *(*) name
      real*8 val
*/    
   

extern int findval_(char * name,double * val,int len);
/*  
     integer function findVal(name,val)
     character *(*) name
     real*8 val
*/ 

extern void   assignvalw_(char * name, double * val, int len);
/*
     subroutine assignValW(name,val)
     character *(*) name
     real*8 val
*/
            
extern double findvalw_(char * name, int len);
/*   
     real*8 function findValW(name)
     character *(*) name
*/


extern int readvar_(char *fname, int len);
/*
     integer  function readVar(fname)
     character *(*) fname
*/

  
/*===========================================
   Checking of parameters 
  ===========================================*/ 
extern void printvar_(int * N);
/*
     subroutine printVar(N)
     integer N
*/
     
extern void  printmasses_(int *N, int*sort);
/*
     subroutine printMasses(N, sort)
     integer N,sort

*/

extern double lopmass_(void);

/*=====================================
    Relic density evaluation 
  =====================================*/ 

extern double darkomega_(double *Xf,int * Fast, double* Beps);
/*
     real*8 function darkOmega(Xf,Fast,Beps)
     real*8 Xf,Fast,Beps
*/


extern double printchannels_(double *Xf,double*cut,double* Beps,int *prcnt,int*N );
/*
     real*8 function printChannels(Xf,cut,Beps,prcn,N)
     real*8 Xf,cut,Beps
     integer prcn,N
*/
/*===============================================
    Annihilation spectra
=================================================*/
extern double calcspectrum_(double *v,int*outP,double *tab,int * err);
/* 
   real*8 function calcSpectrum(v,outP,tab,err)
   real*8 v
   integer outP
   interger err
   real*8 tab(250)
*/
     
extern void spectrinfo_(double*Xmin,double*tab,double*Ntot,double*Etot);
/* 
    subroutine spectrInfo(Xmin,tab,Ntot,Etot)
    real*8 Xmin,tab(250),Ntot,Etot
*/
int  spectrtable_(double*tab, char*fname,char*mess,double *Xmin,int *N,
                   int len1,int len2);
/*
    integer function  spectrTable(tab,fname,mess,Xmin,N)
    real*8 tab(250)
    character*(*) fname
    character*(*) mess
    real*8 Xmin
    integer N                   
*/

extern int displayspectrum_(double*tab, char*fmess,double *Xmin,int len);
/*  integer function displaySpectrum( tab, mess, Xmin)
    real*8 tab(250)
    character*(*) mess
    real*8 Xmin
*/
extern double halofactor_(double *fi,double (*rhoQ)(double*));
/*
  real*8 function HaloFactor( fi, rhoQ)
  real*8  fi, rho 
  external rhoQ
*/ 

extern double rhoqisotermic_(double *x);
/* 
  real*8 function RhoQisotermic(x)
  real*8 x 

*/ 

extern double zinterp_(double*x, double*tab);
/*
   real*8 function zinterp_(x,tab)
   real*8 x,tab(250)
*/


extern double findparam_(char *name, int *err, int len);
/* 
   real*8 function findParam(name,err)
   character*(*) name
   integer err
*/   

extern double findparamw_(char * name, int len);
/*
    real*8 function findParamW(name)
        character*(*) name
*/
        

/*===============================
  Direct Detection 
===================================*/


extern void setprotonff_(double *scalar, double *ps_vector, double * sigma);
/*
    subroutine setProtonFF(scalar,ps_vector,sigma)
    real*8 scalar(3),ps_vector(3),sigma(3)
*/

extern  void setneutronff_(double *scalar, double *ps_vector, double * sigma);
/*
    subroutine setNeutronFF(scalar,ps_vector,sigma)
    real*8 scalar(3),ps_vector(3),sigma(3)
*/
extern void getscalarff_(double *muDmd,double *msDmd,double *sigmaPiN,double *sigma0,
                        double * protonFF, double * neutronFF);
/*
    subroutine getScalarFF(muDmd,msDmd,sigmaPiN,sigma0,protonFF,neutronFF);
    real*8 muDmd,msDmd,sigmaPiN,sigma0,protonFF(3),neutronFF(3)                        
*/


extern double noloop_(double*,double*,double*,double*);
/*  real *8 function NoLoop(sgn, mq,msq,mne)
    real*8 sgn, mq,msq,mne
*/     

extern int nucleonamplitudes_(double (*LF)(double*,double*,double*,double*),
                         double*pA0,double*pA5,double*nA0,double*nA5);
/*
    integer function  nucleonAmplitudes(LF,pA0,pA5,nA0,nA5)
    real*8 LF,pA0(2),pA5(2),nA0(2),nA5(2) 
*/

extern void setfermi_(double *C, double *B,double *a);
/*
     subroutine  SetFermi(C,B,a)
     real*8 C,B,a
*/ 
extern double  fermiff_(int *A, double * Qfermi);
/*
    real*8 function FermiFF(A,Qfermi)
    integer A
    real*8 Qfermi       
*/
extern void setfmaxwell_(double *DV,double *v1,double *vmax);
/*
    subroutine SetfMaxwell(DV,v1,vmax)
    real*8 DV,v1,vmax
*/

extern double fdvmaxwell_(double *v);
/*
       real*8 function fDvMaxwell(v)
       real*8  v

*/

extern void setfdelta_(double *v);
/*     subroutine SetfDelta(v)
       real*8 v
*/        

extern double fdvdelta_(double *v);
/*     real*8 function fDvDelta(v)
       real*8  v  
*/

extern double nucleusrecoil_(
     double *rho, double(*fDv)(double*),int*A, int*Z, double*J, 
     double(*S00)(double*),double(*S01)(double*),double(*S11)(double*),
     double (*LF)(double*,double*,double*,double*), double * dNdE );
/*
     real*8 function nucleusRecoil(rho, fDv,A,Z,J,S00,S01,S11,LF,dNdE )
     real*8 rho, fDv,J,S00,S01,S11,LF,dNdE(200)
     integer A,Z
*/     

extern double nucleusrecoil0_(double *rho, double (*fDv)(double*),
 int*A,int*Z,double*J,double*Sp,double*Sn,
 double (*LF)(double*,double*,double*,double*),double*dNdE);
/*
    real*8 nucleusRecoil0(rho, fDv,A,Z,J,Sp,Sn,LF,dNdE)
    real*8 rho, fDv,J,Sp,Sn,LF,dNdE(200)
    integer A,Z
*/

extern int displayrecoilplot_(double * tab, char * text, int *E1, int *E2,int len);
/*
    integer function displayRecoilPlot(tab, text, E1, E2)
    real*8  tab(*)
    character *(*) text
    integer E1,E2   
*/

extern double cutrecoilresult_(double *tab, int *E1, int *E2);
/*
     real*8 function cutRecoilResult(tab, E1, E2)
     dimension tab
     real*8 tab
     integer E1,E2

*/

typedef double (doubleFORTdouble)(double*);

extern doubleFORTdouble 
s00f19_   ,s11f19_   ,s01f19_   ,s00na23_  ,s11na23_  ,s01na23_  ,s00al27_  ,s11al27_  ,s01al27_,  
s00si29_  ,s11si29_  ,s01si29_  ,s00K39_   ,s11K39_   ,s01K39_   ,s00ge73_  ,s11ge73_  ,s01ge73_,  
s00nb93_  ,s11nb93_  ,s01nb93_  ,s00te125_ ,s11te125_ ,s01te125_ ,s00i127_  ,s11i127_  ,s01i127_,  
s00xe129_ ,s11xe129_ ,s01xe129_ ,s00xe131_ ,s11xe131_ ,s01xe131_ ,s00pb207_ ,s11pb207_ ,s01pb207_, 
s00na23a_ ,s11na23a_ ,s01na23a_ , s00si29a_,s11si29a_,s01si29a_ ,s00te125a_,s11te125a_,s01te125a_,
s00i127a_ ,s11i127a_ ,s01i127a_ ,s00xe129a_,s11xe129a_,s01xe129a_,s00xe131a_,s11xe131a_,s01xe131a_,
s00ge73a_ ,s11ge73a_ ,s01ge73a_ ,s00xe131b_,s11xe131b_,s01xe131b_;                                                                                

#endif
