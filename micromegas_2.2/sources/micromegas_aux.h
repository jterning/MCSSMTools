#ifndef __MICRO_AUX__
#define __MICRO_AUX__

#include"micromegas.h"
#include"../CalcHEP_src/include/V_and_P.h"

#ifdef __cplusplus
extern "C" {
#endif 

/*======================
   Lock/UnLock service
========================*/
extern int checkLockFile(int *delay);
extern void removeLock(int fd);

extern char*WORK;

/*=======================
      Odd particles
========================*/
extern int Nodd;
extern  ModelPrtclsStr * OddPrtcls;
extern int createTableOddPrtcls(void);
extern char * txtListOddParticles(void);     
/*=============================================
C->F and F->C    string convertation     
===============================================*/

extern void  cName2f(char*c_name,char*f_name,int len);
extern void  fName2c(char*f_name,char*c_name,int len);

/*=============================================
Fortran output
===============================================*/
extern void fortreread_(int* N, char * fname, int len);

/*====================================
   Tools for integration 
=====================================*/

extern double simpson(double (*func)(double), double a, double b, double eps);
extern double gauss(double (*func)(double), double a, double b, int n);
extern int  odeint(double*ystart,int nvar,double x1,double x2, double eps,
                   double h1, void(*derivs)(double,double*,double *));

/*==== Tool  for interpolation  ====*/
extern double  polint3(double x, int n,  double *xa, double *ya);

extern double K2pol(double x); /*bessk1(1/x)*exp(1/x)*sqrt(2/M_PI/x);*/
extern double K1pol(double x); /*bessk2(1/x)*exp(1/x)*sqrt(2/M_PI/x);*/

/* Hidden interface with CalcHEP */
extern int FError;
extern int OnlyTEQ0;

extern void pname2lib(char*pname, char * libname);
extern numout*newProcess_(int twidth, int model,char*Process,
          char * excludeVirtual,char*excludeOut,char*lib,int usr);

extern double decayPcm(double am0,  double  am1,  double  am2);
extern int  kin22(double PcmIn,double * pmass);
extern double  dSigma_dCos(double  cos_f);
extern double (*sqme)(int nsub,double *pvect, int * err_code);
extern int _nsub_;

#define NTOF(X) extern forCalchep1 X; double X(double ok){return findValW(#X);}

typedef  double (forCalchep1)(double);
typedef  double (forCalchep2)(double,double);

/*  Loop integrals I1 ... I5 */

extern double   LintIk(int I,double MSQ,double MQ,double MNE);

extern int readVarSpecial(char *fname, int nVar, char ** names);

extern double parton_x( int pNum, double  Q);
extern double parton_alpha(double q);
extern double GluWeight(void);
extern double parton_distr(int pNum, double x, double q);



extern double convStrFun2(double x, double q, int pc1, int pc2, int ppFlag ); 
                           /* result of convolution of structure functions 
                              of pc1 and pc2 particles  */

extern int  pTabPos(char *name);
extern int decodeProcess(char *txt,int*inList,int*outList);

typedef struct{ double width; int dim; txtList pdList[2];}  decayTableStr;

extern  decayTableStr* decayTable;
extern void cleanDecayTable(void);

extern double 
/*scalar */
   wS0P__[6],wS0N__[6],
/* pseudo-vertor */
   wV5P__[3],wV5N__[3];
#ifdef __cplusplus
}
#endif 


#endif
