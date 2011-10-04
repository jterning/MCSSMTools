/*
 Copyright (C) 1997,2006, Alexander Pukhov 
*/
#include <math.h>
#include"err_code.h"
#include"simpson.h"



static int  NMAX;

static void r_gauss( double(*func)(double),double a,double b, 
double eps, double * aEps, double * ans, double * aAns,int* N)
{
  int i,n;
double X3[3]={0.11270166537925831147, 0.50000000000000000000, 0.88729833462074168853};
double F3[3]={0.27777777777777777778, 0.44444444444444444444, 0.27777777777777777778};
double X4[4]={0.06943184420297371240, 0.33000947820757186760, 0.66999052179242813240, 0.93056815579702628757};
double F4[4]={0.17392742256872692869, 0.32607257743127307132, 0.32607257743127307132, 0.17392742256872692869};
double X5[5]={0.046910077030668003594,0.230765344947158454491,0.50000000000000000000, 0.769234655052841545509, 0.953089922969331996378};
double F5[5]={0.118463442528094543757,0.239314335249683234015,0.284444444444444444444,0.239314335249683234015, 0.118463442528094543756};

  double s1,s2,s3,e_err,d=b-a;

  for(n=0,s1=0;n<3;n++) s1+=F3[n]*func(a+ d*X3[n]); s1*=d; *N+=3;
  for(n=0,s2=0;n<4;n++) s2+=F4[n]*func(a+ d*X4[n]); s2*=d; *N+=4;
 
  i=0;
  e_err=eps*fabs(s2);
 
  if( fabs(s1-s2) < 30*e_err) i=1;
  else if( fabs(s1-s2) < 1.6*(*aEps)) i=2; 
  if(i)   for(n=0,s3=0;n<5;n++) s3+=F5[n]*func(a+ d*X5[n]); s3*=d; *N+=5;

  if(i==1 && fabs(s3-s2) < e_err)  i=3; 
  if(i==2 && fabs(s3-s2) <  0.1*(*aEps)){i=3;  *aEps -= fabs(s3-s2);}

  
  if(i==3|| *N>NMAX)
  {*ans+=s3;
   *aAns+=fabs(s3);
    return;
  }  
  r_gauss(func,a,(a+b)/2,eps,aEps,ans,aAns,N);
  r_gauss(func,(a+b)/2,b,eps,aEps,ans,aAns,N);
}   

double gauss345( double (*func)(double),double a,double b, double eps)
{
  double aEps; /* absolute error  */
  int i,k;	
  double X4[4]={6.943185E-02,3.300095E-01 ,6.699905E-01 ,9.305682E-01 };
  double F4[4]={1.739274E-01,3.260726E-01 ,3.260726E-01 ,1.739274E-01 };

  if(a==b) return 0;
  for(i=0,aEps=0;i<4;i++) aEps+=F4[i]*fabs(func(a+ (b-a)*X4[i]));
  if(err_code) return 0;

  if(aEps==0.)  return 0;

  eps=eps/2;
  aEps = eps*aEps*fabs(b-a);
  NMAX=50000*pow(2., (-log10(eps)-2)/2.); 
  for(k=0;;k++)
  {  double ans=0., aAns=0., aEps0=aEps;
     int N=4;
     r_gauss(func,a,b,eps,&aEps,&ans,&aAns,&N);
     if(N> NMAX)  err_code=3; 
     if(err_code) return 0;
     if( aEps0-aEps < eps*aAns)  return ans;
     if(k>5) { err_code=1; return ans;}
     aEps=aAns*eps;
  }
}
