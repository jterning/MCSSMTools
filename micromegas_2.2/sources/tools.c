#include"micromegas.h"
#include"micromegas_aux.h"

static void r_simpson( double(*func)(double),double * f,double a,double b, 
double eps, double * aEps, double * ans, double * aAns, int *deepness)
{
  double f1[5];
  int i;
  int d1=*deepness+1,d2=*deepness+1;
  double s1,s2,s3,e_err;

/*printf("a=%E b=%E d=%d\n",a,b,*deepness);*/

  s1=(f[0]+4*f[4]+f[8])/6;
  s2=(f[0]+4*f[2]+2*f[4]+4*f[6]+f[8])/12;
  s3=(f[0]+4*f[1]+2*f[2]+4*f[3]+2*f[4]+4*f[5]+2*f[6]+4*f[7]+f[8])/24;


  e_err=eps*fabs(s3);
  i=0;
  if( ( fabs(s3-s2) < e_err && fabs(s3-s1) < 16*e_err)) i=1; else
  if( fabs(s3-s2)*(b-a) < 0.1*(*aEps) && fabs(s3-s1)*(b-a) < 1.6*(*aEps)) 
  { i=1;  *aEps -= fabs((s3-s2)*(b-a));}
  
  if(i || *deepness>20)
  { *ans+=s3*(b-a);
    *aAns+=(fabs(f[0])+4*fabs(f[2])+2*fabs(f[4])+4*fabs(f[6])+fabs(f[8]))
          *fabs(b-a)/12;
    return ;
  }
  
  for(i=0;i<5;i++) f1[i]=f[4+i];
  for(i=8;i>0;i-=2)f[i]=f[i/2];

  for(i=1;i<8;i+=2) f[i]=(*func)(a+i*(b-a)/16);

  r_simpson(func,f,a,(a+b)/2,eps,aEps,ans,aAns,&d1);
  for(i=0;i<5;i++) f[2*i]=f1[i];
  for(i=1;i<8;i+=2) f[i]=(*func)((a+b)/2+i*(b-a)/16);
  r_simpson(func, f,(a+b)/2,b,eps,aEps,ans,aAns,&d2);
  if(d1>d2) *deepness=d1; else *deepness=d2;   
}

double simpson( double (*func)(double),double a,double b, double  eps)
{
  double f[9];
  double aEps; /* absolute error  */
  int i;	

  aEps=0;
  if(a==b) return 0;
  for(i=0;i<9;i++) 
  { f[i]=(*func)(a+i*(b-a)/8); aEps +=fabs(f[i]);}
  if(aEps==0.)  return 0;
  eps=eps/2;
  aEps = eps*aEps*fabs(b-a)/9;

  for(;;)
  {  double ans=0., aAns=0.; 
     int deepness=1;
     r_simpson(func,f,a,b,eps,&aEps,&ans,&aAns,&deepness);
     if(5*aAns*eps > aEps) return ans;
     for(i=0;i<9;i++)  f[i]=(*func)(a+i*(b-a)/8);
     aEps=aAns*eps;
  }  
}

double gauss( double (*func)(double),double a,double b, int n)
{
  double X2[2]={2.113249E-01,7.886751E-01 };
  double F2[2]={5.000000E-01,5.000000E-01 };
  double X3[3]={1.127017E-01,5.000000E-01 ,8.872983E-01 };
  double F3[3]={2.777778E-01,4.444444E-01 ,2.777778E-01 };
  double X4[4]={6.943185E-02,3.300095E-01 ,6.699905E-01 ,9.305682E-01 };
  double F4[4]={1.739274E-01,3.260726E-01 ,3.260726E-01 ,1.739274E-01 };
  double X5[5]={4.691008E-02,2.307653E-01 ,5.000000E-01 ,7.692347E-01 ,9.530899E-01 };
  double F5[5]={1.184634E-01,2.393143E-01 ,2.844445E-01 ,2.393143E-01 ,1.184634E-01 };
  double X6[6]={3.376523E-02,1.693953E-01 ,3.806904E-01 ,6.193096E-01 ,8.306047E-01 ,9.662348E-01 };
  double F6[6]={8.566223E-02,1.803808E-01 ,2.339570E-01 ,2.339570E-01 ,1.803808E-01 ,8.566225E-02 };
  double X7[7]={2.544604E-02,1.292344E-01 ,2.970774E-01 ,5.000000E-01 ,7.029226E-01 ,8.707656E-01 ,9.745540E-01 };
  double F7[7]={6.474248E-02,1.398527E-01 ,1.909150E-01 ,2.089796E-01 ,1.909150E-01 ,1.398527E-01 ,6.474248E-02 };
        
        
        
  double ans=0;
 
  switch(n)
  {  int i;
    case 1: ans=(b-a)*func((a+b)/2);  break;
    case 2:
      for(i=0;i<n;i++) ans+=F2[i]*func(a+ (b-a)*X2[i]); break;
    case 3: 
      for(i=0;i<n;i++) ans+=F3[i]*func(a+ (b-a)*X3[i]); break;
    case 4:
      for(i=0;i<n;i++) ans+=F4[i]*func(a+ (b-a)*X4[i]); break;
    case 5:
      for(i=0;i<n;i++) ans+=F5[i]*func(a+ (b-a)*X5[i]); break;
    case 6:
      for(i=0;i<n;i++) ans+=F6[i]*func(a+ (b-a)*X6[i]); break;
    case 7:
      for(i=0;i<n;i++) ans+=F7[i]*func(a+ (b-a)*X7[i]); break;      
    default: 
      return 0;
  }
  return ans*(b-a);                       
 }

/* Numerical recipes codes */

static double*dym,*dyt,*yt,*dysav,*ysav,*ytemp;
static int RKQCprnFlag=1;

static void rk4(double*y, double*dydx, int n, double x,double h,double * yout,
    void (*derivs)(double,double*,double*))
{
   int i;
   double hh=h/2, h6=h/6, xh=x+hh;

   for (i=0;i<n;i++) yt[i]=y[i]+hh*dydx[i];
   (*derivs)(xh,yt,dyt);
   for (i=0;i<n;i++) yt[i]=y[i]+hh*dyt[i];
   (*derivs)(xh,yt,dym);
   for (i=0;i<n;i++) { yt[i]=y[i]+h*dym[i]; dym[i] += dyt[i]; }
   (*derivs)(x+h,yt,dyt);
   for (i=0;i<n;i++) yout[i]=(y[i]+h6*(dydx[i]+dyt[i]+2.0*dym[i]));
}


#define PGROW -0.20
#define PSHRNK -0.25
#define FCOR 0.06666666		/* 1.0/15.0 */
#define SAFETY 0.9
#define ERRCON 6.0e-4


static int rkqc(double * y, double * dydx, int n, double * x, double htry, 
     double eps,  double * yscal,  double* hdid, double* hnext,
     void (*derivs)(double,double *,double *))
{
   int i;
   double xsav=(*x),h=htry;

   for(i=0;i<n;i++) {ysav[i]=y[i]; dysav[i]=dydx[i];}

   for (;;) 
   {  double hh= 0.5*h, errmax=0;
      rk4(ysav,dysav,n,xsav,hh,ytemp,derivs);
      *x=xsav+hh;
      (*derivs)(*x,ytemp,dydx);
      rk4(ytemp,dydx,n,*x,hh,y,derivs);
      *x=xsav+h;
      if (*x == xsav && RKQCprnFlag) 
      { /*printf("Step size too small in routine RKQC\n");*/
        RKQCprnFlag=0;
        return 1;
      }
      rk4(ysav,dysav,n,xsav,h,ytemp,derivs);
      for (i=0;i<n;i++) 
      {  double temp;
         ytemp[i]=y[i]-ytemp[i];
         if(!finite( ytemp[i])) { errmax=ytemp[i]; break;} 
         temp= fabs(ytemp[i]/yscal[i]);
	 if (errmax < temp) errmax=temp;
      }
      if(!finite(errmax)){ h=h/10; continue;}
      errmax /= eps;

      if (errmax <= 1.0) 
      {
	 *hdid=h;
	 *hnext=((errmax > ERRCON ? SAFETY*h*exp(PGROW*log(errmax)) : 4*h));
	  break;
      }
      {  double h_=(SAFETY*h*exp(PSHRNK*log(errmax)));
              if(h_>10*h ) h=10*h;
         else if(h_<0.1*h) h=0.1*h; 
         else  h=h_;
      }  
   }
   for (i=0;i<n;i++) y[i] += (double) (ytemp[i]*FCOR);
   return 0;
}

#define MAXSTP 10000
#define TINY 1.0e-30


int  odeint(double * ystart, int nvar, double x1, double x2, double eps, 
         double h1, void (*derivs)(double,double *,double *))
{
   int nstp,i;
   double x,hnext,hdid,h;

   double *yscal,*y,*dydx;


   double ** allAlloc[9]={NULL,NULL,NULL,&dym,&dyt,&yt,&dysav,&ysav,&ytemp};
   allAlloc[0]=&yscal;allAlloc[1]=&y; allAlloc[2]=&dydx;
   for(i=0;i<9;i++) *allAlloc[i]=(double*)malloc(nvar*sizeof(double)); 

   RKQCprnFlag=1; 
   x=x1;
   h=((x2 > x1) ? fabs(h1) : -fabs(h1));
   for (i=0;i<nvar;i++) y[i]=ystart[i];
   for (nstp=1;nstp<=MAXSTP;nstp++) 
   {
      (*derivs)(x,y,dydx);
      for (i=0;i<nvar;i++) yscal[i]=(fabs(y[i])+fabs(dydx[i]*h)+TINY);

      if ((x+h-x2)*(x+h-x1) > 0.0) h=x2-x;
      if(rkqc(y,dydx,nvar,&x,h,eps,yscal,&hdid,&hnext,derivs))break;
      
      if ((x-x2)*(x2-x1) >= 0.) 
      {
         for (i=0;i<nvar;i++) ystart[i]=y[i];
         for(i=0;i<9;i++) free(*allAlloc[i]);
         return 0;
      }
      h=hnext;
   }
   for(i=0;i<9;i++) free(*allAlloc[i]);
   return 1;
}


static double bessi0(double  x)
{
	double ax,ans;
	double y;

	if ((ax=fabs(x)) < 3.75) {
		y=x/3.75;
		y*=y;
		ans= (1.0+y*(3.5156229+y*(3.0899424+y*(1.2067492
			+y*(0.2659732+y*(0.360768e-1+y*0.45813e-2))))));
	} else {
		y=3.75/ax;
		ans= ((exp(ax)/sqrt(ax))*(0.39894228+y*(0.1328592e-1
			+y*(0.225319e-2+y*(-0.157565e-2+y*(0.916281e-2
			+y*(-0.2057706e-1+y*(0.2635537e-1+y*(-0.1647633e-1
			+y*0.392377e-2)))))))));
	}
	return ans;
}

static double bessi1(double x)
{
	double ax,ans;
	double y;

	if ((ax=fabs(x)) < 3.75) {
		y=x/3.75;
		y*=y;
		ans=(ax*(0.5+y*(0.87890594+y*(0.51498869+y*(0.15084934
			+y*(0.2658733e-1+y*(0.301532e-2+y*0.32411e-3)))))));
	} else {
		y=3.75/ax;
		ans=(0.2282967e-1+y*(-0.2895312e-1+y*(0.1787654e-1
			-y*0.420059e-2)));
		ans= (0.39894228+y*(-0.3988024e-1+y*(-0.362018e-2
			+y*(0.163801e-2+y*(-0.1031555e-1+y*ans)))));
		ans *= ((exp(ax)/sqrt(ax)));
	}
	return x < 0.0 ? -ans : ans;
}

static double bessk0(double x)
{
	double y,ans;

	if (x <= 2.0) {
		y=x*x/4.0;
		ans=(-log(x/2.0)*bessi0(x))+(-0.57721566+y*(0.42278420
			+y*(0.23069756+y*(0.3488590e-1+y*(0.262698e-2
			+y*(0.10750e-3+y*0.74e-5))))));
	} else {
		y=2.0/x;
		ans=(exp(-x)/sqrt(x))*(1.25331414+y*(-0.7832358e-1
			+y*(0.2189568e-1+y*(-0.1062446e-1+y*(0.587872e-2
			+y*(-0.251540e-2+y*0.53208e-3))))));
	}
	return  ans;
}


static double bessk1(double x)
{
	double y,ans;

	if (x <= 2.0) {
		y=x*x/4.0;
		ans=(log(x/2.0)*bessi1(x))+(1.0/x)*(1.0+y*(0.15443144
			+y*(-0.67278579+y*(-0.18156897+y*(-0.1919402e-1
			+y*(-0.110404e-2+y*(-0.4686e-4)))))));
	} else {
		y=2.0/x;
		ans=(exp(-x)/sqrt(x))*(1.25331414+y*(0.23498619
			+y*(-0.3655620e-1+y*(0.1504268e-1+y*(-0.780353e-2
			+y*(0.325614e-2+y*(-0.68245e-3)))))));
	}
	return  ans;
}

static double bessk2(double x)
{
	double bk,bkm,bkp,tox;

	tox= 2.0/x;
	bkm=bessk0(x);
	bk=bessk1(x);
	bkp=bkm+tox*bk;
	bkm=bk;
	bk=bkp;
	return bk;
}

double K2pol(double x)
{
   if(x<0.1) return 1+ 1.875*x*(1+0.4375*x*(1-0.375*x));
   else      return bessk2(1/x)*exp(1/x)*sqrt(2/M_PI/x);
}

double K1pol(double x)
{
  if(x<0.1) return 1+ 0.375*x*(1-0.3125*x*(1+0.875*x));
  else      return bessk1(1/x)*exp(1/x)*sqrt(2/M_PI/x); 
}

static double  polintN(double x, int n,  double *xa, double *ya)
{  double z[10];
   int i,m;
                                                                                
   for(i=0;i<n;i++) z[i]=ya[i];
                                                                                
   for(m=1;m<n;m++) for(i=0;i<n-m;i++)
   z[i]=(z[i]*(xa[i+m]-x) - z[i+1]*(xa[i]-x))/(xa[i+m]-xa[i]);
   return z[0];
}

static int  leftX(int dim, double * xa, double x)
{  int k1,k2,k3;
                                                                                
   if(x<=xa[0]) return 0;
   if(x>=xa[dim-3]) return dim-3;
                                                                                
   k1=0;
   k2=dim-3;
                                                                                
   while(k2-k1>1)
   { k3=(k1+k2)/2;
     if(xa[k3]>x)k2=k3; else k1=k3;
   }
   return k1;
}


double polint3(double x, int n,  double *xa, double *ya)
{ int shift=leftX(n, xa, x);
   return polintN(x,3,xa+shift, ya+shift);
}

#define SQ(x)  ((x)*(x))

double   LintIk(int I,double MSQ,double MQ,double MNE)
{  
  double LAM,SPPM,SPMM,R1,R2,R3,del,CMD;
  double msq2=MSQ*MSQ, mq2=MQ*MQ, mne2=MNE*MNE;

  SPPM=  msq2+mq2-mne2, SPMM= msq2-mq2-mne2;
  R1  =(msq2-mq2)/mne2;
  R2  =(mq2-mne2)/msq2;
  R3  =(msq2-mne2)/mq2;

  del =2.*mne2*(mq2+msq2)-mne2*mne2-SQ(msq2-mq2);
  
  if(del>0) LAM=2.*atan(sqrt(del)/SPPM)/sqrt(del);
  if(del<0) LAM=log((SPPM+sqrt(-del))/(SPPM-sqrt(-del)))/sqrt(-del);

  switch(I)
  { 
    case 1:  
      CMD=1./del*(R2/3-2/3.*R3-5/3.+ (2*msq2-2/3.*mne2)*LAM);	
    break; 
    case 2: 
      CMD=(log(msq2/mq2)-SPMM*LAM)/2./mne2/mne2+
         ( ((mq2*mq2-mq2*msq2)/mne2-7/3.*mq2+2/3.*(mne2-msq2))*LAM+R2/3+R1+2/3.
         )/del;
    break;
    case 3:
      CMD=-3/SQ(del)*SPPM+LAM/del*(-1+6*mq2*msq2/del);
    break;	
    case 4:
      CMD=((log(msq2/mq2) - SPMM*LAM)/2/mne2-1/msq2 -mq2*SPMM/del*LAM)/mne2/mne2
      
         +( mq2/mne2/mne2-SQ(1-mq2/mne2)/msq2+0.5/mne2
            +3*mq2/del*(1 +  R1 + (-R1*mq2-2*mq2-msq2+mne2)*LAM)
          )/del;
    break;
    case 5:
     CMD=(log(msq2/mq2)-SPMM*LAM)/(2*mne2*mne2)-(LAM*(2*(msq2-mne2)+3*mq2+R1*mq2)-3-R1)/del;
     break;
    default: CMD=0.; 
  }
  return CMD;
}

double MaxGapLim(double x, double mu) 
/* S.Yellin, Phys.Rev. D66,032005(2002)

   There is a theoretical model which predicts homogenious event distribution 
   with everage number of events mu. Let experiment gets a gap bitween points 
   where according to theory x point are expected. Then the theoretical model 
   is non-confirmed with probability MaxGap   
*/  
{
  int k;
  double C0,kf;
  for(k=0,C0=0,kf=1;k<mu/x; k++,kf*=k) {C0+= pow(k*x-mu,k)*exp(-k*x)*(1+k/(mu-k*x))/kf;}
  return C0;   
}

/*
int main(void)
{ double x;
  for(x=0.00001; x<1; x*=1.5)
  printf("x=%e bessk2=%e\n",x,  bessk2(x)*x*x); 

}
*/
