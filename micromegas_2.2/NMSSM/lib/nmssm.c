#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <ctype.h>
#include <string.h>

#include"../../sources/micromegas.h"
#include"../../sources/micromegas_aux.h"
#include"pmodel.h"


extern forCalchep1 nMass,vev,aLambda,mu,aKappa,lass;

extern int FError;
	
static struct { double  tb, mu,Lambda,Kappa,aLambda, aKappa,g2,vev,mw; 
              } extpar;

static int nCall=0;
static double La[10], Las[10],Lass;
extern forCalchep1 la,las;

static double Zn_[5][5],NMassM[5][5],Zh_[3][3],Pa_[2][2];

extern forCalchep2  Zn,Zh,Za,Zu,Zv,Pa,NMM;

static double khi2(void)
{
  double tb=extpar.tb;
  double hl=extpar.Lambda;
  double hk=extpar.Kappa;
  double hls=extpar.aLambda; 
  double hks=extpar.aKappa;
  double xvev=extpar.mu/extpar.Lambda;  
  double g2=extpar.g2;
  double cb=1.0/sqrt(1.0+tb*tb);
  double sb=tb*cb;
  double vev=extpar.vev; 
  double sum;

  double neum[3][3],orimm[3][3],mh[3];
  double mw=extpar.mw;
  int i1,i2,k;
  
  for(i1=0;i1<10;i1++) {La[i1]=0;Las[i1]=0;}
  Lass=0;

  mh[0]=findValW("Mha");
  mh[1]=findValW("Mhb");
  for(i1=0;i1<2;i1++) for(i2=0;i2<2;i2++)  for(k=0,orimm[i1][i2]=0 ;k<2;k++) 
  orimm[i1][i2]+=Pa_[k][i1]*Pa_[k][i2]*mh[k]*mh[k];

  xvev=-(orimm[0][1] -hl*vev*hls)/(2*hk*hl*vev); 
  extpar.mu=xvev*hl;
  

  neum[0][0]= hl*xvev*(hls+hk*xvev)/(sb*cb);

  La[5] =(orimm[0][0] -neum[0][0])/(-2.0*vev*vev);
  
  hks=-(orimm[1][1]-4*vev*vev*hk*hl*cb*sb-vev*vev*hl*hls/xvev*cb*sb)/(3*hk*xvev);
  extpar.aKappa=hks;    

  mh[0]=findValW("MHc");
  La[4]=-(mh[0]*mh[0]-mw*mw-xvev*hl*(hls+hk*xvev)/cb/sb+vev*vev*(hl*hl+La[5])
  )/(vev*vev);

  neum[0][0]= g2*vev*vev*cb*cb+hl*hls*tb*xvev+hk*hl*tb*xvev*xvev;  
  neum[1][1]= g2*vev*vev*sb*sb+hl*hls/tb*xvev+hk*hl/tb*xvev*xvev;              
  neum[2][2]= 4*hk*hk*xvev*xvev+vev*vev*hl*hls/xvev*cb*sb+hk*hks*xvev;
  neum[0][1]=neum[1][0]= -g2*vev*vev*cb*sb+2.0*vev*vev*cb*sb*hl*hl
               -xvev*hl*(hls+xvev*hk) +2*vev*vev*cb*sb*(La[4]+La[5]);
  neum[0][2]=neum[2][0]=vev*(2.0*hl*hl*xvev*cb-2*hk*hl*xvev*sb-hl*hls*sb);
  neum[1][2]=neum[2][1]=vev*(2.0*hl*hl*xvev*sb-2*hk*hl*xvev*cb-hl*hls*cb);

  mh[0]=findValW("Mh1"); 
  mh[1]=findValW("Mh2"); 
  mh[2]=findValW("Mh3");              
  for(i1=0;i1<3;i1++) for(i2=0;i2<3;i2++)  for(k=0,orimm[i1][i2]=0;k<3;k++)
  orimm[i1][i2]+= Zh_[k][i1]*Zh_[k][i2]*mh[k]*mh[k];


  La[1] =(orimm[0][0]-neum[0][0])/(2.0*vev*vev*cb*cb);
  La[2] =(orimm[1][1]-neum[1][1])/(2.0*vev*vev*sb*sb); 
  La[3] =(orimm[0][1]-neum[0][1])/(2.0*vev*vev*sb*cb);
  Las[1]=(orimm[0][2]-neum[0][2])/(2.0*vev*xvev*cb);
  Las[2]=(orimm[1][2]-neum[1][2])/(2.0*vev*xvev*sb);
  Lass  =(orimm[2][2]-neum[2][2])/(2.0*xvev*xvev);

  sum=   (La[1])*(La[1]) + (La[2])*(La[2]) 
          + La[3]*La[3] +La[4]*La[4] + La[5]*La[5] 
          + Las[1]*Las[1] + Las[2]*Las[2] + Lass*Lass;

  /*printf("Al = %E",hls);
  printf("sum = %E\n",sum);*/

  return sum;
}

/* SASHA */
static double  min1(double xx1, double xx2, double (*FQ)(double))
{  double x1=xx1,x2=xx2;
   double x3,x4,f1,f2,f3,f4;

   nCall=0; 
   f1=FQ(x1);
   f2=FQ(x2); 
 
   if(f1<f2) { x3=x1;f3=f1;f1=f2;x1=x2;f2=f3;x2=x3;}
   for(;;) 
   {  x3=x2+ (x2-x1);
      f3=FQ(x3); if(nCall>1000) return x2;
      if(f3<f2) { x1=x2;x2=x3; f1=f2;f2=f3;}  else  break;
   }

   for(;;)
   {
      double d1,d2,d3,ff4,xx4;
      if(f1>f3) {f4=f3;x4=x3;f3=f1;x3=x1;x1=x4;f1=f4;}
            
      if(x1==x2|| x2==x3) return x2;
      d1=f1/(x2-x1)/(x3-x1);
      d2=f2/(x1-x2)/(x3-x2);   
      d3=f3/(x1-x3)/(x2-x3); 
      xx4=(d1*(x2+x3)+d2*(x1+x3)+d3*(x1+x2))/(2*(d1+d2+d3));
      ff4=d1*(x2-xx4)*(x3-xx4)+d2*(x1-xx4)*(x3-xx4)+d3*(x1-xx4)*(x2-xx4);
      
           if(fabs(x1-x2)>10*fabs(x2-x3)) x4=(x1+9*x2)/10;
      else if(fabs(x2-x3)>10*fabs(x1-x2)) x4=(x3+9*x2)/10; 
      else x4=xx4;
         
      if(x4==x2) return x2;
      f4=FQ(x4);  if(nCall>1000) return x2;
      if( (x4-x2)*(x1-x2) >0)
      { if(f4<f2) {x3=x2; f3=f2; x2=x4;f2=f4;} else {x1=x4;f1=f4;}}       
      else 
      { if(f4<f2) {x1=x2; f1=f2; x2=x4;f2=f4;} else {x3=x4;f3=f4;}} 
         
      if( (x1==x2|| x2==x3) ||
          (fabs(x1-x3)<0.01)||
          (fabs(f2-ff4) < 0.001*fabs(ff4)) ||
          (fabs(x1-x3)<0.1 &&  fabs(f2-ff4) < 0.01*fabs(ff4))
        ) return x2;
   } 
}

/* CYRIL */

/*static double  min1(double xx1, double xx2, double (*FQ)(double))
{  double x1=xx1,x2=xx2;
   double x3,x4,f1,f2,f3,f4;

   nCall=0; 
   f1=FQ(x1);
   f2=FQ(x2); 

for(;;)
{  
   if(f1<f2) { x3=x1;f3=f1;f1=f2;x1=x2;f2=f3;x2=x3;}
   for(;;) 
   {  x3=x2+(x2-x1);
      if(fabs(x1-x2)<0.01) return x2;
      f3=FQ(x3);
      if(nCall>1000) return x2;
      if(f3<f2) { x1=x2;x2=x3; f1=f2;f2=f3;} else break;
   }

   for(;;)
   {
      double d1,d2,d3;
      if(f1>f3) {f4=f3;x4=x3;f3=f1;x3=x1;x1=x4;f1=f4;}
      if(fabs(x1-x2)<0.01||fabs(x1-x3)<0.01||fabs(x2-x3)<0.01) return x2;
      
      d1=f1/(x2-x1)/(x3-x1);
      d2=f2/(x1-x2)/(x3-x2);   
      d3=f3/(x1-x3)/(x2-x3); 
      x4=(d1*(x2+x3)+d2*(x1+x3)+d3*(x1+x2))/(2*(d1+d2+d3));
      if(fabs(x2-x4)<0.01) return x2;
      f4=FQ(x4);
      if(nCall>1000) return x2;
      if( (x4-x2)*(x1-x2) >0)
      {if(f4<f2) {x3=x2; f3=f2; x2=x4;f2=f4;} else {x1=x4;f1=f4;}}      
      else 
      {if(f4<f2) {x1=x2; f1=f2; x2=x4;f2=f4;} else {x3=x4;f3=f4;}}

      if(100*fabs(f1-f2) < fabs(f2-f3)) break; 
      if(100*fabs(f3-f2) < fabs(f2-f1)) {f1=f3; x1=x3; break;}
      
   }
}   
   return x2;
}
*/
static double  varAlam(double x){ extpar.aLambda=x;  nCall++;  return  khi2();}

static void CharginoZM(void)
{
  double M2=NMassM[1][1];
  double mu=-NMassM[2][3];
  double g2v1= -NMassM[1][3]*sqrt(2.); 
  double g2v2=  NMassM[1][2]*sqrt(2.);

  double offQ=g2v1*g2v1+g2v2*g2v2;
  double TrX2=offQ +M2*M2+mu*mu;
  double detX=mu*M2 - g2v1*g2v2;
  double D=TrX2*TrX2 - 4.*detX*detX;

  double tU=(g2v2*g2v2-g2v1*g2v1-M2*M2+mu*mu-sqrt(D))/2./(M2*g2v2+mu*g2v1);
  double tV=(g2v1*g2v1-g2v2*g2v2-M2*M2+mu*mu-sqrt(D))/2./(M2*g2v1+mu*g2v2);
  
  double CU=cos(atan(tU));
  double SU=sin(atan(tU));
  double CV=cos(atan(tV));
  double SV=sin(atan(tV));
  double mc[2],Zv_[2][2],Zu_[2][2];
  char name[10];
  int i,j;
  
  Zu_[0][0]=CU;
  Zu_[0][1]=SU;
  Zu_[1][0]=-SU;
  Zu_[1][1]=CU;
  Zv_[0][0]=CV;
  Zv_[0][1]=SV;
  Zv_[1][0]=-SV;
  Zv_[1][1]=CV;

  for(i=0;i<2;i++) mc[i]=g2v1*Zu_[i][0]*Zv_[i][1]
                        +g2v2*Zu_[i][1]*Zv_[i][0]
                        +  M2*Zu_[i][0]*Zv_[i][0]
                        +  mu*Zu_[i][1]*Zv_[i][1];
 
 for(i=1;i<=2;i++)
 { sprintf(name,"MC%d",i);
   assignValW(name,mc[i-1]);
   for(j=1;j<=2;j++)
   { 
      sprintf(name,"Zu%d%d",i,j);
      assignValW(name,Zu_[i-1][j-1]);
      sprintf(name,"Zv%d%d",i,j);
      assignValW(name,Zv_[i-1][j-1]);
   }    
 }

}

static void fillNeutralinoMassMatrix(void)
{ int i,j,k;
  double  M[5];
  
  M[0]=findValW("MNE1");
  M[1]=findValW("MNE2");
  M[2]=findValW("MNE3");
  M[3]=findValW("MNE4");
  M[4]=findValW("MNE5");
 
  for(i=0;i<5;i++) for(j=0;j<5;j++) for(k=0,NMassM[i][j]=0;k<5;k++)
  NMassM[i][j]+=Zn_[k][i]*M[k]*Zn_[k][j];

   assignValW("NMM55", NMassM[4][4]);                                                                                                                                                                                                                                                                                                                                                                                                   
   assignValW("NMM34", NMassM[2][3]);                                                                                                                                                                                                                                                                                                                                                                                                                                                  
   assignValW("NMM35", NMassM[2][4]);                                                                                                                                                                                                                                                                                                                                                                                                                                                  
   assignValW("NMM45", NMassM[3][4]);                                                                                                                                                                                                                                                                                                                                                                                                                                                  
   assignValW("NMM13", NMassM[0][2]);                                                                                                                                                                                                                                                                                                                                                                                                                                                  
   assignValW("NMM14", NMassM[0][3]);                                                                                                                                                                                                                                                                                                                                                                                                                                                  
   assignValW("NMM23", NMassM[1][2]);                                                                                                                                                                                                                                                                                                                                                                                                                                                  
   assignValW("NMM24", NMassM[1][3]);                                                                                                                                                                                                                                                                                                                                                                                                                                                  

 CharginoZM(); 

}

double higgspotent_(void)
{ int i,j;
  double ee,sw,cw,lmax;
  double MZ,alfEMZ,tb,At,Ab;
  char name[10]; 

  MZ=findValW("MZ");
  alfEMZ=findValW("alfEMZ");
  At=findValW("At");
  Ab=findValW("Ab");
  tb=findValW("tb");
  
  for(i=1;i<=2;i++) for(j=1;j<=2;j++) 
  { sprintf(name,"Pa%d%d",i,j); Pa_[i-1][j-1]=findValW(name);}

  for(i=1;i<=3;i++) for(j=1;j<=3;j++) 
  { sprintf(name,"Zh%d%d",i,j); Zh_[i-1][j-1]=findValW(name);}
  
  for(i=1;i<=5;i++) for(j=1;j<=5;j++) 
  { sprintf(name,"Zn%d%d",i,j); Zn_[i-1][j-1]=findValW(name);}

  extpar.tb      = findValW("tb");
  extpar.Lambda  = findValW("Lambda");
  extpar.Kappa   = findValW("Kappa");
  extpar.aLambda = findValW("aLmbd0");
  
  extpar.vev=174.102741;
  assignValW("vev",extpar.vev);

  ee=sqrt(4*M_PI*alfEMZ);
  sw=sin(asin(extpar.vev*sqrt(2.0)*ee/MZ)/2);
  cw=sqrt(1.0-sw*sw);
  extpar.mw=MZ*cw; 
  extpar.g2=ee*ee/(sw*cw)/(sw*cw)/2; 

  extpar.aLambda=min1(extpar.aLambda-1, extpar.aLambda+1, varAlam);
  varAlam(extpar.aLambda);
  assignValW("mu",extpar.mu);
  assignValW("aLmbda",extpar.aLambda);
  assignValW("aKappa",extpar.aKappa);

  fillNeutralinoMassMatrix();

  for(i=1;i<=5;i++) {sprintf(name,"la%d",i);  assignValW(name, La[i]);}
  for(i=1;i<=2;i++) {sprintf(name,"la%ds",i); assignValW(name,Las[i]);}
  assignValW("lass", Lass);
  lmax=0;
  for(i=0;i<10;i++)
  {  if(lmax<fabs(La[i]))lmax=fabs(La[i]);
     if(lmax<fabs(Las[i]))lmax=fabs(Las[i]);
  }
  if(lmax<fabs(Lass))lmax=fabs(Lass);
  if(nCall>2000) return -1;
  else return lmax;
}
