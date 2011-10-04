/*  Direct Detection */

#include "micromegas.h"
#include "micromegas_aux.h"
#include "micromegas_f.h"

void setprotonff_(double *scalar, double *ps_vector, double * sigma)
{ int i; 
  double * p[3];
  p[0]=scalar,p[1]=ps_vector,p[2]=sigma;
  for(i=0;i<3;i++) if(*(p[i])<-12340) p[i]=NULL;
  setProtonFF(p[0], p[1],p[2]);
} 

 void setneutronff_(double *scalar, double *ps_vector, double * sigma)
{  int i; 
   double * p[3];
   p[0]=scalar,p[1]=ps_vector,p[2]=sigma;
   for(i=0;i<3;i++) if(*p[i]<-12340.) p[i]=NULL;
   setNeutronFF(p[0], p[1],p[2]);
}

void getscalarff_(double *muDmd,double *msDmd,double *sigmaPiN,double *sigma0,
                        double * protonFF, double * neutronFF)
{ getScalarFF(*muDmd,*msDmd,*sigmaPiN,*sigma0,protonFF, neutronFF);}


double fescloop_(double * sing, double * mq, double * msq, double *mne)
{ return FeScLoop(*sing,*mq,*msq,*mne);}

static double(*_XYloop)(double*,double*,double*,double*);

double XYloop_(double sing, double mq, double msq, double mne)
{
   return (*_XYloop)(&sing, &mq, &msq, &mne);
}

double noloop_(double *sing, double *mq, double *msq, double *mne)
{
   return 0;
}


int nucleonamplitudes_(double (*LF)(double*,double*,double*,double*),double*pA0,double*pA5,double*nA0,double*nA5)
{ double (*LF_)(double,double,double,double);
  if(LF==&noloop_) LF_=NULL; else
  if(LF== fescloop_) LF_=FeScLoop;else {_XYloop=LF; LF_=&XYloop_;}
/*printf("OK  LF=%p nullloop=%p &nullloop=%p \n",LF, noloop_,&noloop_);*/
  return nucleonAmplitudes(*LF_, pA0,pA5,nA0,nA5);
}

void setfermi_(double *C,double *B, double *a) { SetFermi(*C,*B, *a);}

double fermiff_(int *A, double * Qfermi)
{ return  FermiFF(*A, *Qfermi);}
 

void setfmaxwell_(double *DV,double *v1,double *vmax)
{ SetfMaxwell(*DV,*v1,*vmax);}

double fdvmaxwell_(double *v){return fDvMaxwell(*v);}

double fdvdelta_(double *v){return 1;}
void   setfdelta_(double *v){SetfDelta(*v);}



static double(*_fDv)(double*);
static double(*_S00)(double*);
static double(*_S01)(double*);
static double(*_S11)(double*); 

static double fDv_(double v){ return (*_fDv)(&v);}
static double S00_(double p){ return (*_S00)(&p);}
static double S01_(double p){ return (*_S01)(&p);}
static double S11_(double p){ return (*_S11)(&p);}



double nucleusrecoil_(
     double *rho, double(*fDv)(double*),int*A, int*Z, double*J, 
     double(*S00)(double*),double(*S01)(double*),double(*S11)(double*),
     double (*LF)(double*,double*,double*,double*), double * dNdE )
{  
  double (*LF_)(double,double,double,double);
  if(LF==noloop_) LF_=NULL; else
  if(LF== fescloop_) LF_=FeScLoop;else {_XYloop=LF; LF_=&XYloop_;}
  
  _S00=S00;_S01=S01,_S11=S11;
  if(fDv  == fdvmaxwell_)
    return  nucleusRecoil(*rho,fDvMaxwell, *A,*Z,*J,S00_,S01_,S11_,LF_,dNdE);
  else  if(fDv=fdvdelta_)  
    return  nucleusRecoil(*rho,fDvDelta, *A,*Z,*J,S00_,S01_,S11_,LF_,dNdE);
  else  
  { _fDv=fDv;
    return  nucleusRecoil(*rho,fDv_, *A,*Z,*J,S00_,S01_,S11_,LF_,dNdE);
  }    
}

double nucleusrecoilaux_(
     double *rho, double(*fDv)(double*),int*A, int*Z, double*J, 
     double(*S00)(double*),double(*S01)(double*),double(*S11)(double*),
     double *Mwimp, 
     double *LmbdP, double*XiP, double *LmbdN, double*XiN,  double * dNdE )
{  
   double (*c_fDv)(double);
  _S00=S00;_S01=S01,_S11=S11;
  
  if(fDv  == fdvmaxwell_) c_fDv=fDvMaxwell;
  else  if(fDv=fdvdelta_) c_fDv=fDvDelta;
  else { _fDv=fDv; c_fDv=fDv_;}
   
    return  nucleusRecoilAux(*rho,c_fDv, *A,*Z,*J,S00_,S01_,S11_,
    *Mwimp,*LmbdP,*XiP,*LmbdN,*XiN, dNdE);
}



double nucleusrecoil0_(double *rho, double (*fDv)(double*),
 int*A,int*Z,double*J,double*Sp,double*Sn,
 double (*LF)(double*,double*,double*,double*),double*dNdE)
{
  double (*LF_)(double,double,double,double);
  if(LF==noloop_) LF_=NULL; else
  if(LF== fescloop_) LF_=FeScLoop;else {_XYloop=LF; LF_=&XYloop_;}

  if(fDv  == fdvmaxwell_)
   return nucleusRecoil0(*rho,fDvMaxwell,*A,*Z,*J,*Sp,*Sn,LF_,dNdE);               
  else 
  { _fDv=fDv;  
    return nucleusRecoil0(*rho, fDv_,*A,*Z,*J,*Sp,*Sn,LF_,dNdE);
  }
}

double nucleusrecoil0aux_(double *rho, double (*fDv)(double*),
 int*A,int*Z,double*J,double*Sp,double*Sn,
 double *Mwimp, double *LmbdP, double*XiP, double *LmbdN, double*XiN,
 double*dNdE)
{
  double (*c_fDv)(double);

  if(fDv  == fdvmaxwell_) c_fDv=fDvMaxwell;
  else  if(fDv=fdvdelta_) c_fDv=fDvDelta;
  else { _fDv=fDv; c_fDv=fDv_;}

  return nucleusRecoil0Aux(*rho, c_fDv,*A,*Z,*J,*Sp,*Sn,*Mwimp,
   *LmbdP, *XiP, *LmbdN, *XiN, dNdE);
}

int displayrecoilplot_(double * tab, char * text, int *E1, int *E2,int len)
{
   char c_name[200];
   fName2c(text,c_name,len);
   return  displayRecoilPlot(tab, c_name, *E1, *E2);
}


double cutrecoilresult_(double *tab, int *E1, int *E2)
{
   return  cutRecoilResult(tab, *E1, *E2);
}

 double s00f19_   (double *p){ return    S00F19(*p);}
 double s11f19_   (double *p){ return    S11F19(*p);}
 double s01f19_   (double *p){ return    S01F19(*p);}
 double s00na23_  (double *p){ return    S00Na23(*p);}
 double s11na23_  (double *p){ return    S11Na23(*p);}
 double s01na23_  (double *p){ return    S01Na23(*p);}
 double s00al27_  (double *p){ return    S00Al27(*p);}
 double s11al27_  (double *p){ return    S11Al27(*p);}
 double s01al27_  (double *p){ return    S01Al27(*p);}
 double s00si29_  (double *p){ return    S00Si29(*p);}
 double s11si29_  (double *p){ return    S11Si29(*p);}
 double s01si29_  (double *p){ return    S01Si29(*p);}
 double s00K39_   (double *p){ return    S00K39(*p);}
 double s11K39_   (double *p){ return    S11K39(*p);}
 double s01K39_   (double *p){ return    S01K39(*p);}
 double s00ge73_  (double *p){ return    S00Ge73(*p);}
 double s11ge73_  (double *p){ return    S11Ge73(*p);}
 double s01ge73_  (double *p){ return    S01Ge73(*p);}
 double s00nb93_  (double *p){ return    S00Nb93(*p);}
 double s11nb93_  (double *p){ return    S11Nb93(*p);}
 double s01nb93_  (double *p){ return    S01Nb93(*p);}
 double s00te125_ (double *p){ return    S00Te125(*p);}
 double s11te125_ (double *p){ return    S11Te125(*p);}
 double s01te125_ (double *p){ return    S01Te125(*p);}
 double s00i127_  (double *p){ return    S00I127(*p);}
 double s11i127_  (double *p){ return    S11I127(*p);}
 double s01i127_  (double *p){ return    S01I127(*p);}
 double s00xe129_ (double *p){ return    S00Xe129(*p);}
 double s11xe129_ (double *p){ return    S11Xe129(*p);}
 double s01xe129_ (double *p){ return    S01Xe129(*p);}
 double s00xe131_ (double *p){ return    S00Xe131(*p);}
 double s11xe131_ (double *p){ return    S11Xe131(*p);}
 double s01xe131_ (double *p){ return    S01Xe131(*p);}
 double s00pb207_ (double *p){ return    S00Pb207(*p);}
 double s11pb207_ (double *p){ return    S11Pb207(*p);}
 double s01pb207_ (double *p){ return    S01Pb207(*p);}
 double s00na23a_ (double *p){ return    S00Na23A(*p);}
 double s11na23a_ (double *p){ return    S11Na23A(*p);}
 double s01na23a_ (double *p){ return    S01Na23A(*p);}
 double  s00si29a_(double *p){ return    S00Si29A(*p);}
 double s11si29a_ (double *p){ return    S11Si29A(*p);}
 double s01si29a_ (double *p){ return    S01Si29A(*p);}
 double s00te125a_(double *p){ return    S00Te125A(*p);}
 double s11te125a_(double *p){ return    S11Te125A(*p);}
 double s01te125a_(double *p){ return    S01Te125A(*p);}
 double s00i127a_ (double *p){ return    S00I127A(*p);}
 double s11i127a_ (double *p){ return    S11I127A(*p);}
 double s01i127a_ (double *p){ return    S01I127A(*p);}
 double s00xe129a_(double *p){ return    S00Xe129A(*p);}
 double s11xe129a_(double *p){ return    S11Xe129A(*p);}
 double s01xe129a_(double *p){ return    S01Xe129A(*p);}
 double s00xe131a_(double *p){ return    S00Xe131A(*p);}
 double s11xe131a_(double *p){ return    S11Xe131A(*p);}
 double s01xe131a_(double *p){ return    S01Xe131A(*p);}
 double s00ge73a_ (double *p){ return    S00Ge73A(*p);}
 double s11ge73a_ (double *p){ return    S11Ge73A(*p);}
 double s01ge73a_ (double *p){ return    S01Ge73A(*p);}
 double s00xe131b_(double *p){ return    S00Xe131B(*p);}
 double s11xe131b_(double *p){ return    S11Xe131B(*p);}
 double s01xe131b_(double *p){ return    S01Xe131B(*p);}
