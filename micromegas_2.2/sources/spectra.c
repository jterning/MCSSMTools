#include"micromegas.h"
#include"micromegas_aux.h"
#include"micromegas_f.h"

#define NEn 18  /* Number of Energy points for simulation */
#define NZ 250  /* part of energy */
#define Nin 8 
#define Nout 7
#define LnXmin (-10*M_LN10)

static float phidiff[Nin][Nout][NEn][NZ];

static int readSpectra(void)
{ static int rdOk=0;
  int k,l,i,n;
  FILE *f;
  char * buff;
  char *fnames[8]={"cc.dat","bb.dat","tt.dat","ll.dat","ww.dat","zz.dat","gg.dat","mm.dat"};

  if(rdOk) return 0;
  buff=malloc(strlen(WORK)+50);

  for(n=0;n<8;n++)
  {  sprintf(buff,"%s/../../sources/data/%s",WORK,fnames[n]);
     f=  fopen(buff,"r"); 
     if(f==NULL) { free(buff);return 1;}    
     for(i=0;i<6;i++)
     { fscanf(f,"%*s"); 
       for(k=0;k<18;k++)
       { for(l=0;l<NZ;l++) if(1!=fscanf(f,"%f",phidiff[n][i][k]+l)) break;
         fscanf(f,"%*s"); 
         for(;l<NZ;l++)phidiff[n][i][k][l]=0;
       }
     }
     fclose(f);
  }
   
  free(buff);
  rdOk=1;   
  return 0;
}



double zInterp(double zz, double * tab)  /* zz=log(E/Nmass) */
{  
   double z=zz/LnXmin,dz,r;
   int j0;   
   
   if(zz>=0) return 0;
   z=zz/LnXmin;  /*  z_j=0.5/NZ +j/NZ  */
   if(z >= 1-0.5000001/NZ) return tab[NZ-1];
   j0=NZ*z-0.5;
   if(j0<0) j0=0;
   dz=z*NZ-0.5-j0;
   r=(1-dz)*tab[j0]+dz*tab[j0+1];
   if(r<0)r=0;
   return r; 
}

void mInterp(double Nmass,  int  CHin,int  CHout, double*tab)
{  
   float mi[NEn]={10,25,50,80.3,91.2,100,150,176,200,250,350,500,750,1000,1500,2000,3000,5000};
   int l,i0;
   double c0,c1;
   float *p0,*p1;

   for(i0=0; i0<NEn && Nmass>=mi[i0] ;i0++);
   if(i0) i0--;    
   p0=phidiff[CHin][CHout][i0];
   p1=phidiff[CHin][CHout][i0+1];
  
   c1=(Nmass*Nmass -mi[i0]*mi[i0])/(mi[i0+1]*mi[i0+1] -mi[i0]*mi[i0]);
   c0=1-c1;
   for(l=0;l<NZ;l++) tab[l]= c0*p0[l]+c1*p1[l];
}

static double*tabStat;
static double FUNN(double x){return  zInterp(x,tabStat);}
static double FUNE(double x){return  zInterp(x,tabStat)*exp(x);}

void spectrInfo(double Xmin,double*tab, double * Ntot,double*Etot)
{
/*
  if(Xm)
  {  int j,j_m=-1;
    double m=0,jm;
    for(j=0;j<NZ;j++) if(tab[j]>m) {m=tab[j];j_m=j;}
    jm=j_m+(tab[j_m-1]-tab[j_m+1])/(tab[j_m-1]+tab[j_m+1]-2*tab[j_m]);
    *Xm=exp((0.5/NZ+jm/NZ)*LnXmin);
  }
*/  
  tabStat=tab;
  if(Ntot)*Ntot=simpson(FUNN, log(Xmin), 0.,1.E-4);
  if(Etot)*Etot=simpson(FUNE, log(Xmin), 0.,1.E-4)/2;
}

static double m_2;

static double FUNB(double p)
{ 
  return zInterp(log(p),tabStat)/(p*p);
}

static double boost(double Y, double m, double*tab)
{ double chY=cosh(Y), shY=sinh(Y);
  int l;
  double E=chY+sqrt(1-m*m)*shY;
  double tab_out[NZ];  

  m_2=m*m; 
  tabStat=tab;

  for(l=0;l<NZ;l++)
  { double p=E*exp(LnXmin*(l+0.5)/NZ);
    double e=sqrt(p*p+m_2); 
    double p1=fabs(chY*p-shY*e);
    double p2=chY*p+shY*e;
    if(p2>1) p2=1;
    if (p1>=p2) tab_out[l]=0; else tab_out[l]=p*simpson(FUNB,p1,p2,1.E-4)/2/shY;
  }
  
  for(l=0;l<NZ;l++) tab[l]=tab_out[l];

  return E;
}

static double outMass[6]= {0.,0.,1.,0.,0.,0.};
char* outNames[6]={"gamma","e+","p-","nu_e","nu_mu","nu_tau"};

static int _err_=0;

static void getSpectrum(double M, double m1,double m2,char*n1,char*n2, long N1, long N2,int outP, double *tab)
{
  int i;
  int inP=-1;
  long N;
  for(i=0;i<NZ;i++) tab[i]=0; 
  if(N1+N2==0)
  { N=labs(N1);
    switch(N)
    { case 1: case 2: case 3: case 4:   inP=0; break; /* d,u,s,c*/
      case 5: inP=1; break; /* b */ 
      case 6: inP=2; break; /* t */
      case 11:case 13: inP=7; break; /*e,m*/  
      case 15: inP=3; break; /*l*/
      case 24:inP=4; break;  /*w*/
      case 23:inP=5; break;  /*z*/
      case 21:inP=6; break;  /*glu*/ 
      case 12: case 14: case 16: return; /*neutrinos*/
    }
  }
  if(inP>=0)mInterp(M/2,inP, outP,tab); else
  { char* nn[2];
    double mm[2];
    double E[2]; 
    double p2=(M*M-(m1+m2)*(m1+m2))*(M*M-(m1-m2)*(m1-m2))/(4*M*M);
    int k;
    
    nn[0]=n1;
    nn[1]=n2;
    mm[0]=m1;
    mm[1]=m2;
    E[0]=sqrt(m1*m1+p2);
    E[1]=sqrt(m2*m2+p2);
    
    for(i=0;i<NZ;i++) tab[i]=0;
    
    for(k=0;k<2;k++)
    { double tabAux[NZ];
      N=(k==0?labs(N1):labs(N2));
      inP=-1; 
      switch(N)
      {  
         case  1: case 2: case 3: case 4:   inP=0; break; /* d,u,s,c*/
         case  5: inP=1; break; /* b */ 
         case  6: inP=2; break; /* t */
         case 11:case 13: inP=7; break; /*e,m*/
         case 15: inP=3; break; /*l*/
         case 24:inP=4; break;  /*w*/
         case 23:inP=5; break;  /*z*/
         case 21:inP=6; break;  /*glu*/ 
         case 12: case 14: case 16: continue; /*neutrinos*/
      }
      if(inP>=0)
      { double dY=log(M/E[k]/2);
        mInterp(E[k]/2,inP, outP,tabAux);   
        for(i=0;i<NZ;i++)tab[i]+=0.5*zInterp((0.5+i)/NZ*LnXmin+dY,tabAux);
      }else
      { double w=0;
        numout * d2Proc;
        int l; 
        char* n[4];
        double m[4];
        double M_,Y,dY;
        double tab_p[NZ];
        char process[40],plib[40];
        int ntot;
        
        strcpy(plib,"2width_");
        sprintf(process,"%s->2*x",nn[k]);
        pname2lib(nn[k],plib+7);
        for(i=0;i<NZ;i++) tabAux[i]=0;
        if(mm[k]==0) { fprintf(stderr,"Can not hadronize zero mass %s\n",nn[k]);
                       _err_=_err_|1;
                       continue;
                     }
                     
        d2Proc=newProcess_(0,1,process,NULL,NULL,plib,0);
        if(!d2Proc) continue;
        procInfo1(d2Proc,&ntot,NULL,NULL);      
        for(l=1;l<=ntot ;l++)
        {    
          double wP=pWidth2(d2Proc,l);
            
          if(wP>0)
          { long N2=d2Proc->interface->pinfAux(l,2,NULL,NULL,NULL);
            long N3=d2Proc->interface->pinfAux(l,3,NULL,NULL,NULL); 
            procInfo2(d2Proc,l,n,m);    
            getSpectrum(m[0],m[1],m[2],n[1],n[2],N2,N3,outP, tab_p);
            for(i=0;i<NZ;i++) tabAux[i]+=wP*tab_p[i];
            w+=wP;
          }
        }
         
        if(w==0) { fprintf(stderr,"Can find decays for  %s\n",nn[k]);
                   _err_=_err_||2;      
                               continue;
                 }
        Y=acosh(E[k]/mm[k]);
        M_=mm[k]*boost(Y,2*outMass[outP]/mm[k],tabAux);
        dY=log(M/M_);
        for(i=0;i<NZ;i++)tab[i]+=zInterp((0.5+i)/NZ*LnXmin+dY,tabAux)/w; 
      }
    } 
  }
}


static double calcSpectrum0(double v, numout *libPtr,int outP,double *tab)
{
  int k,i;
  double vcsSum=0,M; 
  int ntot;
  for(i=0;i<NZ;i++) tab[i]=0;  

  procInfo1(libPtr,&ntot,NULL,NULL); 
  for(k=1;k<=ntot ;k++) 
  { char * n[4];
    double m[4];

    int err=procInfo2(libPtr,k,n,m);
    if(err==0)
    { long N3=libPtr->interface->pinfAux(k,3,NULL,NULL,NULL);
      long N4=libPtr->interface->pinfAux(k,4,NULL,NULL,NULL);
      double vcs= v*cs22(libPtr,k,v*m[0]/2,-1.,1.,&err);

      M=m[0];
      if(vcs>1.E-8) 
      {  double tab2[NZ]; 
         getSpectrum(m[0]+m[1],m[2],m[3],n[2],n[3],N3,N4,outP,tab2);
         for(i=0;i<NZ;i++) tab[i]+=tab2[i]*vcs;
         vcsSum+=vcs;
      } 
    } else break;
  }
  if(vcsSum)  for(i=0;i<NZ;i++) tab[i]/=vcsSum;
  return 2.9979E-26*vcsSum;
}  

double calcSpectrum(double v,  int outP,double *tab, int *errcode)
{ int n,i,err;
  double vcs, vcs_,tab_[NZ];
  char nameL[10],anameL[10], lib[20],process[20];
  numout* cc;
  char  lop[20];
  _err_=0;
    
  for(i=0;i<NZ;i++) tab[i]=0;  
  vcs=0;
  
  err=readSpectra(); 
  if(err) { printf("calcSpectrum: Can not read data files for spectra\n");
            if(errcode) *errcode=-1; 
            return 0;
          }
  err=sortOddParticles(lop); 
  if(err) { printf("calcSpectrum: Can not calculate %s\n",lop);
            if(errcode) *errcode=-1;
            return 0;
          }
                                            
  for(n=0;n<Nodd;n++) if(strcmp(lop,OddPrtcls[n].name)==0 ||
                         strcmp(lop,OddPrtcls[n].aname)==0  ) break;
  pname2lib(OddPrtcls[n].name,nameL);
  pname2lib(OddPrtcls[n].aname,anameL);
  sprintf(lib,"omg_%s%s",nameL,nameL);
  sprintf(process,"%s,%s->2*x",OddPrtcls[n].name,OddPrtcls[n].name);
  cc=newProcess_(1,1,process,NULL,txtListOddParticles(),lib,0); 
  if(cc) vcs=calcSpectrum0(v,cc, outP,tab);

  if(strcmp(nameL,anameL))
  {  sprintf(lib,"omg_%s%s",nameL,anameL);
     sprintf(process,"%s,%s->2*x",OddPrtcls[n].name,OddPrtcls[n].aname);
     cc=newProcess_(1,1,process,NULL,txtListOddParticles(),lib,0);                
     if(cc)
     {  vcs_=calcSpectrum0(v,cc, outP,tab_);
        vcs=(vcs+vcs_)/2; 
        if(vcs) for(i=0;i<NZ;i++)tab[i]=(tab[i]*vcs+tab_[i]*vcs_)/(vcs+vcs_);
        else    for(i=0;i<NZ;i++)tab[i]=0;
     } 
  }
  return vcs; 
}

int spectrTable(double*tab, char*fname,char*mess,double Xmin,int N)
{ int i;
  FILE *f;
  
  f=fopen(fname,"w");
  if(f==NULL) return 1;  
  fprintf(f,"#type 0 %%curve\n");
  fprintf(f,"#title %s\n",mess);
  fprintf(f,"#yName dN/dx\n");
  fprintf(f,"#xMin %E\n",log10(Xmin));
  fprintf(f,"#xMax %d\n",0);
  fprintf(f,"#xDim %d\n",N);
  fprintf(f,"#xName  x=log10(E/M)\n");
  
  for(i=0;i<N;i++)
  { double x= log10(Xmin) - i*log10(Xmin)/(N-1); 
    fprintf(f,"%E\n",zInterp(log(10.)*x ,tab)*log(10.));
  } 
  fclose(f);
  return 0;
}

int displaySpectrum(double*tab, char*mess,double Xmin)
{
  int err=spectrTable(tab, "_plot_" ,mess,Xmin,250);
  if(!err)system(" trap \"rm -f _plot_\" 0 1 2 3 9; ../CalcHEP_src/bin/plot_view < _plot_ &");
  return err;     
}

#define  SunGCD 8.5   /* Sun distance for galactic center [kpc] */

double rhoQisotermic(double x) 
{ 
/* modified isotermic distribution */
double al=2,be=2,ga=0;
double r0=SunGCD, a0=3.5;  /* kpc           */ 
double rho0=0.3;        /* GeV/c/sm^3    */ 
double rho;
  rho = rho0*pow(r0/x,ga)* pow( (1+ pow(x/a0,al))/(1+pow(r0/a0,al)) , (ga-be)/al);   
  return rho*rho;
}


static double cosfi_stat;
static double (*rhoQ_stat)(double);    


#define rMin 0.01

static double Halo_integrand(double x)
{
  double r,rc,rho;
  if(x==0.) return 0;
  r=1/x-SunGCD;
  rc=sqrt(r*r+SunGCD*SunGCD-2*cosfi_stat*SunGCD*r);
  if(rc<rMin) return 0;
  rho=rhoQ_stat(rc)/x/x;
  if(rc<1.5*rMin) return  2*(rc/rMin-1)*rho; else return rho;
}



double HaloFactor(double fi,double (*rhoQ)(double))
{
  double sm_in_kpc=3.0856775807E21;
  cosfi_stat=cos(fi);
  rhoQ_stat=rhoQ;

  return 1/(8*M_PI)*sm_in_kpc*simpson(Halo_integrand,0.,1/SunGCD,0.0001);      
}




/*=================== Fortran version ============= */

double calcspectrum_(double *v, int*outP,double *tab, int *err)
{  return calcSpectrum(*v, *outP,tab,err); }

void spectrinfo_(double *Xmin,double*tab, double * Ntot,double*Etot)
{ spectrInfo(*Xmin,tab, Ntot,Etot);}


int  spectrtable_(double*tab, char*fname,char*mess,double *Xmin,int *N,
                   int len1,int len2)
{  
   char cfname[200];
   char cmess[200];
   fName2c(fname,cfname, len1);
   fName2c(mess,cmess,len2);

   return spectrTable(tab, cfname,cmess,*Xmin,*N);
}

double rhoqisotermic_(double *x){ return  rhoQisotermic(*x);}
static double (*rhoQ_stat_F)(double*);

static double Halo_integrand_F(double x)
{
  double r,R;
  if(x==0.) return 0;
  r=1/x-SunGCD;
  R=sqrt(r*r+SunGCD*SunGCD-2*cosfi_stat*SunGCD*r);
  return rhoQ_stat_F(&R)/x/x; 
}

double halofactor_(double *fi,double (*rhoQ)(double*))
{
  double sm_in_kpc=3.0856775807E21;
  cosfi_stat=cos(*fi);
  rhoQ_stat_F=rhoQ;

  return 1/(8*M_PI)*sm_in_kpc*simpson(Halo_integrand_F,0.,1/SunGCD,0.0001);      
}

double zinterp_(double*x, double*tab)
{  return zInterp(*x, tab);}


int displayspectrum_(double*tab, char*fmess,double *Xmin,int len)
{
   char cmess[200];
   fName2c(fmess,cmess, len);
   return  displaySpectrum(tab, cmess,*Xmin);
} 
         
