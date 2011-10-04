#include "micromegas.h"
#include "micromegas_aux.h"
#include "micromegas_f.h"

/*===========================================================*/
static double PcmOut, totcoef;
static double pvect[16]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

double (*sqme)(int nsub,double *pvect, int * err_code)=NULL;
static double eps=0.001;
int  _nsub_=0;

double  decayPcm(double am0,  double  am1,  double  am2)
{
  double  summ, diffm, pout;
  summ = am1 + am2;
  diffm = am1 - am2;
  if(am0<summ) return 0;
  return sqrt((am0-summ)*(am0+ summ)*(am0-diffm)*(am0+diffm))/(am0*2);
}
          

int  kin22(double PcmIn,double * pmass)
{  
   double ms,md,sqrtS;

   sqrtS=sqrt(pmass[0]*pmass[0]+PcmIn*PcmIn)+sqrt(pmass[1]*pmass[1]+PcmIn*PcmIn);
   PcmOut = decayPcm(sqrtS,pmass[2],pmass[3]);
   if(PcmOut<sqrtS*1.E-4) return 1;
   totcoef = PcmOut /(32.0*M_PI*PcmIn*sqrtS*sqrtS);
   pvect[3] = PcmIn;
   pvect[7] =-PcmIn;
   pvect[0] = sqrt(PcmIn*PcmIn   + pmass[0]*pmass[0]);
   pvect[4] = sqrt(PcmIn*PcmIn   + pmass[1]*pmass[1]);
   pvect[8] = sqrt(PcmOut*PcmOut + pmass[2]*pmass[2]);
   pvect[12]= sqrt(PcmOut*PcmOut + pmass[3]*pmass[3]);

   return 0;
}

double  dSigma_dCos(double  cos_f)
{
   double  r;
   double sin_f=sqrt(fabs((1-cos_f)*(1+cos_f)));
   int err_code=0;
   pvect[11]=PcmOut*cos_f;
   pvect[15]=-pvect[11];
   pvect[10]=PcmOut*sin_f;
   pvect[14]=-pvect[10];
   r = (*sqme)(_nsub_,pvect,&err_code);
   err_code=0;
   return r * totcoef;
}


double cs22(numout * cc, int nsub, double P, double cos1, double cos2 , int * err) 
{
  int i;
  double pmass[4];
  for(i=1;i<=cc->interface->nvar;i++) if(cc->link[i]) cc->interface->va[i]=*(cc->link[i]);

  if( cc->interface->calcFunc()>0 ) {*err=4; return 0;}
  *(cc->interface->gtwidth)=0;
  *(cc->interface->twidth)=0;
  *(cc->interface->gswidth)=0;
  for(i=0;i<4;i++) cc->interface->pinf(nsub,1+i,pmass+i,NULL); 
  *err=0;
  sqme=cc->interface->sqme;
  _nsub_=nsub; 
  if(kin22(P,pmass)) return 0; else return 3.8937966E8*simpson(dSigma_dCos,cos1,cos2,0.3*eps);
}
/*===================  Collider production ==========*/
static double * GGaddress;
static double * Qaddress;
static double sMin,sMax,pcmOut;
static double pmass[4];
static double x0;
static int pc1_,pc2_;
static int ppFlag;

static numout * colliderProduction(char * name1,char *name2)
{ 
  char libname[100], process[100], lName1[20], lName2[20];
  numout *cc;
  int i,first;
    
  pname2lib(name1,lName1);
  pname2lib(name2,lName2);
  if(strcmp(lName1,lName2)>0)sprintf(libname,"PP_%s%s",lName1,lName2);
  else                       sprintf(libname,"PP_%s%s",lName2,lName1); 
  sprintf(process,"proton,proton->%s,%s{",name1,name2);
  
  for(i=0,first=1;i<nModelParticles;i++)
  { switch(abs(ModelPrtcls[i].NPDG))
    { case 1: case 2: case 3: case 4: case 5: case 81: case 83:
       if(!first) strcat(process,","); else  first=0;
       sprintf(process+strlen(process),"%s,%s", 
                     ModelPrtcls[i].name,ModelPrtcls[i].aname);  
       break;
       case 21:
       if(!first) strcat(process,",");else  first=0;
       sprintf(process+strlen(process),"%s",
                             ModelPrtcls[i].name);
    }                              
  }
  
  cc=newProcess(process,libname);

  return cc;
}


static double  cos_integrand(double xcos)
{ int err;
  double xsin=sqrt(1-xcos*xcos);
  double Q,s,t,u;
  pvect[9]=pcmOut*xcos;
  pvect[10]=pcmOut*xsin;
  pvect[13]=-pvect[9];
  pvect[14]=-pvect[10];
  s=(pvect[0]+pvect[4]); s=s*s;
  t=pmass[0]*pmass[0]+pmass[2]*pmass[2]-2*(pvect[0]*pvect[8] -pvect[1]*pvect[9]);
  u=pmass[0]*pmass[0]+pmass[3]*pmass[3]-2*(pvect[0]*pvect[12]-pvect[1]*pvect[13]);
  Q=sqrt(2*s*fabs(t*u)/(s*s+t*t+u*u));
  
  if(GGaddress) *GGaddress=sqrt(4*M_PI*parton_alpha(Q));  
  return  sqme(_nsub_,pvect,&err)*convStrFun2(x0, Q, pc1_,pc2_,ppFlag);  
}


static double  s_integrand(double y)
{  double r;
   double pp=6;
   double s,pcmIn,Qstat;
   s=sMin+pow(y,pp)*(sMax-sMin);
   
   pcmIn=decayPcm(sqrt(s),pmass[0], pmass[1]);
   if(pcmIn==0) return 0;
   pvect[0]=sqrt(pmass[0]*pmass[0]+pcmIn*pcmIn);
   pvect[1]=pcmIn; pvect[2]=0; pvect[3]=0;
   pvect[4]=sqrt(pmass[1]*pmass[1]+pcmIn*pcmIn);
   pvect[5]=-pcmIn; pvect[6]=0; pvect[7]=0;
   pcmOut=decayPcm(sqrt(s),pmass[2], pmass[3]);
   pvect[8]=sqrt(pmass[2]*pmass[2]+pcmOut*pcmOut);
   pvect[11]=0;
   pvect[12]=sqrt(pmass[3]*pmass[3]+pcmOut*pcmOut);
   pvect[15]=0;
   x0=s/sMax;
 
   r=  3.8937966E8*pcmOut/(32*M_PI*pcmIn*s)*simpson(cos_integrand,-1.,1.,1.E-3);
   r*=pp*pow(y,pp-1)*(sMax-sMin);
   return r; 
}

double hCollider(double Pcm, int pp, char * name1,char *name2)
{ 
  double  sigma_tot=0, Qstat;
  int nsub,i;
  numout *cc;
 
  sMax=4*Pcm*Pcm; 
  sMin=pMass(name1)+pMass(name2); sMin*=sMin; sMin+=1;
  ppFlag=pp;   
  cc=colliderProduction( name1,name2);

  GGaddress=NULL; 
  Qaddress=NULL;
  for(i=1;i<=cc->interface->nvar;i++) 
  { if(cc->link[i]) cc->interface->va[i]=*(cc->link[i]);
    { if(strcmp(cc->interface->varName[i],"GG")==0) GGaddress=cc->interface->va+i;
      if(strcmp(cc->interface->varName[i],"Q")==0) Qaddress=cc->interface->va+i;
    }  
  } 
  if(Qaddress)
  { Qstat=*Qaddress;
    *Qaddress=sqrt(sMin);
  }  
  if( cc->interface->calcFunc()>0 )  return 0;
  *(cc->interface->gtwidth)=0;
  *(cc->interface->twidth)=0;
  *(cc->interface->gswidth)=0;
  sqme=cc->interface->sqme;
      
  sigma_tot=0;   
  for(nsub=1;nsub<=cc->interface->nprc; nsub++) 
  { long pc[4];
    double tmp;

    for(i=0;i<4;i++) cc->interface->pinf(nsub,i+1,pmass+i,pc+i);

    if(pc[0]<=pc[1])
    { pc1_=pc[0];
      pc2_=pc[1];
      _nsub_=nsub;
      tmp=simpson(s_integrand,0.,1.,1.E-2)/sMax;     
      sigma_tot+=tmp;
    }  
  }
  if(Qaddress){ *Qaddress=Qstat;}
            
  return sigma_tot;
}

/* ======================decayTable ===================*/

 decayTableStr* decayTable=NULL;
 
 void cleanDecayTable(void)
 { int i,j;
   if(decayTable)
   {  
     for(i=0;i<nModelParticles;i++) for(j=0;j<2;j++) 
     if(decayTable[i].pdList[j]) cleanTxtList(decayTable[i].pdList[j]);
   }  
   else  decayTable=malloc(nModelParticles*sizeof(decayTableStr));
   for(i=0;i<nModelParticles;i++)
   { for(j=0;j<2;j++) decayTable[i].pdList[j]=NULL;
     decayTable[i].width=0;
     decayTable[i].dim=0; 
   }
 }

/* 
static char *  trim(char* p)
{  int n1=0, n2, k=-1;
   n2=(int)strlen(p)-1;
   while(!isgraph(p[n1]) && n1 <= n2) n1++;
   while(!isgraph(p[n2]) && n1 <  n2) n2--;
   while(++k < n2-n1+1)  p[k] = p[k+n1];
   p[k] = '\0';
   return p;
}
*/ 
 

/*======================  1->2 decay ==================*/

double pWidth2(numout * cc, int nsub)
{
  double pvect[12];
  double width=0.;
  double m1,m2,m3; 
  int i,ntot,nin,nout;

  procInfo1(cc,&ntot,&nin,&nout);
  if(nsub<1 ||  nsub>ntot|| nin!=1||nout !=2)  return 0;
     
  for(i=1;i<=cc->interface->nvar;i++) if(cc->link[i]) cc->interface->va[i]=*(cc->link[i]);

  if( cc->interface->calcFunc()>0 ) { return -1;}
  if(cc->Q) 
  {
    cc->interface->pinf(nsub,1,&m1,NULL);
    *(cc->Q)=m1;
    if( cc->interface->calcFunc()>0 ) { return -1;}
    if(cc->SC && cc->GG)  cc->GG=cc->SC;           
  }  
  

  cc->interface->pinf(nsub,1,&m1,NULL);
  cc->interface->pinf(nsub,2,&m2,NULL); 
  cc->interface->pinf(nsub,3,&m3,NULL);
  if(m1 >m2 + m3)
  {   int i,err_code; 
      double md=m2-m3;
      double ms=m2+m3;
      double pRestOut=sqrt((m1*m1 - ms*ms)*(m1*m1-md*md))/(2*m1);
      double totcoef= pRestOut/(8. * M_PI * m1*m1);
           
      for(i=1;i<12;i++) pvect[i]=0;
      pvect[0]=m1;
      pvect[7]=pRestOut;
      pvect[4]=sqrt(pRestOut*pRestOut+m2*m2);
      pvect[11]=-pRestOut;
      pvect[8]=sqrt(pRestOut*pRestOut+m3*m3);
      width = totcoef * (cc->interface->sqme)(nsub,pvect,&err_code);
  }
  return width;
}

 
double decay2Info(char * pname, FILE* f)
{ int i,j,ntot;
  numout * cc;
  double wtot;
  char pname2[20],process[20],plib[20];
  char * dname[8];

  for(i=0,j=0;pname[i];i++)
  if(pname[i]!=' ') pname2[j++]=pname[i];
  pname2[j]=0;
  strcpy(plib,"2width_");
  pname2lib(pname2,plib+7);
  sprintf(process,"%s->2*x",pname2);
  cc=newProcess_(0,1,process,NULL,"",plib,0);
  if(!cc) return -1; 
  procInfo1(cc,&ntot,NULL,NULL);
  if(f) fprintf(f,"\n Partial width for %s->2x decays in GeV\n",pname2); 
  for(wtot=0,i=1;i<=ntot;i++)
  { double w;
    procInfo2(cc,i,dname,NULL);
    w=pWidth2(cc,i);
    if(w!=0)
    { wtot+=w;
      if(f) fprintf(f,"%3.3s %3.3s  %.2E\n",dname[1],dname[2],w); 
    }
  }
  if(f) fprintf(f," Total width %.2E GeV\n",wtot);
  return  wtot;
}

static double decay22List(char * pname, txtList *LL)
{ int i,j,ntot;
  numout * cc;
  double wtot,w;
  char pname2[20],process[20],plib[20];
  char * dname[8];
  txtList L=NULL,l;
  char buff[100];

  for(i=0,j=0;pname[i];i++)
  if(pname[i]!=' ') pname2[j++]=pname[i];
  pname2[j]=0;
  strcpy(plib,"2width_");
  pname2lib(pname2,plib+7);
  sprintf(process,"%s->2*x",pname2);
  cc=newProcess_(0,1,process,NULL,"",plib,0);
  if(!cc) { if(LL) *LL=NULL; return -1;} 
  procInfo1(cc,&ntot,NULL,NULL);
  for(wtot=0,i=1;i<=ntot;i++)
  { 
    procInfo2(cc,i,dname,NULL);
    w=pWidth2(cc,i);
    if(w!=0)
    {  if(LL)
       { l=malloc(sizeof(txtList*));
         l->next=L;
         L=l;
         sprintf(buff,"%E  %s -> %s,%s",w,pname2,dname[1],dname[2]);
         l->txt=malloc(20+strlen(buff));
         strcpy(l->txt,buff);
       } 
       wtot+=w; 
    }
  }

  if(LL)
  {  for(l=L;l;l=l->next)
     { 
       sscanf(l->txt,"%lf %[^\n]",&w,buff);
       sprintf(l->txt,"%E %s",w/wtot,buff);
     }   
    *LL=L;
  }
  return  wtot;
}


/*==================  1->3 decay ===================*/

static double _x_;

static double kimematic_1_3(double *pmass, int kin, double xm2, double xcos, double * P)
{ 
  double factor,pout,mQmin,mQmax,mQ,m12,chY,shY,xsin;
  double p4,p8,p5,p9;
  double m0=pmass[0],m1=pmass[1],m2=pmass[2],m3=pmass[3];
  int i;
  
  P[0]=m0; P[1]=P[2]=0;
  for(i=0;i<4;i++) P[3+i*4]=0;
  
  mQmin=(m1+m2)*(m1+m2);
  mQmax=(m0-m3)*(m0-m3);
  mQ=mQmin*(1-xm2)+mQmax*xm2;
  m12=sqrt(mQ);
   
  factor=(mQmax-mQmin)/(128*M_PI*M_PI*M_PI*m0*m0*m12);
  
  pout=decayPcm(m0,m12,m3);
  P[12]=sqrt(pout*pout+m3*m3); P[13]=-pout; P[14]=0;

  factor*=pout;  
  
  chY=sqrt(1+pout*pout/mQ);
  shY=sqrt(pout*pout/mQ);
  
  pout=decayPcm(m12,m1,m2);
  factor*=pout;
  xsin=sqrt(1-xcos*xcos);
  p4=sqrt(m1*m1+pout*pout);    p8=sqrt(m2*m2+pout*pout);
  p5=xcos*pout;                p9=-p5;
  P[6]=xsin*pout;              P[10]=-P[6];
  
  P[4]=chY*p4 + shY*p5;    P[8]=chY*p8 + shY*p9;
  P[5]=shY*p4 + chY*p5;    P[9]=shY*p8 + chY*p9;
    
  return factor;
}


static double dWidthdCos(double xcos)
{
  double factor;
  double P[16];
  int err_code;

  factor=kimematic_1_3(pmass,1,_x_,xcos, P);
  return factor*(*sqme)(_nsub_,P,&err_code);
}

static double dWidthdM(double x)
{ _x_=x; return simpson(dWidthdCos,-1.,1.,1.E-4);}



static double width13(numout * cc, int nsub, int * err) 
{
  int i;
  for(i=1;i<=cc->interface->nvar;i++) if(cc->link[i]) cc->interface->va[i]=*(cc->link[i]);

  if( cc->interface->calcFunc()>0 ) {*err=4; return 0;}
  *(cc->interface->gtwidth)=0;
  *(cc->interface->twidth)=0;
  *(cc->interface->gswidth)=0;
  for(i=0;i<4;i++) cc->interface->pinf(nsub,1+i,pmass+i,NULL); 
  *err=0;
  sqme=cc->interface->sqme;
  _nsub_=nsub; 
  return simpson(dWidthdM,0.0001,0.9999,1.E-2);
}

static txtList conBrList(txtList BrList)
{ txtList out=NULL;
  char buff[100];
  double br;
  int inCode[10], outCode[10],i;
  for(;BrList;BrList=BrList->next)
  { txtList new=malloc(sizeof(txtList*));
    new->next=out;out=new;
    sscanf(BrList->txt,"%lf %[^\n]",&br,buff); 
    decodeProcess(buff,inCode,outCode);
    if(inCode[0]>0) sprintf(buff,"%E  %s -> ",br,ModelPrtcls[inCode[0]-1].aname);   
    else            sprintf(buff,"%E  %s -> ",br,ModelPrtcls[-inCode[0]-1].name);
    for(i=0;outCode[i];i++)
    { if(i) strcat(buff,",");
      if(outCode[i]>0) strcat(buff,ModelPrtcls[outCode[i]-1].aname);   
      else             strcat(buff,ModelPrtcls[-outCode[i]-1].name);  
    }
    new->txt=malloc(strlen(buff)+1);
    strcpy(new->txt,buff);
  }
  return out;  
}


double pWidth(char *name, txtList * LL,int *dim)
{
  txtList L,l,Lout;
  char libName[100];
  double sum=0,width;
  int i,i0,j,j0;
  
  for(i=0;i<nModelParticles;i++)
  { char *pnames[2]={ModelPrtcls[i].name,ModelPrtcls[i].aname};
    for(j=0;j<2;j++)
    if(strcmp(name,pnames[j])==0) 
    { if(decayTable[i].pdList[j])
      { if(dim)*dim=decayTable[i].dim;
        if(LL) *LL=decayTable[i].pdList[j];
        return decayTable[i].width;
      } else break;
    } if(j!=2) break;    
  }    
  i0=i,j0=j;
  if(i0==nModelParticles)
  { printf("%s out of model particles\n",name);
    if(LL) *LL=NULL;
    if(dim)*dim=0;
    return 0;
  }  
  
  width=decay22List(name,LL);
  
  if(width) 
  { if(dim) *dim=2;
    decayTable[i0].pdList[j0]=*LL;
    if(strcmp(ModelPrtcls[i0].name,ModelPrtcls[i0].aname)) 
                 decayTable[i0].pdList[1-j0]=conBrList(*LL);
    
    decayTable[i0].dim=2;
    decayTable[i0].width=width;  
    return width;
  }

  L= makeDecayList(name,3);
  massFilter(pMass(name),&L);
  gammaGluFilter(&L);
  Lout=NULL;

  for(sum=0,l=L;l;l=l->next)  
  { numout* cc;
    int err;
    process2Lib(l->txt ,libName);
    cc=newProcess(l->txt,libName);
    width=width13(cc, 1, &err);
    sum+=width;
    if(LL)
    {  txtList new=malloc(sizeof(txtListStr));
       new->next=Lout;
       Lout=new;
       new->txt=malloc(strlen(l->txt)+20);
       sprintf(new->txt,"%E  %s",width,l->txt);
    }    
  }
  cleanTxtList(L);
  if(Lout)
  for(L=Lout;L;L=L->next)
  { char buff[100];
    sscanf(L->txt,"%lf %[^\n]",&width,buff);
    sprintf(L->txt,"%E %s",width/sum,buff);  
  }   
  if(LL) *LL=Lout;
  if(dim) *dim=3;
  decayTable[i0].pdList[j0]=Lout;
  if(strcmp(ModelPrtcls[i0].name,ModelPrtcls[i0].aname)) 
               decayTable[i0].pdList[1-j0]=conBrList(Lout);
  decayTable[i0].dim=3;
  decayTable[i0].width=sum;  
  return sum;
}

static int pListEq(char * txt1, char * txt2)  
{  char buff[100];
   char rd1[10][10];
   char rd2[10][10];
   int n1,n2,i1,i2;
   char *ch;
    
   strcpy(buff,txt1); while(ch=strchr(buff,',')) ch[0]=' ';
   
   n1=sscanf(buff,"%s %s %s %s %s %s %s %s %s %s",
   rd1[0],rd1[1],rd1[2],rd1[3],rd1[4],rd1[5],rd1[6],rd1[7],rd1[8],rd1[9]); 
   
   strcpy(buff,txt2); while(ch=strchr(buff,',')) ch[0]=' ';
   
   n2=sscanf(buff,"%s %s %s %s %s %s %s %s %s %s",
   rd2[0],rd2[1],rd2[2],rd2[3],rd2[4],rd2[5],rd2[6],rd2[7],rd2[8],rd2[9]); 
   
   if(n1!=n2) return 0;
   for(i1=0;i1<n1;i1++)
   { for(i2=0;i2<n2;i2++) if(strcmp(rd1[i1],rd2[i2])==0){rd2[i2][0]=0; break;}
     if(i2==n2) return 0;
   } 
   return 1;
}      

double findBr(txtList L, char * pattern)
{ char buff[100];
  char *ch;
  double width;
  
  for(;L;L=L->next)
  { 
     sscanf(L->txt,"%lf %[^\n]",&width,buff);
     ch=strstr(buff,"->");
     ch+=2;
     if( pListEq(ch,pattern)) return width;
  }
  return 0;   
}

        
/*==============End of C-codes ======================*/
/*============ Fortran ==========*/

double cs22_(int*ccf,int*nsub,double*P,double*cos1,double*cos2,int*err)
{ numout*cc;
  memcpy(&cc,ccf,sizeof(cc));
  return cs22(cc,*nsub,*P,*cos1,*cos2 ,err);
} 

double pwidth2_(int * ccf, int *nsub)
{  numout*cc;
   memcpy(&cc,ccf,sizeof(cc));
   return pWidth2(cc,*nsub);
}

void  sethelicities_(double *h1,double *h2) { Helicity[0]=*h1; Helicity[1]=*h2;} 
