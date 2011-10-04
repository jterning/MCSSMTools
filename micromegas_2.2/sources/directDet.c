#include"micromegas.h"
#include"micromegas_aux.h"
#include"micromegas_f.h"
#include"../CalcHEP_src/c_source/model_aux/include/alpha_s.h"
#define  NEWMODEL
#define  GG_NLO

int QCDcorrections=1;
int Twist2On=1;

#define INMAX 100

static int  create2_th_model(char * pname,char *ProcNameSI, char* ProcNameSD)
{  
   char *fname,*command;
   char buff[100];
   char name1[20],name2[20],mass[20],aux[20];
   long num;
   int spin,color,line,i,k,err;
   int neutr;
   char aP_P[20];
   char * auxVert[10];
   int nAuxVerts;  
   int K;
   FILE*f;
   
   struct pProp
   { int num;
     char name[10];
     char aname[10];
     char mass[10];
     int spin2;
     int color;
     int line;
   } pInvolved[INMAX];
       
   for(i=0;i<INMAX;i++) pInvolved[i].num=0;   

   fname=malloc(strlen(WORK)+40);
   sprintf(fname,"%s/models/prtcls1.mdl",WORK);
  
   f=fopen(fname,"r");
   free(fname);
   if(f==NULL) return 1;
   K=1;
   fscanf(f,"%*[^\n]\n%*[^\n]\n%*[^\n]\n");
   for(line=1;7==fscanf(f,"%*[^|]|%[^|]|%[^|]|%ld |%d |%[^|]|%*[^|]|%d |%[^|]|%*[^\n]\n",
                             name1,name2,&num,&spin,mass,        &color,aux);line++)
   { 
      sscanf(name1,"%s",buff);  strcpy(name1,buff);
      sscanf(name2,"%s",buff);  strcpy(name2,buff);
      sscanf(mass,"%s",buff);   strcpy(mass,buff);
      sscanf(aux,"%s",buff);    strcpy(aux,buff);
   
      if((strcmp(pname,name1)==0||strcmp(pname,name2)==0||
      (abs(color)!=1 && abs(num)!=21))&& strcmp(aux,"*"))
      { 
        if(color!=1) k=K++; else k=0;
         pInvolved[k].num=num;
         strcpy(pInvolved[k].name,name1);
         strcpy(pInvolved[k].aname,name2);
         strcpy(pInvolved[k].mass,mass);
         pInvolved[k].spin2=spin;
         pInvolved[k].color=color;
         pInvolved[k].line=line;
      }
   }
   fclose(f);

   if(pInvolved[0].num==0)
   { printf("Direct detection module can not work because the model contains a particle\n");
     printf("with zero PDG code %s\n",  pInvolved[i].name);
     return 10;     
   }
   
   for(i=1;i<K;i++) switch(pInvolved[i].num)
   { case  81: pInvolved[i].num= 1; break;
     case  83: pInvolved[i].num= 3; break;
     case -81: pInvolved[i].num=-1; break;
     case -83: pInvolved[i].num=-3; break;
   }
   
   for(i=1;i<K;i++)
   if(pInvolved[i].num>=-6 && pInvolved[i].num<=6 && strcmp(pInvolved[i].mass,"0")==0)
   {  printf("\n Error! Direct detection module can not work because the model contains\n");
      printf("zero mass quark  '%s'\n",  pInvolved[i].name); 
      return 1;
   }
#ifdef NEWMODEL  
 
   command=malloc(strlen(WORK)+6000);
   sprintf(command," cd %s/models; cp vars1.mdl vars2.mdl;"
   "cp prtcls1.mdl prtcls2.mdl; cp func1.mdl func2.mdl;"
   "cp lgrng1.mdl lgrng2.mdl; cp extlib1.mdl extlib2.mdl\n",WORK);

   err=system(command); 
   if(err) {free(command);return 2;} 
   
   sprintf(command,"cd %s; rm -f tmp/safe; ../../CalcHEP_src/bin/s_calchep -blind \"[{[[{[[[[[{Auxiliary_Model{}0\"",WORK);
   err=system(command);
   if(err) {free(command);return 2;}
      
   sprintf(command,"cd %s; rm -f tmp/safe; ../../CalcHEP_src/bin/s_calchep -blind "
   "\"[{[[{{{MXXX\\091\\09\\09"

   "{sqrt6D\\092.449489742783178\\09\\09"
   
   "}}0\"",WORK);
   err=system( command);
   if(err) {free(command);return 3;}

   sprintf(command,"cd %s;  ../../CalcHEP_src/bin/s_calchep -blind \"[[{[[{"
   "\\19%d{\\04\\19\\${{}}0\";",WORK,pInvolved[0].line);
   err=system( command);
   if(err) {free(command);return 3;}
      
 
    sprintf(command,"cd %s;  ../../CalcHEP_src/bin/s_calchep -blind \"[[{[[{"   
   "{S0\\09S0\\09S0\\090\\090\\09MXXX\\090\\091\\09*\\09S0\\09S0\\09"
   "{V5\\09V5\\09V5\\090\\092\\09MXXX\\090\\091\\09*\\09V5\\09V5\\09"
      
   "}"
   
   "}0\";",WORK);
   err=system( command);
   if(err) {free(command);return 3;}

   sprintf(command,"cd %s; ../../CalcHEP_src/bin/s_calchep -blind \"[[{[[[{",WORK);
   for(i=1;i<K;i++) 
   { sprintf(aP_P,"%s\\09%s",pInvolved[i].aname,pInvolved[i].name); 
 
     if(pInvolved[i].spin2==1)
     {      sprintf(command+strlen(command),     
       "{%s\\09S0\\09\\09MXXX\\091\\09"
       "{%s\\09V5\\09\\09-MXXX\\09G5*G(m3)\\09"
       ,aP_P,aP_P);        
     }
     else  if(pInvolved[i].spin2==0)
     {
       sprintf(command+strlen(command),
       "{%s\\09S0\\09\\09MXXX\\091\\09",aP_P);
     }
   }
   sprintf(command+strlen(command),"}}0\""); 
   err=system( command);  
   
   if(err) {free(command);return 4;}

   sprintf(command,"cd %s; ../../CalcHEP_src/bin/s_calchep -blind \"[[{[[[{",WORK);

   neutr=strcmp(pInvolved[0].name,pInvolved[0].aname)==0; 
   sprintf(aP_P,"%s\\09%s", pInvolved[0].aname, pInvolved[0].name); 
   if (pInvolved[0].spin2==0)
   {  sprintf(buff,"S0\\09\\092*MXXX*2*%s\\091\\09",pInvolved[0].mass);
      nAuxVerts=1;
      auxVert[0]=buff;
   }else 
   if(pInvolved[0].spin2==1)
   {  
      nAuxVerts=2;
      auxVert[0]="S0\\09\\092*MXXX\\091\\09";               
      auxVert[1]="V5\\09\\092*MXXX\\09G5*G(m3)\\09";
   }else 
   if (pInvolved[0].spin2==2) 
   { 
      nAuxVerts=2;
      sprintf(buff,"S0\\09\\092*MXXX*2*%s\\09m1.m2\\09",pInvolved[0].mass);
      auxVert[0]=buff;
      auxVert[1]="V5\\09\\092*i*MXXX*sqrt6D\\09eps(m1,m2,(p1-p2),m3)\\09";
   }   
   for(i=0;i<nAuxVerts;i++) 
   sprintf(command+strlen(command),"{%s\\09%s",aP_P,auxVert[i]);
    
   strcat(command, "}}0\"");

   err=system(command);
   free(command); 
    
   if(err) return 5; 
#endif   
   sprintf(ProcNameSI,"}}[{[{]{QUARKS,%s->QUARKS,%s{",pname, pname);
   strcpy(ProcNameSD,ProcNameSI);
   { char lquarks[50], allcol[300];
      for(i=1,lquarks[0]=0, allcol[0]=0;i<K;i++)
      { int forlQ=(abs(pInvolved[i].num)<=3);
        sprintf(allcol+strlen(allcol),"%s,",pInvolved[i].name);
        if(forlQ)sprintf(lquarks+strlen(lquarks),"%s,",pInvolved[i].name);
        if(strcmp(pInvolved[i].name,pInvolved[i].aname))
        {
          sprintf(allcol+strlen(allcol),"%s,",pInvolved[i].aname);
          if(forlQ)sprintf(lquarks+strlen(lquarks),"%s,",pInvolved[i].aname);
        }        
      }
      
      if(strlen(lquarks))lquarks[strlen(lquarks)-1]=0;
      if(strlen(allcol))  allcol[strlen( allcol)-1]=0; 
      if( allcol[0]) strcat(ProcNameSI, allcol); else ProcNameSI[0]=0; 
      if(lquarks[0]) strcat(ProcNameSD,lquarks); else ProcNameSD[0]=0;
      if(pInvolved[0].spin2==0) ProcNameSD[0]=0;
   }    
   return 0;
}     

static int getAuxCodesForDD(char*pname,numout**ccSI,numout**ccSD)
{ 
  char libnameSD[40]="dirSD";
  char libnameSI[40]="dirSI";
  
  char ProcNameSI[500], ProcNameSD[500];
  int err=0;
  pname2lib(pname,libnameSD+5);
  pname2lib(pname,libnameSI+5);

  if(ccSI)*ccSI=newProcess_(1,1,"",NULL,NULL,libnameSI,0);
  if(ccSD)*ccSD=newProcess_(1,1,"",NULL,NULL,libnameSD,0);
  
  if(  ccSI&&(!*ccSI) || ccSD&&(!*ccSD) ) 
  { int LFid,delay;  
    LFid=checkLockFile(&delay);
    err=create2_th_model(pname,ProcNameSI,ProcNameSD);
    removeLock(LFid);
    if(err) return err;
  }  
   
  if(ccSI && !*ccSI && ProcNameSI[0])
      *ccSI=newProcess_(1,1,ProcNameSI,"S0!=1,V5",NULL,libnameSI,0);
 
  if(ccSD && !*ccSD && ProcNameSD[0])
      *ccSD=newProcess_(1,1,ProcNameSD,"V5!=1,S0",NULL,libnameSD,0);

  return 0;
}


double
/*scalar */
   wS0P__[6]={0.03302,0.02348,0.2594,0.05067,0.05067,0.05067},
   wS0N__[6]={0.04241,0.01823,0.2594,0.05037,0.05037,0.05037},
/* pseudo-vertor */ 
   wV5P__[3]={-0.427, 0.842,-0.085},    
   wV5N__[3]={0.842 ,-0.427,-0.085};
static double 
   
/* sigma */
   wSM0P[3]={-0.23, 0.84,-0.046},   
   wSM0N[3]={ 0.84,-0.23,-0.046};

double GluWeight(void)
{ double s=0;
  int k;
  for(k=0;k<3;k++) s+=wS0P__[k];
  return 1-s;
}  
void setProtonFF(double *scalar, double *ps_vector, double * sigma)
{  int i;
   double s;
   if(scalar){  for(s=0,i=0;i<3;i++) {wS0P__[i]=scalar[i];s+=scalar[i];}
                s=2./27.*(1-s);
                for(i=3;i<6;i++)wS0P__[i]=s;
             }
   if(ps_vector) for(i=0;i<3;i++) wV5P__[i]=ps_vector[i];
   if(sigma)     for(i=0;i<3;i++) wSM0P[i]=sigma[i];
}

void setNeutronFF(double *scalar, double *ps_vector, double * sigma)
{  int i;
   double s;
   if(scalar) {   for(s=0,i=0;i<3;i++) {wS0N__[i]=scalar[i];s+=scalar[i];}
                  s=2./27.*(1-s);
                  for(i=3;i<6;i++) wS0N__[i]=s;
              }
   if(ps_vector) for(i=0;i<3;i++) wV5N__[i]=ps_vector[i];
   if(sigma)     for(i=0;i<3;i++) wSM0N[i]=sigma[i];
}

void getScalarFF(double muDmd, double msDmd, double sigmaPiN, double sigma0,
double * protonFF, double * neutronFF)
{ double Mp=0.9383,Mn=0.9396;
  double y=1-sigma0/sigmaPiN, z=1.49,alpha=(2*z-(z-1)*y)/(2+(z-1)*y);
  protonFF[0]=2*sigmaPiN/((1+ muDmd)*Mp*(1+alpha))/1000.;
  protonFF[1]=muDmd*alpha*protonFF[0];
  protonFF[2]=sigmaPiN*y*msDmd/(1+ muDmd)/Mp/1000.;
  neutronFF[0]=2*sigmaPiN*alpha/((1+ muDmd)*Mn*(1+alpha))/1000.;
  neutronFF[1]=muDmd*neutronFF[0]/alpha;
  neutronFF[2]=sigmaPiN*y*msDmd/(1+ muDmd)/Mn/1000.;   
}


double FeScLoop(double sgn, double mq,double msq,double mne)
{ 
/*  return   1/(msq*msq-mne*mne); */
return  1.5*mq*( 
       mq*(LintIk(1,msq,mq,mne)
  -2./3.*mne*mne*LintIk(3,msq,mq,mne))

 -sgn*mne*(LintIk(2,msq,mq,mne)-1./3.*LintIk(5,msq,mq,mne)-2./3.*mne*mne*LintIk(4,msq,mq,mne))
                ); 
}

static double pdfQnum,p1_p2;

double twist2FF__(double sgn, double mq,double msq,double mne)
{ 
  double D=(msq*msq-mne*mne -mq*mq);
  double D2=D*D-4*mq*mq*mne*mne;
  
  return  parton_x(pdfQnum,1.2*(msq-mne*mne/msq))/D2*(D +sgn*2*p1_p2);
}


double twist2_subtractionFF__(double sgn, double mq,double msq,double mne)
{ double D=(msq*msq-mne*mne -mq*mq);
  double D2=D*D-4*mq*mq*mne*mne;
  return   1/D2*(D +sgn*2*p1_p2);
}


double zeroloopFactor(double sgn, double mq,double msq,double mne)
{ 
   if(msq>0.5)return 0;
  return   1/(msq*msq-mne*mne-mq*mq -2*sgn*mne*mq);
}



extern double (*loopFF__)(double,double,double,double);

int nucleonAmplitudes(double(*LF)(double,double,double,double), double*pA0,double*pA5,double*nA0,double*nA5) 
{
  double wimpMass; 
  int wimpN,Qnum,aQnum,i,I,sgn,ntot,n;
  double  MN=0.939; 
  numout *ccSI,*ccSD,*cc; 
  double wGluP=wS0P__[4]*27./2., wGluN=wS0N__[4]*27./2.;
  double  wV0P[3]={1,2,0}, wV0N[3]={2,1,0}; /* vector current */
  double pvect[16]; 
  char WIMP[20];

  for(i=0;i<2;i++) {pA0[i]=0; pA5[i]=0; nA0[i]=0; nA5[i]=0;}
  
  if(sortOddParticles(WIMP))return 1;
  wimpMass= lopmass_();
  
  if(sWidth>0) for(i=0;i< Nodd;i++) assignVal(OddPrtcls[i].width,0*wimpMass*sWidth);
  for(wimpN=0; wimpN<Nodd;wimpN++) if(strcmp(OddPrtcls[wimpN].name,WIMP)==0 ||
                           strcmp(OddPrtcls[wimpN].aname,WIMP)==0) break;

  if(OddPrtcls[wimpN].spin2) getAuxCodesForDD(WIMP,&ccSI,&ccSD);
  else { getAuxCodesForDD(WIMP,&ccSI,NULL); ccSD=NULL;}

  if(!ccSI && !ccSD)  return -1;
  
  for(I=0;I<2;I++)
  { if(I) cc=ccSD; else cc=ccSI;
    if(!cc) continue;
    procInfo1(cc,&ntot,NULL,NULL);

    for(i=1;i<=cc->interface->nvar;i++) if(cc->link[i]) cc->interface->va[i]=*(cc->link[i]);
    if(cc->Q) *(cc->Q)=2*wimpMass;

    cc->interface->calcFunc();
    
    *(cc->interface->gtwidth)=0;
    *(cc->interface->twidth)=1;
    *(cc->interface->gswidth)=0;

    for(n=1;n<=ntot;n++)
    { double cs0, ampl, qMass;
      char * names[4];
      double masses[4];
      long pdg[4];
      int l,err_code;
      int spin2,color,neutral,neutralWIMP;
      double alphaMq, qcdNLO;
      for(l=0;l<4;l++)   names[l]= cc->interface->pinf(n,l+1,masses+l,pdg+l);
      if(strcmp(names[0],names[2]) ) continue; 

      cc->interface->pinfAux(n,2,NULL,NULL,&neutralWIMP); 
      cc->interface->pinfAux(n,1,&spin2,&color,&neutral);
    
      switch(pdg[0])
      { case  81: Qnum=1;break;
        case  83: Qnum=3;break;
        case -81: Qnum=-1;break;
        case -83: Qnum=-3;break;
         default: Qnum=pdg[0];
      }

      qMass=masses[0];
      if(!qMass) continue;
      sgn=Qnum>0?1:-1;
      aQnum=abs(Qnum);

      if(QCDcorrections && aQnum>3 )
      {  switch(aQnum)
         { case 4: alphaMq=0.39; break; 
           case 5: alphaMq=0.22; break; 
           default:alphaMq=parton_alpha(qMass);
         }  
      } else alphaMq=0;

      loopFF__=NULL;
      if(I==0)
      {
        if(LF) { if(aQnum>3 && aQnum<7 )      loopFF__=LF; 
                   else if(aQnum>6&&abs(color)==3) loopFF__= zeroloopFactor;
               }
      }
      
      for(i=0;i<16;i++) pvect[i]=0;     
      for(i=0;i<4;i++) pvect[4*i]=masses[i];
      cs0= (*cc->interface->sqme)(n,pvect,&err_code);
      loopFF__=NULL; 
      if(spin2==1)  ampl= cs0/(128*masses[0]*masses[0]*masses[1]*masses[1]);
      else if(spin2==0) ampl= cs0/(32*masses[1]*masses[1]);
/* Comment: 128=2*64.  2-because of we have calculated interference term. */       
/* if(ampl) printf("I=%d %E for %s %s -> %s %s\n",I,ampl, names[0],names[1],names[2],names[3]); */     
      if(I)  /* Spin dependent case */
      { if(aQnum<=3 ) 
        {  /* here (ampl/3) because of normalization of SD amplitude */ 
         pA5[0]+=(ampl/3)*wV5P__[aQnum-1]/2;
         nA5[0]+=(ampl/3)*wV5N__[aQnum-1]/2; 
         pA5[1]+=sgn*(ampl/3)*wSM0P[aQnum-1]/2;
         nA5[1]+=sgn*(ampl/3)*wSM0N[aQnum-1]/2;                                     
        }
      }
      else /* Spin independent case */
      { 
          /* scalar SI amplitude */
        if(aQnum>3)
        { int color,spin2,neutral;
          cc->interface->pinfAux(n,1,&spin2,&color,&neutral);
        
          switch(color)
          { case  3:
            case -3:
              switch(spin2)
              { case 0:
                  qcdNLO=1+(25./6.-16./9.) *alphaMq/M_PI;            
                  pA0[0]+=ampl*qcdNLO*wGluP*(1./108.)*MN/qMass/qMass/2;
                  nA0[0]+=ampl*qcdNLO*wGluN*(1./108.)*MN/qMass/qMass/2;
                break;
                case 1:
                  qcdNLO=1+(11./4. -16./9.)*alphaMq/M_PI;
                  pA0[0]+=ampl*qcdNLO*wGluP*(2./27.)*MN/qMass/2 ;
                  nA0[0]+=ampl*qcdNLO*wGluN*(2./27.)*MN/qMass/2;
                break;           
                default:
                printf("2*spin=%d, color=3 - not implemented \n",spin2);
              }
              break;
            case  8:
              switch(spin2)
              { case 0:
                  pA0[0]+=ampl*wGluP*(6./108.)*MN/qMass/qMass/2;
                  nA0[0]+=ampl*wGluN*(6./108.)*MN/qMass/qMass/2;
                break;
                case 1:

                  pA0[0]+=ampl*wGluP*(2.*6./27.)*MN/qMass/2 ;
                  nA0[0]+=ampl*wGluN*(2.*6./27.)*MN/qMass/2;
                break;           
                default:
                printf("2*spin=%d, color=8 - not implemented \n",spin2);
              }
              break;  
           }    
        }else
        {
            pA0[0]+=ampl*wS0P__[aQnum-1]*MN/qMass/2 ; 
            nA0[0]+=ampl*wS0N__[aQnum-1]*MN/qMass/2;
            pA0[1]+=sgn*ampl*wV0P[aQnum-1]/2;
            nA0[1]+=sgn*ampl*wV0N[aQnum-1]/2;
        }        
        if(aQnum<=6 && Twist2On) /* twist-2  SI amplitude  */   
        { int k;
          double g=0,gp,gn;
          double h=masses[0]/100;
          double d4[3]={0.5,-1,0.5};
          h=0.001;
          if(aQnum<=3 || !LF)  /* subtraction of trace */
          { loopFF__=twist2_subtractionFF__; 
            for(k=0;k<3;k++)
            { double E=masses[0]+h*k;
              double P= (k==0)?0: sqrt(E*E-masses[0]*masses[0]);
              pvect[0]=pvect[8]=E;
              pvect[3]=pvect[11]=P;
              p1_p2=E*masses[1];
              g+= (*cc->interface->sqme)(n,pvect,&err_code)*d4[k];
/* printf("E=%E P=%E for %s %s -> %s %s\n",E,P, names[0],names[1],names[2],names[3]) ; */     

            }
            g/=-pow(h*masses[1],2)*128*masses[0]*masses[1]  * 2 ;
/* printf("aQnum=%d g=%E\n",aQnum,g); */  
            qcdNLO=1+(11./4. -16./9.)*alphaMq/M_PI;          
            pA0[0]+=1.5*g*masses[1]*MN*wS0P__[aQnum-1]/2*qcdNLO ;
            nA0[0]+=1.5*g*masses[1]*MN*wS0N__[aQnum-1]/2*qcdNLO ;
            loopFF__=NULL; 
          }
          if(aQnum<6)
          {          
            gp=0;
            gn=0;          
            loopFF__=twist2FF__;  /* twist=spin=2 contribution  */
            for(k=0;k<3;k++)
            { double E=masses[0]+h*k;
              double P= (k==0)?0: sqrt(E*E-masses[0]*masses[0]);
              pvect[0]=pvect[8]=E;
              pvect[3]=pvect[11]=P;
              p1_p2=E*masses[1];
              pdfQnum=aQnum;
              gp+= (*cc->interface->sqme)(n,pvect,&err_code)*d4[k];
              if(aQnum==1) 
              { pdfQnum=2;
                gn+= (*cc->interface->sqme)(n,pvect,&err_code)*d4[k];
              }else if(pdfQnum==2)
              {
                pdfQnum=1;
                gn+= (*cc->interface->sqme)(n,pvect,&err_code)*d4[k];
              }else gn=gp;
          }          
          gp/=pow(h*masses[1],2)*128*masses[0]*masses[1]  * 2 ;
          gn/=pow(h*masses[1],2)*128*masses[0]*masses[1]  * 2 ;
 
          pA0[0]+=1.5*gp*masses[1]*MN/2;
          nA0[0]+=1.5*gn*masses[1]*MN/2; 
          loopFF__=NULL;
          }
        }
      }              
    }    
  }
/*  OnlyTEQ0=0; */
  { double mem;  
    mem=pA0[1]; pA0[1]=pA0[0]-mem; pA0[0]+=mem;
    mem=pA5[1]; pA5[1]=pA5[0]-mem; pA5[0]+=mem;

    mem=nA0[1]; nA0[1]=nA0[0]-mem; nA0[0]+=mem;
    mem=nA5[1]; nA5[1]=nA5[0]-mem; nA5[0]+=mem;     
  }   
  return 0; 
}

/*===== Intergrands for  Fermi nucleus density =====*/

static double R_, p_, a_, C_=1.23, B_=-0.6, A_=0.52;

static double FermiND0(double r){ return r*r/(1+exp((r-R_)/a_));}
static double FermiNDP(double r){ return sin(p_*r)*r/(1+exp((r-R_)/a_));}

void SetFermi(double C,double B, double a){ C_=C; B_=B, A_=a; }

double FermiFF(int A, double Qfermi)
{
  R_=C_*pow(A,1./3.)+B_;
  a_=A_;
  p_=Qfermi;
  return simpson(FermiNDP,0., R_+5, 1.E-4)/simpson(FermiND0,0., R_+5, 1.E-4)/p_;
}

/*===== Maxwell velosity distribution =====*/ 
static double DV_=220,v1_=225.2,vmax_=700,norm_=1.000151;
static double fMaxwell(double v){ return v*fDvMaxwell(v);}


void SetfMaxwell(double DV, double v1, double vmax)
{ DV_=DV, v1_=v1,vmax_=vmax;
  norm_=1/simpson(fMaxwell,0,vmax_+v1_,1.E-5);
}

double fDvMaxwell(double v) 
{  
   double res,vsum,vdif,DV2;

   if(v>vmax_+v1_) return 0;
   DV2=DV_*DV_;
   if(v1_*v<0.001*DV2) res= 4*v*exp((-v*v-v1_*v1_)/(DV2))/(DV2);
   else 
   {   
     vsum=v1_+v;
     if(vsum>vmax_) vsum=vmax_;
     vdif=v1_-v;
     res=(exp(-vdif*vdif/DV2)-exp(-vsum*vsum/DV2))/v1_;
   }
   return res*norm_/(DV_*sqrt(M_PI));
}


/*===== Delta-function velosity distribution =====*/ 
static double deltaV_=220;


void SetfDelta(double V0){ deltaV_=V0;}

double fDvDelta(double v) {  return 1;}


/*==== nucleusRecoil: main functions ================*/
static double eStep=1;
static int    eGrid=200;

void setRecoilEnergyGrid(double step, int nSteps)
{ eStep=step; eGrid=nSteps;} 

static double nucleusRecoil_stat( double rho, double(*vfv)(double),
      int A, int Z, double J,
      double(*S00)(double),double(*S01)(double),double(*S11)(double),
      double  Mwimp, double css,double csv00, double csv01,double csv11,
      double * dNdE)
{
  const double  Kg=1./1.782662E-27;          /* GeV  */
  const double  vC=299792;                   /* km/s */
  const double lDay=60.*60.*24.;             /* sec  */  

  const double step=1.E-6*eStep;              /* 1KeV */

  int i;
  double MA,ffs,Rs,sum; 
  double FFs=1, s00,s01,s11;
  double E0,vmin,vmax;

  s00=S00(0.),s01=S01(0.),s11=S11(0.);

  if(vfv==fDvMaxwell) vmax=vmax_+v1_; 
  else if(vfv==fDvDelta) vmax=deltaV_;
  else vmax=1200.;
  
  MA=0.94*A;


  { double pA0[2],pA5[2],nA0[2],nA5[2];
    double Mr,SCcoeff;

    E0=Mwimp/(Mwimp+MA),E0=2*MA*E0*E0;
    Mr=MA*Mwimp/(MA+Mwimp);

    SCcoeff=4/M_PI*3.8937966E-28*Mr*Mr;

    css*=SCcoeff;
    csv00*=4*M_PI*SCcoeff/(2*J+1);    
    csv11*=4*M_PI*SCcoeff/(2*J+1);
    csv01*=4*M_PI*SCcoeff/(2*J+1);
  }
/*printf("css=%E, csv00=%E, csv11=%E, csv01=%E\n",css,csv00,csv11,csv01);*/

  if(A>1)
  {
    Rs=C_*pow(A,1./3.)+B_;
    R_=Rs; a_=A_; 
    ffs=simpson(FermiND0,0., R_+10*0.5, 1.E-4);
  }
  
  for(i=0,sum=0;i<eGrid;i++)
  { double E=i*step;

    vmin=sqrt(E/E0)*vC;
    
    if(vmin>=vmax)  { dNdE[i]=0;continue;}
    if(i && A>1)
    { double p=sqrt(E*MA*2)/0.197327;
      if(J){s00=S00(p); s01=S01(p); s11=S11(p);}
      p_=p;    
      R_=Rs; a_=A_; 
      FFs=simpson(FermiNDP,0., R_+10*0.5, 1.E-4)/p_/ffs;      
    }  
    
    dNdE[i]=FFs*FFs*css;
    if(J) dNdE[i]+=(csv00*s00+csv01*s01+csv11*s11);
    if(vfv==fDvDelta) dNdE[i]/=deltaV_;
    else dNdE[i]*=simpson(vfv,vmin,vmax,1.E-4);
    dNdE[i]*=1/E0*1.E5*vC*vC*(rho/Mwimp)*lDay*(Kg/MA)*1.E-6;
    if(i==0 ||i==eGrid-1) sum+=dNdE[i]/2; else sum+=dNdE[i];
  }
  return sum*eStep;
}

static double S00_0,S01_0,S11_0, RS_;


static double S00_(double p)
{ double r;
  r=exp(-p*p*RS_*RS_/4);
  return S00_0*r;
}

static double S01_(double p)
{ double r;
  r=exp(-p*p*RS_*RS_/4);
  return S01_0*r;
}

static double S11_(double p)
{ double r;
  r=exp(-p*p*RS_*RS_/4);
  return S11_0*r;
}


static double nucleusRecoil0_stat(double rho, double(*vfv)(double),
int A, int Z, double J,double Sp,double Sn,
double Mwimp, double css, double csv00, double csv01, double csv11,
double * dNdE)
{
  if(J)
  { double A3=pow(A,1./3.);
    RS_= 1.7*A3 -0.28 - 0.78*(A3-3.8 + sqrt((A3-3.8)*(A3-3.8) +0.2) );
      
    S00_0=  (Sp+Sn)*(Sp+Sn)*(2*J+1)*(J+1)/(4*M_PI*J);
    S11_0=  (Sp-Sn)*(Sp-Sn)*(2*J+1)*(J+1)/(4*M_PI*J);
    S01_0=2*(Sp+Sn)*(Sp-Sn)*(2*J+1)*(J+1)/(4*M_PI*J);
  }  
  return nucleusRecoil_stat(rho,vfv,A,Z,J,S00_,S01_,S11_,
                        Mwimp,css,csv00,csv01,csv11,dNdE);
}


 
double nucleusRecoil( double rho, double(*vfv)(double),
      int A, int Z, double J,
      double(*S00)(double),double(*S01)(double),double(*S11)(double),
      double(*LF)(double,double,double,double),
      double * dNdE)
{
  double Mwimp,css,csv00,csv01,csv11; 
  double pA0[2],pA5[2],nA0[2],nA5[2];
  int i;
  
  nucleonAmplitudes(LF,pA0,pA5,nA0,nA5);
  Mwimp=lopmass_();
  for(i=0,css=0,csv00=0,csv01=0,csv11=0;i<2;i++)
  { double AS=Z*pA0[i]+(A-Z)*nA0[i];
    double AVplus =pA5[i]+nA5[i];
    double AVminus=pA5[i]-nA5[i];
    css+=AS*AS/2;
    csv00+=AVplus*AVplus/2;
    csv11+=AVminus*AVminus/2;
    csv01+=AVplus*AVminus/2;
  }  
  return nucleusRecoil_stat(rho,vfv,A,Z,J,S00,  S01,  S11,
                                Mwimp,css,csv00,csv01,csv11,dNdE);
} 

double nucleusRecoilAux( 
      double rho, double(*vfv)(double),
      int A, int Z, double J,
      double(*S00)(double),double(*S01)(double),double(*S11)(double),
      double Mwimp, 
      double  siP, double siN, double sdP,  double sdN,
      double * dNdE)
{

  double AS,AVplus, AVminus, css,csv00,csv01,csv11; 
  
  double LmbdP, XiP,LmbdN, XiN;
  double MN=0.939;
  double Mr=MN*Mwimp/(MN+Mwimp);
  double sCoeff= 2*sqrt(3.8937966E8/M_PI)*Mr;
   
  LmbdP=sqrt(fabs(siP))/sCoeff; if(siP<0) LmbdP*=-1;
  XiP  =sqrt(fabs(sdP)/3)/sCoeff; if(sdP<0) XiP*=-1;
  LmbdN=sqrt(fabs(siN))/sCoeff; if(siN<0) LmbdN*=-1;
  XiN  =sqrt(fabs(sdN)/3)/sCoeff; if(sdN<0) XiN*=-1;
   
  AS=Z*LmbdP+(A-Z)*LmbdN;
  AVplus =XiP+XiN;
  AVminus=XiP-XiN;
  css=AS*AS;
  csv00=AVplus*AVplus;
  csv11=AVminus*AVminus;
  csv01=AVplus*AVminus;

  return nucleusRecoil_stat(rho,vfv,A,Z,J,S00,  S01,  S11,
                                Mwimp,css,csv00,csv01,csv11,dNdE);
} 

double nucleusRecoil0Aux( 
      double rho, double(*vfv)(double),
      int A, int Z, double J,
      double Sp,double Sn,
      double Mwimp, 
      double siP, double siN, double sdP,   double sdN,
      double * dNdE)
{
  double AS,AVplus, AVminus, css,csv00,csv01,csv11; 
  double LmbdP,XiP,LmbdN,XiN;
  double MN=0.939;
  double Mr=MN*Mwimp/(MN+Mwimp);
  double sCoeff= 2*sqrt(3.8937966E8/M_PI)*Mr;
 
  LmbdP=sqrt(fabs(siP))/sCoeff; if(siP<0) LmbdP*=-1;
  XiP  =sqrt(fabs(sdP)/3)/sCoeff; if(sdP<0) XiP*=-1;
  LmbdN=sqrt(fabs(siN))/sCoeff; if(siN<0) LmbdN*=-1;
  XiN  =sqrt(fabs(sdN)/3)/sCoeff; if(sdN<0) XiN*=-1;

  AS=Z*LmbdP+(A-Z)*LmbdN;
  AVplus =XiP+XiN;
  AVminus=XiP-XiN;
  css=AS*AS;
  csv00=AVplus*AVplus;
  csv11=AVminus*AVminus;
  csv01=AVplus*AVminus;

  return nucleusRecoil0_stat(rho,vfv,A,Z,J,Sp,Sn,
                                Mwimp,css,csv00,csv01,csv11,dNdE);
} 


double nucleusRecoil0(double rho, double(*vfv)(double),
int A, int Z, double J,double Sp,double Sn,double (*LF)(double,double,double,double),
double * dNdE)
{
  if(J)
  { double A3=pow(A,1./3.);
    RS_= 1.7*A3 -0.28 - 0.78*(A3-3.8 + sqrt((A3-3.8)*(A3-3.8) +0.2) );
      
    S00_0=  (Sp+Sn)*(Sp+Sn)*(2*J+1)*(J+1)/(4*M_PI*J);
    S11_0=  (Sp-Sn)*(Sp-Sn)*(2*J+1)*(J+1)/(4*M_PI*J);
    S01_0=2*(Sp+Sn)*(Sp-Sn)*(2*J+1)*(J+1)/(4*M_PI*J);
  }  
  return nucleusRecoil(rho,vfv,A,Z,J,S00_,S01_,S11_,LF,dNdE);
}



/*====== Auxilarry service functions ======*/
int displayRecoilPlot(double * tab, char * text, double  E1, double E2)
{ 
  FILE*f;
  int i;
  int i1=(E1/eStep), i2=(E2/eStep);

  if(E1<0 || E1>=E2-eStep|| i2>eGrid-1|| i2-i1>299  ) return 1;
  f=fopen("_plot_","w");
  if(!f) return 2;
  fprintf(f,"#type 0 %%curve\n");
  fprintf(f,"#title %s\n",text);
  fprintf(f,"#yName dN/dE\n");
  fprintf(f,"#xMin %f\n", i1*eStep);
  fprintf(f,"#xMax %f\n", i2*eStep);
  fprintf(f,"#xDim %d\n",i2-i1+1);
  fprintf(f,"#xName  E[KeV]\n");

  for(i=i1;i<=i2;i++) fprintf(f,"%E\n",tab[i]);  
  fclose(f);
  system("trap \"rm -f _plot_\" 0 1 2 3 9; ../CalcHEP_src/bin/plot_view < _plot_ &");
  return 0;
}

double cutRecoilResult(double *tab, double E1, double E2)
{ int i,i1,i2;
  double sum,dx;
  if(E2>eStep*(eGrid-1)) 
  { 
    printf("cutRecoilResult:  Maximal bound %.3E is larger than data limit %.3E\n",
       E2,eStep*(eGrid-1));
    E2=eStep*(eGrid-1); 
  } 
  if(E1<0)   E1=0;
  if(E1>=E2) return 0;
  
  i1=E1/eStep;  i2=E2/eStep;
  if(i1<E1/eStep) i1++;
  
  for(i=i1,sum=0;i<=i2;i++) sum+=tab[i];
  sum-=(tab[i1]+tab[i2])/2;
  dx=i1 -E1/eStep;
  if(dx>0) sum+=tab[i1]*dx + (tab[i1-1]-tab[i1])*dx*dx/2  ;
  dx=E2/eStep-i2;
  if(dx) sum+=tab[i2]*dx  +(tab[i2+1]-tab[i2])*dx*dx/2;
  
  return sum*eStep;
}

int PlotSS(double (*fun)(double), int A, char * title, double Emax)
{ 
  FILE*f=fopen("_plot_","w");
  int i;
  double M=0.940*A;

  if(!f) return 2;
  fprintf(f,"#type 0 %%curve\n");
  fprintf(f,"#title  %s  A=%d \n",title,A);
  fprintf(f,"#yName Str.Func\n");
  fprintf(f,"#xMin 0\n");
  fprintf(f,"#xMax %E\n",Emax);
  fprintf(f,"#xDim %d\n",(int)Emax+1);
  fprintf(f,"#xName  E[KeV]\n");

  for(i=0;i<(int)Emax+1;i++) fprintf(f,"%E\n",fun(1.E-3*sqrt(2*M*i)/0.197327));
  fclose(f);
  system(" trap \"rm -f _plot_\" 0 1 2 3 9; ../CalcHEP_src/bin/plot_view < _plot_ &");
 
  return 0;
}

int Plot3SS(double xiP,double xiN,
  double (*S00)(double), double (*S01)(double),double (*S11)(double),
    int A, double J, char* title, double Emax)
{ 
  FILE*f=fopen("_plot_","w");
  int i;
  double M=0.940*A;
  double Sp,Sn;
  
  Sp= sqrt( M_PI*J/(2*J+1)/(J+1)*(S00(0.)+S11(0.)+S01(0.)) );
  Sn= sqrt( M_PI*J/(2*J+1)/(J+1)*(S00(0.)+S11(0.)-S01(0.)) );
  if(S00(0.)-S11(0.) < 0) { if(Sp>Sn)Sn*=-1; else Sp*=-1;}
  
  if(!f) return 2;
  fprintf(f,"#type 0 %%curve\n");
  fprintf(f,"#title xiP=%.2f xiN=%.2f %s (A=%d J=%.1f) Sp=%.2E Sn=%.2E\n",
                  xiP,xiN,title,A,J,Sp,Sn);
  fprintf(f,"#yName Str.Func\n");
  fprintf(f,"#xMin 0\n");
  fprintf(f,"#xMax %E\n",Emax);
  fprintf(f,"#xDim %d\n",(int)Emax+1);
  fprintf(f,"#xName  E[KeV]\n");

  for(i=0;i<(int)Emax+1;i++)
  { double q=1.E-3*sqrt(2*M*i)/0.197327;
    fprintf(f,"%E\n", (xiP+xiN)*(xiP+xiN)*S00(q)
                     +(xiP+xiN)*(xiP-xiN)*S01(q)
                     +(xiP-xiN)*(xiP-xiN)*S11(q));  
  }   
  fclose(f);
  system(" trap \"rm -f _plot_\" 0 1 2 3 9; ../CalcHEP_src/bin/plot_view < _plot_ &");
 
  return 0;
}

int Plot3SS0(double xiP,double xiN, double Sp, double Sn,
    int A, double J, char* title, double Emax)
{
  { double A3=pow(A,1./3.);

    RS_= 1.7*A3 -0.28 - 0.78*(A3-3.8 + sqrt((A3-3.8)*(A3-3.8) +0.2) );      
    S00_0=  (Sp+Sn)*(Sp+Sn)*(2*J+1)*(J+1)/(4*M_PI*J);
    S11_0=  (Sp-Sn)*(Sp-Sn)*(2*J+1)*(J+1)/(4*M_PI*J);
    S01_0=2*(Sp+Sn)*(Sp-Sn)*(2*J+1)*(J+1)/(4*M_PI*J);
  }  
  return Plot3SS(xiP,xiN, S00_,S01_,S11_,A,J,title,Emax);
}

/* =====  Nuclear Form Factors for vector currents =====*/

#define beta0(u,upi) ( ((upi)*(upi)/((u)+(upi))/((u)+(upi))-1.)/3.)
#define pimass (0.135/0.197327)
#define P01(u)  ( 0.1145*(u)*(u) - 0.6667*(u) + 1.)
#define P21(u)  (-0.0026*(u)*(u) + 0.0100*(u) )
#define Q01(u)  ( 0.1088*(u)*(u) - 0.6667*(u) + 1.)
#define Q21(u)  ( 0.0006*(u)*(u) + 0.0041*(u) )

#define J 0.5
#define A 19

double S00F19(double q)
{ double q_fm=q;
  double bsq=pow(A,1./3.), u=0.5*bsq*q_fm*q_fm;
  double upi=0.5*bsq*pimass*pimass;
  double p01=P01(u), p21=P21(u); 
  double be0=beta0(u,upi), be2=2*be0,be02=-sqrt(2.)*be0;
  
  return ((2*J+1)/16/M_PI)*2.610*exp(-u)*(p01*p01*(1+be0)+p21*p21*(1+be2)
  -2*be02*p01*p21);
} 

double S11F19(double q)
{ double q_fm=q;
  double bsq=pow(A,1./3.), u=0.5*bsq*q_fm*q_fm;
  double upi=0.5*bsq*pimass*pimass;
  double q01=Q01(u), q21=Q21(u); 
  double be0=beta0(u,upi), be2=2*be0,be02=-sqrt(2.)*be0;
  
  return ((2*J+1)/16/M_PI)*2.807*exp(-u)*(q01*q01*(1+be0)+q21*q21*(1+be2)
  -2*be02*q01*q21);
} 

double S01F19(double q)
{ double q_fm=q;
  double bsq=pow(A,1./3.), u=0.5*bsq*q_fm*q_fm;
  double upi=0.5*bsq*pimass*pimass;
  double p01=P01(u), p21=P21(u); 
  double q01=Q01(u), q21=Q21(u); 
  double be0=beta0(u,upi), be2=2*be0,be02=-sqrt(2.)*be0;
  
  return ((2*J+1)/8/M_PI)*2.707*exp(-u)*(p01*q01*(1+be0)+p21*q21*(1+be2)
  -be02*(p01*q21+p21*q01));
} 


#undef P01
#undef P21
#undef Q01
#undef Q21
#undef J 
#undef A 

#define P01(u)  ( 0.2843*(u)*(u) - 0.6667*(u)  + 1.)
#define P21(u)  ( -0.0567*(u)*(u) + 0.4566*(u))
#define Q01(u)  ( 0.2710*(u)*(u) - 0.6667*(u) + 1.)
#define Q21(u)  ( -0.0621*(u)*(u) + 0.4680*(u) )

#define J 0.5
#define A 29

double S00Si29(double q)
{ double q_fm=q;
  double bsq=pow(A,1./3.), u=0.5*bsq*q_fm*q_fm;
  double upi=0.5*bsq*pimass*pimass;
  double p01=P01(u),p21=P21(u); 
  double be0=beta0(u,upi), be2=2*be0,be02=-sqrt(2.)*be0;
  
  return ((2*J+1)/16/M_PI)*0.208*exp(-u)*(p01*p01*(1+be0)+p21*p21*(1+be2)
  -2*be02*p01*p21);
} 

double S11Si29(double q)
{ double q_fm=q;
  double bsq=pow(A,1./3.), u=0.5*bsq*q_fm*q_fm;
  double upi=0.5*bsq*pimass*pimass;
  double q01=Q01(u), q21=Q21(u); 
  double be0=beta0(u,upi), be2=2*be0,be02=-sqrt(2.)*be0;
  
  return ((2*J+1)/16/M_PI)*0.220*exp(-u)*(q01*q01*(1+be0)+q21*q21*(1+be2)
  -2*be02*q01*q21);
} 

double S01Si29(double q)
{ double q_fm=q;
  double bsq=pow(A,1./3.), u=0.5*bsq*q_fm*q_fm;
  double upi=0.5*bsq*pimass*pimass;
  double p01=P01(u), p21=P21(u); 
  double q01=Q01(u), q21=Q21(u); 
  double be0=beta0(u,upi), be2=2*be0,be02=-sqrt(2.)*be0;
    
  return ((2*J+1)/8/M_PI)*(-sqrt(0.208*0.220))*exp(-u)*(p01*q01*(1+be0)+p21*q21*(1+be2)
  -be02*(p01*q21+p21*q01));
} 

#undef P01
#undef P21
#undef Q01
#undef Q21
#undef J 
#undef A 

#define P01(u)  ( 0.0477*u*u - 0.6667*u + 1)
#define P21(u)  (-0.0177*u*u + 0.1048*u)
#define P23(u)  (-0.0767*u*u + 0.6092*u)
#define P43(u)  ( 0.0221*u*u)

#define Q01(u)  ( 0.0465*u*u - 0.6667*u + 1)
#define Q21(u)  (-0.0349*u*u + 0.1494*u)
#define Q23(u)  (-0.0894*u*u + 0.7405*u)
#define Q43(u)  ( 0.0287*u*u)

#define J 1.5
#define A 23
/*#define b 1.69*/

double S00Na23A(double q)
{ double q_fm=q;
  double bsq=pow(A,1./3.), u=0.5*bsq*q_fm*q_fm;
/*  double bsq=b*b, u=0.5*bsq*q_fm*q_fm; */
  double p01=P01(u),p21=P21(u),p23=P23(u),p43=P43(u); 
  return ((2*J+1)/16/M_PI)*(0.478)*exp(-u)*(p01*p01+p21*p21+ p23*p23+ p43*p43);
} 

double S11Na23A(double q)
{ double q_fm=q;
  double bsq=pow(A,1./3.), u=0.5*bsq*q_fm*q_fm; 
/*  double bsq=b*b, u=0.5*bsq*q_fm*q_fm;*/ 
  double q01=Q01(u),q21=Q21(u),q23=Q23(u),q43=Q43(u); 
  return ((2*J+1)/16/M_PI)*(0.346)*exp(-u)*(q01*q01+q21*q21+q23*q23+q43*q43);
} 

double S01Na23A(double q)
{ double q_fm=q;
  double bsq=pow(A,1./3.), u=0.5*bsq*q_fm*q_fm; 
/*  double bsq=b*b, u=0.5*bsq*q_fm*q_fm;*/
  double p01=P01(u),p21=P21(u),p23=P23(u),p43=P43(u);
  double q01=Q01(u),q21=Q21(u),q23=Q23(u),q43=Q43(u);
  return ((2*J+1)/8/M_PI)*(0.406)*exp(-u)*(p01*q01+p21*q21+p23*q23+p43*q43);
} 

#undef P01
#undef P21
#undef P23
#undef P43
#undef Q01
#undef Q21
#undef Q23
#undef Q43
#undef J 
#undef A 
#undef b

/*--------- exponential parametrization --------*/

double static ss_exp_aux(double q, double A, double a,double b)
{   double q_fm=q, y=0.25*pow(A,1./3.)*q_fm*q_fm;
    return a*exp(-b*y);
}    

/*================= Si29A ===================== */ 
double S00Si29A(double q)  {return ss_exp_aux(q, 29,0.00818,4.428);}   
double S11Si29A(double q)  {return ss_exp_aux(q, 29,0.00818*1.06,6.264);}
double S01Si29A(double q)  {return ss_exp_aux(q, 29,0.00818*(-2.06),5.413);}

/*================= Ge73A ===================== */
double S00Ge73A(double q){return ss_exp_aux(q,73,0.20313*1.102,7.468);}
double S11Ge73A(double q){return ss_exp_aux(q,73,0.20313,8.856);}
double S01Ge73A(double q){return ss_exp_aux(q,73,0.20313*(-2.099),8.191);} 
  
/*-------polinomial parametrization ----------------*/
 
static double SSxxYYzz(double q,double * data, int nData,int A,int expKey,double yN,double y_max)
{ int i;
  double q_fm=q, 
/*  y=0.25*pow((double)A,1./3.)*q_fm*q_fm, y_i,res;*/
    y=0.25*pow(A,1./3.)*(41.467/(45-25/pow(A,1./3.)))*q_fm*q_fm, y_i,res;

  if(y>y_max) return 0.;
  for(i=1,res=data[0],y_i=1;i<nData;i++) {y_i*=y; res+=y_i*data[i]; }
  res+=yN/(1+y);
  if(expKey) res*=exp(-2*y);
  return res;
} 

/*========================= Na23 ===================*/
double S00Na23(double q)
{ 
  double data[4]={0.0380,-0.1743,0.3783,-0.3430};
  return SSxxYYzz(q,data,4,23,0,0.,0.2);
}
double S01Na23(double q)
{ 
  double data[4]={0.0647, -0.3503,0.9100,-0.9858};
  return SSxxYYzz(q,data,4,23,0,0.,0.2);
}
double S11Na23(double q)
{ 
  double data[4]={0.0275,-0.1696,0.5077,-0.6180};
  return SSxxYYzz(q,data,4,23,0,0.,0.2);
}                            

/* ========================= Al27 ==================*/
double S00Al27(double q)
{ 
  double data[4]={0.0930,-0.4721,1.0600,-1.0115};
  return SSxxYYzz(q,data,4,27,0,0.,0.3);
}
double S11Al27(double q)
{ 
  double data[4]={0.0657,-0.4498,1.3504,-1.6851};
  return SSxxYYzz(q,data,4,27,0,0.,0.3);
}
double S01Al27(double q)
{ 
  double data[4]={0.1563,-0.9360,2.4578,-2.7262};
  return SSxxYYzz(q,data,4,27,0,0,0.3);
}                            
/*============================== K39 =================*/

double S00K39(double q)
{ 
  double data[5]={0.0094999,-0.0619718,0.162844,-0.194282,0.0891054};
  return SSxxYYzz(q,data,5,39,0,0.,0.6);
}
double S11K39(double q)
{ 
  double data[5]={0.0298127,-0.2176360,0.623646,-0.814418,0.4050270};
  return SSxxYYzz(q,data,5,39,0,0.,0.6);
}
double S01K39(double q)
{ 
  double data[5]={0.0332044,-0.2319430,0.638528,-0.798523,0.3809750};
  return SSxxYYzz(q,data,5,39,0,0.,0.6);
}                            


/*=========================== Ge73 =====================*/

double S00Ge73(double q)
{ 
  double data[7]={0.1606, -1.1052,3.2320,-4.9245,4.1229,-1.8016,0.3211};
  return SSxxYYzz(q,data,7,73,0,0.,1.3);
}
double S11Ge73(double q)
{ 
  double data[7]={0.1164,-0.9228,2.9753,-4.8709,4.3099,-1.9661,0.3624 };
  return SSxxYYzz(q,data,7,73,0,0.,1.3);
}
double S01Ge73(double q)
{ 
  double data[7]={-0.271006,2.018922,-6.226466,9.860608,-8.502157,3.800620,-0.689352};
  return SSxxYYzz(q,data,7,73,0,0.,1.3);
}

/* ========================= Nb93 ==================== */
double S00Nb93(double q)
{
 double  q2[9]= {0,    0.002,0.004,0.006,0.008,0.01, 0.014,0.02, 0.03};
 double  s00[9]={0.284,0.19, 0.122,0.082,0.058,0.039,0.024,0.021,0.02};
 double x=q*0.197327;
 x*=x;
 if(x>0.03) return 0;
 return polint3(x,9,q2,s00);
}
double S11Nb93(double q)
{
 double q2[9]= {0,   0.002,0.004,0.006,0.008,0.01, 0.014, 0.02, 0.03};
 double s11[9]={0.14,0.085,0.053,0.034,0.02, 0.012,0.0052,0.003,0.0025};
 double x=q*0.197327;
 x*=x;
  if(x>0.03) return 0;
 return  polint3(x,9,q2,s11);
}
double S01Nb93(double q)
{
 double q2[9]= {0,  0.002,0.004,0.006,0.008,0.01, 0.014,0.02, 0.03};
 double s01[9]={0.4,0.26, 0.161,0.105,0.07, 0.046,0.029,0.026,0.024};
 double x=q*0.197327;
 x*=x;
 if(x>0.03) return 0;
 return polint3(x,9,q2,s01);
}

/*============================ Xe131B =======================*/
double S00Xe131B(double q)
{
 double q2[11]= {0,0.0025,0.005,0.01,0.015,0.02,0.025,0.03,0.04,0.05,0.06};
 double  s00[11]={0.04,0.0215,0.014,0.01,0.009,0.008,0.0075,0.0066,0.005,0.0035,0.0017};
 double x=q*0.197327;
 x*=x;
 if(x>0.06) return 0;
 return polint3(x,11,q2,s00);
}
double S11Xe131B(double q)
{
 double q2[11]= {0,0.0025,0.005,0.01,0.015,0.02,0.025,0.03,0.04,0.05,0.06};
 double s11[11]={0.020,0.009,0.006,0.004,0.003,0.0027,0.0025,0.0023,0.0019,0.0015,0.001};
 double x=q*0.197327;
 x*=x;
 if(x>0.06) return 0;
 return polint3(x,11,q2,s11);
}
double S01Xe131B(double q)
{
 double q2[11]= {0,0.0025,0.005,0.01,0.015,0.02,0.025,0.03,0.04,0.05,0.06};
 double s01[11]={ -0.056,-0.028,-0.019,-0.013,-0.01,-0.009,-0.008,-0.007,-0.005,-0.003,-0.001 };
 double x=q*0.197327;
 x*=x;
 if(x>0.06) return 0;
 return polint3(x,11,q2,s01);
}


/*============================ Te125 =======================*/
double S00Te125(double q)
{
 double data[9]={0.03971,-0.19610, 0.47265,-0.65023, 0.54193,-0.26456, 0.07489,-0.01146,0.00075};
 return SSxxYYzz(q,data,9,125,1,0.,10.);
}
double S11Te125(double q)
{
 double data[9]={ 0.03922,-0.22938,0.62215,-0.92253,0.78465,-0.38245,0.1057,-0.01542,0.00093 };
 return SSxxYYzz(q,data,9,125,1,0.,10.);
}
double S01Te125(double q)
{
 double data[9]={-0.07894, 0.42738,-1.09331, 1.55324,-1.28933, 0.61844,-0.16964, 0.02481,-0.00152  }; 
 return SSxxYYzz(q,data,9,125,1,0.,10.);
}

/*============================ I127 =======================*/
double S00I127(double q)
{
 double data[9]={0.0983,-0.4891,1.1402,-1.4717,1.1717,-0.5646,0.1583,-0.0239,0.0015};
 return SSxxYYzz(q,data,9,127,1,0.,10.);
}
double S11I127(double q)
{
 double data[9]={0.0366,-0.1950,0.5049,-0.7475,0.7043,-0.3930,0.1219,-0.0192,0.0012};
 return SSxxYYzz(q,data,9,127,1,0.,10.);
}
double S01I127(double q)
{
 double data[9]={0.1199,-0.6184,1.5089,-2.0737,1.7731,-0.9036,0.2600,-0.0387,0.0024}; 
 return SSxxYYzz(q,data,9,127,1,0.,10.);
}

/*============================ Xe129 =======================*/
double S00Xe129(double q)
{
 double data[9]={0.07132,-0.34478, 0.75590,-0.93345, 0.69006,-0.30248, 0.07653,-0.01032, 0.00057};
 return SSxxYYzz(q,data,9,129,1,0.,10.);
}
double S11Xe129(double q)
{
 double data[9]={-2.05825, 1.80756,-1.27746, 0.65459,-0.22197, 0.04546,-0.00427,-0.00014, 0.00004};
 return SSxxYYzz(q,data,9,129,1,2.11016,10.);
}
double S01Xe129(double q)
{
 double data[9]={-0.12166, 0.64435,-1.52732, 2.02061,-1.57689, 0.72398,-0.19040, 0.02638,-0.00149}; 
 return SSxxYYzz(q,data,9,129,1,0.,10.);
}

/*============================ Xe131 =======================*/
double S00Xe131(double q)
{
 double data[9]={ 0.02964,-0.13343, 0.37799,-0.57961, 0.57890,-0.34556, 0.11595,-0.02012, 0.00142};
 return SSxxYYzz(q,data,9,131,1,0.,10.);
}
double S11Xe131(double q)
{
 double data[9]={ 0.02510,-0.13772,0.36661,-0.53851, 0.49255,-0.26990, 0.08369,-0.01340, 0.00087};
 return SSxxYYzz(q,data,9,131,1,0.,10.);
}
double S01Xe131(double q)
{
 double data[9]={-0.05455, 0.27176,-0.72302, 1.05450,-0.97133, 0.53842,-0.16899, 0.02742,-0.00181}; 
 return SSxxYYzz(q,data,9,131,1,0.,10.);
}

/*============================ Te125A =======================*/
 
double S00Te125A(double q)
{double data[9]={0.04960,-0.24777,0.54766,-0.66553,0.47462,-0.19944,0.04819,-0.00616,0.00032};
 return SSxxYYzz(q,data,9,125,1,0.,10.);
}
double S11Te125A(double q)
{
 double data[9]={ -1.92941,1.68075,-1.16336,0.58650,-0.20730,0.05141,-0.00870,0.00087,0.00004};
 return SSxxYYzz(q,data,9,125,1,1.97923,10.);
}
double S01Te125A(double q)
{
 double data[9]={-0.09939,0.54303,-1.28816,1.67206,-1.26883,0.56728,-0.14545,0.01959,-0.00107}; 
 return SSxxYYzz(q,data,9,125,1,0.,10.);
}

/*============================ I127 =======================*/
double S00I127A(double q)
{
 double data[9]={0.1166,-0.5721,1.3380,-1.7252,1.3774,-0.6700,0.1905,-0.0292,0.0019};
 return SSxxYYzz(q,data,9,127,1,0.,10.);
}
double S11I127A(double q)
{
 double data[9]={0.0563,-0.3038,0.7948,-1.1703,1.0637,-0.5713,0.1722,-0.0266,0.0017};
 return SSxxYYzz(q,data,9,127,1,0.,10.);
}
double S01I127A(double q)
{
 double data[9]={0.1621,-0.8363,2.0594,-2.8319,2.3973,-1.2121,0.3486,-0.0522,0.0032};
 return SSxxYYzz(q,data,9,127,1,0.,10.);
}

/*============================ Xe129A =======================*/

double S00Xe129A(double q)
{
 double data[9]={0.04649,-0.22551,0.49905,-0.62244,0.46361,-0.20375,0.05109,-0.00671,0.00036};
 return SSxxYYzz(q,data,9,129,1,0.,10.);
}
double S11Xe129A(double q)
{
 double data[9]={-1.28214,1.09276,-0.71295,0.31489,-0.08351,0.01059,0.00023,-0.00024,0.00002};
 return SSxxYYzz(q,data,9,129,1,1.32136,10.);
}
double S01Xe129A(double q)
{
 double data[9]={-0.08538,0.45343,-1.06546,1.38670,-1.05940,0.47576,-0.12208,0.01643,-0.00089};
 return SSxxYYzz(q,data,9,129,1,0.,10.);
}

/*============================ Xe131A =======================*/
double S00Xe131A(double q)
{
 double data[9]={0.02773,-0.12449,0.32829,-0.48140,0.47565,-0.28518,0.09682,-0.01710,0.00124};
 return SSxxYYzz(q,data,9,131,1,0.,10.);
}
double S11Xe131A(double q)
{
 double data[9]={0.02234,-0.12206,0.31949,-0.46695,0.42877,-0.23679,0.07408,-0.01197,0.00079};
 return SSxxYYzz(q,data,9,131,1,0.,10.);
}
double S01Xe131A(double q)
{
 double data[9]={-0.04978,0.24725,-0.63231,0.89642,-0.81645,0.45235,-0.14267,0.02335,-0.00156};
 return SSxxYYzz(q,data,9,131,1,0.,10.);
}

/*============================ Pb207 =======================*/

static double _u_;
static double  _fPb_integrand(double x)
{ double u=_u_*x;
  return exp(-u/2)*(1.-u*(5./3.-u*(17./15.-u*(31./105.-u*( 9./280.-u*1./840.)))));
}

double S00Pb207(double q)
{
 _u_=0.5*pow(207.,1./3.)*(q)*(q);
 return   ((2*0.5+1)/16/M_PI)*0.305*simpson(_fPb_integrand,0.,1.,1.E-4);
}
double S11Pb207(double q){return S00Pb207(q)*0.231/0.305;}
double S01Pb207(double q){return  -2*S00Pb207(q)*0.266/0.305;}

