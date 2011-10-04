/*
 Copyright (C) 1997, Alexander Pukhov 
*/
#include"interface.h"
#include"sf_epa.h"
#include"sf_isr.h"
#include"sf_lsr.h"
#include"sf_pdt.h"
#include"sf_pdf.h"
#include"sf_epp.h"
#include"strfun.h"
#include"subproc.h"
#include"crt_util.h"
#include"syst.h"
#include"alphas2.h"

#define MAXFUN 6
#define FUNLEN 40  

int sf_num[2]={0,0}; /* switched off */
double sf_mass[2]={0,0};
double sf_be[2]={0,0};


static struct 
{  int (*myParticle)(long);
   void(*fullName)(int,char *); 
   int (*readName)(int,char*);
   int (*menu)(int);
   int (*mc_code)(int);
   int(*init)(int,double*,double*);
   double(*val)(int,double,double);
}  strFun [MAXFUN]={ {p_epa__,n_epa__,r_epa__,m_epa__,mc_epa__, i_epa__,c_epa__},
                     {p_lsr__,n_lsr__,r_lsr__,m_lsr__,mc_lsr__, i_lsr__,c_lsr__},              
                     {p_isr__,n_isr__,r_isr__,m_isr__,mc_isr__, i_isr__,c_isr__},
                     {p_epp__,n_epp__,r_epp__,m_epp__,mc_epp__ ,i_epp__,c_epp__},
                     {p_pdt  ,n_pdt  ,r_pdt  ,m_pdt  ,mc_pdt  ,be_pdt  ,c_pdt},
                     {p_pdf  ,n_pdf  ,r_pdf  ,m_pdf , mc_pdf  ,be_pdf  ,c_pdf}
                   };

int initStrFun(int mode )
{
  int l,i;
  int returnCode=0;

  if(nin_int!=2) return returnCode;
   
  if(mode<=0||mode>2)
  {
    sf_alpha=NULL;
    if(initStrFun(1)||initStrFun(2)) returnCode=2; 
    return returnCode;   
  }

  i=mode-1;
  l=sf_num[i];  
  if(l)
  {  long N;
     pinf_int(Nsub,i+1,NULL,&N); 
     l--;
     if(!strFun[l].myParticle(N)
         ||!strFun[l].init(i+1,sf_be+i,sf_mass+i))
     {  char txt[60];
        sprintf(txt,"%d-th Stucture function is switched OFF",i+1);
        sf_num[i]=0;
        messanykey(10,15,txt); 
        returnCode=2;
      }
  } else sf_be[i]=0;
    
  return returnCode;  
}


void strFunName(int i, char * mess)
{
  if(sf_num[i-1]) strFun[sf_num[i-1]-1].fullName(i,mess); 
  else strcpy(mess,"OFF");
}

int sf_menu(int i)
{
    int  k;
    char name[STRSIZ];
    int  nfun[MAXFUN];
    long N;
    char strmen[2+MAXFUN*(FUNLEN+1)];

    void * pscr =NULL;
    int mode,l;

    strmen[0]=FUNLEN+1;

    pinf_int(Nsub,i,NULL,&N);
    
    sprintf(strmen+1," %-*.*s",FUNLEN,FUNLEN,"OFF");

    k = 0;          
    for(l=0;l<MAXFUN;l++)
    {    
       if ( strFun[l].myParticle(N)  ) 
       {  
          nfun[k++] = l;
          strFun[l].fullName(i, name); 
	  sprintf(strmen+1+(FUNLEN+1)*k," %-*.*s",FUNLEN,FUNLEN,name);
       }
    }
    if(!k) 
    { messanykey(15,15,"Structure functions for this particle\n"
                       "are not known\n");
      return 0;
    }

    menu1(77-FUNLEN,7,"",strmen, "n_strfun", &pscr, &mode);
    if (mode == 0) return 0;
    put_text(&pscr);
    if (mode == 1) sf_num[i-1]=0;
    else
    {  int ok=strFun[nfun[mode - 2]].menu(i);
       if(ok) sf_num[i-1]=nfun[mode -2]+1; else{ sf_num[i-1]=0; return 0;}
    }
    return 1;
} /* sf_menu__ */

double strfun_(int i, double x,double q)
{ return strFun[sf_num[i-1]-1].val(i,x,q);}

int loadStrFun(char *  name1, char*name2)
{ 
  int l, i;
  char *sf_txt[2];

  sf_txt[0]=name1;
  sf_txt[1]=name2;
  
  for(i=0;i<2;i++) sf_num[i]=0;
  for(i=0;i<2;i++)
  { long N; 
    pinf_int(Nsub,i+1,NULL,&N);
    for(l=0;l<MAXFUN;l++)
    if(strFun[l].myParticle(N)  && strFun[l].readName(i+1,sf_txt[i]))
    { sf_num[i] = l+1; 
      break; 
    }
  }
  initStrFun(0);
  return 0;                                                                       
}

int rd_sf__(FILE *mode)
{ char sf_txt1[500],sf_txt2[500];
  int  i;

  for(i=0;i<2;i++) sf_num[i]=0;
 
  if(1 != fscanf(mode," StrFun1=\"%[^\"]%*c", sf_txt1)) return 1;
  if(strcmp(sf_txt1,"OFF")) fscanf(mode," %*ld"); 

  if(1 != fscanf(mode," StrFun2=\"%[^\"]%*c", sf_txt2)) return 1;
  if(strcmp(sf_txt2,"OFF")) fscanf(mode," %*ld"); 
    
  loadStrFun(sf_txt1,sf_txt2);
  return 0;                                                                       
}


int wrt_sf__(FILE *mode)
{  
  int i;

  for(i=0;i<2;i++) 
  { fprintf(mode,"  StrFun%d=",i+1);
    if(sf_num[i])
    { char sf_txt[500];
      strFun[sf_num[i]-1].fullName(i+1,sf_txt);
      fprintf(mode,"\"%s\" %d\n",sf_txt,strFun[sf_num[i]-1].mc_code(i+1));
    }
    else fprintf(mode,"\"OFF\"\n");
  }  
  return 0; 
}
