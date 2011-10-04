#include"micromegas.h"
#include"micromegas_aux.h"
#include"micromegas_f.h"

#include <sys/types.h>
#include <unistd.h>
              
int readvar_(char *fname, int len)
{ int err;
  char * cname=malloc(len+1);
  fName2c(fname,cname,len);
  err=readVar(cname);
  free(cname);
  return err;
}

void printvar_(int * Nch) 
{ char fname[20];
  FILE*f;
  
  sprintf(fname,"%d.tmptxt",getpid());
  f=fopen(fname,"w");
  printVar(f);
  fclose(f);
  fortreread_(Nch,fname,strlen(fname));
  unlink(fname); 
}
                                                                                

void printmasses_(int * Nu, int* sort)
{ 
  char fname[20];
  FILE*f;

  sprintf(fname,"%d.tmptxt",getpid());
  f=fopen(fname,"w");
  printMasses(f,*sort);
  fclose(f);
  fortreread_(Nu,fname,strlen(fname));
  unlink(fname);
}

int sortoddparticles_(char * f_name, int len) 
{ 
  char c_name[20];
  int err=sortOddParticles(c_name);
  cName2f(c_name, f_name,len);
  return err;  
}
                                                                                                   
double darkomega_(double * Xf_, int * Fast, double *Beps){return darkOmega(Xf_, *Fast, *Beps);}

double printchannels_(double*Xf,double*cut,double*Beps,int* prcnt,int *Nu)
{
  char fname[20];
  FILE*f;
  double res;
                                                                                   
  sprintf(fname,"%d.tmptxt",getpid());
  f=fopen(fname,"w");
  res=printChannels(*Xf,*cut,*Beps,* prcnt,f);
  fclose(f);
  fortreread_(Nu,fname,strlen(fname));
  unlink(fname);
 
  return res;
}

double decay2info_(char * pname, int *Nch, int len)
{ double res;
  char cname[20]; 
  char fname[20];
  FILE*f;
  
  sprintf(fname,"%d.tmptxt",getpid());
  f=fopen(fname,"w");
  fName2c(pname,cname,len);    
  res=decay2Info(cname,f);
  fclose(f);
  fortreread_(Nch,fname,strlen(fname));
  unlink(fname); 
  return res;
}
   
