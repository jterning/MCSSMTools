#include <stdio.h>
#include "histogram.h"
#include "interface.h"
/*
#define const
#include"num_out.h"
#undef const


int  nvar_int=0;
char * varName_int[1];
double  va_int[1];
int nfunc_int=0; 
*/
static char* pinf_ext(int nsub,int num , double * mass, long * N)
{ static char * pname="e";
  return pname;
}

 
int main(int np, char ** par)
{
  int n;
  char process[200];
  FILE*f;

  nin_int=2;
  nout_int=2;
  pinf_int=&pinf_ext;     
  nvar_int=0;
  nfunc_int=0;

  if(np<3) 
  { printf(" This  routine is intended to sum the  distributions produced\n" 
           " in  calcHEP numerical sessions (files dist_#).\n"
           " The names of files must be submitted as parameters\n"
           " Resulting distribution is directed to stdout(screen)\n");
    return 1;
  }

  for(n=1;n<np;n++)
  { 
    f=fopen(par[n],"r"); 
    if(!f) return 2;
    if(add_hist(f,process)) {fclose(f);return 3;}
    fclose(f);
  } 
  wrt_hist2(stdout,process);

  return 0;
}
