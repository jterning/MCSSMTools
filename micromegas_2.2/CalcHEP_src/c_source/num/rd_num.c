/*
 Copyright (C) 1997, Alexander Pukhov 
*/
#include <stdio.h>
#include <ctype.h>
#include "interface.h"
#include "rd_num.h"

int  rd_num(char* s, double *p)
{
  int n;   

  if (isdigit(*s)) { sscanf(s,"%lf",p); return 1;}

  for(n=1;n<=nvar_int+nfunc_int;n++) if(strcmp(varName_int[n],s)==0) 
  {*p=va_int[n];return 1;}
  if(strcmp("PI",s)==0) { *p=M_PI; return 1;}  
  
  return 0;   
}
