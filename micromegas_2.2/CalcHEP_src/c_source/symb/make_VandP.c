#include<stdio.h>
#include <math.h>
#include "model.h"
#include "read_mdl.h"
#include "procvar.h"
#include "reader_c.h"
#include "parser.h"

int main(int argv, char**argc)
{
  char fname[200];
  int i,nLn,L;
  FILE*f;  
  int nVar=0,nFunc=0,first; 
  char path[200];
  
  if(argv!=3) { printf("Arguments expected: 1)path to model files; 2) model number.\n"); return 1;}

  sscanf(argc[1],"%s",path);
  if(sscanf(argc[2],"%d",&L)!=1) { printf("Second argument should be a number\n"); return 1;}
  { 
     char ext[3][8]={"vars","func","prtcls"};
 
     for(i=0;i<3;i++)
     { 
        cleartab(modelTab+i);
        sprintf(fname,"%s/%s%d.mdl",path,ext[i],L);
        if(readtable(modelTab+i,fname)) 
            {printf("Problem in reading of Table %s\n", fname); return 2; }
     }  

     if(read2VarsParticles()) {printf("Error in model\n"); return 3;}
  }
  
  f=fopen("VandP.c","w");
  fprintf(f,"#include <stdio.h>\n");
  fprintf(f,"#include <math.h>\n");
  fprintf(f,"#include \"V_and_P.h\"\n");
  fprintf(f,"extern int  FError;\n");  


  for(nLn=0,i=0;i<nparticles;i++)
     if(!strchr("*fcCtT",prtclbase[i].hlp)&&i+1<=prtclbase[i].anti)nLn++;

  fprintf(f,"int nModelParticles=%d;\n",nLn);
  fprintf(f,"ModelPrtclsStr ModelPrtcls[%d]=\n{\n",nLn);

  for(nLn=0,i=0;i<nparticles;i++) if(!strchr("*fcCtT",prtclbase[i].hlp))
  { int anti=prtclbase[i].anti; 
    if(i+1>anti)  continue;
    if(nLn)fprintf(f,","); else fprintf(f," "); nLn++; 
      fprintf(f," {\"%s\",",prtclbase[i].name);
      if(i+1==anti)   fprintf(f,"\"%s\", ",prtclbase[i].name);
           else       fprintf(f,"\"%s\", ",prtclbase[anti-1].name);
     fprintf(f,"%ld, \"%s\",\"%s\",%d,%d}\n",
       prtclbase[i].N,  prtclbase[i].massidnt, prtclbase[i].imassidnt,  
       prtclbase[i].spin, prtclbase[i].cdim);
    
  }
  fprintf(f,"};\n");

  if (vararr) free(vararr);
  vararr = (singlevardescription*)m_alloc((nCommonVars+1)
                                            * sizeof(singlevardescription));

  for(i=0;i<=nCommonVars;i++)
  {  if(i==0 || strcmp(modelvars[i].varname,"i")==0) 
     { sprintf(vararr[i].alias,"XXX");
       vararr[i].tmpvalue=0;
       vararr[i].num=0;
       vararr[i].used = 0;
     }else
     { sprintf(vararr[i].alias,"V[%d]",nVar+nFunc);
       vararr[i].tmpvalue=modelvars[i].varvalue;
       vararr[i].num=nVar+nFunc;
       vararr[i].used = 1;
       if(modelvars[i].func) nFunc++; else nVar++;
     }
  }

  fprintf(f,"int nModelVars=%d;\n",nVar);
  fprintf(f,"int nModelFunc=%d;\n",nFunc);
  fprintf(f,"char*varNames[%d]={\n ", nVar+nFunc);
 
  for(first=1,i=0;i<=nCommonVars;i++) if(vararr[i].used) 
  { if(first)  first=0;else fprintf(f,",");
    fprintf(f,"\"%s\"\n",modelvars[i].varname);
  }
  fprintf(f,"};\n");
  fprintf(f,"double varValues[%d]={\n ", nVar+nFunc);
  for(first=1,i=0;i<=nCommonVars;i++) if(vararr[i].used) 
  { if(first)  first=0;else fprintf(f,",");
    fprintf(f,"%E\n",modelvars[i].varvalue);
  }
  fprintf(f,"};\n");

/*static void  writesubroutineinit(void) */

  fprintf(f,"#include\"extern.h\"\n");
  ext_h=fopen("extern.h","w");
  fprintf(f,"int calcMainFunc(void)\n{\n");
  fprintf(f," double *V=varValues;\n");
  fprintf(f," FError=0;\n");
  for(i=1;i<=nCommonVars;i++)
  {
     if (vararr[i].used &&  modelvars[i].func)
     { checkNaN=0;
        { char * ss=(char *)readExpression(modelvars[i].func,rd_c,act_c,free);
           fprintf(f,"   %s=%s;\n",vararr[i].alias,ss+3);
           free(ss);
        }
        if(checkNaN)
        fprintf(f,"   if(!finite(%s) || FError) return %d;\n",vararr[i].alias,vararr[i].num);
     }
  }

  fprintf(f,"return 0;\n}\n");
  fclose(ext_h);

  fclose(f);
  return 0;
}

