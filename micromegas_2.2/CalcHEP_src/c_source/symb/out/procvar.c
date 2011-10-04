/*
 Copyright (C) 1997, Alexander Pukhov 
*/
#include "syst.h"
#include "syst2.h"
#include "physics.h"
#include "parser.h"
#include "pvars.h"
#include "s_files.h"
#include "out_serv.h"
#include "procvar.h"
#include "saveres.h"
#include "process.h"

int nProcessVar;

singlevardescription *vararr = NULL;

 
static void * PP=(void *)"PP";

static void *  rd_hiddenVars(char* s)
{ int l;
  if (isdigit(*s)) return PP;
  for(l=0;l<=nmodelvar;l++)
     {if (!strcmp(modelvars[l].varname,s)) {vararr[l].used=1; break;}}
  return PP;
}

static void * act_hiddenVars(char * ch,int n,void ** mm1) {  return PP; }

int initvararray(int nsub, char key, int width)
{  int i,j,k,kk,l; 
   catrec    cr;
   FILE * catalog_;   

   polyvars  allVars={0,NULL}; 
   int nvar,nfunc;        
   int ArcNum=0;

   nProcessVar = nmodelvar  +PPSHIFT +((MAXINOUT*(MAXINOUT-1)/2)); 
   if (vararr) free(vararr); 
   vararr = (singlevardescription*)m_alloc(nProcessVar 
                                            * sizeof(singlevardescription)); 

   for(k=0; k<nProcessVar;k++)
   {  sprintf(vararr[k].alias,"#%d",k); 
      vararr[k].tmpvalue=0;
      vararr[k].num=0;
      vararr[k].used = 0;
   }
   sprintf(vararr[0].alias,"0");
   
   for(k=1;k<=nmodelvar;k++) vararr[k].used=modelvars[k].need;
      
   vardef=&allVars;
   
   if(width &1) catalog_=fopen(CATALOG_NAME,"rb"); else catalog_=NULL;
   
   if(catalog_)while (FREAD1(cr,catalog_))
   {  
      if(!nsub || cr.nsub_ == nsub)
      {
         whichArchive(cr.nFile,'r',0); 
/*         if(ArcNum!=0 && ArcNum!= cr.nFile) fclose(archiv);
         if(ArcNum==0 || ArcNum!= cr.nFile)
         {  char archivName[40];
            ArcNum=cr.nFile; 
            sprintf(archivName,ARCHIV_NAME,ArcNum);
            archiv=fopen(archivName,"rb");
          }
*/          
         seekArchiv(cr.factpos);
         readvardef(archiv);
         for(l=0;l<vardef->nvar;l++) vararr[vardef->vars[l].num].used=1;
         clearvardef();
         seekArchiv(cr.rnumpos);
         readvardef(archiv);
         for(l=0;l<vardef->nvar;l++) vararr[vardef->vars[l].num].used=1;
         clearvardef();
         seekArchiv(cr.denompos);
         readDenominators();
         for(l=0;l< denrno;l++) 
         { if( denom[l].mass)   vararr[denom[l].mass].used=1;
           if( denom[l].width)  vararr[denom[l].width].used=1;
         }
         clearvardef();
      }
   }
/*   if(ArcNum) fclose(archiv); */
   whichArchive(0,0,0);
   if(catalog_) fclose(catalog_);
   ArcNum=0;

   if(width&2)catalog_=fopen(CATALOG_NAME "2","rb"); else catalog_=NULL;
   if(catalog_) while(FREAD1(cr,catalog_))
   {
     if(!nsub || cr.nsub_ == nsub)
     {  
        whichArchive(cr.nFile,'r',1);
/*
        if(ArcNum!=0 && ArcNum!= cr.nFile) fclose(archiv);
        if(ArcNum==0 || ArcNum!= cr.nFile)
        {  char archivName[40];
           ArcNum=cr.nFile; 
           sprintf(archivName,ARCHIV_NAME "2",ArcNum);
           archiv=fopen(archivName,"rb");
         } 
*/                    
        fseek(archiv,cr.factpos,SEEK_SET);
        readvardef(archiv);
        for(l=0;l<vardef->nvar;l++) vararr[vardef->vars[l].num].used=1;
        clearvardef();
        fseek(archiv,cr.rnumpos,SEEK_SET);
        readvardef(archiv);
        for(l=0;l<vardef->nvar;l++) vararr[vardef->vars[l].num].used=1;
        clearvardef();
     }  
   }
/*   if(ArcNum) fclose(archiv);  */
   whichArchive(0,0,0);
   if(catalog_)fclose(catalog_);
  
   for (k = nmodelvar ; k >=0; k--) 
   if( vararr[k].used && modelvars[k].func)
       readExpression(modelvars[k].func,rd_hiddenVars,act_hiddenVars,NULL); 

   kk=0;
   for (i = 2; i <= nin+nout; i++)
   for (j = 1; j <= i - 1; j++)
   { k=scalarProductPos(i,j); 
     switch(key)
     { case 'R':
       case 'F': sprintf( vararr[k].alias, "p%d.p%d", j, i);    break;
       case 'M': sprintf( vararr[k].alias, "SC[p%d,p%d]", j, i);break;
       case 'c': sprintf( vararr[k].alias, "DP[%d]",kk); break;
       case 'f': { char c=kk+1;
                    if(c<10) sprintf( vararr[k].alias, "P%c",'0'+c);
                    else     sprintf( vararr[k].alias, "P%c",'A'+c-10);
                 } break;
     } 
      vararr[k].used = 1;
      kk++; 
   }
   
   for(i=1;i<=2;i++)
   { k=i+nmodelvar+1; 
     switch(key)
     { case 'R':
       case 'F':
       case 'M': sprintf( vararr[k].alias, "Helicity%d", i);  break;
       case 'c': sprintf( vararr[k].alias, "Helicity[%d]",i-1); break;
       case 'f': sprintf( vararr[k].alias, "Helicity(%d)",i); break;
     } 
     vararr[k].used = 1;
   }
   for(i=1;i<=2;i++)
   { k=i+nmodelvar+3; 
     switch(key)
     { case 'R':
       case 'F':
       case 'M': sprintf( vararr[k].alias, "HelicityN%d", i);  break;
       case 'c': sprintf( vararr[k].alias, "HelicityN[%d]",i-1); break;
       case 'f': sprintf( vararr[k].alias, "HelicityN(%d)",i); break;
     } 
     vararr[k].used = 1;
   }

      
   nvar=0;
   for(k=0;k<=nmodelvar;k++) 
   if (vararr[k].used && (!modelvars[k].func && !modelvars[k].pwidth))
   {  vararr[k].tmpvalue = modelvars[k].varvalue;
   
      switch(key)
      {
          case 'R':
          case 'F':
          case 'M':  strcpy( vararr[k].alias, modelvars[k].varname);break;
          case 'c':  sprintf(vararr[k].alias,"V[%d]",++nvar); 
                     vararr[k].num=nvar; break;
          case 'f':  sprintf(vararr[k].alias,"A(%d)",++nvar); break;
      }
   }
   nfunc=0;
   for(k=0;k<=nmodelvar;k++) 
   if (vararr[k].used && (modelvars[k].func ||modelvars[k].pwidth))
   {  vararr[k].tmpvalue = modelvars[k].varvalue;
      switch(key)
      {
          case 'R':
          case 'F':
          case 'M':  strcpy( vararr[k].alias, modelvars[k].varname);break;
          case 'c':  sprintf(vararr[k].alias,"V[%d]",++nfunc+nvar); 
                     vararr[k].num=nfunc+nvar; break;
          case 'f':  sprintf(vararr[k].alias,"A(%d)", ++nfunc+nvar); break;
      }
   }

   return 1;
} 
