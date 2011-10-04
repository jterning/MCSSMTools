/*
 Copyright (C) 1997-2006, Alexander Pukhov
*/

#include <limits.h>
#include <math.h>
#include <ctype.h>
#include <stdarg.h>

#include "chep_crt.h"
#include "syst2.h"
#include "physics.h"
#include "parser.h"
#include "s_files.h"
#include "file_scr.h"
#include "pre_read.h"
#include "pvars.h"
#include "getmem.h"
#include "read_func.h"

#define minvarmem  ((unsigned)sizeof(struct varrec)  + 1 - STRSIZ)
#define minlagrmem ((unsigned)sizeof(struct algvert) + 1 - STRSIZ)

#include "read_mdl.h"

int  nCommonVars=0;
table  modelTab[5] =
{{"","","",NULL},
 {"","","",NULL},
 {"","","",NULL},
 {"","","",NULL},
 {"","","",NULL}};

static int lastModel=0;
static char Sqrt2FuncTxt[10]="sqrt(2)";
static char * tabName;
static int nLine;

static void errorMessage( char * fieldName, char * format, ...)
{ 
   char  dump[200];
   va_list args;
   va_start(args,format);
   vsprintf(dump,format,args);
   va_end(args);

   if (strcmp(dump,"*") == 0)
    sprintf(errorText,"Error in table '%s' line %d field '%s'\nposition %u: %s",
	  tabName,nLine,fieldName,rderrpos,errmesstxt(rderrcode) );
   else
	sprintf(errorText,"Error in table '%s' line %d field '%s' \n %s",
		  tabName,nLine,fieldName,dump);
   if(blind) printf("ERROR:%s\n",errorText); else messanykey(2,10,errorText);
}


static int tabCharPos(char *  str, int  n)
{
 int k=0;
 int nn=0;
 if (n==0) return 0;
 while (str[k] != 0)
 { if (str[k]=='|')
	{ nn++;
	  if (nn==n) return (k+1);
	}
	k++;
 }
 return k;
}

static int  isVarName(char*  s)
{
  int i=1;

  if ( !isalpha(s[0]) ) return 0;
  while  (s[i] !=0)
  { 
    if (!(isalnum(s[i])||s[i]=='_' ) )  return  0;
    i++;
  }
  
  if(strlen(s)>6) return 0;

  if(s[0]=='p'||s[0]=='m'||s[0]=='M')
  {  i=1;
     while  (isdigit(s[i])) i++;
     if (s[i]==0) return 0;
  }

  if(!strcmp(s,"G5")) return 0;
  if(!strcmp(s,"g5")) return 0;
  if(!strcmp(s,"I")) return 0;
  if(!strcmp(s,"i")) return 0;  
  return 1;
}

static int  isOriginName(char* s)
{
  char newName[10], oldName[10];
  int i,k;

  strcpy(newName,s);
  i=0; while (newName[i] !=0) {newName[i]=toupper(newName[i]);i++;}

  for(k=1;k<=nmodelvar;k++)
  {
     strcpy(oldName,modelvars[k].varname);
      i=0; while (oldName[i] !=0) {oldName[i]=toupper(oldName[i]);i++;}
      if (strcmp(newName,oldName) == 0) return 0;
  }  
  return 1;
}

static int Number,String;

static void * rd_33(char*s)
{  int i;

   if(isdigit(s[0])) return &Number;
   if(s[0]=='"')  return  &String;
   if (strcmp(s,strongconst) == 0 || strcmp(s,"i")==0  ){rderrcode=unknownidentifier; return NULL;} 
   for(i=1;i<=nmodelvar;i++) if(strcmp(modelvars[i].varname,s) == 0)  return &Number;
   rderrcode=unknownidentifier; return NULL;   
}  
                                
static void*  act_33(char* ch,int n, void**args)
{ 
   if(ch[0]=='.') rderrcode=typemismatch;
   
   switch (ch[0])
   { 
     case '-': 
     case '+': 
     case '*': 
     case '/': 
     case '^':  if(args[0]==&String || args[1]==&String) { rderrcode=typemismatch; return NULL;}
                else return &Number;
     case '.': rderrcode=typemismatch; return NULL;
   }
   if(n==1 && args[0]==&String) if(!strcmp(ch,"sqrt")||
                !strcmp(ch,"sin")  ||
                !strcmp(ch,"cos")  ||
                !strcmp(ch,"tan")  ||
                !strcmp(ch,"asin") ||
                !strcmp(ch,"acos") ||
                !strcmp(ch,"atan") ||
                !strcmp(ch,"exp")  ||
                !strcmp(ch,"log")  ||
                !strcmp(ch,"fabs")    )  {rderrcode=typemismatch; return NULL;}
                
   if(n==2 && (args[0]==&String || args[1]==&String)  ) 
   if(!strcmp(ch,"atan2") || !strcmp(ch,"if")) {rderrcode=typemismatch; return NULL;} 
   return &Number;
}


void  readEXTLIB(void)
{ linelist  ln;
  if(EXTLIB) free(EXTLIB);
  EXTLIB=malloc(2); EXTLIB[0]=0;
  for(ln=modelTab[4].strings ; ln; ln=ln->next)
  { char buff[100];
    if(sscanf(ln->line,"%[^%\n|]", buff)!=1) continue;
    trim(buff);
    if(strlen(buff))
    { EXTLIB=realloc(EXTLIB, strlen(EXTLIB)+3+strlen(buff));
      sprintf(EXTLIB+strlen(EXTLIB)," %s",buff); 
    }  
  }
}

static int  readvars(int  check, varrec * recForGG)
{
  char      numtxt[60];
  char      name[60];
  char *    ss, * endstr;
  int       commentShift,funcShift;
  varlist   mvars;
  double    varvalue_tmp;
  linelist  ln;
  int       ggOn=0;
  
  int nv=3; /* Plus 0, i,Sqrt2 */
  recForGG->need=-1;
  ln=vars_tab.strings; while (ln) { ln=ln->next; nv++;} 
  ln=func_tab.strings; while (ln) { ln=ln->next; nv++;} 
  ln=prtcls_tab.strings; while (ln)
  { char imassname[20];
    char *p=ln->line; 
    int i;
    for(i=0;i<6;i++) p=strchr(p,'|')+1;  
    sscanf(p,"%[^|]", imassname);
    trim(imassname);
    if(imassname[0]=='!') nv++; 
    ln=ln->next;
  }  

  if(nv>SHRT_MAX){errorMessage("Name","too many parameters"); goto errExi1;}


  tabName=vars_tab.headln;
  commentShift=tabCharPos(vars_tab.format,2);

  if(modelvars) free(modelvars);

   modelvars = m_alloc(nv*sizeof(*modelvars));
   nmodelvar=0;

   strcpy(modelvars[0].varname,"0");
   modelvars[0].varvalue=0.;
   modelvars[0].func=NULL;   
   modelvars[0].need=0;
   modelvars[0].hidden=0;
   modelvars[0].pwidth=0;

   for(nLine=1,ln=vars_tab.strings; ln; nLine++, ln=ln->next)
   {  int need=0;
      ss=ln->line;

      sscanf(ss,"%[^|]%*c%[^|]", name,numtxt);

      trim(name); trim(numtxt);		
      if(name[0]=='*') {need=1; name[0]=' '; trim(name);}
      else if(name[0]=='%') continue;
       
      if (check && (!isVarName(name)) )
      {  errorMessage("Name","incorrect name '%s'",name);
         goto errExi1;
      }


      if (check &&  (! isOriginName(name)) )
      {  errorMessage("Name","duplicate name '%s'",name);   
         goto errExi1;
      }

      varvalue_tmp=strtod(trim(numtxt), &endstr);
      
      if ( ( check )&&( endstr !=numtxt+strlen(numtxt)) )
      {  errorMessage("Value","wrong number '%s'",numtxt);
	 goto errExi1;
      }


      if(strcmp(name,strongconst) == 0)  
      {
         if(ggOn)
         {  errorMessage("Name","duplicate name '%s'",name);
            goto errExi1;
         }
         ggOn=1;
         strcpy(recForGG->varname,strongconst);
         recForGG->func = NULL;
         recForGG->varvalue = varvalue_tmp;
         recForGG->need=need;
         recForGG->hidden=0;
         recForGG->pwidth=0;
         continue; 
      } 
      nmodelvar++;
      mvars=modelvars+nmodelvar; 
      strcpy(mvars->varname,name);
      mvars->varvalue=varvalue_tmp;
      mvars->func= NULL;
      mvars->need=need;
      mvars->hidden=0;
      mvars->pwidth=0;
   }

   nmodelvar++;
   mvars=modelvars+nmodelvar;
   strcpy(mvars->varname,"i");
   mvars->func = NULL;
   mvars->varvalue = 0.;
   mvars->need=0;
   mvars->hidden=0;
   mvars->pwidth=0;

   mvars++; 
   nmodelvar++;
   strcpy(mvars->varname,"Sqrt2");
   mvars->func = Sqrt2FuncTxt;
   mvars->varvalue = sqrt(2.0);
   mvars->need=0;
   mvars->hidden=1;
   mvars->pwidth=0;

   nCommonVars=0; 
  
   tabName=func_tab.headln;
   commentShift=tabCharPos(func_tab.format,2);
   funcShift=tabCharPos(func_tab.format,1);

   for(nLine = 1,ln=func_tab.strings; ln; nLine++,ln=ln->next)
   {  
      int need=0,hidden=0;
      ss=ln->line;
      sscanf(ss,"%[^|]",name);
      trim(name);
      if(name[0]=='*') {name[0]=' '; trim(name); need=1;} 
      else if(name[0]=='#') {name[0]=' '; trim(name); hidden=1;}
      else if(name[0]=='%') { if(strcmp(name,"%Local!")==0) nCommonVars=nmodelvar;  continue;}
      
      if (! isVarName(name))
      {  errorMessage("Name","incorrect name '%s'",name);
	 goto errExi1;
      }

      if ( check && (! isOriginName(name)) )
      {  errorMessage("Name","duplicate name '%s'",name);
			goto errExi1;
      }

      if(readExpression(ln->line+funcShift,rd_33, act_33,NULL)==&String) rderrcode=typemismatch;

      if(rderrcode &&(rderrcode!=unknownfunction))
      {    
	  errorMessage("Expression","*");
	  goto errExi1;
      }
      
      nmodelvar++;
      
      mvars=modelvars+nmodelvar;
      if(rderrcode==unknownfunction) varvalue_tmp=0; 
      mvars->varvalue = varvalue_tmp;
      strcpy(mvars->varname,name);
      mvars->func=ln->line+funcShift;
      mvars->need=need;
      mvars->hidden=hidden;
      mvars->pwidth=0;
   }  
   return 1;

errExi1:
   free (modelvars);
   modelvars=NULL; 
   return 0;
}



static int  findvar(char* txt,double* num, int *pos)
{
   int  i;

   trim(txt);
  
   for(i=1;i<=nmodelvar;i++)
   { 
     if(strcmp(txt,modelvars[i].varname)==0)
     {
        if(num) *num = modelvars[i].varvalue;
        if(pos) *pos =i;
        return 0;
      }
   }
   return -1;
}


static int  ghostaddition(void)
{  int  i, nPrim;
   nPrim=nparticles;
   if( (prtclbase[nPrim-1].spin == 2) /*&& (prtclbase[nPrim-1].cdim !=1)*/ )
   { 
      nPrim ++;   
      prtclbase[nPrim -1] = prtclbase[nparticles-1];      
      prtclbase[nparticles -1].hlp = 'T';
      prtclbase[nparticles -1].spin=4;
      nparticles ++; 

      nPrim ++;   
      prtclbase[nPrim -1] = prtclbase[nparticles-1];      
      prtclbase[nparticles -1].hlp = 't';
      prtclbase[nparticles -1].spin=4;
      nparticles ++; 
   }
   
   if (gaugep(nPrim))
   {
      nparticles ++;
      prtclbase[nparticles -1] = prtclbase[nPrim-1];
      prtclbase[nparticles -1].hlp = 'c';
      nparticles ++;
      prtclbase[nparticles -1] = prtclbase[nPrim-1];
      prtclbase[nparticles -1].hlp = 'C';
      
      if (strcmp(prtclbase[nPrim-1].massidnt,"0") != 0)
      {
         nparticles ++;
         prtclbase[nparticles -1] = prtclbase[nPrim-1];
         prtclbase[nparticles -1].hlp = 'f';
      }
      for (i=nPrim+1; i<=nparticles;i++)  prtclbase[i -1].spin=0;
   }
   return nPrim;
}


static void  clearlgrgn(void)
{algvertptr  l1;

   while (lgrgn != NULL)
   {
      l1 = lgrgn;
      lgrgn = lgrgn->next;
      free(l1);
   }
}


static void  cleardecaylist(void)
{decaylink   v1, v2;
 int        j;

   if(!prtclbase) return; 
   for (j = 0; j < nparticles; j++)
   {
      v1 = prtclbase[j].top;
      prtclbase[j].top = NULL;
      while (v1 != NULL)
      {
         v2 = v1;
         v1 = v1->next;
         free(v2);
      }
   }
}

static void  clearLatexNames(void)
{  int  j;
   if(!prtclbase) return; 
   for (j = 0; j < nparticles; j++)
   {  if (!strchr("fcCtT",prtclbase[j].hlp)) free(prtclbase[j].latex);}
}

static int isPrtclName(char*p){ return (strlen(p)<=P_NAME_SIZE);}

static int  readparticles(int  check, int ugForce )
{  char      *ss,*endstr;
   char      fullname[60], massname[60], imassname[60], p1[60], p2[60],numtxt[20];
   char      latex[STRSIZ], latex_[STRSIZ], s[60],c[60], chlp[40];
   int       itmp,i,j, errcode,np1,np2,nparticleLimit =128;

   linelist  ln=prtcls_tab.strings;

   tabName=prtcls_tab.headln;

   if(prtclbase) { cleardecaylist(); clearLatexNames(); free(prtclbase);}
 
   prtclbase=(prtcl_base*) malloc(nparticleLimit*sizeof(prtcl_base));
   prtclbase1=prtclbase-1;	
   nparticles = 0;

   for(i=nparticles;i<nparticleLimit;i++)
           {prtclbase[i].top=NULL;prtclbase[i].latex=NULL;}
   
   nLine=1;
   while (ln != NULL)
   {  ss=ln->line;
      if (nparticles >= nparticleLimit-16)
      {  nparticleLimit+=128;
         prtclbase=re_alloc(prtclbase,nparticleLimit*sizeof(prtcl_base));
         if(!prtclbase) 
         { errorMessage(" P ","too many particles");
	   return 0;
	 }
	 prtclbase1=prtclbase-1;
         for(i=nparticles;i<nparticleLimit;i++)
               {prtclbase[i].top=NULL;prtclbase[i].latex=NULL;} 
      }

      sscanf(ss,"%[^|]%*c%[^|]%*c%[^|]%*c%[^|]%*c%[^|]%*c%[^|]%*c%[^|]%*c%[^|]%*c%[^|]%*c%[^|]%*c%[^|]",
	    fullname,p1,p2,numtxt,s,massname,imassname,c,chlp,latex,latex_);
      trim(p1); trim(p2); trim(latex); trim(latex_); 
      {
         static char fldName[2][5]={" P "," aP"};
         char * pName[2];
         pName[0]=p1;
         pName[1]=p2;

         for ( i=0;i<=1;i++)
         {
            if (check && (! isPrtclName(pName[i])))
            {  errorMessage(fldName[i],"incorrect particle name '%s'",pName[i]);
               return 0;
            }

            if (check )
            {
               locateinbase(pName[i],&j);
               if (j != 0)
               {
                  errorMessage(fldName[i],"duplicate particle name '%s'",pName[i]);
                  return 0;
               }
            }
         }
      }
      nparticles++;
      strcpy(prtclbase[nparticles-1].name,p1);

            
      itmp=strtol(trim(s),&endstr,10);
      if(check)
      {
         if (s+strlen(s) != endstr)
         {  errorMessage("2*spin","number expected");
            return 0 ;
         }
         if((itmp!=0)&&(itmp!=1)&&(itmp!=2)&&(itmp!=3)&&(itmp!=4))
         {  errorMessage("2*spin","value out of range");
            return 0;
         }
      }
      prtclbase[nparticles-1].spin=itmp;
      if( 1!=sscanf(numtxt,"%ld",&prtclbase[nparticles-1].N))
      {    errorMessage("number","can't read MC numeration parameter");
            return 0;
      }
      trim(massname);
      if(strcmp(massname,"0")==0)
      { if(prtclbase[nparticles-1].spin==3 || prtclbase[nparticles-1].spin==4)
         errorMessage("mass","spin 3/2 and spin 2 particles should be massive");
      }
      else 
      {  int pos;
         errcode=findvar(massname,NULL,&pos);
         if (check && (errcode != 0))
         {
            errorMessage("mass","unknown variable %s",massname);
            return 0;
         }else 
         {  if(pos>nCommonVars) nCommonVars=pos;

         } 
      }
      strcpy(prtclbase[nparticles-1].massidnt,massname);

      trim(imassname);
      
      if( check && strcmp(imassname,"0")&& !strcmp(massname,"0")) 
      {  errorMessage("width","non zero width for zero mass particle '%s'", 
            prtclbase[nparticles-1].name);
         return 0;
      }
                                       
      if(strcmp(imassname,"0") != 0)
      { 
         errcode=findvar(imassname,NULL,NULL);
         
         if(errcode)  
         {  if(imassname[0]=='!') 
            {  imassname[0]=' ';  
               trim(imassname);
               if(check && (!isVarName(imassname)) )
               {  errorMessage("width","incorrect name '%s'",imassname);
                  return 0;
               }
               if(check && (!isOriginName(imassname)) )                            
               { errorMessage("width","this identifier  '%s' already was used",imassname);
                 return 0;
               }
               { varlist mvars=modelvars+1+nmodelvar; 
                  nmodelvar++;
                  strcpy(mvars->varname,imassname);
                  mvars->func = NULL;
                  mvars->varvalue = 0.;
                  mvars->need=0;
                  mvars->hidden=0;
                  mvars->pwidth=nparticles;
               }                                                     
            }else 
            { 
               errorMessage("width","unknown variable %s",imassname);
               return 0;
            }
         }
      }
      strcpy(prtclbase[nparticles-1].imassidnt,imassname);

      itmp=strtol(trim(c),&endstr,10);
      if(check)
      {
         if (c+strlen(c) != endstr)
         {  errorMessage("color","number expected");
            return 0;
         }
         if (((itmp!=1)&&(itmp!=3)&&(itmp!=8))||((itmp==3)&&(strcmp(p1,p2)==0))  )
         {  errorMessage("color","value out of range");
            return 0;
         }
      }
      prtclbase[nparticles-1].cdim=itmp;
      trim(chlp);
      if (strcmp(chlp,"") == 0) strcpy(chlp," ");
      prtclbase[nparticles-1].hlp = toupper(chlp[0]);
      if(check)
      {  int ner;
         ner=1;
         switch(prtclbase[nparticles-1].hlp)
         {
            case ' ':if(prtclbase[nparticles-1].spin==2 &&
                        !strcmp(prtclbase[nparticles-1].massidnt,"0")
                       )
              {  errorMessage("aux","Massless vector boson must\n"
                                     " be a gauge particle");                  
                 return 0;
              }
              break;
            case 'L':
            case 'R':
              if ((prtclbase[nparticles-1].spin !=1)
               ||((prtclbase[nparticles-1].massidnt[0])!='0')
               ||(strcmp(p1,p2)==0))  ner=0;
              break;

            case '*':
              if(prtclbase[nparticles-1].massidnt[0]=='0')   ner=0;
              else if(prtclbase[nparticles-1].imassidnt[0]!='0')
              { errorMessage("aux","for aux='*' zero width is  expected"); return 0;}
              break;
            case 'G':
              if(prtclbase[nparticles-1].spin!=2)   ner=0;
              break;
            default: ner=0;
         }
         if(!ner){ errorMessage("aux","unexpeted character"); return 0;}
         if(prtclbase[nparticles-1].N==0 && prtclbase[nparticles-1].hlp!='*')
         { errorMessage("number","Zero PDG code."); return 0;}  
      }
      prtclbase[nparticles-1].latex=malloc(1+strlen(latex));
      strcpy(prtclbase[nparticles-1].latex,latex);
      
      np1 = ghostaddition();
      if (strcmp(p1,p2) == 0) prtclbase[np1-1].anti = np1;
      else
      {
        ++(nparticles);
        prtclbase[nparticles-1] = prtclbase[np1-1];
        prtclbase[nparticles-1].N *=(-1);
        strcpy(prtclbase[nparticles-1].name,p2);
        prtclbase[nparticles-1].latex=malloc(1+strlen(latex_));
        strcpy(prtclbase[nparticles-1].latex,latex_);
        if (prtclbase[np1-1].cdim == 3) prtclbase[nparticles-1].cdim = -3;
        np2=ghostaddition();
        prtclbase[np1-1].anti = np2;
        prtclbase[np2-1].anti = np1; 
      }
      ln=ln->next;
      nLine++;
   }

   for (i = 1; i <= nparticles; i++)
   {  prtcl_base *with1 = &prtclbase[i-1];
      with1->top = NULL;
      
      if (strchr("fcCtT",with1->hlp) != NULL)
      {  
         sprintf(with1->name+strlen(with1->name),".%c",with1->hlp); 
         switch (with1->hlp)
         {
            case 'c':                        
               with1->anti = prtclbase[i-1 - 1].anti +2;
               break;
            case 'C':
               with1->anti = prtclbase[i-1 - 2].anti +1;
               break;
            case 'f':
               with1->anti = prtclbase[i-1 - 3].anti +3;
               break;
            case 't':
               with1->anti = prtclbase[i-1 + 1].anti -1;
               break;               
            case 'T':
               with1->anti = prtclbase[i-1 + 2].anti -2;
               break;               
         }
      }
   }
   
   if(ugForce)for(i=1; i <= nparticles; i++) 
               if(gaugep(i) && (!zeromass(i))) prtclbase[i-1].hlp=' ';

  for(i=1;i<=nmodelvar;i++)
  if(modelvars[i].pwidth) modelvars[i].pwidth=ghostmother(modelvars[i].pwidth);
   
   return 1;
}

static int  testLgrgn(algvertptr lgrgn)
{
  preres  m;
  int n;
/*  goto_xy(1,20); print("%d           ",nLine); */
  m = (preres) readExpression(lgrgn->comcoef,rd_pre, act_preF,NULL);
  if (rderrcode )
  {  errorMessage("Factor","*");
    return 0;
  }
  m->free=1;

  if (m->tp >rationtp)
  {  errorMessage("Factor","scalar expected");
    return 0;
  }

  if (m->maxp>0)
  {  errorMessage("Factor","moments p%d are not permitable here",m->maxp);
    return 0;
  }

  for (n = 0; n < vardef->nvar; n++)
  { int err;
    err=findvar(vardef->vars[n].name,NULL,NULL);
    if (err)
    {  errorMessage("Factor","unknown variable '%s'", vardef->vars[n].name);
      return 0;
    }
  }
  
  clearVars(vardef);

   
  m=(preres) readExpression(lgrgn->description,rd_pre,act_pre,NULL);
  if(rderrcode) {  errorMessage("Lorentz part","*"); return 0; }
  m->free=1;
 
  if (m->tp == rationtp)
  {  errorMessage("Lorentz part","division is not permited here");  return 0; }

  if( (m->tp == spintp) &&( prtclbase1[lgrgn->fields[0]].spin&1 !=1) 
                        &&( prtclbase1[lgrgn->fields[1]].spin&1 !=1)
                        &&( prtclbase1[lgrgn->fields[2]].spin&1 !=1) )
  {
    errorMessage("Lorentz part","Dirac gamma matrix not expected");
    return 0;
  }

  if ((m->maxp == 4)&&(lgrgn->fields[3] == 0))
  {  errorMessage("Lorentz part","p4 are not permited here");
    return 0;
  }


  for (n = 0; n < vardef->nvar; n++)
  {  int err;
    err=findvar (vardef->vars[n].name,NULL,NULL);
    if (err)
    {  errorMessage("Lorentz part","unknown variable '%s'",vardef->vars[n].name);
      return 0;
    }
  }

  clearVars(vardef);


  for (n = 0; n <= 3; n++)
  {  int  ind1,ind2,np ;
     
     ind1=0;
     np  = lgrgn->fields[n];
     if ( np != 0 )   switch  (prtclbase[np-1].spin)
     { case 2: 
       case 3: ind1=1; break;
       case 4: ind1=3;
     }     
     
     ind2=0;
     if( (1<<n)& m->indlist)     ind2 += 1;
     if( (1<<(n+4))& m->indlist) ind2 += 2;
     
     if (ind1 != ind2 )
     {  errorMessage("Lorentz part","index 'm%d'  unbalanced",n+1);
        return 0;
     }
  }
  
  return 1;
}


static int  readlagrangian(int check, int ugForce)
{ 
  algvertptr  lgrgn1,lgrgn2;
  int    i, j, mm;
  char * ss;
  char   pPtr[4][60];
  int  factorShift,lorentzShift;
  arr4byte  f_copy;
  int mLine,totcolor,color,spinorNumb;
  linelist ln;
  static char fName[4][5] = {"P1","P2","P3","P4"};
  polyvars var_testing={0,NULL}; 

  vardef=&(var_testing);


  clearlgrgn();
  factorShift=tabCharPos(lgrng_tab.format,4);
  lorentzShift=tabCharPos(lgrng_tab.format,5);
  tabName=lgrng_tab.headln;

  for(ln=lgrng_tab.strings, nLine=1; ln; ln=ln->next, nLine++)
  { 
    ss=ln->line;
    sscanf(ss,"%[^|]%*c%[^|]%*c%[^|]%*c%[^|]",pPtr[0],pPtr[1],pPtr[2],pPtr[3]);
    for(i=0;i<4;i++) trim(pPtr[i]);
    if(pPtr[0][0]=='%') continue;
    for(i=0;i<4;i++)
    { 
      if(pPtr[i][0]) 
      { locateinbase(pPtr[i],&j);
        if(check && j == 0) 
        { errorMessage( fName[i]," unknown particle %s" ,pPtr[i]);
                          return 0;
        } 
        f_copy[i]=j;
      }
      else if(i==3) f_copy[3]=0; 
      else { errorMessage( fName[i],"particle name is expected");return 0;}
    }                            

    if(ugForce)
    { for(i=0;i<4;i++)
      { j=f_copy[i];
         if(j && ghostp(j) &&(!zeromass(j))) i=10;
      } 
      if(i>=10) continue;
    }

    lgrgn1=(algvertptr)m_alloc( sizeof(*lgrgn1));
    lgrgn1->next = lgrgn;
    lgrgn = lgrgn1;
    lgrgn->comcoef=    ln->line+factorShift;
    lgrgn->description=ln->line+lorentzShift;
    for (i=0;i<4;i++) lgrgn->fields[i] = f_copy[i];

    if(check)
    {
      totcolor=1;
      for (mm=0;((mm<4)&&(lgrgn->fields[mm] !=0));mm++)
      {
        color=prtclbase[lgrgn->fields[mm] -1].cdim;
        if (color==-3) color=5;
        totcolor=totcolor*color;
      }
      if( (totcolor!=1)&&(totcolor!=15)&&(totcolor!=64)&&(totcolor!=120)&&(totcolor!=512) )
      {   errorMessage("Lorentz part","wrong color structure");
         return 0;
      }
      spinorNumb=0;
      for (mm=0;((mm<4)&&(lgrgn->fields[mm] !=0));mm++)
      {
        if( prtclbase1[lgrgn->fields[mm]].spin&1 )  spinorNumb++ ;
      }
      if( (spinorNumb!=0)&&(spinorNumb!=2) )
      {  errorMessage("Lorentz part","wrong spinor  structure");
        return 0;
       }
    }
    if (! testLgrgn(lgrgn) )  { clearVars(vardef); return 0;}
   }

   clearVars(vardef);
   clearpregarbage();

   lgrgn1 = lgrgn;   /*     Sorting    */
   do
   {  lgrgn1->factor=1;
      for(i=0;i<4 && lgrgn1->fields[i];i++)
      { int hlp=prtclbase[lgrgn1->fields[i]-1].hlp;
        if(hlp=='C') break; 
        else if(hlp=='c') {lgrgn1->factor=-1; break;}
      }    
      for (i = 1; i <= 4; i++) lgrgn1->perm[i-1] = i;
      i = 1;
      while (i < 4)
         if (lgrgn1->fields[i-1] >= lgrgn1->fields[i + 1-1]) ++(i);
         else
         {
            mm = lgrgn1->fields[i-1];
            lgrgn1->fields[i-1] = lgrgn1->fields[i + 1-1];
            lgrgn1->fields[i + 1-1] = mm;
            mm = lgrgn1->perm[i-1];
            lgrgn1->perm[i-1] = lgrgn1->perm[i + 1-1];
            lgrgn1->perm[i + 1-1] = mm;
            if (i == 1)
               ++(i);
            else
               --(i);
         }
      lgrgn1 = lgrgn1->next;
  }  while (lgrgn1 != NULL);

  if (check)
  {
    mLine=nLine;
    lgrgn1 = lgrgn;   /*    check1       */
    do
    {
      nLine--;
      lgrgn2=lgrgn1->next;
      while (lgrgn2 != NULL )
      { if( (lgrgn1->fields[0]==lgrgn2->fields[0]) &&
            (lgrgn1->fields[1]==lgrgn2->fields[1]) &&
            (lgrgn1->fields[2]==lgrgn2->fields[2]) &&
            (lgrgn1->fields[3]==lgrgn2->fields[3])
           )
        {  char sss[20]="";
           for(i=0;i<4;i++) if(lgrgn1->fields[i])
           { strcat(sss, prtclbase1[lgrgn1->fields[i]].name); strcat(sss," ");}                                            
           errorMessage("P1,P2,P3,P4","duplicate vertex {%s}",sss);
           return 0;
        }
        lgrgn2=lgrgn2->next;
      }
      lgrgn1= lgrgn1->next;
     }  while (lgrgn1 != NULL);


    nLine=mLine;
    lgrgn1 = lgrgn;   /*    check2       */
    do
    {
      nLine--;
      for (i=0;i<4;i++)
      {  f_copy[i]=lgrgn1->fields[i];
        if (f_copy[i] !=0)
        {
          mm=ghostmother(f_copy[i]);
          f_copy[i]=prtclbase[mm-1].anti  + f_copy[i]-mm   ;
         }
      }

      i = 1;
      while (i < 4)
        if (f_copy[i-1] >= f_copy[i ]) ++(i);
        else
        {
          mm = f_copy[i-1];
          f_copy[i-1] = f_copy[i ];
          f_copy[i ] = mm;
          if (i == 1)
            ++(i);
          else
            --(i);
        }

      lgrgn2=lgrgn;
      while ((lgrgn2 != NULL ) && (  (f_copy[0] !=lgrgn2->fields[0]) ||
                          (f_copy[1] !=lgrgn2->fields[1])  ||
                          (f_copy[2] !=lgrgn2->fields[2])  ||
                          (f_copy[3] !=lgrgn2->fields[3])
                         )
          )
         {
         lgrgn2=lgrgn2->next;
          }
      if (lgrgn2 == NULL)
      {  char sss[10];
        strcpy(sss,"");
        for (i=0;i<3;i++)
        { strcat (sss, prtclbase[lgrgn1->fields[i]-1].name);
          strcat(sss," ");
        }
        if (lgrgn1->fields[3] !=0   )
           strcat(sss,prtclbase[lgrgn1->fields[3]-1].name);

      errorMessage("P1,P2,P3,P4","conjugated vertex for %s not found",sss);
            return 0;
       }
      lgrgn1= lgrgn1->next;
     }  while (lgrgn1 != NULL);
   }
  return 1;
}


static void  filldecaylist(void)
{ algvertptr  lgrgn1;
  int        i, j, k, n;
  particleNumType   pn[5], cc[3];
  decaylink   kk, qq;

   lgrgn1 = lgrgn;
   do
   {
      for (i = 1; i <= 4; i++)
      {
         pn[i-1] = ghostmother(lgrgn1->fields[i-1]);
         if (pn[i-1] != 0) pn[i-1] = prtclbase[pn[i-1]-1].anti;
      }
      pn[4] = 0;

      for (i = 1; i <= 4; i++)
         if (pn[i-1] != pn[i + 1-1] && pn[i-1] != 0)
         {
            j = 1;
            for (k = 1; k <= 4; k++)
               if (k != i)
               {  cc[j-1] = pn[k-1];
                  ++(j);
               }
            n = prtclbase[pn[i-1]-1].anti;

            if (prtclbase[n-1].top == NULL)
            {
               prtclbase[n-1].top = (decaylink)m_alloc(sizeof(modeofdecay));
               prtclbase[n-1].top->next = NULL;
               memcpy(prtclbase[n-1].top->part,cc,3*sizeof(particleNumType));
            }
            else
            {
               qq = prtclbase[n-1].top;
               while (1)
               {
                  k = 1;
                  while (k < 4 && qq->part[k-1] == cc[k-1]) ++(k);
                  if (k == 4) goto exi;
                  if (qq->part[k-1] > cc[k-1])
                  {
                     kk = (decaylink)m_alloc(sizeof(modeofdecay));
                     kk->next = qq->next;
                     qq->next = kk;
                     memcpy(kk->part,qq->part,3*sizeof(particleNumType));
                     memcpy(qq->part,cc,3*sizeof(particleNumType));
                     goto exi;
                  }

                  if (qq->next == NULL)
                  {
                     kk = (decaylink)m_alloc(sizeof(modeofdecay));
                     kk->next = qq->next;
                     qq->next = kk;
                     memcpy(kk->part,cc,3*sizeof(particleNumType));
                     goto exi;
                  }
                  qq = qq->next;
               }
exi:;
            }
         }
      lgrgn1 = lgrgn1->next;
   }  while (lgrgn1 != NULL);
}


int  loadModel(int check,int ugForce)
{ varrec recForGG;
  errorText[0]=0;
  if( (!check)&&(lastModel == n_model) ) return 1;
  if( !readvars(check,&recForGG) )       if(blind) sortie(125); else return 0; 
  if( !readparticles(check,ugForce))     if(blind) sortie(125); else return 0;
  if( recForGG.need!=-1)   
  {
      nmodelvar++;
      *(modelvars+nmodelvar)=recForGG;
  }
  if( !readlagrangian(check,ugForce))    if(blind) sortie(125); else return 0;
  filldecaylist();
  readEXTLIB();
  lastModel = n_model;
  return 1;
}

void readModelFiles(int l)
{  char fname[100];
   char ext[5][10]={"vars","func","prtcls","lgrng","extlib"};
   int i;

   for(i=0;i<5;i++) cleartab(modelTab+i);
   lastModel=0;
   for(i=0;i<5;i++)
   { sprintf(fname,"./models/%s%d.mdl",ext[i],l);
     readtable(modelTab+i,fname);
   }  
}


int read2VarsParticles(void)
{
   varrec recForGG;
   errorText[0]=0;
   blind=1;
  
   if(!readvars(1,&recForGG) ) return 1;
   if(!readparticles(1,0))     return 2;
   return 0;
}
