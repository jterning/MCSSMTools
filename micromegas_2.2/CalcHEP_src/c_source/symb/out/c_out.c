/*
 Copyright (C) 1997, Alexander Pukhov 
*/

#include <unistd.h>
#include "chep_crt.h"
#include "getmem.h"
#include "syst2.h" 
#include "physics.h"
#include "s_files.h"
#include "procvar.h"
#include "pvars.h"
#include "diaprins.h"
#include "optimise.h"
#include "l_string.h"
#include "parser.h"
#include "reader_c.h"
#include "out_serv.h"
#include "saveres.h"
#include "denominators.h"
#include "process.h"
#include "cweight.h"
#include "c_out.h"
#include "writeF.h"

typedef char s6[7];

/* ======================================================== */
typedef char shortname[6];
typedef shortname prtclsarray[MAXINOUT];

int noPict=0;
int noCChain=0;
int tWidths=0;
static int sumDiag=0;
static int flag_12=0;

static void  getprtcls(char* txt1,prtclsarray pnames)
{  int   j; 
   char   txt[STRSIZ]; 
   char*c0,*c1;
   
   strcpy(txt,txt1);
   c0=strstr(txt,"->"); 
   c0[0]=','; 
   c0[1]=' '; 

   for(j=nin+nout; j<MAXINOUT; j++) strcpy(pnames[j],"***"); 

   for(c0=txt,j=0; ;j++)
   { 
      c1=strchr(c0,',');      
      if(c1) c1[0]=0;
      trim(c0);
      strcpy(pnames[j],c0);
      if(!c1) return; else c0=c1+1; 
   }
} 


/*===========================================================*/


#define procinfptr struct procinfrec *
typedef struct procinfrec
   {
      procinfptr     next;
      unsigned       tot;
      unsigned       firstdiagpos;
      prtclsarray    p_name;
      int            p_masspos[MAXINOUT];
      long           p_code[MAXINOUT];
   }  procinfrec;
#undef procinfptr

typedef struct procinfrec *procinfptr;

static procinfptr   inf, inftmp;  /*  information about subProcess  */

       /*  statictics  */
static unsigned  ndiagrtot,  diagrcount;
static int  nvars,  nfunc;

static marktp   heapbeg;/*  for RELEASE  */

static int   nden_w, nden_s, nden_t, nden_0,  nsub1; /* from writesubprocess */


static int cBasisPower;
static int nC, *cChains=NULL;
static long  *cCoefN, *cCoefD;


static void clearstatistic(void)
{int  i; for (i = 17; i < 24; i++) { goto_xy(1,i); clr_eol();} }

static void init_stat(void)
{
   goto_xy(1,17);
   scrcolor(Yellow,Blue);
   print(" C Source Codes \n");
   scrcolor(Red,BGmain);
   print(" Process..........\n");
   print(" Total diagrams...\n");
   print(" Processed........\n");
   print(" Current..........\n");
   scrcolor(Yellow,Blue);
   print(" Press Esc to stop    ");
   scrcolor(Black,BGmain);
   goto_xy(20,18); print("%s",processch);
   goto_xy(20,19); print("%4u",ndiagrtot);
   goto_xy(20,20); print("   0");
   scrcolor(Yellow ,BGmain);
   goto_xy(20,21); print("   1");
   scrcolor(Yellow,BGmain);
}


static void writestatistic(void)
{
   scrcolor(Black ,BGmain);
   goto_xy(20,19); print("%4u",ndiagrtot);
   goto_xy(20,20);
   print("%2u (%%)",(((diagrcount - 1) * 100) / ndiagrtot));
   goto_xy(20,21); print("%4u",diagrcount);
}


static void writpict(unsigned ndiagr)
{  vcsect vcs;
   csdiagram  csdiagr;
   fseek(diagrq,ndiagr*sizeof(csdiagr),SEEK_SET);
   FREAD1(csdiagr,diagrq);
   transfdiagr(&csdiagr,&vcs);
   writeF("/*\n");
   DiagramToOutFile(&vcs,1,' ');
   writeF("*/\n");
}


static void labl(void)
{
  writeF("/*******************************\n");
  writeF("*    %s*\n",version);
  writeF("*******************************/\n");
}

  /* =========== Preliminary  calculations ================ */

static void calc_nvars_nfunc(void)
{ int   k;
  for(nvars=0, nfunc=0, k=1; k<=nmodelvar; k++) if(vararr[k].used)
  { if( modelvars[k].func || modelvars[k].pwidth) nfunc++; else  nvars++; }
} 


static void prepareprocinform(void)
{int ndel, ncalc, nrest;
 long recpos;
 char        txt[STRSIZ];
 int     i, k;
 csdiagram   csd;
 s6          mass;
 int        nn;
 int nsubs;

   inf = NULL;
   fseek(menuq,0,SEEK_SET);
   for (nsubs=1;nsubs<=subproc_sq;nsubs++)
   {
      inftmp = inf;
      inf = (procinfptr)getmem_((unsigned)sizeof(procinfrec));
      inf->next = inftmp;
      rd_menu(2,nsubs,txt,&ndel,&ncalc,&nrest,&recpos);
      inf->firstdiagpos = recpos;
      getprtcls(txt,inf->p_name);
      for (i = 0; i < nin + nout; i++)
      {
         locateinbase(inf->p_name[i],&nn);
         strcpy(mass,prtclbase[nn-1].massidnt);
         if (strcmp(mass,"0")) for(k=1; strcmp(modelvars[k].varname,mass); k++); 
         else k=0;
         inf->p_masspos[i] = k;
         inf->p_code[i]=prtclbase1[nn].N; 
      }
      for (i = nin + nout; i < MAXINOUT; i++)
      {
         strcpy(inf->p_name[i],"***");
         inf->p_masspos[i] = 0;
      }
      fseek(diagrq,recpos*sizeof(csdiagram),SEEK_SET);
      inf->tot = 0;
      for (i = 1; i <= ndel + ncalc + nrest; i++)
      {
         FREAD1(csd,diagrq);
         if (csd.status == 1) ++(inf->tot);
      }
      if(inf->tot==0) for (i = 0; i < nin + nout; i++) inf->p_masspos[i]=0;
   }
   nsubs--;
   revers((void **)&inf);
}


static void scanvars(int mode)
{
 int k;

   for(k=1;k<=nmodelvar;k++)
   {  varlist  modl=modelvars+k;  
      if (vararr[k].used && !modl->func && !modl->pwidth)
      {
         switch(mode)
         {
	    case 2: writeF("\n,%E",modl->varvalue);    break;
            case 3: writeF("\n,\"%s\"",modl->varname); break;
         }
      }
   }

   
   for(k=1;k<=nmodelvar;k++)
   {  varlist  modl=modelvars+k;  
      if (vararr[k].used && (modl->func || modl->pwidth ))
      {
         switch(mode)
         {
	    case 2: if(modl->hidden) writeF("\n,1"); else  writeF("\n,0");    break;
            case 3: writeF("\n,\"%s\"",modl->varname); break;
         }
      }
   }
}

/* ======= Information functions =========== */

static void geninf(char* name,int value)
{  if(flag_12) writeF("static int %s = %d;\n\n",name, value);
   else        writeF("const int %s = %d;\n\n",name, value);
}


static void  writesubroutineinit(void)
{
   int         l;
   char        *ss;
   
   if(!flag_12) writeF("#include\"extern.h\"\n");
   ext_h=fopen("results/extern.h","w");   
   fprintf(ext_h,"extern double width2_ext(char*,long );\n");
   writeF("static int calcFunc_stat(void)\n{\n");
   writeF(" double * V=va_ext;\n");
   writeF(" FError=0;\n");
   for(l=1;l<=nmodelvar;l++)
   {
      if (vararr[l].used &&  (modelvars[l].func||modelvars[l].pwidth))
      {  int num;
         checkNaN=0;
         if(modelvars[l].pwidth)
         {  writeF("   %s=width2_ext(\"%s\",%d);\n",vararr[l].alias,
             prtclbase1[modelvars[l].pwidth].name, prtclbase1[modelvars[l].pwidth].N); 
             checkNaN=1;
         } else
         {
           ss=(char *)readExpression(modelvars[l].func,rd_c,act_c,free);
/*	   writeF("   %s=%s;\n",vararr[l].alias,ss+3);*/
	   fprintf(outFile,"   %s=%s;\n",vararr[l].alias,ss+3);
	   free(ss);
	 }
	 if(checkNaN)
	 {
	 sscanf(vararr[l].alias,"va_ext[%d]",&num); 
         writeF("   if(!finite(%s) || FError) return %d;\n",vararr[l].alias,num);
         }
      }
   }

   writeF("return 0;\n}\n"); 
   fclose(ext_h); 
}

static void  writeDenominators(deninforec* dendescript)
{   
int i,k;
int nden_t2=nden_t, nden_s2=nden_s;

   for (i = 1; i <= dendescript->tot_den; i++)
   {   int numm = dendescript->denarr[i-1].order_num; 
      if(dendescript->denarr[i-1].width==0)  numm += nden_w; else  
      { if(dendescript->denarr[i-1].stype)  nden_s2--; 
        else { numm += nden_s; nden_t2--;}
      } 
      if(dendescript->denarr[i-1].power == 1) writeF("*Q1[%d]",numm);
                                         else writeF("*Q2[%d]",numm);
   }

   writeF(";\n");

   if(nden_s2+nden_t2)
   {   
      if(nden_s2) 
      { 
        writeF("if(gswidth_ext )");    
        writeF(" R=R");  
        for (k = 1; k <= nden_s; k++)					  	
        { int addpr = 1;							       
          for (i = 1; i <= dendescript->tot_den; i++)			       
          if (dendescript->denarr[i-1].width &&	dendescript->denarr[i-1].stype&&			       
   	     k == dendescript->denarr[i-1].order_num)  addpr = 0;	  	
          if (addpr)  writeF("*Q0[%d]",k);
        }
      writeF(";\n");
      }  
      if(nden_t2) 
      {  writeF("if(gtwidth_ext) R=R"); 
 
        for (k = nden_s+1; k <= nden_w; k++)					  	
        { int addpr = 1;							       
          for (i = 1; i <= dendescript->tot_den; i++)			       
          if (dendescript->denarr[i-1].width &&	(!dendescript->denarr[i-1].stype) &&			       
   	     k == dendescript->denarr[i-1].order_num+nden_s)  addpr = 0;	  	
          if (addpr)  writeF("*Q0[%d]",k);
        }
       writeF(";\n");
      }      					       
   }									       
}

static void calcColor(long diag)
{
   csdiagram  csdiagr;

   fseek(diagrq, (diag-1)*sizeof(csdiagr),SEEK_SET);
   FREAD1(csdiagr,diagrq);

   if(cBasisPower&&generateColorWeights(&csdiagr,cBasisPower,nC,cChains,cCoefN,cCoefD))
   { int k;
      writeF(" if(cb_coeff_ext)\n {\n");
     
      for(k=0; k<cBasisPower; k++) 
      if(cCoefN[k]) writeF("  cb_coeff_ext[%d] += (R*%d)/(%d);\n",k,
      cCoefN[k],cCoefD[k]);

      writeF(" }\n");
   }
}



static void  onediagram(deninforec* dendescript)
{  catrec      cr;
   marktp      bh;
   varptr      totnum, totdenum, rnum;
   long pos_c;
   int deg1,nConst;

   mark_(&bh);
   tmpNameMax=0;
   initinfo();
   initdegnames();
   
   fseek(catalog,dendescript->cr_pos,SEEK_SET);
   FREAD1(cr,catalog);
   ++(diagrcount);
   whichArchive(cr.nFile,'r',flag_12);

   fseek(archiv,cr.factpos,SEEK_SET);
    
   readvardef(archiv);
   readpolynom(&totnum);
   readpolynom(&totdenum);
   clearvardef();
   
   fseek(archiv,cr.rnumpos,SEEK_SET);

   readvardef(archiv);
   readpolynom(&rnum);
   clearvardef();

   if(!flag_12)
   {  outFileOpen("%sresults%cf%d.c",pathtouser,f_slash,diagrcount);
      labl();
      writeF("#include\"num_out.h\"\n");
      writeF("#include\"num_in.h\"\n");   
   }
   writeF("extern double F%d_ext(void);\n",diagrcount);
   writeF("static void C%d(double * C)\n{\n",diagrcount);
   writeF("double* V=va_ext;\n");
   pos_c= ftell(outFile); writeF("%80s\n",""); 
   nConst=write_const();
   deg1=cleardegnames();       
   writeF("}\n"); 

   fseek(outFile,pos_c,SEEK_SET);
   if(deg1) writeF("double S[%d];",deg1);
   if(tmpNameMax) writeF("double tmp[%d];",tmpNameMax );                                

   fseek(outFile,0,SEEK_END);
   tmpNameMax=0;
   initdegnames();

   writeF("double F%d_ext(void)\n{\n",diagrcount);

   if(!noPict) writpict(cr.ndiagr_ + inftmp->firstdiagpos - 1);

   writeF("double N,D,R;\n");
   pos_c= ftell(outFile); writeF("%80s\n","");
  
   writeF("if(CalcConst) C%d(C);\n",diagrcount);
   
   fortwriter("N",totnum);
   fortwriter("D",totdenum);
   fortwriter("R",rnum);
   
   writeF("R*=(N/D)");
   writeDenominators(dendescript);
 
   if(!noCChain)calcColor(cr.ndiagr_+inftmp->firstdiagpos);

   writeF(" return R;\n");  
   writeF("}\n");

   deg1=cleardegnames();
   if(nConst==0) nConst=1;
   fseek(outFile,pos_c,SEEK_SET);
   writeF("static double C[%d];",nConst);
   if(deg1) writeF("double S[%d];",deg1);
   if(tmpNameMax) writeF("double tmp[%d];",tmpNameMax );
   fseek(outFile,0,SEEK_END);
   if(!flag_12) outFileClose();
   release_(&bh);
}


static int  alldiagrams(FILE * fd,  int nsub)
{  
   marktp     bh;
   varptr     totnum, totdenum, rnum;
   long       pos_c1,pos_c2; int deg1,deg2,tmpn1,tmpn2, nC;
   catrec     cr;
   deninforec dendescript;

   mark_(&bh); tmpNameMax=0; initinfo(); initdegnames();

   writeF("{\n");
   pos_c1= ftell(outFile); writeF("%70s\n","");   
   writeF("double N,D,R;\n");
   writeF("if(CalcConst) C%d(C);\n",nsub);

   while(FREAD1(dendescript,fd) == 1)
   {
      fseek(catalog,dendescript.cr_pos,SEEK_SET);
      FREAD1(cr,catalog); ++(diagrcount);
      if(!noPict)writpict(cr.ndiagr_ + inftmp->firstdiagpos - 1);
      whichArchive(cr.nFile,'r',flag_12);
      fseek(archiv,cr.factpos,SEEK_SET);

      readvardef(archiv);
      readpolynom(&totnum);
      readpolynom(&totdenum);
      clearvardef();
   
      fseek(archiv,cr.rnumpos,SEEK_SET);

      readvardef(archiv);
      readpolynom(&rnum);
      clearvardef();
   
      fortwriter("N",totnum);
      fortwriter("D",totdenum);

      fortwriter("R",rnum);

      writeF("R*=(N/D)");
      if(nin+nout>3)
      {
         writeDenominators(&dendescript);
         if(!flag_12) writeF(" if(R>Fmax) Fmax=R; else if(R<-Fmax) Fmax=-R;\n");
      } else  writeF(";\n");
      if(!flag_12) if(!noCChain)calcColor(cr.ndiagr_+inftmp->firstdiagpos);
      writeF("ans+=R;\n");
   }   
   whichArchive(0,0,0);
   writeF("\n}\nreturn ans;\n}\n");

   deg1=cleardegnames();
   tmpn1=tmpNameMax;
   tmpNameMax=0;
   initdegnames();

   writeF("\nstatic void C%d(double*C)\n{\n",nsub); 
   writeF("  double* V=va_ext;\n");
   pos_c2= ftell(outFile); writeF("%70s\n","");   

   nC=write_const(); 
   if(nC==0) nC=1; 
   writeF("}\n");

   fseek(outFile,pos_c1,SEEK_SET); 
   writeF("static double C[%d];",nC); 
   if(deg1) writeF("double S[%d];",deg1);
   if(tmpn1) writeF("double tmp[%d];",tmpn1); 

 
   fseek(outFile,pos_c2,SEEK_SET);
   deg2=cleardegnames();
   tmpn2=tmpNameMax;
   if(deg2) writeF("double S[%d];",deg2) ;
   if(tmpn2) writeF("double tmp[%d];",tmpn2 );
   fseek(outFile,0,SEEK_END);

   release_(&bh);
   if( escpressed()) return 1; else return 0;
}



static void  writesubprocess(int nsub,long firstDiag,long totDiag,int* breaker)
{  denlist    den_;
   long      i;
    
   deninforec   dendescript;
   FILE * fd;                /* file of (deninforec)  */
   char fd_name[STRSIZ];
   marktp mem_start;

   nsub1 = nsub;
   if(!flag_12)
   { outFileOpen("%sresults%cd%d.c",pathtouser,f_slash,nsub);
     labl();
     writeF("#include\"num_in.h\"\n");
     writeF("#include\"num_out.h\"\n");
   }

   if(totDiag==0) 
   { if(flag_12) writeF("static ");else writeF("extern DNN S%d_ext;\n",nsub); 
     writeF("double S%d_ext(double * momenta,int * err)\n{",nsub); 
     writeF("  return 0;\n}\n");
     if(!flag_12) outFileClose(); 
     return;
   }

   if(sumDiag) writeF("static void C%d(double *);\n",nsub); else 
   {  writeF("extern FNN F%d_ext",firstDiag);
      for(i=1;i<totDiag;i++) writeF(",F%d_ext",i+firstDiag);
      writeF(";\n");
      writeF("static FNN *Farr[%d]={&F%d_ext",totDiag,firstDiag);
      for(i=1;i<totDiag;i++) writeF(",&F%d_ext",i+firstDiag);
      writeF("};\n");
   } 
   if(flag_12) writeF("static "); else writeF("extern DNN S%d_ext;\n",nsub);
   writeF("double S%d_ext(double * momenta,int * err)\n{",nsub);
   writeF("double  ans=0;\n");
   if(nin+nout>3) 
   {  
/*    writeF("double s0max=0;\n");*/
      writeF("int i;\n");
   }
   sprintf(fd_name,"%stmp%cden.inf",pathtouser,f_slash);
   fd=fopen(fd_name,"wb"); 

   mark_(&mem_start);
   denominatorStatistic(nsub, &nden_s, &nden_t, &nden_0, &den_, fd,flag_12); 
   fclose(fd);
   nden_w=nden_s+nden_t;

   if(nin+nout>3)
   {  int nden= nden_w+nden_0+1; 
      writeF("double mass[%d],width[%d];\n",nden,nden);
      writeF("double* V=va_ext;\n");
      writeF("char * Qtxt[%d];\n",nden);
/*      writeF(" for(i=0;i<nin_ext;i++) s0max+=momenta[4*i];\n"); */
   
      for(;den_;den_ = den_->next)
      {  int m=0;
         i=den_->order_num;
         if(den_->width) 
         {
           if(den_->stype) fprintf(outFile,"width[%d]=%s; ",i,vararr[den_->width].alias);
           else 
           {  i+=nden_s;
             fprintf(outFile,"width[%d]=(twidth_ext)? %s : 0.; ",i,vararr[den_->width].alias); 
           }
         }else 
         { i+=nden_w;
           writeF("width[%d]=0.; ",i);
         }
         fprintf(outFile,"mass[%d]=%s; ",i,vararr[den_->mass].alias);
         fprintf(outFile," Qtxt[%d]=\"",i);       
/*         fprintf(outFile," Q[%d]=mass[%d]*mass[%d]-sqrMom(nin_ext,\"",i,i,i);*/

         while(den_->momStr[m]) fprintf(outFile,"\\%o",den_->momStr[m++]);
         fprintf(outFile,"\";\n");
/*         fprintf(outFile,"\",momenta);\n");   */ 
      }  
      writeF("*err=*err|prepDen(%d,nin_ext, BWrange_ext*BWrange_ext,mass,width, Qtxt,momenta);\n",
      nden_w+nden_0);
   }
   release_(&mem_start);
   fd=fopen(fd_name,"rb"); 
   if(sumDiag)       
   {   
     *breaker = alldiagrams(fd,nsub); 
     writestatistic();   
     if(!flag_12) outFileClose();    
   } else 
   {
      writeF("for(i=0;i<%d;i++) \n",totDiag);
      writeF(
      "{ double r=Farr[i]();\n"
      "  if(r>Fmax) Fmax=r;\n"
      "  ans+=r;\n"
      "}\n"
      "return ans;\n}\n"
            );

      if(!flag_12) outFileClose();

      *breaker = 0;
      while(FREAD1(dendescript,fd) == 1)
      {
         if (escpressed())
         {  *breaker = 1;
            break;
         }
         onediagram(&dendescript);
         writestatistic();
      } 
   }
   fclose(fd);
   unlink(fd_name);
}  /*  WriteSubprocess  */



static void  make_pinf(void)
{
   int    i;
   if(flag_12) writeF("static ");
   writeF("char * pinf_ext(int nsub,int nprtcl,double* pmass,long * num)\n{\n");
   writeF("int n;\n");

   writeF(" static char *names[%d][%d] ={\n",subproc_sq,nin + nout);
   inftmp = inf;
   for (nsub = 1; nsub <= subproc_sq; nsub++)
   {  writeF("{");
      for (i = 1; i <= nin + nout; i++)
      {  if(i!=1) writeF(",");
         writeF("\"%s\"",inftmp->p_name[i-1]);
      }
      if (nsub== subproc_sq) writeF("}};\n"); else  writeF("},\n");
      inftmp = inftmp->next;
   }
   
   writeF("int const nvalue[%d][%d]={\n",subproc_sq,nin + nout);
   inftmp = inf;
   for (nsub = 1; nsub <= subproc_sq; nsub++)
   {  writeF("{");
      for (i = 1; i <= nin + nout; i++)
      {  int k=inftmp->p_masspos[i-1];
         if(k) sscanf(vararr[k].alias,"V[%d]",&k);
         if(i!=1) writeF(",");  
	 writeF("%d",k);
      }
      if (nsub== subproc_sq) writeF("}};\n"); else  writeF("},\n"); 
      inftmp = inftmp->next;
   }

   writeF("long const pcode[%d][%d]={\n",subproc_sq,nin + nout);
   inftmp = inf;
   for (nsub = 1; nsub <= subproc_sq; nsub++)
   {  writeF("{");
      for (i = 0; i < nin + nout; i++)
      { 
         if(i) writeF(",");  
	 writeF("%d",inftmp->p_code[i]);
      }
      if (nsub== subproc_sq) writeF("}};\n"); else  writeF("},\n"); 
      inftmp = inftmp->next;
   }

   writeF("if  (nsub<0 ||nsub>%d||nprtcl<0||nprtcl>%d) return NULL;\n",
   subproc_sq,nin + nout);
   writeF("if(pmass)\n{\n");
   writeF("  n=nvalue[nsub-1][nprtcl-1];\n");
   writeF("  if (n==0) *pmass=0; else *pmass=va_ext[n];\n"); 
   writeF("  if (*pmass<0) (*pmass)=-(*pmass);\n");  
   writeF("}\n");  
   writeF("if(num)*num=pcode[nsub-1][nprtcl-1];\n");
   writeF("return names[nsub-1][nprtcl-1];\n}\n");

   if(nin==1) writeF("char * polarized_ext[3]={\"\",\"\",\"\"};\n");
   else 
   {  writeF("char * polarized_ext[3]={\"\",\",");       
      for (i = 1; i <= nparticles; i++)
      if(polarized(1,i)) writeF("%s,",prtclbase1[i].name);
      writeF("\",\",");
      for (i = 1; i <= nparticles; i++)
      if(polarized(2,i)) writeF("%s,",prtclbase1[i].name);
      writeF("\"};\n");
   }

   if(flag_12) writeF("static ");
   writeF("long pinfAux_ext(int nsub,int nprtcl,int*spin2,int*color,int*neutral)\n{\n");
   writeF("int n;\n");

   writeF("long const pcode[%d][%d]={\n",subproc_sq,nin + nout);

   inftmp = inf;
   for (nsub = 1; nsub <= subproc_sq; nsub++)
   {  writeF("{");
      for (i = 0; i < nin + nout; i++)
      { 
         if(i) writeF(",");  
	 writeF("%d",inftmp->p_code[i]);
      }
      if (nsub== subproc_sq) writeF("}};\n"); else  writeF("},\n"); 
      inftmp = inftmp->next;
   }

   writeF("int const Spin2[%d][%d]={\n",subproc_sq,nin + nout);

   inftmp = inf;
   for (nsub = 1; nsub <= subproc_sq; nsub++)
   {  writeF("{");
      for (i = 0; i < nin + nout; i++)
      { int pos;
        if(i) writeF(",");  
	locateinbase(inftmp->p_name[i], &pos);
         writeF("%d",prtclbase1[pos].spin);	
      }
      if (nsub== subproc_sq) writeF("}};\n"); else  writeF("},\n"); 
      inftmp = inftmp->next;
   }

   writeF("int const Color[%d][%d]={\n",subproc_sq,nin + nout);

   inftmp = inf;
   for (nsub = 1; nsub <= subproc_sq; nsub++)
   {  writeF("{");
      for (i = 0; i < nin + nout; i++)
      { int pos;
        if(i) writeF(",");  
	locateinbase(inftmp->p_name[i], &pos);
         writeF("%d",prtclbase1[pos].cdim);	
      }
      if (nsub== subproc_sq) writeF("}};\n"); else  writeF("},\n"); 
      inftmp = inftmp->next;
   }

 writeF("int const Neutral[%d][%d]={\n",subproc_sq,nin + nout);

   inftmp = inf;
   for (nsub = 1; nsub <= subproc_sq; nsub++)
   {  writeF("{");
      for (i = 0; i < nin + nout; i++)
      { int pos;
        if(i) writeF(",");  
	locateinbase(inftmp->p_name[i], &pos);
	if(pos==prtclbase1[pos].anti)
         writeF("1");	else  writeF("0"); 
      }
      if (nsub== subproc_sq) writeF("}};\n"); else  writeF("},\n"); 
      inftmp = inftmp->next;
   }
   
   writeF("if(nsub<0 ||nsub>%d||nprtcl<0||nprtcl>%d) return 0;\n",
   subproc_sq,nin + nout);
   writeF("if(spin2) *spin2=Spin2[nsub-1][nprtcl-1];\n");
   writeF("if(color) *color=Color[nsub-1][nprtcl-1];\n");
   writeF("if(neutral) *neutral=Neutral[nsub-1][nprtcl-1];\n");
   writeF("return pcode[nsub-1][nprtcl-1];\n}\n");   
}

static void  make_den_info(void)
{  int nden_s,nden_t,nden_0;

   writeF("\n char * den_info_ext(int nsub,int n, int * mass, int * width)\n{\n");
   writeF(" switch(nsub){\n");

   for (nsub = 1; nsub <= subproc_sq; nsub++)
   {  int n=1;
      marktp mem_start; 
      denlist    den_;
 
      writeF(" case %d: switch(n){",nsub);
      mark_(&mem_start);
      denominatorStatistic(nsub, &nden_s, &nden_t, &nden_0, &den_, NULL,flag_12); 
      for(n=1 ;den_;den_ = den_->next,n++)
      { int m=0;
          writeF("\n    case %d: *mass=%d; *width=%d; return \"",
          n, vararr[den_->mass].num, vararr[den_->width].num);
         while(den_->momStr[m]) fprintf(outFile,"\\%o",den_->momStr[m++]);
         fprintf(outFile,"\";");    
      }  
      writeF("\n    default:*mass=0; *width=0; return NULL;\n                  }\n");
 
      release_(&mem_start); 
   }
   writeF("   default: *mass=0; *width=0; return NULL;\n            }\n}\n");
}


static void  make_infbasis(void)
{
   int    i,j;
   int pcolor[MAXINOUT];

   writeF("static void cStrings(int nsub,int *nC, int * power, int **  chains)\n{\n");
   writeF("   switch(nsub)\n   {\n");
   
   
   if(!noCChain)
   for (nsub = 1, inftmp = inf; nsub <= subproc_sq; nsub++)
   {  
      writeF("   case %d : ",nsub);

      for (i = 0; i < nin + nout; i++)
      {  int l;
         locateinbase(inftmp->p_name[i], &l);
         pcolor[i]=prtclbase[l-1].cdim;
         if(i<nin)
         { if(pcolor[i]==3) pcolor[i]=-3; else if(pcolor[i]==-3) pcolor[i]=3;} 
      }

      infCbases(nin+nout,pcolor,&nC,&cBasisPower,&cChains);
      if(cBasisPower)
      {   writeF("\n     { static int cc[%d]=\n       {\n",2*nC*cBasisPower);
          for(i=0;i<cBasisPower;i++)
          {  writeF("       ");
             for(j=0; j<nC; j++)  
             { writeF(" %d,%d",cChains[2*i*nC+2*j]+1,cChains[2*i*nC+2*j+1]+1);
               if(i==cBasisPower-1 && j==nC-1) writeF("\n       };");
               else  writeF(",");
             }
             writeF("\n");
          }  
              
          writeF("       *nC=%d; *power=%d; *chains=cc;\n     }\n     break;\n",
                      nC,cBasisPower);
          free(cChains); cChains=NULL;                  
      } else  writeF("   *nC=0; *power=0; *chains=NULL; break;\n"); 
          
      inftmp = inftmp->next;
      
   }
   writeF("   default: *nC=0; *power=0; *chains=NULL;\n");
   writeF("   }\n");

   writeF("}\n\n");    
}

static void  make_vinf(void)
{  
  if(flag_12) writeF("static ");
  writeF("char * varName_ext[%d]={\"P(cms)\"",nvars+nfunc+1);
  scanvars(3);
  writeF("};\n");
}

static void zeroHeep(void)
{ goto_xy(1,1);print("Heep is empty!!!");inkey();
  sortie(70);
}


static int c_prog_int(int for12)
{
   int breaker;
   int i;
   long dfirst;
      
   if(nin+nout<=4) sumDiag=1; else sumDiag=0;   
   flag_12=for12;   

   memerror=zeroHeep;
   mark_(&heapbeg);

   if(for12)  initvararray(0,'c',2 ); else initvararray(0,'c',3);
  /* ======= Initialisation parth ======= */

   firstVar=nmodelvar;
   if(!strcmp( modelvars[firstVar].varname,strongconst))  firstVar--;
   prepareprocinform();
   calc_nvars_nfunc();
  /* ======= End of Initialisation ====== */
   if(!flag_12) 
   {  outFileOpen("%sresults%cservice.c",pathtouser,f_slash); 
      labl();
      writeF("#include<math.h>\n");                 
      writeF("#include\"num_out.h\"\n");
      writeF("#include\"num_in.h\"\n");

      writeF("double BWrange_ext=2.7;\n");
      writeF("int twidth_ext=0;\n");
      writeF("int gtwidth_ext=0;\n");
      writeF("int gswidth_ext=0;\n");
      writeF(" double va_ext[%d] ={0",nvars+nfunc+1);
      scanvars(2);
      writeF("};\n"); 
   } else writeF(" static double va_ext[%d];\n",nvars+nfunc+1);  
   geninf("nin_ext",nin);
   geninf("nout_ext",nout);
   geninf("nprc_ext",subproc_sq);
   make_pinf();
   geninf("nvar_ext",nvars);
   geninf("nfunc_ext",nfunc);
    
   make_vinf();
   
   if(!flag_12)
   { 
      make_den_info();
      writeF("\nCalcHEP_interface interface_ext={ %d, %d, varName_ext,va_ext,"
          "%d, %d, %d, &pinf_ext, &pinfAux_ext, polarized_ext, &calcFunc_ext, &BWrange_ext,&twidth_ext,"
          "&gtwidth_ext,&gswidth_ext, &sqme_ext,&den_info_ext,&build_cb_ext, &cb_pow_ext,"
          "&cb_nc_ext, &cb_chains_ext, &cb_coeff_ext, &destroy_cb_ext};\n", 
      nvars, nfunc, nin,nout,subproc_sq);

      writeF("\nCalcHEP_interface * PtrInterface_ext=&interface_ext;\n");

      outFileClose();
      outFileOpen("%sresults%csqme.c",pathtouser,f_slash); 
      labl();
      writeF("#include<stdio.h>\n");
      writeF("#include\"num_out.h\"\n");
      writeF("#include\"num_in.h\"\n");
   }
   writeF("static int calcall[%d];\n",subproc_sq+1);
   if(!flag_12)
   {
   writeF("static int particles[%d]={0",1+nin+nout); 
   for(i=0;i<nin+nout;i++) writeF(",0");
   writeF("};\n");
   }
   if(flag_12) writeF("static DNN "); else writeF("extern DNN ");
   for(i=1;i<subproc_sq;i++)  writeF("S%d_ext,",i); 
   writeF("S%d_ext;\n",subproc_sq); 
   
   writeF("static  DNN * darr[%d]={",subproc_sq);
   for(i=1;i<subproc_sq;i++)  writeF("&S%d_ext,",i);
   writeF("&S%d_ext};\n",subproc_sq);
   fseek(catalog,0,SEEK_END);
   ndiagrtot = ftell(catalog)/sizeof(catrec);
   
   writesubroutineinit();
   
   if(!flag_12)
   {  make_infbasis();
      writeF("#include\"sqme.inc\"\n");
      outFileClose();
   }
   diagrcount = 0;
   inftmp = inf;
   init_stat();
   for (nsub = 1,dfirst=1; nsub <= subproc_sq; nsub++)
   {  int colors[MAXINOUT];

      if (inftmp->tot != 0)   /*  this subprocess IN archive  */
      {

         for(i=0;i<nin+nout;i++) 
         {  int l;
            locateinbase(inftmp->p_name[i], &l);
            colors[i]=prtclbase[l-1].cdim;
         }
         for(i=0;i<nin; i++) 
         if(colors[i]==3) colors[i]=-3; else if(colors[i]==-3) colors[i]=3;
         if(noCChain) for(i=0;i<nin+nout; i++) colors[i]=1; 
         infCbases(nin+nout,colors,&nC,&cBasisPower,&cChains);
         if(cBasisPower)
         { 
            cCoefN=malloc(cBasisPower*sizeof(long));
            cCoefD=malloc(cBasisPower*sizeof(long));
         }

         writesubprocess(nsub,dfirst,inftmp->tot, &breaker);
         dfirst+=inftmp->tot;
         if (breaker) goto exi;
      
         if(cBasisPower)
         {
            if(cChains){free(cChains); cChains=NULL;} 
            free(cCoefN); free(cCoefD);
         }

      } else writesubprocess(nsub,dfirst,0, NULL);
      inftmp = inftmp->next;
   }
   
exi:
   clearstatistic();
   release_(&heapbeg);
   return !breaker;
}


int  c_prog(void)
{  
 int result;
    
 catalog=fopen(CATALOG_NAME "2","rb");
 if(catalog)
 { 
   int nin_=nin,nout_=nout, subproc_sq_=subproc_sq;
/*   archiv=fopen(ARCHIV_NAME "2","rb"); */
   diagrq=fopen(DIAGRQ_NAME "2","rb");
   menuq=fopen(MENUQ_NAME "2","rb");
   nin=1; nout=2;
   { char txt[100]; int ndel, ncalc, nrest;
     long recpos;
     for(subproc_sq=0;rd_menu(2,subproc_sq+1,txt,&ndel,&ncalc,&nrest,&recpos); subproc_sq++);
   }   
   outFileOpen("%sresults%cwidth_2.c",pathtouser,f_slash);
   labl();
   writeF("#include<math.h>\n");
   writeF("#include<stdio.h>\n");
   writeF("#include\"extern_2.h\"\n");                  
   writeF("#include\"num_out.h\"\n");
   writeF("#include\"num_in.h\"\n");

   result=c_prog_int(1);
   outFileClose();
   system("cd results; sed  s/_ext/_2/g width_2.c > width2.c; "
          "cat width2.c $CALCHEP/include/width_2+ > width_2.c; rm width2.c; "
          "mv extern.h extern_2.h");
            
   fclose(catalog);
   whichArchive(0,0,0);
   fclose(diagrq);
   fclose(menuq);
   nin=nin_; nout=nout_;
   subproc_sq=subproc_sq_;
 } else if(nin+nout>3) 
 { outFileOpen("%sresults%cwidth_2.c",pathtouser,f_slash);
   writeF("extern double width2_ext(char * pname);\n");
   writeF("double width2_ext(char * pname){return 0.;}\n");
   outFileClose();
 }
 catalog=fopen(CATALOG_NAME,"rb"); 
 diagrq=fopen(DIAGRQ_NAME,"rb");
 menuq=fopen(MENUQ_NAME,"rb");

 result=c_prog_int(0);
 fclose(catalog);
 whichArchive(0,0,0);
 
 fclose(diagrq);
 fclose(menuq);

 return result; 
}
