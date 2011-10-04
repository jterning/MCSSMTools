/*
 Copyright (C) 1997, Alexander Pukhov, e-mail: pukhov@theory.npi.msu.su
*/
#include <sys/types.h>
#include <unistd.h>
#include <sys/wait.h>

#include "chep_crt.h"
#include "syst2.h"
#include "physics.h"
#include "screen.h"
#include "lbl.h"
#include "constr.h"
#include "read_mdl.h"
#include "sos.h"
#include "process.h"
#include "s_files.h"
#include "squar.h"
#include "r_code.h"
#include "c_out.h"
#include "red_out.h"
#include "math_out.h"
#include "form_out.h"
#include "symbolic.h"
#include "m_utils.h"
#include "amplitudes.h"
#include "viewdir.h"

#include "saveres.h"
#include "out_serv.h"


static int errorcode=0;
static void*sv_screen=NULL;

static void xw_error(void) { sortie(80);}

static int     width2calc=0;
static char *  inkey2=NULL;
static int     blind2=0;

static char* Plus12code(void)
{
   int l,cont=0;
   char * p=malloc(nparticles);
   char * str=NULL;
   catrec cr;
   int ArcNum=0;
   
   for(l=0; l<nparticles; l++) p[l]=0;  
   catalog=fopen(CATALOG_NAME,"rb");
   while(FREAD1(cr,catalog))
   { char archiveName[40];
     if(ArcNum!=cr.nFile)
     { if(ArcNum!=0) fclose(archiv);
       ArcNum=cr.nFile;
       sprintf(archiveName,ARCHIV_NAME,ArcNum);
       archiv=fopen(archiveName,"rb");
     }  
     fseek(archiv,cr.denompos,SEEK_SET);
     
     readDenominators();
     for(l=0;l< denrno;l++) if(denom[l].width)
     { int np= modelvars[denom[l].width].pwidth;
       if(!np) continue;
       np=ghostmother(np);
       if(p[np-1]==0) {p[np-1]=1; cont++;}
     }
   }
 
   fclose(catalog); 
   if(ArcNum)fclose(archiv);
   if(cont) 
   { str=malloc(30+cont*4);
     strcpy(str,"}}{wildPrt->2*x{");
     for(l=0;l<nparticles;l++) 
          if(p[l])sprintf(str+strlen(str),"%s,",prtclbase[l].name);
     sprintf(str+strlen(str)-1,"{ { {[{[{"); 
   }
   free(p);
   return str;
}


int main(int argc,char** argv)
{   

/*===================================
0 - First start.
1 - Model menu.
2 - Enter process menu.
3 - Feynman diagrams menu; squaring.
4 - Squared diagram menu; symbolic calculation.                        
5 - Write results; new process.
10 -Restart symbolic calculations
==========================================================*/

/* 0-Start; 1-Restart; 2-Heap Error,3-Edit Model,4-UserBreak */

  void *pscr1=NULL,*pscr2=NULL,*pscr3=NULL,*pscr4=NULL,*pscr5=NULL;
  int   k1=1,k2=1,k3=1,k4=1,k5=1;

  int n;
  int pid=0;
  char LOCKtxt[]="Directory 'results/' contains the .lock file created by n_calchep.\n"
                 "To continue you may a)Close the n_calchep session, or b)Rename 'results/',\n";
  blind=0;
  strcpy(pathtocomphep,argv[0]);
  for(n=strlen(pathtocomphep)-1; n && pathtocomphep[n]!=f_slash; n--);
  pathtocomphep[n-3]=0;                                         
     
   for ( n=1;n<argc;n++) 
   { if (strcmp(argv[n],"-blind")==0 && n<argc-1 )
     {  blind=1;
        inkeyString=argv[++n];
     }
     if (strcmp(argv[n],"+blind")==0 )  blind=2;                                     
   }    

   if(!writeLockFile(".lock")) 
   { fprintf(stderr,"locked by other s_calchep. See .lock file\n");
      exit(100);
   }
   strcpy(pathtouser,"");  
   sprintf(pathtohelp,"%shelp%c",pathtocomphep,f_slash);
   outputDir="results/";        
   { char * icon=(char *) malloc(strlen(pathtocomphep)+10);
     sprintf(icon,"%sicon",pathtocomphep);
     start1("CalcHEP/symb",icon,"calchep.ini",&xw_error);
     free(icon);
   }
   fillModelMenu();
   
   f3_key[0]=f3_key_prog;   f3_mess[0]="Model"; 
   f3_key[1]=f4_key_prog;   f3_mess[1]="Diagrams"; 
   f3_key[2]=f5_key_prog;   f3_mess[2]="Switches";
   f3_key[3]=f6_key_prog;   f3_mess[3]="Results"; 
                            f3_mess[4]="Del"; 
                            f3_mess[5]="UnDel";
   f3_key[6]=f9_key_prog;   f3_mess[6]="Ref";    
   f3_key[7]=f10_key_prog;  f3_mess[7]="Quit";

   restoreent(&menulevel);


   if(!blind && menulevel<2) cheplabel(); 

   switch (menulevel)
   {
      case 10: 
      case 6:
      case 5: 
      case 4:
      case 3: readModelFiles(n_model); 
              modelinfo();
              loadModel(0,forceUG);
              processinfo();
              diagramsinfo();
      case 2: k1=n_model;

      case 1: break;
   }

   switch (menulevel)
   {
      case 2:  menuhelp();
               goto label_20;
      case 3:  goto label_31;
      case 4:  goto label_40;
      case 5:  goto label_50;
      case 6:  
      case 10: goto restart2;
   }
        
label_10:   /*   Menu2(ModelMenu): */
   f3_key[0]=NULL; /*models*/ 
   f3_key[1]=NULL; /*diagrams*/
   menulevel = 1;
   forceUG=0;
   menuhelp();
   for(;;)
   { 
      showheap();
      k1=n_model;
      menu1(56,4,"",modelmenu,"s_1",&pscr1,&k1);
      n_model=k1;
      if(n_model == 0)
      {
	if( mess_y_n(56,4,"Quit session")) {  saveent(menulevel); goto exi; }         
      }
      else  if(n_model == maxmodel+1)
      {
         clrbox(1,4,55,18);
         makenewmodel();
         menuhelp();
      }
      else if (n_model > 0)
      { 
	put_text(&pscr1);
	goto label_20;
      }
   }

label_20:   /*  Menu3:Enter Process  */
   f3_key[0]=NULL; 
   f3_key[1]=NULL; 

   menulevel = 2;
   saveent(menulevel);
   readModelFiles(n_model);
   modelinfo();
   k2 = 1;
   do
   {  char strmen[]="\026"
        " Enter Process        "
        " Force Unit.Gauge OFF "
        " Edit model           "
        " Delete model         ";

      if(forceUG)improveStr(strmen,"OFF","ON");
      menu1(56,4,"",strmen,"s_2_*",&pscr2,&k2);
      switch (k2)
      {
         case 0:  goto_xy(1,1); clr_eol(); goto label_10;
         case 2:  forceUG=!forceUG;   modelinfo(); break;
	 case 3:  editModel(1); break;
         case 4: 
	    if ( deletemodel(n_model) )
            {
               goto_xy(1,1);
               clr_eol();
               n_model=1;
               fillModelMenu();
               goto label_10;
            }
            else   readModelFiles(n_model);
      }
   }  while (k2 != 1);

   loadModel(0,forceUG);
label_21:

   f3_key[0]=NULL; 
   f3_key[1]=NULL; 

   menulevel=2;
   errorcode=enter();   /*  Enter a process  */
   newCodes=0;
   showheap();
   if (errorcode)   /*  'Esc' pressed  */ { menuhelp(); goto label_20;}
   errorcode=construct();          /*  unSquared diagrams  */
   if (errorcode) 
   {  if(blind)
      {  if(width2calc) goto absent12;else
         { printf("Processes of this type are absent\n"); sortie(111);}
      } else 
      { messanykey(5,22,"Processes of this type are absent");  
        clrbox(1,19,80,24); 
        goto label_21;
      }
   }
   else if(!blind)
   { int dirStat=checkDir("results"); 
     if(dirStat!=0)
     {  messanykey( 10,10,"There are files in directory 'results/'.\n"  
                          "To continue you has to clean or rename this directory.");
        viewresults();
        if(checkDir("results")!=0)  goto label_21;                   
     }
     clr_scr(FGmain,BGmain);
     modelinfo();
     processinfo();
     diagramsinfo();
     goto label_31;	
   }
   
label_30: /*  Menu4: Squaring,...*/
   clr_scr(FGmain,BGmain);
   modelinfo();
   processinfo();
   diagramsinfo();
label_31: 

   f3_key[0]=f3_key_prog; 
   f3_key[1]=NULL; 

   menulevel=3;  
   k3 = 1;
   do
   {
      menu1(56,4,"","\026"
         " View diagrams        "
/*         " Amplitude calculation" */
         " Squaring technique   "
         " Write down processes ","s_squa_*",&pscr3,&k3);
      switch (k3)
      {
         case 0: clrbox(1,2,55,11); menuhelp(); goto label_20;
         case 1: viewfeyndiag(1);   break;
         case 3: { FILE*f=fopen("results/list_prc.txt","w");
                   int k,ndel,ncalc,nrest;
                   char process[100];
                   long recpos; 
                   menup=fopen(MENUP_NAME,"r");
                   for(k=1;;k++) 
                   { int err=rd_menu(1,k,process,&ndel,&ncalc,&nrest,&recpos);
                     if(!err) break;
                     trim(process);
                     fprintf(f,"%s\n",process);
                   } 
                   fclose(f);
                   fclose(menup);
                    messanykey(20,14,"See file 'results/list_prc.txt'");
 
                 } break; 

/*       case 2: messanykey(10,10,"Not implemented yet"); Amplitudes(); */
      }
   }  while (k3 != 2);      

   if (!squaring()) goto label_30;  /*  process is absent  */

   clear_tmp();

   saveent(menulevel);
   restoreent(&menulevel);  

label_40:   /*  Menu5: View squared diagrams.....   */

   f3_key[0]=f3_key_prog; 
   f3_key[1]=f4_key_prog; 

   menulevel=4;
   clr_scr(FGmain,BGmain);
   modelinfo();
   processinfo();
   diagramsinfo();
   sq_diagramsinfo(); /*   ????????   */

   k4=1;
   saveent(menulevel);
   pscr4=NULL;
   for(;;)   
   {  int res;
      static char* key2str=NULL;
      menu1(56,4,"","\026"
         " View squared diagrams"
         " Symbolic calculations"
         " Make&Launch n_calchep"        
         " Make n_calchep       "
         " REDUCE program       "
         ,"s_calc_*",&pscr4,&k4);

      res=checkDir("results");
      if(res==1)
      {
        int  n_calchep_id=setLockFile("results/.lock");
        if(n_calchep_id) unLockFile(n_calchep_id); else res=2;                          
      }
      switch (k4)
      {  case 0:
            if (mess_y_n(50,3,"Return to previous menu?"))goto label_30;
         break;

         case 1:
            viewsqdiagr();  break;

         case 2:     /*  Compute all diagrams   */
restart2:
            f3_key[0]=f3_key_prog; 
            f3_key[1]=f4_key_prog; 

            menulevel=4;
            calcallproc();
            sq_diagramsinfo();

            if(!continuetest()) break;            
            if(width2calc==0 &&  nin+nout>3 )
            { 
               key2str=Plus12code();  
               if(key2str==NULL) break;  
               saveent(menulevel);
               system("cd tmp; rm -f *2; all=`ls`; mkdir  tmp; mv $all tmp;");
               inkey2=inkeyString;
               blind2=blind;
               inkeyString=key2str;
               width2calc=1;
               blind=1;                  
               get_text(1,1, maxCol(),maxRow(),&sv_screen);
               continue;          
            } else  if(width2calc)
            { 
               system("for FILE in tmp/* ;do if(test ! -d $FILE) then mv  $FILE $FILE'2' ;fi;done;");
absent12:      k4=2;         
               system("mv tmp/tmp/* tmp; rmdir tmp/tmp;");
               restoreent(&menulevel);
               if(key2str){ free(key2str);key2str=NULL;}
               width2calc=0;
               blind=blind2;
               inkeyString=inkey2;
               clr_scr(FGmain,BGmain);
               if(sv_screen)put_text(&sv_screen);
            }
            showheap();  
            put_text(&pscr4);
         break;
         case 3:
            { static char keystr[12]="]{{[{";
              if(res==2) messanykey(3,10,LOCKtxt); else inkeyString=keystr;
            }
            break;
         case 4:
            if(res==2) messanykey(3,10,LOCKtxt); else
            { 
              FILE *f=fopen("results/EXTLIB","w");
              fprintf(f,"EXTLIB=\"%s\"\nexport EXTLIB\n",EXTLIB);
              fclose(f); 
              saveent(menulevel);
              finish();
              sortie(22);
            }
         case 5:  mk_reduceprograms();  break;
      }
      if(k4==2 && continuetest()) { put_text(&pscr4); break; }
   }

label_50:
   k5=1;
   pscr5=NULL;
   menulevel=5;
   saveent(menulevel);
   
   for(;;)  
   {  int n_calchep_id;   
      menu1(56,4,"","\026"
         " C code               "
         "     C-compiler       "
         "     Edit Linker      "
         " REDUCE code          "
         " MATHEMATICA code     "
         " FORM code            "
         " Enter new process    " ,"s_out_*",&pscr5,&k5);
         
      if((k5==1||k5==2)&&pid) 
      { int epid=waitpid(pid,NULL,WNOHANG);
        if(epid) pid=0; else
        { 
          messanykey(3,10,"This option is frozen while n_calchep runs or is compiled");
          continue;
        }
      }  

      switch (k5)
      {  case 0: goto label_40;  break;
         case 1:
           c_prog(); newCodes=0; saveent(menulevel);
           break;
         case 2:
           if(newCodes) { c_prog(); newCodes=0; saveent(menulevel);} 
           pid=fork();
           if(pid==0)
           { 
             char * command;
             n_calchep_id=setLockFile("results/.lock"); 
             if(n_calchep_id)
             {  
               command=(char*)malloc(strlen(pathtocomphep)+strlen(EXTLIB)+100);
               sprintf(command,"EXTLIB=\"%s\"\n"
                               "export EXTLIB\n"
                               "cd results\n"
                               "xterm -e %s/make__n_calchep",
                               EXTLIB,pathtocomphep);                                 
               system(command);
               sprintf(command,"cd results; ./n_calchep"); 
               unLockFile(n_calchep_id);
               system(command); 
               free(command);
             }
             exit(0);
           }
           break;
         case 3: if(edittable(1,1,&modelTab[4],1," ",0))
                 {  char fName[STRSIZ];
                    readEXTLIB();
                    sprintf(fName,"%smodels%c%s%d.mdl",pathtouser,f_slash,mdFls[4],n_model);
                    writetable( &modelTab[4],fName); 
                 }
                 break;
         case 4: makeReduceOutput(); break;
         case 5: makeMathOutput();   break; 
	 case 6: makeFormOutput();   break;
   
         case 7:
            put_text(&pscr5);
            clrbox(1,2,55,11);
            menuhelp();
            goto label_20;
      }
   }

exi:
   finish();
   sortie(0);
   return 0;
}
