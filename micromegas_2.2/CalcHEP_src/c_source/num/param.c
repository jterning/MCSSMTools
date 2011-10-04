/*
 Copyright (C) 2000, Alexander Pukhov, e-mail pukhov@theory.npi.msu.su
*/

#include"chep_crt.h"
#include"interface.h"
#include"plot.h"
#include"param.h"
#include"rw_sess.h"
#include"subproc.h"
#include"num_in.h"

/* ************************************************* */
/* Physics model parameters menu                     */
/* ************************************************* */

int checkParam(void)
{ int err=calcFunc_int();
  if(err>0)
  {  char mess[100];
     sprintf(mess,"Can not evaluate constraned parameter '%s' ",
             varName_int[err]);
    
     if(blind) { printf("%s\n",mess); sortie(122);} 
     else     messanykey(10,10, mess);
     return err; 
  }
  return 0;
}

static void param_menu(int sqtrS_on, int vars_on, int func_on, int * polar, char ** strmen)
{
  int k,pos,npos;
  npos=0;
  if(sqtrS_on)  
  {   npos++; for(k=0;k<2;k++) if(is_polarized(k+1,Nsub))
      {polar[k]=1; npos++;} else polar[k]=0;
  }  else {polar[0]=0;polar[1]=0;}    
  if(vars_on) npos += nvar_int;
  if(func_on) { for(k=0;k<nfunc_int;k++) if(!hiddenf[k])npos++;}    
 
  
  if(npos==0) {*strmen=NULL; return;}
  *strmen=malloc(24*npos+2);
  (*strmen)[0]=24;

  pos=1;    
  for(k=0;k<2;k++) if(polar[k])
  { sprintf((*strmen)+pos,"Helicity%d=%-14.5g",k+1,Helicity[k]);
    pos+=24;               
  }  
  for (k = !sqtrS_on; k <= nvar_int + nfunc_int; ++k)
  { 
     if(k==0 || (vars_on && k<=nvar_int) || (func_on && k>nvar_int &&!hiddenf[k-nvar_int-1]))
     { char c=' '; 
       if(k>nvar_int && vars_on  ) c='*';
       sprintf((*strmen)+pos,"%c%7s= %-14.5g",c,varName_int[k],va_int[k]); 
       pos+=24;
     }
  }
  (*strmen)[pos]=0;    
}

int selectParam(int x, int y, char * mess, void ** pscrPtr,  
    int sqtrS_on, int vars_on, int func_on,
     double ** varPos, char * varName,int*mPos)
{
    char* strmen;
    void * pscr;
    int polar[2];
    int position=1;
    
    if(pscrPtr) pscr=*pscrPtr; else pscr=NULL; 

    param_menu( sqtrS_on,vars_on,func_on, polar,  &strmen);
    if(strmen)
    { int i, shift=0;
      menu1(x,y,mess,strmen,"",&pscr,mPos);
      position=*mPos;
      if(position)
      { sscanf(strmen+1+(position-1)*strmen[0],"%[^=]",varName); 
        trim(varName);
      }
      free(strmen);
      if (position == 0){ if(pscrPtr) *pscrPtr=NULL;  return 0;}
      if(pscrPtr)   *pscrPtr=pscr; else  put_text(&pscr);
      if(sqtrS_on)
      { shift=1; for(i=0;i<2;i++) if(polar[i]) shift++;   
        if(position<=shift) switch(position)
        { case 1: if(polar[0]) {*varPos=&Helicity[0]; return 1;}
          case 2: if(polar[1]) {*varPos=&Helicity[1]; return 1;}
          case 3:               *varPos=va_int;     return 1;
        }
           
        position -= shift;       
      } 
      if(!vars_on)  position+=nvar_int;
      for(i=0; i< (position)-nvar_int ; i++) if(hiddenf[i])(position)++;
      *varPos= va_int+position;
      return 1;
    } else return 0;
} 

int change_parameter(int x,int y, int for22)
{ 
  double val;
  char name[20];
  double * vPos;
  int i,err,mPos=1;
  int returnCode=0; 
  void * pscr=NULL;
  double *va_mem=(double*)malloc(sizeof(double)*(nvar_int+1));

  for(i=1;i<=nvar_int;i++) va_mem[i]=va_int[i];
  
  for(err=1; err;)
  {  
    for(;selectParam(x,y,"Change parameter",&pscr,for22,1,0,&vPos,name,&mPos);)
    {  
          strcat(name," = ");
          val=*vPos;
          if(correctDouble(x,y+4,name,&val,1)) { *vPos=val; returnCode=1;}  
    }
    if(returnCode)
    { 
      void * pscr;
      get_text(x,y+1,x+25, y+3,&pscr);
      scrcolor(FGmain,BGmain);
      goto_xy(x+5,y+1); print("Be patient:");
      goto_xy(x,y+3); print("Calculation of constraints");
      escpressed();
      err=checkParam();
      put_text(&pscr);
    } else  err=0;
    if(err && mess_y_n(10,10, "Restore previous parameter set?"))   
       for(i=1;i<=nvar_int;i++) va_int[i]=va_mem[i];
  }
  free(va_mem);
  return returnCode;
}

void show_depend(int x, int y)
{ void *pscr1=NULL;
  int i,mPos=1; 
  double*allfunc=(double*) malloc(sizeof(double)*nfunc_int);
  for(i=0;i<nfunc_int;i++) allfunc[i]=va_int[1+nvar_int+i];

  for(;;) 
  { void *pscr2=NULL;
    char name1[20];
    double * fPos;
    
    if(!selectParam(x,y+1,"Display dependence",&pscr1,0,0,1,&fPos,name1,&mPos)) return;
    for(;;)
    { char name2[20];
      double val2;
      void *pscr3=NULL;
      double xMin, xMax;
      int  nPoints=100;
      int k3;
      double *vPos;
      int mPos_=1;
      
      if(!selectParam(x,y+5,"on parameter",&pscr2,0,1,0,&vPos,name2,&mPos_))break;
   
      val2=*vPos; 
      xMin=val2 - fabs(val2)/10;
      xMax=val2 + fabs(val2)/10; 
      
      for(;;)
      {   
         char strmen[]="\030 "
            " x-Min = XXX            "
            " x-Max = YYY            "
            " Npoints = NNN          "
            " Display                ";
     
         improveStr(strmen,"XXX","%G",xMin);
         improveStr(strmen,"YYY","%G",xMax);
         improveStr(strmen,"NNN","%d",nPoints);

         
         menu1(x,y+9,"Plot",strmen,"",&pscr3,&k3);
         if(!k3) break;
         switch(k3)
         {  case 1: correctDouble(x,y+12,"xMin = ",&xMin,1); break;
            case 2: correctDouble(x,y+12,"xMax = ",&xMax,1); break;
            case 3: correctInt(x,y+12,"nPoints = ",&nPoints,1); break;
            case 4:
            if( xMax>xMin && nPoints>=3 && nPoints<=150)
            {  double dx=(xMax-xMin)/(nPoints-1);
               double f[150];
               int i, NaN=0,Esc=0;
         
               informline(0,nPoints);               
               for(i=0;i<nPoints;i++)
               {  double x=xMin+i*dx;
                  *vPos=x;
                  NaN=checkParam();
                  if(NaN) 
                  {  char mess[100];
                     sprintf(mess,"Can not evaluate constraints for %s=%G",name2, x);
                     messanykey(16,5,mess);        
                     break;
                  }
                  f[i]=*fPos;
                  Esc=informline(i,nPoints);
                  if(Esc) break;  
               }
                  
               *vPos=val2; 
               for(i=0;i<nfunc_int;i++) va_int[1+nvar_int+i]=allfunc[i];

               if(!(NaN||Esc)) plot_1(xMin,xMax,nPoints,f,NULL,"Plot",
                               name2,name1);
                               
            } else messanykey(16,5," Correct input is \n"
                                   "  xMin<xMax,\n"
                                   " 3<=nPoints<=150");
            break;
         }
       }
     }
  }
  free(allfunc);
}
