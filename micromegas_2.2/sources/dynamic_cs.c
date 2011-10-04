#ifdef __hpux
#include<dl.h>
#else
#include <dlfcn.h>
#endif

#include <unistd.h>
#include <ctype.h>
#include <string.h>

#include <sys/wait.h>

#include <sys/types.h>
#include <sys/utsname.h>


#include <sys/stat.h>
#include <fcntl.h>
#include <sys/file.h>


#include "micromegas.h"
#include "micromegas_aux.h"
#include "micromegas_f.h"

int checkLockFile(int * delay)
{ 
  int fd,ret;
  static char *lockTxt=NULL; 
  char*fname;
  
#if defined (FLOCK) || defined (LOCKF)|| defined (FCNTL)
  if(lockTxt==NULL)
  {
     struct utsname buff;
     uname(&buff);
     lockTxt=(char*)malloc(100);
     sprintf(lockTxt," process %d on %s\n",getpid(),buff.nodename);
  }

  fname=malloc(strlen(WORK)+20);
  sprintf(fname,"%s/lock_",WORK);
  fd=open(fname, O_WRONLY|O_CREAT,0666);
  free(fname);
  if(fd<0) return -1;
#else 
  return -1;
#endif

#ifdef FCNTL
{ 
  struct flock myLock;
  myLock.l_type= F_WRLCK;       /*  F_RDLCK ||  F_WRLCK || F_UNLCK */
  myLock.l_whence=SEEK_SET;
  myLock.l_start=0;
  myLock.l_len=10; 
  ret=fcntl(fd, F_SETLK, &myLock);
}
#endif

#ifdef FLOCK  
  ret=flock(fd,LOCK_EX|LOCK_NB);   
#endif

#ifdef LOCKF
  ret = lockf(fd, F_TLOCK, 10);
#endif  
  if(ret)
  { printf("Library generation is temporary locked by other program\n"
           "see file %s/lock_\n",WORK);
    *delay=1;
#ifdef FCNTL
{ 
  struct flock myLock;
  myLock.l_type= F_WRLCK;       /*  F_RDLCK ||  F_WRLCK || F_UNLCK */
  myLock.l_whence=SEEK_SET;
  myLock.l_start=0;
  myLock.l_len=10; 
  ret=fcntl(fd, F_SETLKW, &myLock);
}
#endif

#ifdef FLOCK
    flock(fd,LOCK_EX);
#endif

#ifdef LOCKF
    lockf(fd, F_LOCK, 10);
#endif    
  } else *delay=0;
  write(fd,lockTxt,strlen(lockTxt));
  return fd;
}

void removeLock(int fd)
{
if(fd<=0) return;

#if defined (FLOCK) || defined (LOCKF)||defined (FCNTL)
#ifdef FCNTL
{ 
  struct flock myLock;
  myLock.l_type= F_UNLCK;       /*  F_RDLCK ||  F_WRLCK || F_UNLCK */
  myLock.l_whence=SEEK_SET;
  myLock.l_start=0;
  myLock.l_len=10; 
  fcntl(fd, F_SETLK, &myLock);
}
#endif

#ifdef FLOCK
   flock(fd,LOCK_UN);   
#endif

#ifdef LOCKF
   lockf(fd, F_ULOCK, 10);
#endif
   close(fd);
#endif
}


static void* newSymbol(void*handle,char *name)
{
#ifdef __hpux
void * addr;
     if(shl_findsym((shl_t*)&handle,name,TYPE_UNDEFINED,&addr)) return NULL;
       else return addr;
#else
      return dlsym(handle, name);
#endif
}

static void * dLoad(char * libName)
{
void *q;
if(access(libName,R_OK)) return NULL;

#ifdef __hpux
   return  shl_load(libName,0,0L);
#else
   q= dlopen(libName, RTLD_NOW);
   if(!q) printf("%s\n",dlerror()); 
   return q;
#endif
}


static void dClose(void * handle)
{
#ifdef __hpux
       shl_unload(handle);
#else
       dlclose(handle);
#endif
}

static numout* loadLib(void* handle, char * lib)
{ numout * cc=malloc(sizeof(numout));
  char name[100];
  if(!handle) {free(cc); return NULL;}   
  cc->handle=handle;
  sprintf(name,"interface_%s",lib);
  cc->interface=newSymbol(handle, name);
  if(!cc->interface || cc->interface->nprc==0){free(cc); return NULL;}
  else
  {  int i;
     cc->init=0;
     cc->Q=NULL, cc->SC=NULL, cc->GG=NULL;
     cc->link=malloc(sizeof(double*)*(1+cc->interface->nvar));
     cc->link[0]=NULL;
     for(i=1;i<=cc->interface->nvar;i++) 
     { char *name=cc->interface->varName[i];
       cc->link[i]=varAddress(name);
       if(strcmp(name,"Q")==0) cc->Q=cc->interface->va+i; 
       else if(strcmp(name,"SC")==0) cc->SC=cc->interface->va+i;
       else if(strcmp(name,"GG")==0) cc->GG=cc->interface->va+i;
     }  
  }
  return cc;
}

void pname2lib(char*pname, char * libname)
{
  int n;
  n=pTabPos(pname);
  if(!n) {printf("Wrong particle name %s\n",pname); libname[0]=0; return;}
  if(n>0) sprintf(libname,"p%d",n); else sprintf(libname,"a%d",-n);
}


typedef struct  procRec 
{ struct procRec  * next;
  char * libname;
  numout * cc;
}  procRec;   

static  procRec* allProc=NULL;

numout*newProcess_(int twidth,int model, char*Process, char*excludeVirtual, char*excludeOut,char*lib,int usr)
{
   char *proclibf,*command;
   void * handle=NULL;
   int new=0;
   numout * cc;
   procRec*test;
   int Len;

   for(test=allProc;test; test=test->next)
   { if(strcmp(lib,test->libname)==0) return test->cc;}

   Len=strlen(WORK)+strlen(lib)+200;
   proclibf=malloc(Len);

   if(Process) Len+=strlen(Process);
   if(excludeVirtual) Len+=strlen(excludeVirtual);
   if(excludeOut) Len+=strlen(excludeOut);
   command=malloc(Len);
 
   sprintf(command,"cd %s;" MAKE " -s -f make_check  lib=%s.so model=%d",WORK,lib,model); 
   system(command);
   
   sprintf(proclibf,"%s/so_generated/%s.so",WORK,lib);
   handle=dLoad(proclibf);
   if(!handle)
   {  int delay,fd,i;
   
      for(i=0;Process[i]==' ';i++); if(Process[i]==0)
      {    
        free(command); free(proclibf); 
        return NULL;    
      }      
  
   
      fd=checkLockFile(&delay);       
      if(delay)  handle=dLoad(proclibf);
      if(!handle)
      {
        char options[20];
        char modelCh[50];
        int ret;  
        if(usr)
        {  char * reserved[3]={"omg", "pwidth", "dir"};
           for(i=0;i<3;i++)
           {
             char * c=strstr(lib,reserved[i]);
             if(c==lib)
             {  printf("Library names started from '%s' are reserved for internal use\n",
                           reserved[i]);
                if(i==0) printf("The libraries of WIMP annihilation  become available automatically\n"
                              "after calculation of relic density.\n");              
                free(command); free(proclibf); removeLock(fd);
                
                return NULL;
              }
           }
        }
      
        if(twidth) strcpy(options,"5[[{[{}");else strcpy(options,"");
        for(modelCh[0]=0; model>1;model--) strcat(modelCh,"[");
         
        sprintf(command,"cd %s; ./newProcess %s \"%s\" \"%s\" \"%s\"",
                            WORK, lib,options,modelCh,Process);
  
        if(excludeVirtual) sprintf(command+strlen(command)," \"%s\"",excludeVirtual);
        else  sprintf(command+strlen(command)," \"\"");            
        if(excludeOut) sprintf(command+strlen(command)," \"%s\"",excludeOut);       
        ret=system(command);
    
  
        if(ret<0 || WIFSIGNALED(ret)>0 ) exit(10);
        
        
        if(ret==0) handle=dLoad(proclibf); else 
        { printf(" Can not compile %s \n", Process);
          free(command); free(proclibf);
          removeLock(fd);
          return NULL;
        } 
        if(!handle)
        { printf(" Can not load the compiled library %s \n",proclibf);
           free(command); free(proclibf);
           removeLock(fd);
          return NULL;
        }         
        new=1;   
      }
      removeLock(fd);
      
   }
   cc=loadLib(handle,lib);
   if(!cc && new) dClose(handle);
   if(cc)
   {  test=(procRec*)malloc(sizeof(procRec));
      test->next=allProc; allProc=test;
      test->libname=(char*) malloc(strlen(lib)+1);
      strcpy(test->libname,lib);
      test->cc=cc;
   } else if(new) dClose(handle);  
    free(command); free(proclibf);
    return cc; 
}


int  procInfo1(numout*cc, int *nsub, int * nin, int *nout)
{
  if(nin) *nin=cc->interface->nin;
  if(nout)*nout=cc->interface->nout;
  if(nsub)*nsub=cc->interface->nprc;
  return 0;
}

int procInfo2(numout*cc,int nsub,char**name,double*mass)
{
  int i;
  int ntot=cc->interface->nin+cc->interface->nout;
    
  if(nsub<1 || nsub> cc->interface->nprc) return 2;

  if(name)for(i=0;i<ntot ;i++) 
  name[i]=(cc->interface->pinf)(nsub,i+1,NULL,NULL);

  if(mass)
  {  
    for(i=1;i<=cc->interface->nvar;i++) if(cc->link[i]) cc->interface->va[i]=*(cc->link[i]);  
    if(cc->Q) *(cc->Q)=100;
    if(cc->interface->calcFunc()>0) { printf("cannot calculate constr\n"); return 4;}
    if(cc->Q)
    {
      cc->interface->pinf(nsub,1,mass,NULL);
      *(cc->Q)=mass[0];
      if(cc->SC && cc->GG) *(cc->GG)=*(cc->SC);
      cc->interface->pinf(nsub,1,mass,NULL);
      *(cc->Q)=mass[0];
      if(cc->SC && cc->GG) *(cc->GG)=*(cc->SC);
    }      
       
    for(i=0;i<ntot ;i++) cc->interface->pinf(nsub,i+1,mass+i,NULL);     
  }
  return 0;
}

numout*newProcess(char*Process,char*lib){return newProcess_(0,1,Process,NULL,"",lib,1);}


void newprocess_(char*Process,char*lib, int * address, int len1,int len2)
{ char cProcess[100], clib[100];
  numout*cc;
  fName2c(Process,cProcess,len1);
  fName2c(lib,clib,len2);
  cc=newProcess(cProcess,clib);
  if(!cc) *address=0;else memcpy(address,&cc,sizeof(cc)); 
}
  
     
void  procinfo1_(int*ccf, int *ntot, int * nin, int *nout)
{  numout*cc;
   memcpy(&cc,ccf,sizeof(cc));
   procInfo1(cc, ntot, nin, nout);}

void procinfo2_(int*ccf,int*nsub,char*name,double*mass,int len)
{ int ntot, nin, nout,i;
  numout*cc;
  char ** cname;  
  memcpy(&cc,ccf,sizeof(cc));
  procInfo1(cc, &ntot, &nin, &nout);
  cname=malloc((nin+nout)*sizeof(char*));
  
  procInfo2(cc,*nsub,cname, mass);
  for(i=0;i<nin+nout;i++) cName2f(cname[i],name+i*len,len);
}

txtList  makeDecayList(char * pname, int nx)
{ 
  char command[200],fnameG[200],fnameL[50],lname[20],buff[200];
  FILE *f;
  txtList List=NULL;

  pname2lib(pname,lname);
  sprintf(fnameL,"dList_%s_%dx",lname,nx);
  sprintf(fnameG,"%s/so_generated/%s",WORK,fnameL);

   sprintf(command,"cd %s;" MAKE " -s -f make_check  lib=%s model=1",WORK,fnameL); 
   system(command);
  
    
  if(access(fnameG, R_OK))
  {
     sprintf(command,"cd %s; rm -f tmp/* results/*;"
                     " ../../CalcHEP_src/bin/s_calchep -blind \"{{%s->%d*x{{{[[{0\";"
                     "mv results/list_prc.txt so_generated/%s",WORK,pname,nx, fnameL);
     system(command);
  }
  f=fopen(fnameG,"r");
  if(!f) return NULL;
  for(; 1==fscanf(f,"%[^\n]\n",buff); )
  { txtList l=malloc(sizeof(txtListStr));
    l->next=List; List=l;
    l->txt=malloc(strlen(buff)+1);
    strcpy(l->txt,buff);
  }
  fclose(f);
  return List;
}

int decodeProcess(char *txt,int*inList,int*outList)
{ char name[20];
   char *ch_,*ch;
   int i,p;   
   ch_=strstr(txt,"->");
   if(!ch_) { inList[0]=0; ch=txt;}
   else
   { 
     for(p=0,ch=txt;; )
     { sscanf(ch," %[^,]",name);
       ch_=strstr(name,"->");
       if(ch_) *ch_=0;
       for(i=strlen(name)-1; i>=0 && name[i]==' '; i--) name[i]=0;
       inList[p]=pTabPos(name);
       if(!inList[p]) return -(p+1);
       p++;
       if(ch_) break;     
       ch=strchr(ch,',');
       if(!ch) break; else ch++;
     }  
     inList[p]=0; 
     ch=strstr(txt,"->")+2; 
   }  

   for(p=0;ch; )
   { sscanf(ch," %[^,]",name);
     for(i=strlen(name)-1;i>=0 && name[i]==' '; i--) name[i]=0;
     outList[p]=pTabPos(name);
     if(!outList[p]) return p+1;
     p++;    
     ch=strchr(ch,',');
     if(!ch) break;
     ch++;
   }
   outList[p]=0;
   return 0;
}


 
void massFilter(double M, txtList * List)
{
  txtList lold=*List, lnew=NULL,lnext;

  while(lold)
  {  double Msum=0;
     char *ch;
     ch=strstr(lold->txt,"->")+2;
     for( ; ch; ch=strchr(ch,','))
     { char buff[10];
       ch++;
       sscanf(ch,"%[^,]",buff);
       Msum+=fabs(pMass(buff));
     } 
     lnext=lold->next;
     if(M>Msum) {lold->next=lnew; lnew=lold;}
     else {free(lold->txt); free(lold);}
     lold=lnext;
 } 
 *List=lnew;
}

void gammaGluFilter(txtList * List)
{
  txtList lold=*List, lnew=NULL,lnext;

  while(lold)
  {  int del=0,code;
     char *ch;
     ch=strstr(lold->txt,"->")+2;
     for( ; !del && ch; ch=strchr(ch,','))
     { char buff[10];
       ch++;
       sscanf(ch,"%[^,]",buff); 
       code=pNum(buff);
       if(code==22 || code ==21) { del=1;}
     } 
     lnext=lold->next;
     if(del) {free(lold->txt); free(lold);}
     else    {lold->next=lnew; lnew=lold;}
  
     lold=lnext;
 } 
 *List=lnew;
}


void process2Lib(char * process,char * lib)
{ 
  char * ch, *ch_;
  char buff[20],bufflib[20];
  int in,out;

  ch= strstr(process,"->");
  ch_=strchr(process,',');
  if(ch_>ch) in=1;else in=2;

  for(out=1,ch_=strchr(ch,','); ch_; ch_=strchr(ch_+1,',')) out++;

  if(in==1) 
  { sprintf(lib,"d_");  
    sscanf(process,"%s",buff);
    pname2lib(buff,bufflib);
    strcat(lib,bufflib);
  }
  else
  { sprintf(lib,"c_");
    sscanf(process,"%[^,]",buff);
    pname2lib(buff,bufflib);
    strcat(lib,bufflib);
    ch=strchr(process,',');
    sscanf(ch+1,"%s",buff);
    pname2lib(buff,bufflib);
    strcat(lib,bufflib);
  }
  ch=strstr(process,"->");
  ch+=2;
  for( ;  ch; ch=strchr(ch,','))
  {
     ch++;
     sscanf(ch,"%[^,]",buff);
     pname2lib(buff,bufflib);
     strcat(lib,bufflib);
  } 
}

void cleanTxtList(txtList L)
{ 
  while(L) {txtList l=L; free(L->txt); l=L;  L=L->next; free(l);}  
}

void printTxtList(txtList L, FILE *f)
{ for(;L;L=L->next) fprintf(f,"%s\n",L->txt);}

