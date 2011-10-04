#include <string.h>
#include <ctype.h>
#include <stdio.h>
#include<stdarg.h>
#include<math.h>
#include <stdlib.h>


#include"SLHAreader.h"

extern int  FError;

typedef   struct blockRec
{  struct blockRec* next;
   double val;
   int nkey;
   int  keys[12];
} blockRec;

typedef struct blockStr
{ struct blockStr*  next;
  char name[20]; 
  double scale;
  blockRec* dataList;
} blockStr;

typedef   struct decayRec
{  struct decayRec* next;
   double Br;
   int nkey;
   int  pNum[12];
} decayRec;

typedef struct decayStr
{ struct decayStr*  next;
  int pNum; 
  double pWidth;
  decayRec* dataList;
} decayStr;


static blockStr* blockList=NULL;
static decayStr* decayList=NULL;
static char*Warnings=NULL;
static int nWarnings=0;
static char creator[200],version[200];

static cleanBlockRec(blockRec*List)
{
  while(List){ blockRec*list=List; List=List->next; free(list);}  
} 
static cleanDecayRec(decayRec*List)
{
  while(List){ decayRec*list=List; List=List->next; free(list);}  
} 


static cleanData(void)
{ 
  while(blockList)
  { blockStr*block=blockList; blockList=blockList->next; 
    cleanBlockRec(block->dataList);
    free(block); 
  }
  while(decayList)
  { decayStr*decay=decayList; decayList=decayList->next; 
    cleanDecayRec(decay->dataList);
    free(decay); 
  }
  
  if(Warnings) {free(Warnings); Warnings=NULL;}
  nWarnings=0;
  creator[0]=0;
  version[0]=0;
}

static int nLine=0;

static int readLine(FILE*f,int size, char*buff)
{ int len,i;
  for(;;)
  {
    if(!fgets(buff,size-1,f))
    {  if(blockList==NULL && decayList==NULL)
       { char * mess=" No data in SLHA file\n";
         Warnings=(char*)malloc(strlen(mess)+1);
         sprintf(Warnings,mess);
         nWarnings=1;
         return -3;
       } else return -1;
    }                                                            
    nLine++;
    len=strlen(buff);
    if(len>0 && buff[len-1]=='\n')
    {
      buff[len-1]=0;
      if(len>1 && buff[len-2]==13) buff[len-2]=' ';
    }
    else
    { char ch;
      if(fscanf(f,"%c",ch)==1) 
      {  
         printf("Line %d is too long\n",nLine);  
         cleanData();
         return nLine;
      }else if(len>0 && buff[len-1]==13) buff[len-1]=' ';    
    }
    for(i=0;buff[i]&&buff[i]==' ';i++);
    if(buff[i]!=0 && buff[i]!='#') 
    { if(i)
      { int k;
        for(k=i;buff[k];k++) buff[k-i]=buff[k];
        buff[k-i]=0; 
      }
      return 0;
    }
  }
}

static int readNumber(char *s, int * iN, double *dN)
{ char buff[200];
  if(sscanf(s,"%d%s",iN, buff)==1) return 0;
  if(sscanf(s,"%lf%s",dN,buff)==1) return 1;
  return -1;
}

/*
  mode= 1*m1+2*m1+4*m4

  m1  0 overwrite all;  1 keep old data
  
  m2  0 read DECAY      1: don't read   Decay 
  m4  0 read BLOCK      1: don't read   Blocks

  2 - clean and read only Blocks
  1 - keep old results  and read both BLOCKs and DECAYs  
*/

int slhaRead(char *fname, int mode)
{
  FILE *f;
  char buff[200],name[200],rest[200];
  int n,err,m1,m2,m4;
  double scale;
  
  m1=mode&1;
  m2=mode&2;
  m4=mode&4;
  
  nLine=0;
  FError=0;
  
  if(m1==0)cleanData();
    
  f=fopen(fname,"r");  if(f==NULL) return -1;
  err=readLine(f,200,buff);
  if(err){fclose(f); if(err==-1) return 0; else return err;}
  
  for(;;) 
  { char *c, *block="BLOCK ";
    int L;
    c=strchr(buff,'#'); if(c) c[0]=0;
      
    for(L=0;L<6 && buff[L]&& toupper(buff[L])==block[L] ;L++) continue;
    if(L==6)
    { int i; 
      if( sscanf(buff+6,"%s %s",name , rest)==1) scale=-1;  
      else if(sscanf(buff+6,"%s %*s %lf %s",name , &scale, rest)!=2)
      {  cleanData();
         printf("Unexpected BLOCK specification at line %d\n",nLine);
         fclose(f);
         return nLine;
      }     
      if(strlen(name)>19)
      {  cleanData();
         printf("Too long name of BLOCK  at line %d\n",nLine);
         fclose(f);
         return nLine;
      } 
      for(i=0;name[i];i++) name[i]=toupper(name[i]);
      if(strcmp(name,"SPINFO")==0 ||strcmp(name,"DCINFO")==0   )
      { 
       for(;;)
       { err=readLine(f,200,buff);
         if(err){fclose(f); if(err==-1) return 0; else return err;}
         if(sscanf(buff,"%d",&n)!=1)  break; 
         if(n==1) sscanf(buff,"%*d %[^\n]",creator);
         else if(n==2) sscanf(buff,"%*d %[^\n]",version);
         else if(n==4)
         {  cleanData();  
            Warnings=(char*)malloc(strlen(buff)+2);
            sscanf(buff,"%*d %[^\n]", Warnings);
            nWarnings=1; 
            fclose(f);                       
            return -2;
         }
         else if(n==3)
         {
           if(Warnings==NULL) 
           { Warnings=(char*)malloc(strlen(buff)+2);
             Warnings[0]=0;
           } else   
           Warnings=(char*)realloc(Warnings,strlen(Warnings)+strlen(buff)+4);
           sscanf(buff,"%*d %[^\n]", Warnings+strlen(Warnings)); 
           strcat(Warnings,"\n");
           nWarnings++;
         }     
       }
      }else  
      { blockStr*newBlock; 
        if(m4) for(;;)
        {
          err=readLine(f,200,buff);
          if(err){fclose(f); if(err==-1) return 0; else return err;}
          if(!isdigit(buff[0])) break;
        }
        newBlock=(blockStr*)malloc(sizeof(blockStr));
        newBlock->next=blockList;
        strcpy(newBlock->name,name);         
        newBlock->dataList=NULL;
        newBlock->scale=scale;
        blockList=newBlock;

        for(;;)
        { int err,last,k;
          blockRec*dr;

          err=readLine(f,200,buff);
          if(err){fclose(f); if(err==-1) return 0; else return err;}

          if(isalpha(buff[0])) break;
          c=strchr(buff,'#'); if(c) c[0]=0;
          
          dr=(blockRec*)malloc(sizeof(blockRec));
          dr->next=blockList->dataList;
          dr->nkey=0;
          blockList->dataList=dr;

          for(last=0,k=0; ;k++)
          {  int tp,iVal;
             double dVal;
             
             if(k) c=strtok(NULL," "); else c=strtok(buff," ");
             if(!c) {if(!last) dr->nkey--;     break;}

             tp=readNumber(c,&iVal,&dVal);
             if(last || isalpha(c[0]) || tp==-1 ) 
             { cleanData();
               printf("Unexpected token %d at line %d\n",k+1,nLine);
               fclose(f);
               return nLine;
             }
             
             if(k>10) 
             {
               printf("Too long key sequence at line %d\n",nLine);
               fclose(f);
               return nLine;              
             }
              
             if(tp) 
             { last=1;
               dr->val=dVal;
             } else 
             {  
               dr->keys[dr->nkey]=iVal;
               dr->nkey++;
               dr->val=iVal;
             }      
          }  

        }  
      }
    }else 
    { char  *decay="DECAY ";
      int pNum;
      double pWidth;
      decayStr*newDecay;
        
      for(L=0;L<6 && buff[L]&& toupper(buff[L])==decay[L] ;L++) continue;
      if(L!=6) 
      { 
        printf("Unexpected first word  in  line %d\n",nLine);
        fclose(f);  return nLine;
      }
      
      c=strchr(buff,'#'); if(c) c[0]=0;
      if( sscanf(buff+6,"%d %lf %s",&pNum ,&pWidth, rest)!=2)
      {  printf("Unexpected DECAY specification at line %d\n",nLine);
         printf("buff=%s\n",buff);
         fclose(f);
         return nLine;
      }
      if(m2)
      { 
        for(;;)
        { double x;
          err=readLine(f,200,buff);
          if(err){fclose(f); if(err==-1) return 0; else return err;}
          if(sscanf(buff,"%lf",&x)!=1) break;
        }  
      } else 
      {
        newDecay=(decayStr*)malloc(sizeof(decayStr));
        newDecay->next=decayList;
        newDecay->pNum=pNum;         
        newDecay->dataList=NULL;
        newDecay->pWidth=pWidth;
        decayList=newDecay;
        for(;;)
        { int err,k;
          decayRec*dr;

          err=readLine(f,200,buff);
          if(err){fclose(f); if(err==-1) return 0; else return err;}
          if(isalpha(buff[0])) break;
          c=strchr(buff,'#'); if(c) c[0]=0;
          
          dr=(decayRec*)malloc(sizeof(decayRec));
          dr->next=newDecay->dataList;
          dr->nkey=0;
          newDecay->dataList=dr;
        
          err=sscanf(buff,"%lf %d %d %d %d %d %d %d %d %d",
          &dr->Br, &dr->nkey, dr->pNum, dr->pNum+1, dr->pNum+2,dr->pNum+3,
          dr->pNum+4, dr->pNum+5, dr->pNum+6,dr->pNum+7);
          if(err<2 || err!=dr->nkey+2)
          {  cleanData();
           printf("Wrong decay record  at line %d\n",nLine);
           fclose(f); 
           return nLine;
          } 
        }
      }    
    }     
  }
  return 0;
}         

static double* slhaValAddress(char * Block, int nKey, int *keys)
{
  char BLOCK[40];
  blockStr* blck=blockList;
  blockRec * dr;
  int i;
                            
  for(i=0; Block[i];i++)  BLOCK[i]=toupper(Block[i]);
  BLOCK[i]=0;
  
  while(blck && strcmp(BLOCK,blck->name)) blck=blck->next;
  if(!blck) return NULL;
  dr=blck->dataList;
  for(;dr;dr=dr->next)
  { if(dr->nkey!=nKey) continue;
    for(i=0;i<nKey;i++) if(keys[i]!=dr->keys[i]) break;
    if(i==nKey) return &dr->val; else continue;
  }
  return NULL;  
}

double slhaVal(char * Block, double Q, int nKey, ...)
{ 
  va_list ap; 
  int keys[12];
  int i;
  char BLOCK[40];
  blockStr* blck=blockList;
  blockRec * dr;

  double val[3], scale[3];
  int found[3]={0,0,0};

  if(strlen(Block)>19) {FError=1; return 0;}
      
  va_start(ap,nKey);
  for(i=0;i<nKey;i++)keys[i]=va_arg(ap, int);
  va_end(ap);

  for(i=0; Block[i];i++)  BLOCK[i]=toupper(Block[i]);
  BLOCK[i]=0;

  for(blck=blockList;  blck; blck=blck->next)if(strcmp(BLOCK,blck->name)==0)
  {  
     int pos;
     if(blck->scale < 0)    { if(found[0]) continue; else pos=0;}
     else if(blck->scale<Q) { if(found[1] && scale[1]>blck->scale) continue; else pos=1;}
     else                   { if(found[2] && scale[2]<blck->scale) continue; else pos=2;}

     dr=blck->dataList;
     for(;dr;dr=dr->next)
     { if(dr->nkey!=nKey) continue;
       for(i=0;i<nKey;i++) if(keys[i]!=dr->keys[i]) break;
       if(i==nKey) 
       {
          found[pos]=1;
          scale[pos]=blck->scale;
          val[pos]=dr->val; 
          break;
       }
     }  
  }
  
  if(found[0]==0 && found[1]==0 && found[2]==0) 
  { printf(" Block '%s', key={",BLOCK);
    for(i=0;i<nKey;i++) printf(" %d",keys[i]);
    printf("} - is absent\n"); 
    FError=1;
    return 0;
  }
  if(found[1]==0 && found[2]==0) return val[0];
  if(found[1]==0) return val[2];
  if(found[2]==0) return val[1];
  return  (val[1]*log(scale[2]/Q)+val[2]*log(Q/scale[1]))/log(scale[2]/scale[1]);  
}


int slhaValExists(char * Block, int nKey, ...)
{ 
  va_list ap; 
  int keys[12];
  double * address;
  int i;

  if(nKey>10) return 0;  
  va_start(ap,nKey);
  for(i=0;i<nKey;i++)keys[i]=va_arg(ap, int);
  va_end(ap);
  
  address=slhaValAddress(Block,nKey, keys);
  if(address) return 1; else return 0;
}

int slhaWarnings(FILE*f)
{
  if(f&&Warnings) fprintf(f,Warnings);
  return nWarnings;
}

int slhaWrite(char *fname)
{ 
  blockStr* block;
  blockRec* rec;
  decayStr* decay;
  decayRec* dr;

  FILE*f=fopen(fname,"w");
  if(!f) return 1;
  if(!blockList  && ! decayList) return 2;
  
  if(blockList)
  {
    fprintf(f,"BLOCK SPINFO # General Information\n");
    fprintf(f," 1 %s\n",creator);
    fprintf(f," 2 %s\n",version);
  } else 
  {
    fprintf(f,"BLOCK DCINFO # General Information\n");
    fprintf(f," 1 %s\n",creator);
    fprintf(f," 2 %s\n",version);
  } 
  
  
  
  if(nWarnings) 
  { char buff[100]; 
    char *  c;
    for(c=Warnings;;)
    {  
       sscanf(c,"%[^\n]",buff); 
       fprintf(f," 3 %s \n",buff);
       c=strchr(c,'\n');
       if(!c || c[1]==0) break;
       c++;
    } 
  } 
  
  for(block=blockList;block;block=block->next)
  {
     fprintf(f,"BLOCK %s", block->name);
     if(block->scale>0) fprintf(f," Q= %E", block->scale); 
     fprintf(f," # nc\n");
     for(rec=block->dataList;rec;rec=rec->next)
     {
       int i;
       for(i=0;i<rec->nkey; i++) fprintf(f," %d", rec->keys[i]);
       fprintf(f,"   %E  # nc\n",rec->val);  
     }            
  }

  for(decay=decayList;decay;decay=decay->next)
  {
     fprintf(f," DECAY %d %E # nc " , decay->pNum,decay->pWidth); 

     for(dr=decay->dataList;dr;dr=dr->next)
     {
       int i;
       fprintf(f," %E  %d ", dr->Br,dr->nkey);
       for(i=0;i<dr->nkey; i++) fprintf(f," %d", dr->pNum[i]);
       fprintf(f," # nc\n");  
     }            
  }

  fclose(f);
  return 0;
}

int slhaDecayExists(int pNum)
{ 
   decayStr* decay=decayList;
   for(;decay;decay=decay->next) if( decay->pNum==pNum)
   { decayRec*dr=decay->dataList;
     int n;
     for(n=0; dr; dr=dr->next, n++) continue;
     return n;
   }
   return -1;  
}

double slhaWidth(int pNum)
{  
   decayStr* decay=decayList;
   for(;decay;decay=decay->next) if( decay->pNum==pNum) return decay->pWidth;
   printf("Error: width for particle %d is unknown\n",pNum);
   FError=1;
   return 0;
}

double slhaBranch(int pNum,int N, int * nCh)
{
   decayStr* decay=decayList;
   for(;decay;decay=decay->next) if( decay->pNum==pNum)
   { decayRec*dr=decay->dataList;
     int i;
     for(N--;N,dr;N--,dr=dr->next) continue;
     if(dr==NULL) { FError=1;return 0;}
     for(i=0;i<dr->nkey;i++) nCh[i]=dr->pNum[i];
     nCh[i]=0;
     return dr->Br;
   }
   FError=1;
   return 0;
}

/* ====================== FORTRAN VERSIONS =================== */

#ifdef PLUS_FORT
struct 
{ int error;} ferror_;


static void fName2c(char*f_name,char*c_name,int len)
{ int i; for(i=len-1;i>=0 &&f_name[i]==' ';i--);
  c_name[i+1]=0;
  for(;i>=0;i--) c_name[i]=f_name[i];
}


int slharead_(char * fname, int * mode, int len)
{ char  c_name[200];
  fName2c(fname,c_name,len);
  return slhaRead(c_name, *mode);
}

int slhawrite_(char *fname, int len)
{ char  c_name[200];
  fName2c(fname,c_name,len);
  slhaWrite(c_name);
}

int slhadecayexists_(int *pNum){ return slhaDecayExists(*pNum);}

double slhawidth_(int *pNum){ return slhaWidth(*pNum);}
double slhabranch_(int*pNum,int*N,int*nCh){return slhaBranch(*pNum,*N,nCh);}

double slhaval0_(char * Block, double *Q, int len)
{ char  c_name[200];
  fName2c(Block,c_name,len);
  return slhaVal(c_name, *Q, 0);
} 

double slhaval1_(char * Block, double *Q, int *k1, int len)
{  char  c_name[200];
   fName2c(Block,c_name,len);
   return slhaVal(c_name, *Q, 1,*k1);
} 

double slhaval2_(char * Block, double *Q, int *k1,int*k2, int len)
{  char  c_name[200];
   fName2c(Block,c_name,len);
   return slhaVal(c_name, *Q, 2,*k1,*k2);
} 

double slhaval3_(char * Block, double *Q, int *k1,int*k2,int*k3, int len)
{  char  c_name[200];
   fName2c(Block,c_name,len);
   return slhaVal(c_name, *Q, 3,*k1,*k2,*k3);
} 
 

int slhavalexists0_(char * Block, int len)
{  char  c_name[200];
   fName2c(Block,c_name,len);
   return slhaValExists(Block, 0);
}

int slhavalexists1_(char * Block, int*k1, int len)
{  char  c_name[200];
   fName2c(Block,c_name,len);
   return slhaValExists(Block, 1,*k1);
}

int slhavalexists2_(char * Block, int*k1, int*k2, int len)
{  char  c_name[200];
   fName2c(Block,c_name,len);
   return slhaValExists(Block, 2,*k1,*k2);
}

int slhavalexists3_(char * Block, int*k1, int*k2, int*k3, int len)
{  char  c_name[200];
   fName2c(Block,c_name,len);
   return slhaValExists(Block, 3,*k1,*k2,*k3);
}

#endif
/*int slhaWarnings(FILE*f)*/


#ifdef TEST

int  FError=0;

int main(int nargs, char** vargs)
{ int nv,err;
/*
  if(nargs!=2) return 0;
  err=slhaRead(vargs[1],0);
*/  
  err=slhaRead("decay.dat",0);
  printf("err=%d\n",err);
  
  printf("%E \n", slhaVal("STAuMIX",0., 2,2,2));
  
 printf("0    %E \n", slhaVal("HMIX",0.,   1,2));
 printf("400  %E \n", slhaVal("HMIX",400., 1,2));
 printf("2000 %E \n", slhaVal("HMIX",2000.,1,2));
 printf("4000 %E \n", slhaVal("HMIX",4000.,1,2));
 printf("8000 %E \n", slhaVal("HMIX",8000.,1,2));  

printf("AU 8000 %E \n", slhaVal("AU",8000.,2,3,3));

printf("AU    0 %E \n", slhaVal("AU",0.,2,3,3));
printf("AU    200 %E \n", slhaVal("AU",200.,2,3,3));

  nv = slhaWarnings(NULL);
  
  if(nv)
  { 
    printf("There are %d warnings\n",nv);
    slhaWarnings(stdout);
  }  
  
  slhaWrite("qq");
printf("write qq ok\n");  
  return 0;
}
#endif
