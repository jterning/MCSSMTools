
#include"micromegas_aux.h"

ModelPrtclsStr *OddPrtcls=NULL;
int Nodd=0;

int createTableOddPrtcls(void)
{
   int i;
   if(OddPrtcls) return 0;

   for(i=0,Nodd=0;i<nModelParticles;i++)
   { int o1,o2;
     o1=(ModelPrtcls[i].name[0]=='~');
     o2=(ModelPrtcls[i].aname[0]=='~');
     if(o1&&o2) Nodd++;
     if( (o1&&!o2) || (o2&&!o1) ) return 1;
   }    

   OddPrtcls=( ModelPrtclsStr*)malloc(Nodd*sizeof(ModelPrtclsStr));
   Nodd=0;
   for(i=0,Nodd=0;i<nModelParticles;i++) if(ModelPrtcls[i].name[0]=='~')
         OddPrtcls[Nodd++]=ModelPrtcls[i];
   return 0;
}

int pTabPos(char * name)
{
  int i;
  for(i=0;i<nModelParticles;i++)
  { 
     if(!strcmp(name,ModelPrtcls[i].name )) return   i+1;
     if(!strcmp(name,ModelPrtcls[i].aname)) return -(i+1);
  }
  return 0;
}


double pMass(char * name)
{
  char *nm;
  int n=pTabPos(name);
  if(!n){printf("Wrong particle name %s\n",name); return 0;}
  nm=ModelPrtcls[abs(n)-1].mass;
  if(nm[0]=='0') return 0; else return fabs(findValW(nm));
}

long pNum(char * name)
{
  int n=pTabPos(name);
  if(!n){printf("Wrong particle name %s\n",name); return 0;}
  if(n>0)  return  ModelPrtcls[abs(n)-1].NPDG;
  else     return -ModelPrtcls[abs(n)-1].NPDG;
}
