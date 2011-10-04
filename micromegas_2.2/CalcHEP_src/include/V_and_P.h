
#ifndef __Variables_and_Particles__
#define __Variables_and_Particles__

extern int nModelParticles;
 
typedef struct
{ 
  char* name; char* aname; int NPDG;char* mass; char* width; int spin2; int cdim;
}  ModelPrtclsStr;

extern ModelPrtclsStr ModelPrtcls[];
extern int nModelVars;
extern int nModelFunc;
extern char*varNames[];
extern double varValues[];
extern int calcMainFunc(void);

#endif
