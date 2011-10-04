#ifndef __SLHA_READ__
#define __SLHA_READ__

extern int slhaRead(char *fname,int mode);
extern double slhaVal(char * Block, double Q, int nKey, ...);
extern int slhaValExists(char * Block, int nKey, ...);
extern double slhaWidth(int pNum);
extern int slhaWrite(char *fname);
extern int slhaWarnings(FILE*f);
extern int slhaDecayExists(int pNum);
extern double slhaBranch(int pNum,int N, int * nCh);
#endif
