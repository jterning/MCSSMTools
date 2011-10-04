#ifndef __RDIAG__
#define __RDIAG__
extern int rDiagonal(int id, int nDim,...);
extern int rDiagonal2(int id, int nDim,...);
extern double MassArray(int id,  int i);
extern double MixMatrix(int id, int i,int j);
extern double MixMatrixU(int id, int i,int j);

#endif
