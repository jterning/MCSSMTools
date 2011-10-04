#ifndef __ALPHA_QCD__
#define __ALPHA_QCD__

extern double initQCD(double MZalphaS,double McMc,double MbP,double MtP);
extern double alphaQCD(double Q);
extern double MbRun(double Q);
extern double MbEff(double Q);
extern double MtRun(double Q);
extern double MtEff(double Q);
extern double McRun(double Q);
extern double McEff(double Q);
extern double MbPole;

#endif
