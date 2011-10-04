#ifndef __PROCESS_
#define __PROCESS_

#include"syst2.h"
#include"model.h"
#include"physics.h"

#define whohowMAX 100
typedef struct whohow
{  int    who, how;
}  whohow[whohowMAX];

typedef struct hadron
{
   char        name[8];
   shortstr    contents;
   int         pow;
   int         parton[100];
   int         polarized[100];
}  hadron;

extern whohow     liminsp, LimQ;
extern whohow     limout;
extern void  nilprtcl(whohow      p_list);


extern shortstr processch, limpch, deloutch;
extern int  nin, nout, n_x;   /* Number of X-particles */
extern int  enter(void);
extern  hadron hadrons[MAXINOUT];

extern int polarized(int p, int Prtcl);

#endif
