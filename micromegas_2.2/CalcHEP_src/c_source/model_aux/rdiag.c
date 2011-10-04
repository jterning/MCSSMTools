#include<stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include"rdiag.h"

static double *vector(int nl, int nh)
{
	double *v;

	v=(double *)calloc((unsigned) (nh-nl+1), sizeof(double));
	if (!v) return NULL;
	return v-nl;
}

static double **matrix(int nrl, int nrh, int ncl, int nch)
{
	int i;
	double **m;

	m=(double **) calloc((unsigned) (nrh-nrl+1), sizeof(double*));
	if (!m) return NULL;
	m -= nrl;

	for(i=nrl;i<=nrh;i++) {
		m[i]=(double *) calloc((unsigned) (nch-ncl+1), sizeof(double));
		if (!m[i]) return NULL;
		m[i] -= ncl;
	}
	return m;
}


static void free_vector(double *v, int nl, int nh)
{
	free((char*) (v+nl));
}

static void free_matrix(double** m, int nrl, int nrh, int ncl, int nch)
{
	int i;

	for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
	free((char*) (m+nrl));
}


#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);\
	a[k][l]=h+s*(g-h*tau);

int jacobi(double** a, int n, double d[], double** v, int* nrot)
{
	int j,iq,ip,i;
	double tresh,theta,tau,t,sm,s,h,g,c,*b,*z;

	b=vector(1,n);
	z=vector(1,n);
	for (ip=1;ip<=n;ip++) {
		for (iq=1;iq<=n;iq++) v[ip][iq]=0.0;
		v[ip][ip]=1.0;
	}
	for (ip=1;ip<=n;ip++) {
		b[ip]=d[ip]=a[ip][ip];
		z[ip]=0.0;
	}
	*nrot=0;
	for (i=1;i<=50;i++) {
		sm=0.0;
		for (ip=1;ip<=n-1;ip++) {
			for (iq=ip+1;iq<=n;iq++)
				sm += (double) (fabs(a[ip][iq]));
		}
		if (sm == 0.0) {
			free_vector(z,1,n);
			free_vector(b,1,n);
			return 0;
		}
		if (i < 4)
			tresh=(double) (0.2*sm/(n*n));
		else
			tresh=0.0;
		for (ip=1;ip<=n-1;ip++) {
			for (iq=ip+1;iq<=n;iq++) {
				g=(double) (100.0*fabs(a[ip][iq]));
				if (i > 4 && (double)(fabs(d[ip])+g) == (double)fabs(d[ip])
					&& (double)(fabs(d[iq])+g) == (double)fabs(d[iq]))
					a[ip][iq]=0.0;
				else if (fabs(a[ip][iq]) > tresh) {
					h=d[iq]-d[ip];
					if ((double)(fabs(h)+g) == (double)fabs(h))
						t=(a[ip][iq])/h;
					else {
						theta=(double) (0.5*h/(a[ip][iq]));
						t=(double) (1.0/(fabs(theta)+sqrt(1.0+theta*theta)));
						if (theta < 0.0) t = -t;
					}
					c=(double) (1.0/sqrt(1+t*t));
					s=t*c;
					tau=(double) (s/(1.0+c));
					h=t*a[ip][iq];
					z[ip] -= h;
					z[iq] += h;
					d[ip] -= h;
					d[iq] += h;
					a[ip][iq]=0.0;
					for (j=1;j<=ip-1;j++) {
						ROTATE(a,j,ip,j,iq)
					}
					for (j=ip+1;j<=iq-1;j++) {
						ROTATE(a,ip,j,j,iq)
					}
					for (j=iq+1;j<=n;j++) {
						ROTATE(a,ip,j,iq,j)
					}
					for (j=1;j<=n;j++) {
						ROTATE(v,j,ip,j,iq)
					}
					++(*nrot);
				}
			}
		}
		for (ip=1;ip<=n;ip++) {
			b[ip] += z[ip];
			d[ip]=b[ip];
			z[ip]=0.0;
		}
	}
	return 1;
}

#undef ROTATE
/*==========================================================*/

extern int FError;

static int idMax=0;
static int * dim;
static double *** Rotations=NULL;
static double *** Rotations_=NULL;
static double ** MassV=NULL;

static void extendData(int id, int nDim, int std)
{ int i,j;
  if(idMax==0)
  { Rotations=(double ***) malloc((id+1)*sizeof(double**));
    Rotations_=(double ***) malloc((id+1)*sizeof(double**));
    MassV=(double **) malloc((id+1)*sizeof(double*));
    dim=(int*) malloc((id+1)*sizeof(int));
    for(i=0;i<=id;i++) 
    { Rotations[i]=NULL; Rotations_[i]=NULL; MassV[i]=NULL; dim[id]=0; }
    idMax=id;
  }else if(id>idMax)
  { Rotations=(double ***) realloc( Rotations,(id+1)*sizeof(double**));
    Rotations_=(double ***) realloc( Rotations_,(id+1)*sizeof(double**));
    MassV=(double **) realloc(MassV, (id+1)*sizeof(double*));
    dim=(int*) realloc(dim, (id+1)*sizeof(int));
    for(i=idMax+1;i<=id;i++) 
    { Rotations[i]=NULL; Rotations_[i]=NULL; MassV[i]=NULL; dim[id]=0;}
    idMax=id;
  }

  if(dim[id] && dim[id]!=nDim)
  { if(Rotations[id]){free_matrix(Rotations[id],1,nDim,1,nDim); Rotations[id]=NULL;}
    if(Rotations_[id]){free_matrix(Rotations_[id],1,nDim,1,nDim);Rotations_[id]=NULL;}
    if(MassV[id]) { free_vector( MassV[id],1,nDim);MassV[id]=NULL;}
    dim[id]=0;
  }
  
  if(dim[id]==0) 
  { Rotations[id]=matrix(1,nDim,1,nDim);
    MassV[id]=vector(1,nDim);
    dim[id]=nDim;
  }
  if(Rotations_[id]==NULL && std==2) Rotations_[id]=matrix(1,nDim,1,nDim);   
}

static  sorting(int nDim, double ** Rot, double * MassV)
{ int i,j,fl=0;
  do
  { fl=0;
    for(i=1;i<nDim;i++) if(fabs(MassV[i])>fabs(MassV[i+1]))
    {
      double tmp;
      fl=1;
      tmp=MassV[i];
      MassV[i]=MassV[i+1];
      MassV[i+1]=tmp;
      for(j=1;j<=nDim;j++)
      {
        tmp=Rot[j][i];
        Rot[j][i]=Rot[j][i+1];
        Rot[j][i+1]=tmp;
      }
    }
  }while(fl);
}

int rDiagonal(int id, int nDim,...) 
{ va_list ap;
  double**MassM=matrix(1,nDim,1,nDim);
  int i,j,nrot,fl=0;
   
  va_start(ap,nDim);
  {
    for(i=1;i<=nDim;i++)
    {  MassM[i][i]=va_arg(ap, double);
       for(j=i+1;j<=nDim;j++) MassM[i][j]=MassM[j][i]=va_arg(ap, double);
    }
  } 
  va_end(ap);

  extendData(id,nDim,1);

  FError=FError|jacobi(MassM, nDim, MassV[id], Rotations[id], &nrot); 
  free_matrix(MassM,1,nDim,1,nDim);
  
  if(!FError)sorting(nDim,Rotations[id], MassV[id]);
  
  return id;
}

int rDiagonal2(int id, int nDim,...) 
{ va_list ap;

  double**MassM=matrix(1,nDim,1,nDim), **MassM2=matrix(1,nDim,1,nDim);
  int i,j,nrot,fl=0;
   
  va_start(ap,nDim);
   for(j=1;j<=nDim;j++) for(i=1;i<=nDim;i++)  MassM[i][j]=va_arg(ap, double);
  va_end(ap);

  extendData(id,nDim,2);

  for(i=1;i<=nDim;i++) for(j=1;j<=nDim;j++)
  {
    int k;
    MassM2[i][j]=0.0;
    for(k=1;k<=nDim;k++) MassM2[i][j]+=MassM[k][i]*MassM[k][j];
  }
  FError=FError|jacobi(MassM2, nDim, MassV[id], Rotations[id], &nrot); 
  if(FError){ free_matrix(MassM,1,nDim,1,nDim); free_matrix(MassM2,1,nDim,1,nDim); return id;}
  sorting(nDim,Rotations[id], MassV[id]);
  
  for(i=1;i<=nDim;i++) for(j=1;j<=nDim;j++)
  {
    int k;
    MassM2[i][j]=0.0;
    for(k=1;k<=nDim;k++) MassM2[i][j]+=MassM[i][k]*MassM[j][k];
  }
  FError=FError|jacobi(MassM2, nDim, MassV[id], Rotations_[id], &nrot);
  if(FError){ free_matrix(MassM,1,nDim,1,nDim); free_matrix(MassM2,1,nDim,1,nDim); return id;}
    
  sorting(nDim,Rotations_[id], MassV[id]); 
  
  for(i=1;i<=nDim;i++)
  {
    int k,l;
    MassV[id][i]=0.0;
    for(k=1;k<=nDim;k++)for(l=1;l<=nDim;l++)
    MassV[id][i]+=  Rotations_[id][l][i]*MassM[l][k]*Rotations[id][k][i]  ;
  }

  for(i=1;i<=nDim;i++)
  {
    int l;
    if(MassV[id][i]<0) { MassV[id][i]*=-1; for(l=1;l<nDim;l++) Rotations_[id][l][i]*=1;} 
  }
     
  free_matrix(MassM,1,nDim,1,nDim);
  free_matrix(MassM2,1,nDim,1,nDim);
  
  return id;
}

double MassArray(int id,  int i)
{ if(i<1 ||i>dim[id]) {FError=1; return 0;}
  return  MassV[id][i];  
}

double MixMatrix(int id, int i,int j)
{ if(i<1 || i>dim[id] ||j<1 ||j>dim[id]) {FError=1; return 0;}
  return  Rotations[id][j][i];  
}

double MixMatrixU(int id, int i,int j)
{ if(i<1 || i>dim[id] ||j<1 ||j>dim[id]||!Rotations_[id]) {FError=1; return 0;}
  return  Rotations_[id][j][i];  
}

#ifdef TEST
/*=============== TESTING PROGRAM ===========================*/
int FError;

int main(void)
{ int i,j,k,l;
  double a11,a12,a21,a22;
  
  a11=0.1; a12=3.0; a21=3.0,a22=0; 

printf("=========================\ndiagonazising of symmetry matrix:\n");
printf("%E  %E\n%E  %E\n",a11,a12,a12,a22); 

  rDiagonal(1,2,a11,a12,a22);
  
  printf("masses : %E %E\n",   MassArray(1,1), MassArray(1,2) );
  printf("rotation :\n%E  %E\n%E  %E\n", MixMatrix(1,1,1), MixMatrix(1,1,2),
                                         MixMatrix(1,2,1), MixMatrix(1,2,2) );
printf("restoring of original matrix:\n");  
  for(j=1;j<=2;j++) for(i=1;i<=2;i++) 
  { double s=0;
     for(k=1;k<=2;k++) s+=MixMatrix(1,k,i)*MassArray(1,k)*MixMatrix(1,k,j);
     printf("%E",s); if(i==1)printf("  "); else printf("\n");
  }
printf("=========================\ndiagonazising of asymmetry matrix:\n");  

a11=1; a12=0.01; a21=0.05,a22=-1;
printf("%E  %E\n%E  %E\n",a11,a12,a21,a22); 

  rDiagonal2(2,2,a11,a12,a21,a22);

printf("masses : %E %E\n",   MassArray(2,1), MassArray(2,2) );
  
  printf("U rotation :\n%E  %E\n%E  %E\n", MixMatrixU(2,1,1), MixMatrixU(2,1,2),
                                         MixMatrixU(2,2,1), MixMatrixU(2,2,2) );
  
  printf("V rotation :\n%E  %E\n%E  %E\n", MixMatrix(2,1,1), MixMatrix(2,1,2),
                                           MixMatrix(2,2,1), MixMatrix(2,2,2) );

printf("restoring of original matrix:\n");  
  for(j=1;j<=2;j++) for(i=1;i<=2;i++) 
  { double s=0;
     for(k=1;k<=2;k++) s+=MixMatrixU(2,k,i)*MassArray(2,k)*MixMatrix(2,k,j);
     printf("%E",s); if(i==1)printf("  "); else printf("\n");
  }

}
#endif
