/*
Copyright (C) 1997, Alexander Pukhov
*/
#include"s_files.h"
#include"out_serv.h"
#include"denominators.h"
#include"saveres.h"
#include"prepdiag.h"
#include"pvars.h"
#include"process.h"
#include"c_out.h"


static int  nincount(int v,int l)
{int  i, vv, ll, summ;

   if (IN_PRTCL  & vcs.vertlist[v-1][l-1].prop) return 1;
   if (OUT_PRTCL & vcs.vertlist[v-1][l-1].prop) return 0;
   summ = 0;
   vv = vcs.vertlist[v-1][l-1].nextvert.vno;
   ll = vcs.vertlist[v-1][l-1].nextvert.edno;
   for (i = 1; i <= vcs.valence[vv-1]; i++)
      if (i != ll)
         summ += nincount(vv,i);
   return summ;
}



int  ttypepropag(int v,int l)
{
   if (nin == 1) return 0;
   return (nincount(v,l) == 1);
}
      
      
static int stype(char * momStr)
{ 
  int c,i; 
  if(nin==1) return 1;
  for(i=0,c=0; momStr[i];i++) if (momStr[i]<=2) c++; 
  if(c&1) return 0; else return 1;
}       


/* momdep from prepdiag.h must be calculated before */
void  calcdenominators(vcsect vcs )
{ int  v, l, k;
  char buff[MAXINOUT+1];



   denrno = 0;
   for (v = 1; v <= vcs.sizet; v++)
   for (l = 1; l <= vcs.valence[v-1]; l++)
   {  edgeinvert *ln = &vcs.vertlist[v-1][l-1];
      if(!(ln->moment<0||((IN_PRTCL|OUT_PRTCL)&ln->prop)||pseudop(ln->partcl)))
      {
         for( k=1;k<=momdep[ln->moment-1][0];k++) buff[k-1]=momdep[ln->moment-1][k];
         buff[k-1]=0;

         k=-1; while(buff[++k]) if(buff[k]<0) buff[k]=-buff[k];
         if( (2*strlen(buff) > nin+nout) ||
             (2*strlen(buff)==nin+nout && !strchr(buff,1))
           )
         { int ll=0;
           char buff2[MAXINOUT+1];
           for(k=1;k<=nin+nout;k++){ if(!strchr(buff,k)) buff2[ll++]=k;}
           buff2[ll]=0;
           strcpy(buff,buff2);
         }
         k=0;
         while(buff[k])
         {  if(!k) k++;
            if(buff[k]<buff[k-1])
            { int c=buff[k];
              buff[k]=buff[k-1];
              buff[k-1]=c;
              k--;
            } else k++;
         }

         strcpy(denom[denrno].momStr,buff);
         denom[denrno].power = 1;

         denom[denrno].mass=modelVarPos(prtclbase[ln->partcl-1].massidnt);
         if(ttypepropag(v,l)&&!tWidths) denom[denrno].width = 0; else 
         denom[denrno].width=modelVarPos(prtclbase[ln->partcl-1].imassidnt);

         for (k = 0; k < denrno; k++)
         if ( !strcmp(denom[denrno].momStr,denom[k].momStr) &&
               denom[denrno].mass  ==  denom[k].mass &&
               denom[denrno].width ==  denom[k].width )
         {  ++(denom[k].power); goto label_1;}
         denrno++;
label_1:;
      }
   }
}


void  denominatorStatistic(int nsub, 
   int * n_swidth, int *n_twidth, int * n_0width, denlist * allDenominators, 
   FILE * fd, int for12)
{ 
   int i;
   catrec    cr;
   denlist    den_, den_tmp;
   deninforec   dendescript;
   int ArchNum=0;
    
   (*n_swidth)  = 0;
   (*n_twidth)  = 0;
   (*n_0width) = 0;

   den_ =NULL;
   
   fseek(catalog,0,SEEK_SET);
   while (FREAD1(cr,catalog))
   {
      
      if (cr.nsub_ == nsub)
      { 
         whichArchive(cr.nFile,'r',for12);
/*          if(ArchNum!=0 && ArchNum!=cr.nFile) fclose(archiv);
         if(ArchNum==0 || ArchNum!=cr.nFile) 
         { char archiveName[40];
           ArchNum=cr.nFile;
           sprintf(archiveName,ARCHIV_NAME,ArchNum);
           if(for12) strcat(archiveName,"2");
           archiv=fopen(archiveName,"rb");
         }
            
*/      
         dendescript.cr_pos = ftell(catalog) - sizeof(cr);

         fseek(archiv,cr.denompos,SEEK_SET);
         readDenominators();
         dendescript.tot_den=denrno; 

         for (i = 0; i < dendescript.tot_den; i++)
         {  
            dendescript.denarr[i].power=denom[i].power;
            dendescript.denarr[i].width=denom[i].width;
            den_tmp = den_;  
            while (den_tmp != NULL &&
              (  strcmp(denom[i].momStr,den_tmp->momStr)
              ||  denom[i].mass!=den_tmp->mass 
              ||  denom[i].width!=den_tmp->width ) ) den_tmp=den_tmp->next;
            if(den_tmp == NULL)
            {  
               den_tmp = (denlist)getmem_((unsigned)sizeof(denlistrec));
               den_tmp->next = den_;
               strcpy(den_tmp->momStr,denom[i].momStr);
               den_tmp->mass=denom[i].mass;
               den_tmp->width=denom[i].width;
               den_tmp->stype= stype(denom[i].momStr);
               den_ = den_tmp;
               if(denom[i].width) 
               { if(den_tmp->stype) den_tmp->order_num= ++(*n_swidth);
                         else       den_tmp->order_num= ++(*n_twidth);
               }  else              den_tmp->order_num= ++(*n_0width);
            }
            dendescript.denarr[i].order_num=den_tmp->order_num;
            dendescript.denarr[i].stype=den_tmp->stype;
         }      
         if(fd) FWRITE1(dendescript,fd);
      }  /* if CR.nsub_ =nsub */
   }
   
/*   if(ArchNum) fclose(archiv);  */
   whichArchive(0,0,0);
  *allDenominators=den_;
}
