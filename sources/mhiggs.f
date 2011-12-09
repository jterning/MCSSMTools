      SUBROUTINE MHIGGS(PAR,IFAIL)

***********************************************************************
*
*      This subroutine computes the Higgs masses and couplings in
*       the NMSSM. The relevant parameters are read from
*       COMMON/QGAUGE, /QHIGGS, /QQUARK, /QPAR, /QEXT
*       (computed in RUNPAR) and  /SUSYEXT
*       The squark masses and mixing parameters are read from
*       COMMON/RADCOR (computed in MSFERM)
*
*      On output:
*
*      SMASS(1-3): CP-even masses (ordered)
*
*      SCOMP(1-3,1-3): Mixing angles: if HB(I) are the bare states,
*        HB(I) = Re(H1), Re(H2), Re(S), and HM(I) are the mass eigenstates, 
*        the convention is HB(I) = SUM_(J=1,3) SCOMP(J,I)*HM(J)
*        which is equivalent to HM(I) = SUM_(J=1,3) SCOMP(I,J)*HB(J)
*
*      PMASS(1-2): CP-odd masses (ordered)
*
*      PCOMP(1-2,1-2): Mixing angles: if AB(I) are the bare states,
*        AB(I) = Im(H1), Im(H2), Im(S), and AM(I) are the mass eigenstates, 
*        the convention is 
*        AM(I) = PCOMP(I,1)*(COSBETA*AB(1)+SINBETA*AB(2))
*              + PCOMP(I,2)*AB(3)
*
*      CMASS: Charged Higgs mass
*
*      IFAIL = 0         OK
*              1,3,5,7   m_h1^2 < 0
*              2,3,6,7   m_a1^2 < 0
*              4,5,6,7   m_h+^2 < 0
*
*      The precision in the computation of the lightest Higgs mass is:
*
*      Terms ~ ht^4/hb^4 from (s)top/(s)bottom-loops are computed
*      exactly in the mixing parameters.
*
*      Terms ~ g^2*(ht^2/hb^2) (where g is an electro-weak gauge coupling)
*      are taken into account due to the wave function renormalizations
*      factors, finite self energies (pole masses) and corrections from
*       stop/sbottom D terms
*
*      Leading logs ~ (g, l, k)^4 are added explicitely.
*
*      Leading double logs from two loops ~ ht/hb^6 and ht/hb^4*alpha_s
*      are taken into account.
*
*      For heavy higgses (with masses mhh > mtop) the leading log
*      contributions ~ (ht^2/hb^2)*log(mhh/mtop) to the pole masses
*      are taken into account, but not the corresponding effects on the
*      mixing angles.
*
*      All mixing angles are at the scale m_top.
*
*      The dominant errors come from one loop terms ~ (g,l,k)^4 without
*      large logs, and from two loop terms without large double logs.
*
************************************************************************

      IMPLICIT NONE

      INTEGER IFAIL,I,J,PFLAG

      DOUBLE PRECISION PAR(*),VEC3(3,3),VEC2(2,2),EPS
      DOUBLE PRECISION SMASS(3),SCOMP(3,3),PMASS(2),PCOMP(2,2),CMASS
      DOUBLE PRECISION H(3,3),A(2,2),Alshift,B,X,BT,BB,MGAU
      DOUBLE PRECISION mst1,mst2,s2t,msb1,msb2,s2b,XT,XB,M1,M2,T
      DOUBLE PRECISION emt,fmt,gmt,emb,fmb,gmb
      DOUBLE PRECISION QSTSB,PI,C2TW,sferm,bos
      DOUBLE PRECISION sb,cb,s2,rt,rb,MS2,MP2,M12
      DOUBLE PRECISION LM2,Lmu,Lnu,LM2mu,Lmunu,LQZ,LA,LP,LS,LPP,P1,P2
      DOUBLE PRECISION ALSMZ,ALEMMZ,GF,g1,g2,S2TW,gb,gt,Db,Dt,D1,D2
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW,subdet
      DOUBLE PRECISION G1Q,G2Q,GQ,ALSQ
      DOUBLE PRECISION ZHU,ZHD,ZS,H1Q,H2Q,TANBQ
      DOUBLE PRECISION HTQ,HBQ,MTOPQ,MBOTQ
      DOUBLE PRECISION LQ,KQ,ALQ,AKQ,MUQ,NUQ
      DOUBLE PRECISION MA2,COEF,GAUGE,htau
      DOUBLE PRECISION XIF,XIS,MUP,MSP,M3H,GMCOMB
      DOUBLE PRECISION MH(3),MA(2),MHC
      DOUBLE PRECISION At,Tad,f,mSing2,majS,yy,vv,c2b,sin2b

      LOGICAL SLmom
*      INTEGER ERRSL, INDEXGB
*      DOUBLE PRECISION SLMSS(3),SLMAA(3),SLOS(3,3),SLOP(3,3)
*      DOUBLE PRECISION PCOMPNEW(2,2),PMASSNEW(2),SLCMASS
      DOUBLE PRECISION mZpole,mWpole,mtpole,mbmb,slmtau,asmz,SLGF
      DOUBLE PRECISION MQ3P,MU3P,MD3P,ATP,ABP

      COMMON/GAUGE/ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/HIGGSPEC/SMASS,SCOMP,PMASS,PCOMP,CMASS
      COMMON/RADCOR/mst1,mst2,s2t,msb1,msb2,s2b,XT,XB
      COMMON/STSBSCALE/QSTSB      
      COMMON/QGAUGE/G1Q,G2Q,GQ,ALSQ
      COMMON/QHIGGS/ZHU,ZHD,ZS,H1Q,H2Q,TANBQ
      COMMON/QQUARK/HTQ,HBQ,MTOPQ,MBOTQ
      COMMON/QPAR/LQ,KQ,ALQ,AKQ,MUQ,NUQ
      COMMON/ALSHIFT/ALSHIFT
      COMMON/QEXT/XIF,XIS,MUP,MSP,M3H
      COMMON/PFLAG/PFLAG
      COMMON/SMINPUTS/mZpole,mWpole,mtpole,mbmb,SLmtau,asmz,SLGF
      COMMON/RADCOR2/MQ3P,MU3P,MD3P,ATP,ABP
      COMMON/HMO/MH,MA,MHC
      COMMON/POTENTIAL/Tad,f,mSing2,majS

      IFAIL=0
*   	  write(0,*)"Tad",Tad,f

      EPS=1.D-8
      PI=4d0*DATAN(1d0)
      COEF=1d0/(16d0*PI**2)
      
*   Trig. functions of beta
      cb=1d0/DSQRT(1d0+tanbq**2)
      sb=tanbq*cb
      s2=2d0*sb*cb
*	  write(0,*) "tanbq ",tanbq
      M1=PAR(20)
      M2=PAR(21)
      MA2=PAR(23)**2
      htau=MTAU/H2Q

*   Weak angle theta_W (S2TW= sin(theta_W)**2):
      C2TW=1d0-S2TW

*   Approximate value for the Singlet-like CP even Higgs mass squared:
      MS2=MAX(NUQ*(AKQ+4d0*NUQ+3d0*MUP)-LQ*(XIS+XIF*MUP)/MUQ,MZ**2)

*   Approximate value for the Singlet-like CP odd Higgs mass squared:
      MP2=-3d0*NUQ*AKQ-XIF*(4d0*KQ+LQ*MUP/MUQ)
     .   -2d0*MSP-MUP*NUQ-LQ*XIS/MUQ

*   Approximate value for the off-diag. CP odd mass matrix element:
      M12=LQ*DSQRT(h1q**2+h2q**2)*(ALQ-2d0*NUQ-MUP)

*   One loop functions for stop/sbottom loop corrections

      rt= 3d0/2d0*COEF*htq**2
      IF(mst1.NE.mst2)THEN
       fmt= (mst2*DLOG(mst2/QSTSB)-mst1*DLOG(mst1/QSTSB))/
     .     (mst2-mst1)-1d0
       gmt= s2t**2*((mst2+mst1)/(mst2-mst1)*DLOG(mst2/mst1)-2d0)
       emt= -mtopq*s2t*DLOG(mst2/mst1)
      ELSE
       fmt= DLOG(mst1/QSTSB)
       gmt= 0d0
       emt= 0d0
      ENDIF

      rb= 3d0/2d0*COEF*hbq**2
      IF(msb1.NE.msb2)THEN
       fmb= (msb2*DLOG(msb2/QSTSB)-msb1*DLOG(msb1/QSTSB))/
     .    (msb2-msb1)-1d0
       gmb= s2b**2*((msb2+msb1)/(msb2-msb1)*DLOG(msb2/msb1)-2d0)
       emb= -mbotq*s2b*DLOG(msb2/msb1)
      ELSE
       fmb= DLOG(msb1/QSTSB)
       gmb= 0d0
       emb= 0d0
      ENDIF
       fmt= 0d0
       gmt= 0d0
       emt= 0d0    
       fmb= 0d0
       gmb= 0d0
       emb= 0d0  
*  The subsequent shifts in Alambda simplify the expressions for the
*  one loop rad. corrs. below.
*  The parameter ALQ is defined at the scale QSTSB

      MGAU= (G1q*M1+3d0*G2q*M2)
      LM2mu= DLOG(MAX(M2**2,MUQ**2,MZ**2)/QSTSB)
      LM2MU= 0d0
      Alshift= ALQ+2d0*rt*Atp*fmt+2d0*rb*Abp*fmb
     .      +MGAU*LM2MU*COEF
      B= Alshift+NUQ
      GMCOMB= LQ*XIF+MUP*MUQ+M3H

*   Tree level CP-even Higgs mass matrix
*	  write(0,*)"CP-even Higgs",GQ,h1q, MUQ,B,GMCOMB,tanbq
*	  write(0,*)h2q, NUQ,AKQ,ALQ
*	  H(1,1) = GQ*h1q**2 + (MUQ*B+GMCOMB)/tanbq
*      H(2,2) = GQ*h2q**2 + (MUQ*B+GMCOMB)*tanbq
*      H(3,3) = LQ**2*(Alshift+MUP)*h1q*h2q/MUQ
*     .   + NUQ*(AKQ+4d0*NUQ+3d0*MUP)-LQ/MUQ*(XIS+XIF*MUP)
*      H(1,2) = (2d0*LQ**2-GQ)*h1q*h2q - MUQ*B-GMCOMB
*      H(1,3) = LQ*(2d0*MUQ*h1q - (B+NUQ+MUP)*h2q)
*      H(2,3) = LQ*(2d0*MUQ*h2q - (B+NUQ+MUP)*h1q)


*******************************************************************
*******************************************************************
*******************************************************************
*******************************************************************
*******************************************************************
* Mass matrix for minimal composite supersymmetric standard model
*
*	  write(0,*)"CP even Mass Matrix"
      At=LQ * ALQ
      yy=LQ
      c2b=cb**2-sb**2
      sin2b= 2d0 *sb *cb
*	  c3b=c2b* cb-sin2b *sb

      vv=4 * MZ**2/(g1q+g2q)
*	  vv=246.0d0**2

*	  write(0,*)"MZ,   Tad ,  f, At, LQ,      ALQ"
*	  write(0,*)MZ,Tad ,f, At, LQ, ALQ
*	  write(0,*)"tanbq, v, mSing, s2b, c2b"
*	  write(0,*)tanbq, DSQRT(vv),DSQRT(mSing2),sin2b, c2b


      H(1,1) =  (2*mSing2*MZ**2 + At**2*vv + 
     -    MZ**2*vv*yy**2 + 
     -    c2b*(-2*mSing2*MZ**2 + 
     -       vv*(At**2 - MZ**2*yy**2)) + 
     -    (2*(2*At*Tad + 
     -         f**2*yy**2*(2*mSing2 + vv*yy**2)))/
     -     tanbq)/(4*mSing2 + 2*vv*yy**2)
*      write(0,*)2*mSing2*MZ**2,At**2*vv,MZ**2*vv*yy**2
*      write(0,*)c2b*(-2*mSing2*MZ**2 + vv*(At**2 - MZ**2*yy**2)) 
*      write(0,*)4*mSing2 + 2*vv*yy**2
      H(2,2) = (At**2*sb**2*vv + 
     -    cb**2*MZ**2*(2*mSing2 + vv*yy**2) + 
     -    tanbq*(2*At*Tad + 
     -       f**2*yy**2*(2*mSing2 + vv*yy**2)))/
     -  (2*mSing2 + vv*yy**2)
      H(3,3) =mSing2 + (vv*yy**2)/2.0D0
      H(1,2) = -((2*At*Tad + 2*f**2*mSing2*yy**2 - 
     -      mSing2*sin2b*vv*yy**2 + f**2*vv*yy**4 + 
     -      cb*sb*(2*mSing2*MZ**2 + 
     -         vv*(At**2 + MZ**2*yy**2 - vv*yy**4)))
     -  /(2*mSing2 + vv*yy**2))
      H(1,3) = (DSQRT(vv)*(At*cb - (sb*(4*Tad + At*sin2b*vv)*
     -         yy**2)/(2*mSing2 + vv*yy**2)))/Sqrt(2.0)
      H(2,3) = (DSQRT(vv)*(At*sb - (cb*(4*Tad + At*sin2b*vv)*
     -    yy**2)/(2*mSing2 + vv*yy**2)))/Sqrt(2.0)
      H(2,1) = H(1,2)     
      H(3,1) = H(1,3) 
      H(3,2) = H(2,3) 


*   1-loop radiative corrections

*      H(1,1)= H(1,1) + rt*(4d0*Atp*emt - Atp**2*gmt
*     .      + 4d0*mtopq**2*DLOG(mst1*mst2/mtopq**4))
*      H(2,2)= H(2,2) - rt*MUQ**2*gmt
*      H(3,3)= H(3,3) - rt*LQ**2*h2q**2*gmt
*      H(1,2)= H(1,2) + rt*MUQ*(Atp*gmt - 2d0*emt)      
*      H(1,3)= H(1,3) + rt*LQ*h2q*(Atp*gmt - 2d0*emt)
*      H(2,3)= H(2,3) + rt*(4d0*LQ*h2q*MUQ*fmt - LQ*h2q*MUQ*gmt)

*      H(1,1)= H(1,1) - rb*MUQ**2*gmb
*      H(2,2)= H(2,2) + rb*(4d0*Abp*emb - Abp**2*gmb
*     .      + 4d0*mbotq**2*DLOG(msb1*msb2/mbotq**4))
*      H(3,3)= H(3,3) - rb*LQ**2*h1q**2*gmb      
*      H(1,2)= H(1,2) + rb*MUQ*(Abp*gmb - 2d0*emb)      
*      H(1,3)= H(1,3) + rb*(4d0*LQ*h1q*MUQ*fmb - LQ*h1q*MUQ*gmb)
*      H(2,3)= H(2,3) + rb*LQ*h1q*(Abp*gmb - 2d0*emb)

*   Corrections from higgs/stop/sbottom couplings from D-terms

      gt= GQ/2d0-2d0*g1q/3d0
      Dt= (PAR(7)-PAR(8)+gt*(h2q**2-h1q**2))/2d0
      D1= 3d0*COEF*(-GQ/4d0*emt + gt/2d0*Dt/Xt*gmt)
      D2= 3d0*COEF*(-gt*Dt/Xt*emt
     .    -GQ/2d0*mtopq**2*DLOG(mst1*mst2/QSTSB**2))
      D2=0d0
*      H(1,1)= H(1,1) + 2d0*Atp*D1 + 2d0*D2
*      H(2,2)= H(2,2) + 2d0*MUQ/tanbq*D1
*      H(1,2)= H(1,2) - 1d0/tanbq*D2 - (MUQ+Atp/tanbq)*D1
*      H(1,3)= H(1,3) - LQ*h2q*D1
*      H(2,3)= H(2,3) + LQ*h2q/tanbq*D1

      gb= GQ/2d0-g1q/3d0
      Db= (PAR(7)-PAR(9)+gb*(h1q**2-h2q**2))/2d0
      D1= 3d0*COEF*(-GQ/4d0*emb + gb/2d0*Db/Xb*gmb)
      D2= 3d0*COEF*(-gb*Db/Xb*emb
     .    -GQ/2d0*mbotq**2*DLOG(msb1*msb2/QSTSB**2))
*      H(1,1)= H(1,1) + 2d0*MUQ*tanbq*D1
*      H(2,2)= H(2,2) + 2d0*Abp*D1 + 2d0*D2
*      H(1,2)= H(1,2) - tanbq*D2 - (MUQ+Abp*tanbq)*D1
*      H(1,3)= H(1,3) + LQ*h1q*tanbq*D1
*      H(2,3)= H(2,3) - LQ*h1q*D1

*   2-loop terms

      T= DLOG(QSTSB/MT**2)

*      H(1,1)= H(1,1) + 6d0*COEF**2*htq**4*h1q**2*
*     .       (T**2*(64d0*PI*ALSQ+4d0/3d0*g1q-3d0*htq**2*sb**2
*     .       +3d0*hbq**2*cb**2)+((dlog(MAX(QSTSB,MA2,MT**2)/MT**2))**2
*     .       -(dlog(MAX(MA2,MT**2)/MT**2))**2)*
*     .       (-3d0*htq**2*cb**2-hbq**2*(3d0*cb**2+1d0)))

*      H(2,2)= H(2,2) + 6*COEF**2*hbq**4*h2q**2*
*     .       (T**2*(64d0*PI*ALSQ-2d0/3d0*g1q-3d0*hbq**2*cb**2
*     .       +3d0*htq**2*sb**2)+((dlog(MAX(QSTSB,MA2,MT**2)/MT**2))**2
*     .       -(dlog(MAX(MA2,MT**2)/MT**2))**2)*
*     .       (-3d0*hbq**2*sb**2-htq**2*(3d0*sb**2+1d0)))

*   Leading-log electroweak contributions (1 loop):

*   a) Sfermion contributions

      sferm= 2d0*COEF*GQ*MZ**2*
     .   ((S2TW**2/6d0+3d0/2d0*C2TW**2)*DLOG(MAX(PAR(7),MZ**2)/QSTSB)
     .   +4d0/3d0*S2TW**2*DLOG(MAX(PAR(8),MZ**2)/QSTSB)
     .   +1d0/3d0*S2TW**2*DLOG(MAX(PAR(9),MZ**2)/QSTSB)
     .   +(S2TW**2/3d0+C2TW**2*3d0)*DLOG(MAX(PAR(15),MZ**2)/QSTSB)
     .   +8d0/3d0*S2TW**2*DLOG(MAX(PAR(16),MZ**2)/QSTSB)
     .   +2d0/3d0*S2TW**2*DLOG(MAX(PAR(17),MZ**2)/QSTSB)
     .   +(S2TW**2/2d0+C2TW**2/2d0)*DLOG(MAX(PAR(10),MZ**2)/QSTSB)
     .   +S2TW**2*DLOG(MAX(PAR(11),MZ**2)/QSTSB)
     .   +(S2TW**2+C2TW**2)*DLOG(MAX(PAR(18),MZ**2)/QSTSB)
     .   +2d0*S2TW**2*DLOG(MAX(PAR(19),MZ**2)/QSTSB))

*      H(1,1)= H(1,1) + sferm*sb**2
*      H(2,2)= H(2,2) + sferm*cb**2
*      H(1,2)= H(1,2) - sferm*sb*cb

*    b) Chargino/neutralino contributions

      LQZ= DLOG(QSTSB/MZ**2)
      LM2= DLOG(MAX(M2**2,MZ**2)/QSTSB)
      Lmu= DLOG(MAX(MUQ**2,MZ**2)/QSTSB)
      Lnu= DLOG(MAX((2d0*NUQ+MUP)**2,MZ**2)/QSTSB)
      Lmunu= DLOG(MAX(MUQ**2,(2d0*NUQ+MUP)**2,MZ**2)/QSTSB)
      LQZ= 0
      LM2= 0
      Lmu= 0
      Lnu= 0
      Lmunu= 0

*      H(1,1)= H(1,1) + COEF*
*     .  ((MZ**2*SB**2*GQ*(-20d0+32d0*S2TW-16d0*S2TW**2))*LM2MU
*     .  -4d0*(MUQ*(NUQ+MUP/2d0)*LQ**2/TANBQ+MZ**2*SB**2*LQ**4/GQ)
*     .  *LMUNU)

*      H(2,2)= H(2,2) + COEF*
*     .  ((MZ**2*CB**2*GQ*(-20d0+32d0*S2TW-16d0*S2TW**2))*LM2MU
*     .  -4d0*(MUQ*(NUQ+MUP/2d0)*LQ**2*TANBQ+MZ**2*CB**2*LQ**4/GQ)
*     .  *LMUNU)

*      H(3,3)= H(3,3) + COEF*
*     .  (2d0*KQ*(-16d0*KQ*NUQ**2-12d0*KQ*NUQ*MUP+LQ*MUP**3/MUQ)*LNU
*     .  -8.*LQ**2*MUQ**2*LMU
*     .  +LQ**3*MZ**2*MUP/(MUQ*GQ)*(4d0*KQ-LQ*S2)*Lmunu)

*      H(1,2)= H(1,2) + COEF*
*     .  (4d0*(MUQ*(NUQ+MUP/2d0)*LQ**2-MZ**2*SB*CB*LQ**4/GQ)*LMUNU
*     .  -4d0*GQ*MZ**2*SB*CB*LM2MU)

*      H(1,3)= H(1,3) + COEF*MZ/DSQRT(GQ)*
*     .  (LQ*GQ*MUQ*SB*(-12d0+8d0*S2TW)*LM2MU
*     .  +4d0*LQ**2*(CB*LQ*(2d0*NUQ+MUP/2d0)-SB*(LQ*MUQ
*     .  +4d0*KQ*(NUQ+MUP/2d0)))*LMUNU)

*      H(2,3)= H(2,3) + COEF*MZ/DSQRT(GQ)*
*     .  (LQ*GQ*MUQ*CB*(-12d0+8d0*S2TW)*LM2MU
*     .  +4d0*LQ**2*(SB*LQ*(2d0*NUQ+MUP/2d0)-CB*(LQ*MUQ
*     .  +4d0*KQ*(NUQ+MUP/2d0)))*LMUNU)

*    c) Higgs loop contributions
*    (Only if all masses squared are positive, and only to the lighter
*    CP-even doublet-like state)

      subdet= H(3,3)*(sb**2*H(1,1)+cb**2*H(2,2)+s2*H(1,2))
     .       -(sb*H(1,3)+cb*H(2,3))**2

      IF(subdet.GT.0d0)THEN

       P2= MAX(MA2+MP2,MZ**2)
       P1= MAX((MA2*MP2-M12**2)/P2,MZ**2)
       LA= DLOG(MAX(MA2,MZ**2)/MZ**2)
       LS= DLOG(MS2/MZ**2)
       LP= DLOG(P2/MZ**2)
       LPP= DLOG(P2/P1)
       

       LA= 0
       LS= 0
       LP= 0
       LPP= 0

       bos= COEF*MZ**2/GQ*((GQ**2*(2d0*S2TW**2
     .   -2d0*S2TW*(1d0+s2**2)-11d0/4d0*s2**4+5d0*s2**2+3d0/4d0)
     .   +GQ*LQ**2*(2d0*S2TW*s2**2+11d0/2d0*s2**4-15d0/2d0*s2**2
     .   -1d0)+LQ**4*(-11d0/4d0*s2**4+5d0/2d0*s2**2+1d0))*LA
     .   +(LQ**2*(LQ-KQ*s2)**2+3d0*LQ**2/MS2*(GQ+(LQ**2-GQ)*s2**2)
     .   *(2d0*MUQ-s2*(ALQ+2d0*NUQ+MUP))**2
     .   -LQ**4/MS2**2*(2d0*MUQ-s2*(ALQ+2d0*NUQ+MUP))**4)*LS
     .   +(GQ**2/4d0*(1d0-s2**4)+GQ*LQ**2*(1d0/2d0*s2**4
     .   +1d0/2d0*s2**2-1d0)+LQ**4*(-1d0/4d0*s2**4-1d0/2d0*s2**2
     .   +1d0)+LQ**2*(LQ+KQ*s2)**2)*LP
     .   -((GQ-LQ**2)**2/2d0*MP2/P2*s2**2*(1d0-s2**2)
     .   +(LQ*MA2*(LQ+KQ*s2)-LQ**2*(ALQ-2d0*NUQ-MUP)**2
     .   -MP2/2d0*(GQ*(1d0-s2**2)
     .   -LQ**2*(2d0-s2**2)))**2/P2**2-LQ*MA2*MP2*(LQ+KQ*s2)
     .   *(GQ*(1d0-s2**2)-LQ**2*(2d0-s2**2))/P2**2)*LPP
     .   +(GQ**2*(-4d0+S2**2+2d0*S2TW*(1d0+S2**2)-2d0*S2TW**2)
     .   +GQ*LQ**2*(2d0+S2**2-2d0*S2**2*S2TW)
     .   -LQ**4*(4d0+2d0*S2**2)-2d0*LQ**2*KQ**2*S2**2)*LQZ)

*      H(1,1)= H(1,1) + bos*sb**2
*      H(2,2)= H(2,2) + bos*cb**2
*      H(1,2)= H(1,2) + bos*sb*cb

      ENDIF

*   d) Gauge Loop Corrections

      GAUGE= MZ**2*GQ*COEF*(-9d0+12d0*S2TW-6d0*S2TW**2)*LQZ
      
*      H(1,1)= H(1,1) + GAUGE*sb**2
*      H(2,2)= H(2,2) + GAUGE*cb**2
*      H(1,2)= H(1,2) + GAUGE*sb*cb

*   Take care of the Z factors

*      H(1,1)= H(1,1)/ZHU
*      H(2,2)= H(2,2)/ZHD
*      H(3,3)= H(3,3)/ZS
*      H(1,2)= H(1,2)/DSQRT(ZHU*ZHD)
*      H(1,3)= H(1,3)/DSQRT(ZHU*ZS)
*      H(2,3)= H(2,3)/DSQRT(ZHD*ZS)

*   Diagonalization

      write(0,*)"CP Even Mass Matrix"
      DO I= 1,3
         write(0,*)H(I,1),H(I,2),H(I,3)
      ENDDO
      
      CALL DIAGN(3,H,MH,VEC3,EPS)
      CALL SORTN(3,MH,VEC3)
      DO I= 1,3
       DO J= 1,3
        SCOMP(I,J)= VEC3(J,I)
       ENDDO
      ENDDO

*   CP even pole masses

      DO I= 1,3
      IF(MH(I).GT.0d0)THEN
       IF(MH(I).GT.4d0*MT**2)THEN
        X= DSQRT(1d0-4d0*MT**2/MH(I))
        BT= 2d0-X*DLOG((1d0+X)/(1d0-X))
       ELSE
        X= DSQRT(MH(I)/(4d0*MT**2-MH(I)))
        BT= 2d0*(1d0-DATAN(X)/X)
       ENDIF
       IF(MH(I).GT.4d0*MB**2)THEN
        X= DSQRT(1d0-4d0*MB**2/MH(I))
        BB= 2d0-X*DLOG((1d0+X)/(1d0-X))
       ELSE
        X= DSQRT(MH(I)/(4d0*MB**2-MH(I)))
        BB= 2d0*(1d0-DATAN(X)/X)
       ENDIF
       SMASS(I)= MH(I)
*     .  - 2d0*rt*SCOMP(I,1)**2*(MH(I)-4d0*MT**2)*BT
*     .  - 2d0*rb*SCOMP(I,2)**2*(MH(I)*DLOG(MT**2/MB**2)
*     .   +(MH(I)-4d0*MB**2)*BB)
       MH(I)= DSQRT(MH(I))
      ELSE
       SMASS(I)= MH(I)
      ENDIF
         write(0,*)"SMASS ",I,SMASS(I)
      ENDDO

      IF(SMASS(1).LT.1d0)THEN
      	IFAIL= IFAIL+1
      	IFAIL=0
      ELSE
       DO I= 1,3
        SMASS(I)= DSQRT(SMASS(I))
       ENDDO
       write(0,*)"Higgs mass = ",SMASS(1)
      ENDIF

*   CP-odd Higgs mass matrix including
*   1-loop top/bottom radiative corrections

*      A(1,1)= (MUQ*B+GMCOMB)*(h1q/(ZHD*h2q)+h2q/(ZHU*h1q))
*      A(2,2)= LQ**2*h1q*h2q/MUQ*(B+3d0*NUQ+MUP)-3d0*AKQ*NUQ
*     .       -XIF*(4d0*KQ+LQ*MUP/MUQ)-2d0*MSP-MUP*NUQ-LQ*XIS/MUQ
*      A(1,2)= LQ*(Alshift-2d0*NUQ-MUP)
*     .       *DSQRT(h1q**2/ZHD+h2q**2/ZHU)


*******************************************************************
*******************************************************************
*******************************************************************
*******************************************************************
*******************************************************************
* CP odd Mass matrix for minimal composite supersymmetric standard model
*

      A(1,1)= (4*At*Tad + At**2*sin2b*vv + 
     -    2*f**2*yy**2*(2*mSing2 + vv*yy**2))/
     -  (sin2b*(2*mSing2 + vv*yy**2))
      A(2,2)=mSing2 + (vv*yy**2)/2.0
      A(1,2)= -(At*DSQRT(vv))/DSqrt(2.0d0)
     
     
*   Diagonalization

      CALL DIAGN(2,A,MA,VEC2,EPS)
      CALL SORTN(2,MA,VEC2)
      DO I= 1,2
       DO J= 1,2
        PCOMP(I,J)= VEC2(J,I)
       ENDDO
      ENDDO

*   CP odd pole masses

      DO I= 1,2
      IF(MA(I).GT.0d0)THEN
       IF(MA(I).GT.4d0*MT**2)THEN
        X= DSQRT(1d0-4d0*MT**2/MA(I))
        BT= 2d0-X*DLOG((1d0+X)/(1d0-X))
       ELSE
        X= DSQRT(MA(I)/(4d0*MT**2-MA(I)))
        BT= 2d0*(1d0-DATAN(X)/X)
       ENDIF
       IF(MA(I).GT.4d0*MB**2)THEN
        X= DSQRT(1d0-4d0*MB**2/MA(I))
        BB= 2d0-X*DLOG((1d0+X)/(1d0-X))
       ELSE
        X= DSQRT(MA(I)/(4d0*MB**2-MA(I)))
        BB= 2d0*(1d0-DATAN(X)/X)
       ENDIF
       PMASS(I)= MA(I)
*     .  - 2d0*rt*PCOMP(I,1)**2*cb**2*MA(I)*BT
*     .  - 2d0*rb*PCOMP(I,1)**2*sb**2*(MA(I)*DLOG(MT**2/MB**2)
*     .   +MA(I)*BB)
       MA(I)= DSQRT(MA(I))
      ELSE
       PMASS(I)= MA(I)
      ENDIF
      ENDDO

      IF(PMASS(1).LT.1d0)THEN
        IFAIL= IFAIL+2
        IFAIL=0
      ELSE
       DO I= 1,2
        PMASS(I)= DSQRT(PMASS(I))
       ENDDO
      ENDIF

*   Charged Higgs mass, 1-loop LL + mixing and 2-loop LL

*      MHC=((G2Q/2d0-LQ**2)*h1q*h2q+MUQ*B+GMCOMB)
*     .    *(h1q**2*ZHU+h2q**2*ZHD)/(h1q*h2q*ZHU*ZHD)
*     .    +COEF*(h1q**2*ZHU+h2q**2*ZHD)*(6d0*htQ**2*hbQ**2*T
*     .    +3d0/2d0*htq**2*hbq**2*(Atp+Abp)**2/QSTSB
*     .    -1d0/2d0*htq**2*hbq**2*(muQ**2-Atp*Abp)**2/QSTSB**2
*     .    -3d0/2d0*(htq**2+hbq**2)**2*muQ**2/QSTSB
*     .    +6d0*htq**2*hbq**2*COEF*(32d0*PI*ALSQ-htq**2-hbq**2)
*     .    *((dlog(MAX(QSTSB,MA2,MT**2)/MT**2))**2
*     .    -(dlog(MAX(MA2,MT**2)/MT**2))**2)
*     .    +G2Q**2*(3d0/4d0*DLOG(MAX(PAR(7),MZ**2)/QSTSB)
*     .    +3d0/2d0*DLOG(MAX(PAR(15),MZ**2)/QSTSB)
*     .    +1d0/4d0*DLOG(MAX(PAR(10),MZ**2)/QSTSB)
*     .    +1d0/2d0*DLOG(MAX(PAR(18),MZ**2)/QSTSB))
*     .    +G2Q*(7d0*G1Q-G2Q)/4d0*LQZ+2d0*G2Q*(G1Q-G2Q)*LM2mu)
 
*******************************************************************
* Charged Higgs mass for minimal composite supersymmetric standard model
*

      MHC=(-(vv*(-2*At**2 + 2*mSing2*yy**2 + 
     -         vv*yy**4)) + 
     -    (4*(2*At*Tad + 
     -         f**2*yy**2*(2*mSing2 + vv*yy**2)))/
     -     sin2b)/(4*mSing2 + 2*vv*yy**2)
      write(0,*)"Charged Higgs Mass ",DSQRT(MHC)
 
*   Charged Higgs pole mass (in the LLA only)

*      X=MHC/MT**2
*      IF(X.EQ.0d0)THEN
*       BT=1d0
*      ELSEIF(X.EQ.1d0)THEN
*       BT=0d0
*      ELSE
*       BT=(1d0-1d0/X)*DLOG(DABS(X-1d0))
*      ENDIF
*      CMASS= MHC+3d0*COEF*((htq**2*cb**2+hbq**2*sb**2)
*     .      *(MHC*(BT-2d0) - (MT**2+MB**2)*(BT-1d0))
*     .      - 4d0*HTQ*MT*HBQ*MB*cb*sb*(BT-1d0))

      CMASS= MHC
      IF(MHC.GT.0)THEN
       CMASS= DSQRT(MHC)
      ENDIF

*      IF(CMASS.LT.1d0)THEN
*        IFAIL= IFAIL+4
*        IFAIL=0
*      ELSE
*       CMASS= DSQRT(CMASS)
*      ENDIF

*****************************************************************

*      IF(PFLAG.GE.1)THEN

* Full 1- and 2-Loop corrections from Degrassi/Slavich,
* ``On the radiative corrections to the neutral Higgs boson masses in
*  the NMSSM,'' Nucl. Phys.  B {\bf 825} (2010) 119
* [arXiv:0907.4682 [hep-ph]].
*
*   Parameters for COMMON/SMINPUTS:

*      mzpole=mz
*      mwpole=mw
*      mtpole=mt
*      mbmb=mb
*      SLmtau=mtau
*      asmz=alsmz
*      slgf=gf
            
      IF(PFLAG.EQ.1) THEN
       SLmom  = .false.
      ENDIF
      
      IF(PFLAG.EQ.2) THEN
       SLmom  = .true.
      ENDIF

      IFAIL=0

*     With parameters at QSTSB
*     (except for slepton, gaugino and 1st gener. squark masses)

*      CALL FULLHIG(SLmom,LQ,KQ,MUQ,TANBQ,
*     .   dsqrt(MQ3P),dsqrt(MU3P),dsqrt(MD3P),
*     .   dsqrt(par(15)),dsqrt(par(16)),dsqrt(par(17)),
*     .   dsqrt(par(10)),dsqrt(par(11)),dsqrt(par(18)),dsqrt(par(19)),
*     .   ATP,ABP,par(14),
*     .   ALQ,AKQ,par(20),par(21),par(22),dsqrt(QSTSB),
*     .   SLMSS,SLMAA,SLOS,SLOP,slcmass,errsl)
     
*      write(*,*) "In mhiggs1: cmass, slcmass:",cmass,slcmass
      
*      cmass=slcmass

*      IF(ERRSL.NE.0) THEN
*       IF(errsl.eq.1) IFAIL=IFAIL+1
*       IF(errsl.eq.2) IFAIL=IFAIL+2
*       IF(errsl.eq.3) IFAIL=8
*       IF(errsl.eq.4) IFAIL=7
*       IF(errsl.eq.5) IFAIL=IFAIL+4
*       RETURN
*      ENDIF

*   Search for the Goldstone Boson:

*       INDEXGB=1
*       IF(dabs(SLOP(2,3)).le.dabs(SLOP(1,3))) INDEXGB=2
*       IF((dabs(SLOP(3,3)).le.dabs(SLOP(1,3))).and.
*     .   (dabs(SLOP(3,3)).le.dabs(SLOP(2,3)))) INDEXGB=3

*   New entries in PCOMP(2,2) in PCOMPNEW, try to avoid entries >1:

* *      IF(INDEXGB.EQ.1) THEN
* *      PCOMPNEW(1,1)=
* *    .     SIGN(MIN(DABS(SLOP(2,1))/sb,DABS(SLOP(2,2))/cb),SLOP(2,1))
* *      PCOMPNEW(1,2)=SLOP(2,3)
* *      PCOMPNEW(2,1)=
* *    .     SIGN(MIN(DABS(SLOP(3,1))/sb,DABS(SLOP(3,2))/cb),SLOP(3,1))
* *      PCOMPNEW(2,2)=SLOP(3,3)
* *      PMASSNEW(1)=SLMAA(2)
* *      PMASSNEW(2)=SLMAA(3)
* *      ENDIF
* *     
* *      IF(INDEXGB.EQ.2) THEN
* *      PCOMPNEW(1,1)=
* *    .     SIGN(MIN(DABS(SLOP(1,1))/sb,DABS(SLOP(1,2))/cb),SLOP(1,1))
* *      PCOMPNEW(1,2)=SLOP(1,3)
* *      PCOMPNEW(2,1)=
* *    .     SIGN(MIN(DABS(SLOP(3,1))/sb,DABS(SLOP(3,2))/cb),SLOP(3,1))
* *      PCOMPNEW(2,2)=SLOP(3,3)
* *      PMASSNEW(1)=SLMAA(1)
* *      PMASSNEW(2)=SLMAA(3)
* *      ENDIF

*       IF(INDEXGB.EQ.3) THEN
*       PCOMPNEW(1,1)=
*     .     SIGN(MIN(DABS(SLOP(1,1))/sb,DABS(SLOP(1,2))/cb),SLOP(1,1))
*       PCOMPNEW(1,2)=SLOP(1,3)
*       PCOMPNEW(2,1)=
*     .     SIGN(MIN(DABS(SLOP(2,1))/sb,DABS(SLOP(2,2))/cb),SLOP(2,1))
*       PCOMPNEW(2,2)=SLOP(2,3)
*       PMASSNEW(1)=SLMAA(1)
*       PMASSNEW(2)=SLMAA(2)
*       ENDIF

*      DO I=1,3
*        SMASS(I)=SLMSS(I)
*      ENDDO
      
*      DO I=1,2
*        PMASS(I)=PMASSNEW(I)
*      ENDDO

*      SCOMP(1,1)=SLOS(1,2)
*      SCOMP(1,2)=SLOS(1,1)      
*      SCOMP(1,3)=SLOS(1,3)
*      SCOMP(2,1)=SLOS(2,2)
*      SCOMP(2,2)=SLOS(2,1)
*      SCOMP(2,3)=SLOS(2,3)
*      SCOMP(3,1)=SLOS(3,2)
*      SCOMP(3,2)=SLOS(3,1)
*      SCOMP(3,3)=SLOS(3,3)
      
*      DO I= 1,2
*       DO J= 1,2
*        PCOMP(I,J)= PCOMPNEW(I,J)
*       ENDDO
*      ENDDO
      
*      ENDIF


      END
