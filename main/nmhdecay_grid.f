      PROGRAM MAIN

*  Program to compute the NMSSM Higgs masses, couplings and
*  Branching ratios, with experimental and theoretical constraints
*
*   On input:      
*
*      OMGFLAG      Computation of relic density: 0=no, 1=yes
*      NMSFLAG      Computation of sparticle decays: 0=no, 2=yes
*
*      PAR(1) = lambda
*      PAR(2) = kappa
*      PAR(3) = tan(beta)
*      PAR(4) = mu (effective mu term = lambda*s)
*      PAR(5) = Alambda (if MA is not an input)
*      PAR(6) = Akappa
*      PAR(7) = mQ3**2
*      PAR(8) = mU3**2
*      PAR(9) = mD3**2
*      PAR(10) = mL3**2
*      PAR(11) = mE3**2
*      PAR(12) = AU3
*      PAR(13) = AD3
*      PAR(14) = AE3
*      PAR(15) = mQ2**2
*      PAR(16) = mU2**2
*      PAR(17) = mD2**2
*      PAR(18) = mL2**2
*      PAR(19) = mE2**2
*      PAR(20) = M1
*      PAR(21) = M2
*      PAR(22) = M3
*      PAR(23) = MA (diagonal doublet CP-odd mass matrix element)
*      PAR(24) = MP (diagonal singlet CP-odd mass matrix element)
*      PAR(25) = AE2
*
*      All these parameters are assumed to be defined in DRbar at the scale Q2
*      Q2 is either defined by the user in the input file or computed as Q2 = mQ3*mU3
*
*      The input file contains lower and upper bounds as well as number of steps
*      for the parameters PAR(1..6) on which the scan is performed
*
*   On output:
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
*      CU,CD,CV,CJ,CG(i) Reduced couplings of h1,h2,h3 (i=1,2,3) or
*                        a1,a2 (i=4,5) to up type fermions, down type
*                        fermions, gauge bosons, gluons and photons
*                        Note: CV(4)=CV(5)=0
*
*      WIDTH(i) Total decay width of h1,h2,h3,a1,a2 (i=1..5)
*               with the following branching ratios:
*      BRJJ(i) h1,h2,h3,a1,a2 -> gluon gluon
*      BRMM(i)        "       -> mu mu
*      BRLL(i)        "       -> tau tau
*      BRSS(i)        "       -> ss
*      BRCC(i)        "       -> cc
*      BRBB(i)        "       -> bb
*      BRTT(i)        "       -> tt
*      BRWW(i)        "       -> WW (BRWW(4)=BRWW(5)=0)
*      BRZZ(i)        "       -> ZZ (BRZZ(4)=BRZZ(5)=0)
*      BRGG(i)        "       -> gamma gamma
*      BRZG(i)        "       -> Z gamma
*      BRHIGGS(i)   (i=1..5)  -> other Higgses, including:
*        BRHAA(i,j)   hi      -> a1a1, a1a2, a2a2 (i=1..3, j=1..3)
*        BRHCHC(i)    hi      -> h+h- (i=1..3)
*        BRHAZ(i,j)   hi      -> Zaj  (i=1..3)
*        BRHCW(i)  h1,h2,h3   -> h+W- (i=1..3), a1,a2 -> h+W- (i=4,5)
*        BRHHH(i)     h2      -> h1h1, h3-> h1h1, h1h2, h2h2 (i=1..4)
*        BRAHA(i)     a2      -> a1hi (i=1..3)
*        BRAHZ(i,j)   ai      -> Zhj  (i=1,2, j=1..3)
*      BRSUSY(i)    (i=1..5)  -> susy particles, including:
*        BRNEU(i,j,k)         -> neutralinos j,k (i=1..5, j,k=1..5)
*        BRCHA(i,j)           -> charginos 11, 12, 22 (i=1..5, j=1..3)
*        BRHSQ(i,j)   hi      -> uLuL, uRuR, dLdL, dRdR, t1t1, t2t2,
*                                t1t2, b1b1, b2b2, b1b2 (i=1..3, j=1..10)
*        BRASQ(i,j)   ai      -> t1t2, b1b2 (i=1,2, j=1,2)
*        BRHSL(i,j)   hi      -> lLlL, lRlR, nLnL, l1l1, l2l2, l1l2,
*                                ntnt (i=1..3, j=1..7)
*        BRASL(i)     ai      -> l1l2 (i=1,2)
*
*      HCWIDTH  Total decay width of the charged Higgs
*               with the following branching ratios:
*      HCBRM         h+ -> mu nu_mu
*      HCBRL         "  -> tau nu_tau
*      HCBRSU        "  -> s u
*      HCBRBU        "  -> b u
*      HCBRSC        "  -> s c
*      HCBRBC        "  -> b c
*      HCBRBT        "  -> b t
*      HCBRWHT       "  -> neutral Higgs W+, including:
*        HCBRWH(i)   "  -> H1W+, H2W+, h3W+, a1W+, a2W+ (i=1..5)
*      HCBRSUSY      "  -> susy particles,including
*        HCBRNC(i,j) "  -> neutralino i chargino j (i=1..5, j=1,2)
*        HCBRSQ(i)   "  -> uLdL, t1b1, t1b2, t2b1, t2b2 (i=1..5)
*        HCBRSL(i)   "  -> lLnL, t1nt, t2nt (i=1..3)
*
*      MNEU(i)   Mass of neutralino chi_i (i=1,5, ordered in mass)
*      NEU(i,j)  chi_i components of bino, wino, higgsino u&d, singlino 
*                (i,j=1..5)
*
*      MCHA(i)       Chargino masses
*      U(i,j),V(i,j) Chargino mixing matrices
*
*  Significances for Higgs detection at the LHC:
*
*    At low luminosity (30 fb^-1): in LOWSIG(X,Y), where
*                               
*      X=1: h1
*      X=2: h2
*      X=3: h3
*      X=4: a1
*      X=5: a2
*
*      Y=1: channel bbh/A -> bbtautau
*      Y=2: channel gg -> h/a -> gamma gamma
*      Y=3: channel gg -> h -> ZZ -> 4 leptons
*      Y=4: channel gg -> h -> WW -> 2 leptons 2 neutrinos
*      Y=5: channel WW -> h -> tautau
*      Y=6: channel WW -> h -> WW
*      Y=7: channel WW -> h -> gamma gamma
*
*    At high luminosity (300 fb^-1): in HIGSIG(X,Y), where X as above,
*
*      Y=1: channel h/a -> gamma gamma
*      Y=2: channel h/a -> gamma gamma lepton
*      Y=3: channel tth/a -> bb + X
*      Y=4: channel bbh/a -> bbtautau
*      Y=5: channel gg -> h -> ZZ -> 4 leptons
*      Y=6: channel gg -> h -> WW -> 2 leptons 2 neutrinos
*      Y=7: channel WW -> h -> tautau
*      Y=8: channel WW -> h -> WW
*      Y=9: channel WW -> h -> invisible
*
*  ERRORS: IFAIL = 0..12
*
*  IFAIL = 0         OK
*          1,3,5,7   m_h1^2 < 0
*          2,3,6,7   m_a1^2 < 0
*          4,5,6,7   m_h+^2 < 0
*          8         m_sfermion^2 < 0
*          9         l, k, tan(beta) or mu = 0
*          10        Violation of phenomenological constraint(s)
*          11,12     Problem in integration of RGEs
*
*  Phenomenological constraints:
*
*      PROB(I)  = 0, I = 1..45: OK
*            
*      PROB(1) =/= 0   chargino too light
*      PROB(2) =/= 0   excluded by Z -> neutralinos
*      PROB(3) =/= 0   charged Higgs too light
*      PROB(4) =/= 0   excluded by ee -> hZ 
*      PROB(5) =/= 0   excluded by ee -> hZ, h -> bb
*      PROB(6) =/= 0   excluded by ee -> hZ, h -> tautau
*      PROB(7) =/= 0   excluded by ee -> hZ, h -> invisible 
*      PROB(8) =/= 0   excluded by ee -> hZ, h -> 2jets
*      PROB(9) =/= 0   excluded by ee -> hZ, h -> 2photons
*      PROB(10) =/= 0  excluded by ee -> hZ, h -> AA -> 4bs
*      PROB(11) =/= 0  excluded by ee -> hZ, h -> AA -> 4taus
*      PROB(12) =/= 0  excluded by ee -> hZ, h -> AA -> 2bs 2taus
*      PROB(13) =/= 0  excluded by Z -> hA (Z width)
*      PROB(14) =/= 0  excluded by ee -> hA -> 4bs
*      PROB(15) =/= 0  excluded by ee -> hA -> 4taus
*      PROB(16) =/= 0  excluded by ee -> hA -> 2bs 2taus
*      PROB(17) =/= 0  excluded by ee -> hA -> AAA -> 6bs
*      PROB(18) =/= 0  excluded by ee -> hA -> AAA -> 6taus
*      PROB(19) =/= 0  excluded by ee -> Zh -> ZAA -> Z + light pairs
*      PROB(20) =/= 0  excluded by stop -> b l sneutrino
*      PROB(21) =/= 0  excluded by stop -> neutralino c
*      PROB(22) =/= 0  excluded by sbottom -> neutralino b
*      PROB(23) =/= 0  squark/gluino too light
*      PROB(24) =/= 0  selectron/smuon too light
*      PROB(25) =/= 0  stau too light
*      PROB(26) =/= 0  lightest neutralino is not LSP
*      PROB(27) =/= 0  Landau Pole in l, k, ht, hb below MGUT
*      PROB(28) =/= 0  unphysical global minimum
*      PROB(29) =/= 0  Higgs soft masses >> Msusy
*      PROB(30) =/= 0  excluded by WMAP (checked only if OMGFLAG=1)
*      PROB(31) =/= 0  eff. Higgs self-couplings in Micromegas > 1
*      PROB(32) =/= 0  b->s gamma more than 2 sigma away
*      PROB(33) =/= 0  Delta M_s more than 2 sigma away
*      PROB(34) =/= 0  Delta M_d more than 2 sigma away
*      PROB(35) =/= 0  B_s->mu+mu- more than 2 sigma away
*      PROB(36) =/= 0  B+-> tau+nu_tau more than 2 sigma away
*      PROB(37) =/= 0  (g-2)_muon more than 2 sigma away
*      PROB(38) =/= 0  excluded by Upsilon(1S) -> A gamma
*      PROB(39) =/= 0  excluded by eta_b(1S) mass measurement
*      PROB(40) =/= 0  BR(B-->X_s mu+ mu-) more than 2 sigma away
*      PROB(41) =/= 0  excluded by ee -> hZ, h -> AA -> 4taus (new ALEPH analysis)
*      PROB(42) =/= 0  excluded by top -> b H+, H+ -> c s (CDF, D0)
*      PROB(43) =/= 0  excluded by top -> b H+, H+ -> tau nu_tau (D0)
*      PROB(44) =/= 0  excluded by top -> b H+, H+ -> W+ A1, A1 -> 2taus (CDF)
*      PROB(45) =/= 0  excluded by LHC: A/H -> 2taus
*
************************************************************************

      IMPLICIT NONE

      INTEGER NFL,NPROB,NPAR
      PARAMETER (NFL=12,NPROB=45,NPAR=25)
      INTEGER NFAIL(NFL),IFAIL,I,IDIM
      INTEGER OMGFLAG,MAFLAG,M1FLAG,M3FLAG,TOT,ITOT,NTOT,NMSFLAG,NMSCAN
      INTEGER NL,NK,NTB,NMU,N1,N3,NM1,NM2,NM3
      INTEGER IL,IK,ITB,IMU,I1,I3,IM1,IM2,IM3

      DOUBLE PRECISION PAR(NPAR),PROB(NPROB)
      DOUBLE PRECISION LMIN,LMAX,KMIN,KMAX,TBMIN,TBMAX,MUMIN,MUMAX
      DOUBLE PRECISION ALMIN,ALMAX,AKMIN,AKMAX,MAMIN,MAMAX
      DOUBLE PRECISION M1MIN,M1MAX,M2MIN,M2MAX,M3MIN,M3MAX
      DOUBLE PRECISION LN,LNN,KN,KNN,TBN,TBNN,MUN,MUNN
      DOUBLE PRECISION ALN,ALNN,AKN,AKNN,MAN,MANN
      DOUBLE PRECISION M1N,M1NN,M2N,M2NN,M3N,M3NN
      DOUBLE PRECISION XIF,XIS,MUP,MSP,M3H

      COMMON/BOUNDS/LN,LNN,KN,KNN,TBN,TBNN,MUN,MUNN,
     . ALN,ALNN,AKN,AKNN,MAN,MANN,
     . M1N,M1NN,M2N,M2NN,M3N,M3NN
      COMMON/MINMAX/LMIN,LMAX,KMIN,KMAX,TBMIN,TBMAX,MUMIN,MUMAX,
     . ALMIN,ALMAX,AKMIN,AKMAX,MAMIN,MAMAX,
     . M1MIN,M1MAX,M2MIN,M2MAX,M3MIN,M3MAX
      COMMON/STEPS/NTOT,NL,NK,NTB,NMU,N1,N3,NM1,NM2,NM3
      COMMON/FLAGS/OMGFLAG,MAFLAG
      COMMON/SCANFLAGS/M1FLAG,M3FLAG
      COMMON/SUSYEXT/XIF,XIS,MUP,MSP,M3H
      COMMON/IDIM/IDIM
      COMMON/NMSFLAG/NMSFLAG,NMSCAN

*   Initialization

      CALL INITIALIZE()
      DO I=1,NFL
       NFAIL(I)=0
      ENDDO
      TOT=0
      ITOT=0
      IDIM=0

*   Reading of the input parameters

      CALL INPUT(PAR,NPAR)

*   Initialization of the range of parameters that has passed all tests

      TBN=MAX(TBMAX,TBMIN)
      TBNN=MIN(TBMAX,TBMIN)
      IF(M1FLAG.EQ.0)THEN
       M1N=MAX(M2MAX,M2MIN)/2d0
       M1NN=MIN(M2MAX,M2MIN)/2d0
      ELSE
       M1N=MAX(M1MAX,M1MIN)
       M1NN=MIN(M1MAX,M1MIN)
      ENDIF
      M2N=MAX(M2MAX,M2MIN)
      M2NN=MIN(M2MAX,M2MIN)
      IF(M3FLAG.EQ.0)THEN
       M3N=MAX(M2MAX,M2MIN)*3d0
       M3NN=MIN(M2MAX,M2MIN)*3d0
      ELSE
       M3N=MAX(M3MAX,M3MIN)
       M3NN=MIN(M3MAX,M3MIN)
      ENDIF
      LN=MAX(LMAX,LMIN)
      LNN=MIN(LMAX,LMIN)
      KN=MAX(KMAX,KMIN)
      KNN=MIN(KMAX,KMIN)
      MUN=MAX(MUMAX,MUMIN)
      MUNN=MIN(MUMAX,MUMIN)
      IF(MAFLAG.EQ.0)THEN
       ALN=MAX(ALMAX,ALMIN)
       ALNN=MIN(ALMAX,ALMIN)
       MAN=1d99
       MANN=0d0
      ELSE
       MAN=MAX(MAMAX,MAMIN)
       MANN=MIN(MAMAX,MAMIN)
       ALN=1d99
       ALNN=-1d99
      ENDIF
      AKN=MAX(AKMAX,AKMIN)
      AKNN=MIN(AKMAX,AKMIN)

*   Beginning of the scan

      DO ITB=1,NTB
      IF(NTB.EQ.1)THEN
       PAR(3)=TBMIN
      ELSE
       PAR(3)=TBMIN+(TBMAX-TBMIN)*DFLOAT(ITB-1)/DFLOAT(NTB-1)
      ENDIF

      DO IM2=1,NM2
      IF(NM2.EQ.1)THEN
       PAR(21)=M2MIN
      ELSE
       PAR(21)=M2MIN+(M2MAX-M2MIN)*DFLOAT(IM2-1)/DFLOAT(NM2-1)
      ENDIF

      DO IM1=1,NM1
      IF(M1FLAG.EQ.0)THEN
       PAR(20)=PAR(21)/2d0
      ELSEIF(NM1.EQ.1)THEN
       PAR(20)=M1MIN
      ELSE
       PAR(20)=M1MIN+(M1MAX-M1MIN)*DFLOAT(IM1-1)/DFLOAT(NM1-1)
      ENDIF

      DO IM3=1,NM3
      IF(M3FLAG.EQ.0)THEN
       PAR(22)=PAR(21)*3d0
      ELSEIF(NM3.EQ.1)THEN
       PAR(22)=M3MIN
      ELSE
       PAR(22)=M3MIN+(M3MAX-M3MIN)*DFLOAT(IM3-1)/DFLOAT(NM3-1)
      ENDIF

      DO IL=1,NL
      IF(NL.EQ.1)THEN
       PAR(1)=LMIN
      ELSE
       PAR(1)=LMIN+(LMAX-LMIN)*DFLOAT(IL-1)/DFLOAT(NL-1)
      ENDIF

      DO IK=1,NK
      IF(NK.EQ.1)THEN
       PAR(2)=KMIN
      ELSE
       PAR(2)=KMIN+(KMAX-KMIN)*DFLOAT(IK-1)/DFLOAT(NK-1)
      ENDIF

      DO IMU=1,NMU
      IF(NMU.EQ.1)THEN
       PAR(4)=MUMIN
      ELSE
       PAR(4)=MUMIN+(MUMAX-MUMIN)*DFLOAT(IMU-1)/DFLOAT(NMU-1)
      ENDIF

      DO I1=1,N1
      IF(N1.EQ.1)THEN
       IF(MAFLAG.EQ.0)THEN
        PAR(5)=ALMIN
       ELSE
        PAR(23)=MAMIN
       ENDIF
      ELSE
       IF(MAFLAG.EQ.0)THEN
        PAR(5)=ALMIN+(ALMAX-ALMIN)*DFLOAT(I1-1)/DFLOAT(N1-1)
       ELSE
        PAR(23)=MAMIN+(MAMAX-MAMIN)*DFLOAT(I1-1)/DFLOAT(N1-1)
       ENDIF
      ENDIF
      
      DO I3=1,N3
      IF(N3.EQ.1)THEN
       PAR(6)=AKMIN
      ELSE
       PAR(6)=AKMIN+(AKMAX-AKMIN)*DFLOAT(I3-1)/DFLOAT(N3-1)
      ENDIF

      XIF=0d0
      XIS=0d0
      MUP=0d0
      MSP=0d0
      M3H=0d0

      ITOT=ITOT+1

*   Initialization of PROB and IFAIL

      DO I=1,NPROB
       PROB(I)=0d0
      ENDDO
      IFAIL=0

*   Check for singular parameters l, k, tan(beta) and mu_eff

      IF(PAR(1)*PAR(2)*PAR(3)*PAR(4).EQ.0d0)THEN
       IFAIL=9
       GOTO 11
      ENDIF
      
*   Computation of parameters at QSTSB
      
      CALL RUNPAR(PAR)
      
*   Computation of sfermion masses

      CALL MSFERM(PAR,IFAIL)
      IF(IFAIL.NE.0)GOTO 11
      
*   Computation of Higgs masses

      CALL MHIGGS(PAR,IFAIL)
      IF(IFAIL.NE.0)GOTO 11

*   Computation of gluino mass

      CALL GLUINO(PAR)

*   Computation of chargino masses

      CALL CHARGINO(PAR)

*   Computation of neutralino masses

      CALL NEUTRALINO(PAR)

*   Computation of Higgs + top branching ratios

      CALL DECAY(PAR)
      CALL TDECAY(PAR)

*   Exp. constraints on sparticles/Higgses

      CALL SUBEXP(PAR,PROB)

*   b -> s gamma + B physics

      CALL BSG(PAR,PROB)

*   Anom. magn. moment of the Muon

      CALL MAGNMU(PAR,PROB)

*   Global minimum?

      CALL CHECKMIN(PROB)
      
*   Landau Pole?

      CALL RGES(PAR,PROB,IFAIL)

*   RGEs for the soft terms

      CALL RGESOFT(PAR,IFAIL)

*   Relic density

      CALL RELDEN(PAR,PROB)

*   Computation of the statistical significances at LHC

c      CALL LSIG()
c      CALL HSIG()

*   Phenomenological problem

      DO I=1,NPROB
       IF(PROB(I).NE.0d0.AND.I.NE.31)IFAIL=10
      ENDDO
      
*   Sparticle decays

      IF(NMSFLAG.EQ.1) THEN
        CALL NMSDECAY_INTERFACE(PAR)
      ENDIF

*   Recording of the results

 11   IF(IFAIL.EQ.0)THEN
       TOT=TOT+1
       LN=MIN(PAR(1),LN)
       LNN=MAX(PAR(1),LNN)
       KN=MIN(PAR(2),KN)
       KNN=MAX(PAR(2),KNN)
       TBN=MIN(PAR(3),TBN)
       TBNN=MAX(PAR(3),TBNN)
       MUN=MIN(PAR(4),MUN)
       MUNN=MAX(PAR(4),MUNN)
       ALN=MIN(PAR(5),ALN)
       ALNN=MAX(PAR(5),ALNN)
       AKN=MIN(PAR(6),AKN)
       AKNN=MAX(PAR(6),AKNN)
       M1N=MIN(PAR(20),M1N)
       M1NN=MAX(PAR(20),M1NN)
       M2N=MIN(PAR(21),M2N)
       M2NN=MAX(PAR(21),M2NN)
       M3N=MIN(PAR(22),M3N)
       M3NN=MAX(PAR(22),M3NN)
       MAN=MIN(PAR(23),MAN)
       MANN=MAX(PAR(23),MANN)
      ELSE
       NFAIL(IFAIL)=NFAIL(IFAIL)+1
      ENDIF

      CALL OUTPUT(PAR,PROB,IFAIL,ITOT,NTOT)

      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO

*   Summary of the scanning:
*   Number of points that passed/failed the tests
*   and range for scanned parameters

      CALL ERROR(TOT,NTOT,NFAIL)

      END


      SUBROUTINE INPUT(PAR,NPAR)
      
*******************************************************************
*   This subroutine reads SM and NMSSM parameters from input file   .
*******************************************************************

      IMPLICIT NONE

      CHARACTER CHINL*120,CHBLCK*60,CHDUM*120

      INTEGER I,NLINE,INL,ICH,IX,IVAL,Q2FIX
      INTEGER N0,NLOOP,NBER,NPAR,ERR
      INTEGER OMGFLAG,MAFLAG,M1FLAG,M3FLAG,PFLAG,NMSFLAG,NMSCAN
      INTEGER NTOT,NL,NK,NTB,NMU,N1,N3,NM1,NM2,NM3
      INTEGER NAL,NMA,NAK

      DOUBLE PRECISION PAR(*),VAL
      DOUBLE PRECISION ACC,XITLA,XLAMBDA,MC0,MB0,MT0
      DOUBLE PRECISION ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION VUS,VCB,VUB,Q2,Q2MIN
      DOUBLE PRECISION LMIN,LMAX,KMIN,KMAX,TBMIN,TBMAX,MUMIN,MUMAX
      DOUBLE PRECISION ALMIN,ALMAX,AKMIN,AKMAX,MAMIN,MAMAX
      DOUBLE PRECISION M1MIN,M1MAX,M2MIN,M2MAX,M3MIN,M3MAX

      COMMON/GAUGE/ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/CKM/VUS,VCB,VUB
      COMMON/ALS/XLAMBDA,MC0,MB0,MT0,N0
      COMMON/RENSCALE/Q2
      COMMON/Q2FIX/Q2MIN,Q2FIX
      COMMON/MINMAX/LMIN,LMAX,KMIN,KMAX,TBMIN,TBMAX,MUMIN,MUMAX,
     . ALMIN,ALMAX,AKMIN,AKMAX,MAMIN,MAMAX,
     . M1MIN,M1MAX,M2MIN,M2MAX,M3MIN,M3MAX
      COMMON/STEPS/NTOT,NL,NK,NTB,NMU,N1,N3,NM1,NM2,NM3
      COMMON/SCANFLAGS/M1FLAG,M3FLAG
      COMMON/FLAGS/OMGFLAG,MAFLAG
      COMMON/PFLAG/PFLAG
      COMMON/NMSFLAG/NMSFLAG,NMSCAN

*   INITIALIZATION OF THE SUSY PARAMETERS
      DO I=1,NPAR
       PAR(I)=1d99
      ENDDO

*   INITIALIZATION OF THE SCANNING PARAMETERS
      TBMIN=1d99
      TBMAX=1d99
      M1MIN=1d99
      M1MAX=1d99
      M2MIN=1d99
      M2MAX=1d99
      M3MIN=1d99
      M3MAX=1d99
      LMIN=1d99
      LMAX=1d99
      KMIN=1d99
      KMAX=1d99
      ALMIN=1d99
      ALMAX=1d99
      AKMIN=1d99
      AKMAX=1d99
      MUMIN=1d99
      MUMAX=1d99
      MAMIN=1d99
      MAMAX=1d99
      NM1=0
      NM2=1
      NM3=0
      NL=1
      NK=1
      NTB=1
      NMU=1
      NAL=0
      NMA=0
      NAK=1

*   DEFAULT VALUE FOR OMGFLAG
      OMGFLAG=0

*   DEFAULT VALUE FOR PFLAG
      PFLAG=0

*   DEFAULT VALUE FOR NMSFLAG
      NMSFLAG=0
      NMSCAN=1

*   DEFAULT VALUE FOR SCANFLAGS
      M1FLAG=0
      M3FLAG=0

*   DEFAULT VALUE FOR THE RENSCALE Q2
      Q2=0d0

*   INITIALIZE READ LOOP
      NLINE=0
      CHBLCK=' '

*   START TO READ NEW LINE INTO CHINL
 21   CHINL=' '

*   LINE NUMBER
      NLINE=NLINE+1

      READ(5,'(A120)',END=29,ERR=999) CHINL
      
*   CHECK FOR EMPTY OR COMMENT LINES
      IF(CHINL.EQ.' '.OR.CHINL(1:1).EQ.'#'
     .  .OR.CHINL(1:1).EQ.'*') GOTO 21

*   FORCE UPPER CASE LETTERS IN CHINL (AS REQUIRED BELOW)
      INL=0
 22   INL=INL+1
      IF(CHINL(INL:INL).NE.'#')THEN
       DO ICH=97,122
        IF(CHINL(INL:INL).EQ.CHAR(ICH)) CHINL(INL:INL)=CHAR(ICH-32)
       ENDDO
       IF(INL.LT.120) GOTO 22
      ENDIF

*   CHECK FOR BLOCK STATEMENT
      IF(CHINL(1:1).EQ.'B')THEN
       READ(CHINL,'(A6,A)',ERR=999) CHDUM,CHBLCK
       GOTO 21
      ENDIF

*   CHECK FOR NMSSM MODEL IN MODSEL
*   IF THE RELIC DENSITY SHOULD BE COMPUTED
*   THE BLOCK MODSEL_min MUST CONTAIN THE LINE "  9     1    "
      IF(CHBLCK(1:6).EQ.'MODSEL')THEN
       READ(CHINL,*,ERR=999) IX,IVAL
       IF(IX.EQ.8) PFLAG=IVAL
       IF(IX.EQ.9) OMGFLAG=IVAL
       IF(IX.EQ.13) NMSFLAG=IVAL

*   READ SMINPUTS
      ELSEIF(CHBLCK(1:8).EQ.'SMINPUTS')THEN
       READ(CHINL,*,ERR=999) IX,VAL
       IF(IX.EQ.2) GF=VAL
       IF(IX.EQ.3) ALSMZ=VAL
       IF(IX.EQ.4) MZ=VAL
       IF(IX.EQ.5) MB=VAL
       IF(IX.EQ.6) MT=VAL
       IF(IX.EQ.7) MTAU=VAL

*   READ Q2 AND TANBETA
      ELSEIF(CHBLCK(1:6).EQ.'MINPAR')THEN
       READ(CHINL,*,ERR=999) IX,VAL
       IF(IX.EQ.0) Q2=VAL**2
       IF(IX.EQ.37) TBMIN=VAL
       IF(IX.EQ.38) TBMAX=VAL
 
*   READ EXTPAR
      ELSEIF(CHBLCK(1:6).EQ.'EXTPAR')THEN
       READ(CHINL,*,ERR=999) IX,VAL
       IF(IX.EQ.107) M1MIN=VAL
       IF(IX.EQ.108) M1MAX=VAL
       IF(IX.EQ.207) M2MIN=VAL
       IF(IX.EQ.208) M2MAX=VAL
       IF(IX.EQ.307) M3MIN=VAL
       IF(IX.EQ.308) M3MAX=VAL
       IF(IX.EQ.11) PAR(12)=VAL
       IF(IX.EQ.12) PAR(13)=VAL
       IF(IX.EQ.13) PAR(14)=VAL
       IF(IX.EQ.16) PAR(25)=VAL
       IF(IX.EQ.32) PAR(18)=VAL**2
       IF(IX.EQ.33) PAR(10)=VAL**2
       IF(IX.EQ.35) PAR(19)=VAL**2
       IF(IX.EQ.36) PAR(11)=VAL**2
       IF(IX.EQ.42) PAR(15)=VAL**2
       IF(IX.EQ.43) PAR(7)=VAL**2
       IF(IX.EQ.45) PAR(16)=VAL**2
       IF(IX.EQ.46) PAR(8)=VAL**2
       IF(IX.EQ.48) PAR(17)=VAL**2
       IF(IX.EQ.49) PAR(9)=VAL**2
       IF(IX.EQ.617) LMIN=VAL
       IF(IX.EQ.618) LMAX=VAL
       IF(IX.EQ.627) KMIN=VAL
       IF(IX.EQ.628) KMAX=VAL
       IF(IX.EQ.637) ALMIN=VAL
       IF(IX.EQ.638) ALMAX=VAL
       IF(IX.EQ.647) AKMIN=VAL
       IF(IX.EQ.648) AKMAX=VAL
       IF(IX.EQ.657) MUMIN=VAL
       IF(IX.EQ.658) MUMAX=VAL
       IF(IX.EQ.1247) MAMIN=VAL
       IF(IX.EQ.1248) MAMAX=VAL

*   READ STEPS
       ELSEIF(CHBLCK(1:6).EQ.'STEPS')THEN
       READ(CHINL,*,ERR=999) IX,IVAL
       IF(IX.EQ.109) NM1=IVAL
       IF(IX.EQ.209) NM2=IVAL
       IF(IX.EQ.309) NM3=IVAL
       IF(IX.EQ.39) NTB=IVAL
       IF(IX.EQ.619) NL=IVAL
       IF(IX.EQ.629) NK=IVAL
       IF(IX.EQ.639) NAL=IVAL
       IF(IX.EQ.649) NAK=IVAL
       IF(IX.EQ.659) NMU=IVAL
       IF(IX.EQ.1249) NMA=IVAL

      ENDIF

      GOTO 21

*   END OF READING FROM INPUT FILE

*   Check for errors

 29   ERR=0
      IF(TBMIN.EQ.1d99)THEN
       WRITE(0,1)"TANB_min MUST BE GIVEN IN BLOCK MINPAR"
       ERR=1
      ENDIF
      IF(M1MIN.EQ.1d99 .AND. M1MAX.NE.1d99)THEN
       WRITE(0,1)"M1_min MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(M2MIN.EQ.1d99)THEN
       WRITE(0,1)"M2_min MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(M3MIN.EQ.1d99 .AND. M3MAX.NE.1d99)THEN
       WRITE(0,1)"M3_min MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(PAR(12).EQ.1d99)THEN
       WRITE(0,1)"AU3 MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(PAR(13).EQ.1d99)THEN
       WRITE(0,1)"AD3 MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(PAR(14).EQ.1d99)THEN
       WRITE(0,1)"AE3 MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(PAR(7).EQ.1d99)THEN
       WRITE(0,1)"MQ3 MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(PAR(8).EQ.1d99)THEN
       WRITE(0,1)"MU3 MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(PAR(9).EQ.1d99)THEN
       WRITE(0,1)"MD3 MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(PAR(10).EQ.1d99)THEN
       WRITE(0,1)"ML3 MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(PAR(11).EQ.1d99)THEN
       WRITE(0,1)"ME3 MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(LMIN.EQ.1d99)THEN
       WRITE(0,1)"LAMBDA_min MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(KMIN.EQ.1d99)THEN
       WRITE(0,1)"KAPPA_min MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(MUMIN.EQ.1d99)THEN
       WRITE(0,1)"MUEFF_min MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(MAMIN.EQ.1d99 .AND. MAMAX.NE.1d99)THEN
       WRITE(0,1)"MA_min MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(ALMIN.EQ.1d99 .AND. ALMAX.NE.1d99)THEN
       WRITE(0,1)"AL_min MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(ALMIN.NE.1d99 .AND. MAMIN.NE.1d99)THEN
       WRITE(0,1)"BOTH ALAMBDA AND MA CANNOT BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(ALMIN.EQ.1d99 .AND. MAMIN.EQ.1d99)THEN
       WRITE(0,1)"EITHER ALAMBDA OR MA MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(AKMIN.EQ.1d99)THEN
       WRITE(0,1)"AK_min MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF

*   Set default values

      IF(TBMAX.EQ.1d99)TBMAX=TBMIN
      IF(M1MAX.EQ.1d99)M1MAX=M1MIN
      IF(M2MAX.EQ.1d99)M2MAX=M2MIN
      IF(M3MAX.EQ.1d99)M3MAX=M3MIN
      IF(LMAX.EQ.1d99)LMAX=LMIN
      IF(KMAX.EQ.1d99)KMAX=KMIN
      IF(ALMAX.EQ.1d99)ALMAX=ALMIN
      IF(AKMAX.EQ.1d99)AKMAX=AKMIN
      IF(MUMAX.EQ.1d99)MUMAX=MUMIN
      IF(MAMAX.EQ.1d99)MAMAX=MAMIN
      
      DO I=15,19
       IF(PAR(I).EQ.1d99)PAR(I)=PAR(I-8)
      ENDDO
      IF(PAR(25).EQ.1d99)PAR(25)=PAR(14)

*   Set MAFLAG, SCANFLAGS

      IF(MAMIN.EQ.1d99)MAFLAG=0
      IF(ALMIN.EQ.1d99)MAFLAG=1
      IF(M1MIN.NE.1d99)M1FLAG=1
      IF(M3MIN.NE.1d99)M3FLAG=1

*   Number of points

      IF(MAFLAG.EQ.0)THEN
       IF(NMA.NE.0)THEN
        WRITE(0,1)"NMA CANNOT BE GIVEN IN BLOCK STEPS"
       ERR=1
       ENDIF
       IF(NAL.EQ.0)NAL=1
       N1=NAL
      ENDIF

      IF(MAFLAG.EQ.1)THEN
       IF(NAL.NE.0)THEN
        WRITE(0,1)"NAL CANNOT BE GIVEN IN BLOCK STEPS"
        ERR=1
       ENDIF
       IF(NMA.EQ.0)NMA=1
       N1=NMA
      ENDIF

      N3=NAK

      IF(M1FLAG.EQ.0)THEN
       IF(NM1.NE.0)THEN
        WRITE(0,1)"NM1 CANNOT BE GIVEN IN BLOCK STEPS"
        ERR=1
       ELSE
        NM1=1
       ENDIF
      ELSE
       IF(NM1.EQ.0)NM1=1
      ENDIF

      IF(M3FLAG.EQ.0)THEN
       IF(NM3.NE.0)THEN
        WRITE(0,1)"NM3 CANNOT BE GIVEN IN BLOCK STEPS"
        ERR=1
       ELSE
        NM3=1
       ENDIF
      ELSE
       IF(NM3.EQ.0)NM3=1
      ENDIF

      IF(NM1.LE.0 .OR. NM2.LE.0 .OR. NM3.LE.0 .OR. NTB.LE.0 .OR.
     . NL.LE.0 .OR. NK.LE.0 .OR. N1.LE.0 .OR. N3.LE.0 .OR. NMU.LE.0)THEN
       WRITE(0,1)"WRONG NUMBER OF POINTS IN BLOCK STEPS"
       ERR=1
      ENDIF

*   Stop if error

      IF(ERR.EQ.1)THEN
       WRITE(0,1)"ERROR IN INPUT FILE"
       STOP 1
      ENDIF

*   Warnings

      IF(TBMIN.EQ.TBMAX .AND. NTB.GT.1)THEN
       WRITE(0,1)"WARNING TANB_min = TANB_max => NTB=1"
       NTB=1
      ENDIF
      IF(M1MIN.EQ.M1MAX .AND. NM1.GT.1)THEN
       WRITE(0,1)"WARNING M1_min = M1_max => NM1=1"
       NM1=1
      ENDIF
      IF(M2MIN.EQ.M2MAX .AND. NM2.GT.1)THEN
       WRITE(0,1)"WARNING M2_min = M2_max => NM2=1"
       NM2=1
      ENDIF
      IF(M3MIN.EQ.M3MAX .AND. NM3.GT.1)THEN
       WRITE(0,1)"WARNING M3_min = M3_max => NM3=1"
       NM3=1
      ENDIF
      IF(LMIN.EQ.LMAX .AND. NL.GT.1)THEN
       WRITE(0,1)"WARNING LAMBDA_min = LAMBDA_max => NL=1"
       NL=1
      ENDIF
      IF(KMIN.EQ.KMAX .AND. NK.GT.1)THEN
       WRITE(0,1)"WARNING KAPPA_min = KAPPA_max => NK=1"
       NK=1
      ENDIF
      IF(MUMIN.EQ.MUMAX .AND. NMU.GT.1)THEN
       WRITE(0,1)"WARNING MUEFF_min = MUEFF_max => NMU=1"
       NMU=1
      ENDIF

      IF(TBMIN.NE.TBMAX .AND. NTB.EQ.1)THEN
       WRITE(0,10)"WARNING NTB=1 => TANB_max=TANB_min=",TBMIN
       TBMAX=TBMIN
      ENDIF
      IF(M1MIN.NE.M1MAX .AND. NM1.EQ.1)THEN
       WRITE(0,10)"WARNING NM1=1 => M1_max=M1_min=",M1MIN
       M1MAX=M1MIN
      ENDIF
      IF(M2MIN.NE.M2MAX .AND. NM2.EQ.1)THEN
       WRITE(0,10)"WARNING NM2=1 => M2_max=M2_min=",M2MIN
       M2MAX=M2MIN
      ENDIF
      IF(M3MIN.NE.M3MAX .AND. NM3.EQ.1)THEN
       WRITE(0,10)"WARNING NM3=1 => M3_max=M3_min=",M3MIN
       M3MAX=M3MIN
      ENDIF
      IF(LMIN.NE.LMAX .AND. NL.EQ.1)THEN
       WRITE(0,10)"WARNING NL=1 => LAMBDA_max=LAMBDA_min=",LMIN
       LMAX=LMIN
      ENDIF
      IF(KMIN.NE.KMAX .AND. NK.EQ.1)THEN
       WRITE(0,10)"WARNING NK=1 => KAPPA_max=KAPPA_min=",KMIN
       KMAX=KMIN
      ENDIF
      IF(MUMIN.NE.MUMAX .AND. NMU.EQ.1)THEN
       WRITE(0,10)"WARNING NMU=1 => MUEFF_max=MUEFF_min=",MUMIN
       MUMAX=MUMIN
      ENDIF

      IF(MAFLAG.EQ.0)THEN
       IF(ALMIN.EQ.ALMAX .AND. NAL.GT.1)THEN
        WRITE(0,1)"WARNING ALAMBDA_min=ALAMBDA_max => NAL=1"
        NAL=1
        N1=NAL
       ENDIF
       IF(ALMIN.NE.ALMAX .AND. NAL.EQ.1)THEN
        WRITE(0,10)"WARNING NAL=1 => ALAMBDA_max=ALAMBDA_min=",ALMIN
        ALMAX=ALMIN
       ENDIF
      ENDIF

      IF(MAFLAG.EQ.1)THEN
       IF(MAMIN.EQ.MAMAX .AND. NMA.GT.1)THEN
        WRITE(0,1)"WARNING MA_min=MA_max => NMA=1"
        NMA=1
        N1=NMA
       ENDIF
       IF(MAMIN.NE.MAMAX .AND. NMA.EQ.1)THEN
        WRITE(0,10)"WARNING NMA=1 => MA_max=MA_min=",MAMIN
        MAMAX=MAMIN
       ENDIF
      ENDIF

      IF(AKMIN.EQ.AKMAX .AND. NAK.GT.1)THEN
       WRITE(0,1)"WARNING AKAPPA_min=AKAPPA_max => NAK=1"
       NAK=1
       N3=NAK
      ENDIF
      IF(AKMIN.NE.AKMAX .AND. NAK.EQ.1)THEN
       WRITE(0,10)"WARNING NAK=1 => AKAPPA_max=AKAPPA_min=",AKMIN
       AKMAX=AKMIN
      ENDIF

*   Total number of points

      NTOT=NTB*NM1*NM2*NM3*NL*NMU*NK*N1*N3

*   Set Q2MIN, Q2FIX:
      Q2MIN=100d0**2
      Q2FIX=1
      IF(Q2.LE.Q2MIN)THEN
       Q2FIX=0
      ENDIF

*   Initialization for ALPHAS and RUNM (as in hdecay)
*   The bottom quark pole mass MBP is set in INIT and can be changed
*   only there (changing its running mass MB above has no effect
*   on MBP, since one would have to compute alpha_s(MB) first)

      MC0=MC
      MB0=MBP
      MT0=MT
      N0=5
      NLOOP=2
      NBER=18
      ACC=1.D-8
      XLAMBDA=XITLA(NLOOP,ALSMZ,ACC)
      CALL ALSINI(ACC)
      CALL BERNINI(NBER)

*    g1,g2  and sin(theta)^2 in the on-shell scheme in terms of 
*    GF, MZ(pole) and MW(pole)

      g2=4d0*DSQRT(2d0)*GF*MW**2
      g1=4d0*DSQRT(2d0)*GF*(MZ**2-MW**2)
      S2TW=1d0-(MW/MZ)**2

      RETURN

 999  WRITE(0,1)"READ ERROR ON LINE:", NLINE
      WRITE(0,*)CHINL(1:80)
      STOP 1

 1    FORMAT(2A)
 10   FORMAT(A,E10.3)

      END


      SUBROUTINE OUTPUT(PAR,PROB,IFAIL,ITOT,NTOT) 

*********************************************************************      
*   Subroutine writing all the results in the the output file.
*********************************************************************      
 
      IMPLICIT NONE 

      INTEGER NBIN,IFAIL,ITOT,NTOT,IDIM,I,J,IMAX,JMAX
      PARAMETER (IMAX=1000,JMAX=100)

      DOUBLE PRECISION RES(IMAX,JMAX),PAR(*),PROB(*),S
      DOUBLE PRECISION SMASS(3),SCOMP(3,3),PMASS(2),PCOMP(2,2),CMASS
      DOUBLE PRECISION MGL,MCHA(2),U(2,2),V(2,2),MNEU(5),NEU(5,5)
      DOUBLE PRECISION BRJJ(5),BRMM(5),BRLL(5),BRSS(5),BRCC(5)
      DOUBLE PRECISION BRBB(5),BRTT(5),BRWW(3),BRZZ(3),BRGG(5)
      DOUBLE PRECISION BRZG(5),BRHHH(4),BRHAA(3,3),BRHCHC(3)
      DOUBLE PRECISION BRHAZ(3,2),BRAHA(3),BRAHZ(2,3),BRHCW(5)
      DOUBLE PRECISION BRHIGGS(5),BRNEU(5,5,5),BRCHA(5,3)
      DOUBLE PRECISION BRHSQ(3,10),BRHSL(3,7),BRASQ(2,2),BRASL(2)
      DOUBLE PRECISION BRSUSY(5),WIDTH(5)
      DOUBLE PRECISION HCBRM,HCBRL,HCBRSU,HCBRBU,HCBRSC,HCBRBC
      DOUBLE PRECISION HCBRBT,HCBRWH(5),HCBRWHT,HCBRNC(5,2)
      DOUBLE PRECISION HCBRSQ(5),HCBRSL(3),HCBRSUSY,HCWIDTH
      DOUBLE PRECISION CU(5),CD(5),CV(3),CJ(5),CG(5)
      DOUBLE PRECISION ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION VUS,VCB,VUB
      DOUBLE PRECISION MUR,MUL,MDR,MDL,MLR,MLL,MNL
      DOUBLE PRECISION MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT
      DOUBLE PRECISION CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU,Q2
      DOUBLE PRECISION MGUT,g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      DOUBLE PRECISION MHUS,MHDS,MSS
      DOUBLE PRECISION G1GUT,G2GUT,G3GUT,LGUT,KGUT,HTOPGUT
      DOUBLE PRECISION HBOTGUT,HTAUGUT,M1GUT,M2GUT,M3GUT,ALGUT,AKGUT
      DOUBLE PRECISION ATGUT,ABGUT,ATAUGUT,AMUGUT,MHUGUT,MHDGUT
      DOUBLE PRECISION MSGUT,MQ3GUT,MU3GUT,MD3GUT,MQGUT,MUGUT
      DOUBLE PRECISION MDGUT,ML3GUT,ME3GUT,MLGUT,MEGUT
      DOUBLE PRECISION XIFGUT,XISGUT,MUPGUT,MSPGUT,M3HGUT
      DOUBLE PRECISION XIFSUSY,XISSUSY,MUPSUSY,MSPSUSY,M3HSUSY
      DOUBLE PRECISION OMG,OMGMIN,OMGMAX
      DOUBLE PRECISION Xf,sigmaV,x(100),dNdx(100),EMIN
      DOUBLE PRECISION sigmaPiN,sigma0,csPsi,csNsi,csPsd,csNsd
      DOUBLE PRECISION MHUQ,MHDQ,MSX,LQ,KQ,ALQ,AKQ,MUQ,NUQ
      DOUBLE PRECISION ZHU,ZHD,ZS,H1Q,H2Q,TANBQ,QSTSB
      DOUBLE PRECISION BRSG,BRSGmax,BRSGmin,DMd,DMdmin,DMdmax,DMs,
     . DMsmax,DMsmin,BRBMUMU,BRBMUMUmax,BRBMUMUmin,BRBtaunu,
     . BRBtaunumax,BRBtaunumin
      DOUBLE PRECISION delmagmu,damumin,damumax,amuthmax,amuthmin
      DOUBLE PRECISION LOWSIG(5,7),HIGSIG(5,9)
      DOUBLE PRECISION brtopbw,brtopbh,brtopneutrstop(5,2),toptot
*
      DOUBLE PRECISION chartot2(2),chartot(2),chartotmulti(2)
      DOUBLE PRECISION brcharst1(2),brcharst2(2),brcharsb1(2),
     .          brcharsb2(2),brcharsupl(2),brcharsupr(2),
     .          brcharsdownl(2),
     .          brcharsdownr(2),brcharsnel(2),brcharsn1(2),
     .          brcharsn2(2),brcharsell(2),brcharstau1(2),brcharselr(2),
     .          brcharstau2(2),brcharhcneut(2,5),brcharwneut(2,5),
     .          brcharzchic,brcharHchic(3),brcharAchic(2)
      DOUBLE PRECISION brntaunut(2,5),brnelnue(2,5),brnmunumu(2,5),
     .          brnupdb(2,5),brnchsb(2,5),brntopbb(2,5),
     .          brglupdb(2),brglchsb(2),brgltopbb(2),
     .          brchee,brchmumu,brchtautau,brchnene,
     .          brchnmunmu,brchntauntau,brchupup,brchdodo,
     .          brchchch,brchstst,brchtoptop,brchbotbot
*
      DOUBLE PRECISION neuttot2(5),neuttot(5),neuttot3(5),neuttotrad(5)
      DOUBLE PRECISION brneutst1(5),brneutst2(5),brneutsb1(5),
     .          brneutsb2(5),
     .          brneutsupl(5),brneutsupr(5),brneutsdownl(5),
     .          brneutsdownr(5),brneutsnel(5),brneutsn1(5),
     .          brneutsn2(5),brneutsell(5),brneutselr(5),
     .          brneutstau1(5),brneutstau2(5),brneutwchar(5,2),
     .          brneuthcchar(5,2),brneutzneut(5,5),
     .          brneutHneut(5,5,3),brneutAneut(5,5,2),brnraddec(5,5)
      DOUBLE PRECISION brneutup(5,5),brneutdow(5,5),brneutch(5,5),
     .          brneutst(5,5),brneutbot(5,5),brneuttop(5,5),
     .          brneutel(5,5),brneutmu(5,5),brneuttau(5,5),
     .          brneutnue(5,5),brneutnumu(5,5),brneutnutau(5,5),
     .          brchubd(5,2),brchcbs(5,2),brchtbb(5,2),brchelne(5,2),
     .          brchmunmu(5,2),brchtauntau(5,2),brglup(5),brgldo(5),
     .          brglch(5),brglst(5),brgltop(5),brglbot(5)
*
      DOUBLE PRECISION selltot,selltot2,selltot3,
     . selrtot,selrtot2,selrtot3,stau1tot2,stau2tot2,stau2tot3,stau2tot,
     . sneltot2,sneltot3,sneltot,sntautot2,sntautot3,sntautot
      DOUBLE PRECISION brsellneute(5),brselrneute(5),brsellcharnue(2),
     .          brselrcharnue(2),brsnellneut(5),brsnellchar(5),
     .          brsntauneut(5),brsntauchar(2),brsntau1hcstau(2),
     .          brsntau1wstau(2),brstau1neut(5),brstau2neut(5),
     .          brstau1char(2),
     .          brstau2char(2),brstau1hcsn(2),brstau2hcsn(2),
     .          brstau1wsn(2),brstau2wsn(2),brstau2H(3),brstau2A(2),
     .          brstau2ztau
      DOUBLE PRECISION brselrstau,brselrstaustar,brstau2stau1star,
     .    brstau2stau1,brstau2stau1nn,
     .    brsntaustau1star,brsntaustau1,brsntaustau1nutau,
     .    brsellstau1star,brsellstau1,brsellstau1nutau,
     .    brsnestau1star,brsnestau1,brsnestau1nutau
*
      DOUBLE PRECISION supltot2,suprtot2,sdowltot2,sdowrtot2
      DOUBLE PRECISION brsuplnup(5),brsuplcdow(2),brsuplglui,
     .          brsuprnup(5),brsuprcdow(2),brsuprglui,
     .          brsdowlndow(5),brsdowlchup(2),brsdowlglui,
     .          brsdowrndow(5),brsdowrchup(2),brsdowrglui
*
      DOUBLE PRECISION 
     .          stoptot(2),stoptot2(2),stoptotmulti(2),stoptotrad(2)
      DOUBLE PRECISION brst1neutt(5),brst2neutt(5),brst1charb(2),
     .          brst2charb(2),brst1hcsb(2),brst2hcsb(2),brst1wsb(2),
     .          brst2wsb(2),brst1glui,brst2glui,brst2H(3),brst2A(2),
     .          brst2ztop,brgamma,brgammaup,brgammagluino
      DOUBLE PRECISION brstopw(2,5),brstoph(2,5),brststau(2,2),
     .          brstsntau(2,2),brstsel(2,2),brstsnel(2),
     .          brstbsbst(2,2),brstbbsbt(2,2),brsttausbnu(2,2),
     .          brstelsbnu(2,2),brstupsbdow(2,2),
     .          brst2st1tt,brst2st1startt,brst2st1bb,brst2st1uu,
     .          brst2st1dd,brst2st1ee,brst2st1nunu,brst2st1tautau
*
      DOUBLE PRECISION sbottot(2),sbottot2(2),sbottotmulti(2)
      DOUBLE PRECISION brsb1neutt(5),brsb2neutt(5),brsb1chart(2),
     .          brsb2chart(2),brsb1hcst(2),brsb2hcst(2),
     .          brsb1glui,brsb2glui,brsb1wst(2),
     .          brsb2wst(2),brsb2H(3),brsb2A(2),brsb2zbot
      DOUBLE PRECISION  brsbstau(2,2),brsbsntau(2,2),brsbsel(2,2),
     .          brsbtstsb(2,2),brsbtbstb(2,2),brsbtaustnu(2,2),
     .          brsbelstnu(2,2),brsbupstdow(2,2),brsbsnel(2),
     .          brsb2sb1bb,brsb2sb1starbb,brsb2sb1tt,
     .          brsb2sb1uu,brsb2sb1dd,brsb2sb1ee,brsb2sb1nunu,
     .          brsb2sb1tautau
      DOUBLE PRECISION gluitot,gluitot2,gluitotmulti,gluitotrad
      DOUBLE PRECISION brgst1,brgst2,brgsb1,brgsb2,brgsupl,brgsupr,
     .         brgsdownl,brgsdownr,brglnjgluon(5)
      DOUBLE PRECISION brgoup(5),brgoch(5),brgodn(5),brgost(5),
     .         brgotp(5),
     .         brgobt(5),brgoud(2),brgocs(2),brgotb(2),brhcst1b,brwst1b
*************************************************************
*
      COMMON/CHARGINO_WIDTH/chartot2,chartot,chartotmulti
      COMMON/CHARGINO_BR_2BD/brcharst1,brcharst2,brcharsb1,
     .          brcharsb2,brcharsupl,brcharsupr,brcharsdownl,
     .          brcharsdownr,brcharsnel,brcharsn1,
     .          brcharsn2,brcharsell,brcharstau1,brcharselr,
     .          brcharstau2,brcharhcneut,brcharwneut,
     .          brcharzchic,brcharHchic,brcharAchic
      COMMON/CHARGINO_BR_3BD/brntaunut,brnelnue,brnmunumu,
     .          brnupdb,brnchsb,brntopbb,
     .          brglupdb,brglchsb,brgltopbb,
     .          brchee,brchmumu,brchtautau,brchnene,
     .          brchnmunmu,brchntauntau,brchupup,brchdodo,
     .          brchchch,brchstst,brchtoptop,brchbotbot
      COMMON/NEUTRALINO_WIDTH/neuttot2,neuttot,neuttot3,neuttotrad
      COMMON/NEUTRALINO_BR_2BD/brneutst1,brneutst2,brneutsb1,brneutsb2,
     .         brneutsupl,brneutsupr,brneutsdownl,brneutsdownr,
     .         brneutsnel,brneutsn1,brneutsn2,brneutsell,brneutselr,
     .         brneutstau1,brneutstau2,brneutwchar,brneuthcchar,
     .         brneutzneut,brneutHneut,brneutAneut,brnraddec
      COMMON/NEUTRALINO_BR_3BD/brneutup,brneutdow,brneutch,brneutst,
     .         brneutbot,brneuttop,brneutel,brneutmu,brneuttau,
     .         brneutnue,brneutnumu,brneutnutau,brchubd,brchcbs, 
     .         brchtbb,brchelne,brchmunmu,brchtauntau,brglup,brgldo,
     .         brglch,brglst,brgltop,brglbot
*
      COMMON/SLEPTON_WIDTH/selltot,selltot2,selltot3,
     .selrtot,selrtot2,selrtot3,stau1tot2,stau2tot2,stau2tot3,stau2tot,
     .sneltot2,sneltot3,sneltot,sntautot2,sntautot3,sntautot
      COMMON/SLEPTON_BR_2BD/brsellneute,brselrneute,brsellcharnue,
     .          brselrcharnue,brsnellneut,brsnellchar,
     .          brsntauneut,brsntauchar,brsntau1hcstau,
     .          brsntau1wstau,brstau1neut,brstau2neut,brstau1char,
     .          brstau2char,brstau1hcsn,brstau2hcsn,
     .          brstau1wsn,brstau2wsn,brstau2H,brstau2A,brstau2ztau
      COMMON/SLEPTON_BR_3BD/brselrstau,brselrstaustar,brstau2stau1star,
     .   brstau2stau1,brstau2stau1nn,
     .   brsntaustau1star,brsntaustau1,brsntaustau1nutau,
     .   brsellstau1star,brsellstau1,brsellstau1nutau,
     .   brsnestau1star,brsnestau1,brsnestau1nutau
*
      COMMON/SQUARK_WIDTH/supltot2,suprtot2,sdowltot2,sdowrtot2
      COMMON/SQUARK_BR_2BD/brsuplnup,brsuplcdow,brsuplglui,
     .          brsuprnup,brsuprcdow,brsuprglui,
     .          brsdowlndow,brsdowlchup,brsdowlglui,
     .          brsdowrndow,brsdowrchup,brsdowrglui
*
      COMMON/STOP_WIDTH/stoptot,stoptot2,stoptotmulti,stoptotrad
      COMMON/STOP_BR_2BD/brst1neutt,brst2neutt,brst1charb,
     .          brst2charb,brst1hcsb,brst2hcsb,brst1wsb,
     .          brst2wsb,brst1glui,brst2glui,brst2H,brst2A,
     .          brst2ztop,brgamma,brgammaup,brgammagluino
      COMMON/STOP_BR_3BD/brstopw,brstoph,brststau,
     .          brstsntau,brstsel,brstsnel,
     .          brstbsbst,brstbbsbt,brsttausbnu,
     .          brstelsbnu,brstupsbdow,
     .          brst2st1tt,brst2st1startt,brst2st1bb,brst2st1uu,
     .          brst2st1dd,brst2st1ee,brst2st1nunu,brst2st1tautau
*
      COMMON/SBOTTOM_WIDTH/sbottot,sbottot2,sbottotmulti
      COMMON/SBOTTOM_BR_2BD/brsb1neutt,brsb2neutt,brsb1chart,
     .          brsb2chart,brsb1hcst,brsb2hcst,
     .          brsb1glui,brsb2glui,brsb1wst,
     .          brsb2wst,brsb2H,brsb2A,brsb2zbot
      COMMON/SBOTTOM_BR_3BD/brsbstau,brsbsntau,brsbsel,
     .          brsbtstsb,brsbtbstb,brsbtaustnu,
     .          brsbelstnu,brsbupstdow,brsbsnel,
     .          brsb2sb1bb,brsb2sb1starbb,brsb2sb1tt,
     .          brsb2sb1uu,brsb2sb1dd,brsb2sb1ee,brsb2sb1nunu,
     .          brsb2sb1tautau
*
      COMMON/GLUINO_WIDTH/gluitot,gluitot2,gluitotmulti,gluitotrad
      COMMON/GLUINO_BR_2BD/brgst1,brgst2,brgsb1,brgsb2,brgsupl,brgsupr,
     .         brgsdownl,brgsdownr,brglnjgluon
      COMMON/GLUINO_BR_3BD/brgoup,brgoch,brgodn,brgost,brgotp,
     .         brgobt,brgoud,brgocs,brgotb,brhcst1b,brwst1b
*
      COMMON/BRN/BRJJ,BRMM,BRLL,BRSS,BRCC,BRBB,BRTT,BRWW,BRZZ,
     . BRGG,BRZG,BRHHH,BRHAA,BRHCHC,BRHAZ,BRAHA,BRAHZ,
     . BRHCW,BRHIGGS,BRNEU,BRCHA,BRHSQ,BRHSL,BRASQ,BRASL,
     . BRSUSY,WIDTH
      COMMON/BRC/HCBRM,HCBRL,HCBRSU,HCBRBU,HCBRSC,HCBRBC,
     . HCBRBT,HCBRWH,HCBRWHT,HCBRNC,HCBRSQ,HCBRSL,
     . HCBRSUSY,HCWIDTH
      COMMON/BRSG/BRSG,BRSGmax,BRSGmin,DMd,DMdmin,DMdmax,DMs,
     . DMsmax,DMsmin,BRBMUMU,BRBMUMUmax,BRBMUMUmin,BRBtaunu,
     . BRBtaunumax,BRBtaunumin
      COMMON/MAGMU/delmagmu,damumin,damumax,amuthmax,amuthmin
      COMMON/BR_top2body/brtopbw,brtopbh,brtopneutrstop
      COMMON/topwidth/toptot
      COMMON/REDCOUP/CU,CD,CV,CJ,CG
      COMMON/HIGGSPEC/SMASS,SCOMP,PMASS,PCOMP,CMASS
      COMMON/SUSYSPEC/MGL,MCHA,U,V,MNEU,NEU
      COMMON/GAUGE/ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/CKM/VUS,VCB,VUB
      COMMON/SFSPEC/MUR,MUL,MDR,MDL,MLR,MLL,MNL,
     . MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT,
     . CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU
      COMMON/RENSCALE/Q2
      COMMON/STSBSCALE/QSTSB
      COMMON/MGUT/MGUT
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      COMMON/SUSYMH/MHUS,MHDS,MSS
      COMMON/QMHIGGS/MHUQ,MHDQ,MSX
      COMMON/QPAR/LQ,KQ,ALQ,AKQ,MUQ,NUQ
      COMMON/QHIGGS/ZHU,ZHD,ZS,H1Q,H2Q,TANBQ
      COMMON/GUTCOUP/G1GUT,G2GUT,G3GUT,LGUT,KGUT,HTOPGUT,
     . HBOTGUT,HTAUGUT
      COMMON/GUTPAR/M1GUT,M2GUT,M3GUT,ALGUT,AKGUT,ATGUT,ABGUT,
     . ATAUGUT,AMUGUT,MHUGUT,MHDGUT,MSGUT,MQ3GUT,MU3GUT,MD3GUT,
     . MQGUT,MUGUT,MDGUT,ML3GUT,ME3GUT,MLGUT,MEGUT
      COMMON/GUTEXT/XIFGUT,XISGUT,MUPGUT,MSPGUT,M3HGUT
      COMMON/SUSYEXT/XIFSUSY,XISSUSY,MUPSUSY,MSPSUSY,M3HSUSY
      COMMON/MICROMG/OMG,OMGMIN,OMGMAX,Xf,sigmaV,x,dNdx,EMIN,NBIN
      COMMON/MICROMG2/sigmaPiN,sigma0,csPsi,csNsi,csPsd,csNsd
      COMMON/LOWSIG/LOWSIG
      COMMON/HIGSIG/HIGSIG
      COMMON/IDIM/IDIM
      COMMON/RES/RES

      IDIM=IDIM+1      
      RES(IDIM,1)=PAR(1)
      RES(IDIM,2)=PAR(2)
      RES(IDIM,3)=PAR(3)
      RES(IDIM,4)=PAR(4)
      RES(IDIM,5)=PAR(5)
      RES(IDIM,6)=PAR(6)
      RES(IDIM,7)=PAR(20)
      RES(IDIM,8)=PAR(21)
      RES(IDIM,9)=PAR(22)
      RES(IDIM,10)=PAR(23)
      RES(IDIM,11)=PAR(24)

      IF(IFAIL.EQ.0.OR.IFAIL.GE.10)THEN

       DO I=1,3
        RES(IDIM,10+2*I)=SMASS(I)
        RES(IDIM,11+2*I)=SCOMP(I,3)**2
       ENDDO
       DO I=1,2
        RES(IDIM,16+2*I)=PMASS(I)
        RES(IDIM,17+2*I)=PCOMP(I,2)**2
       ENDDO
       RES(IDIM,22)=CMASS
       RES(IDIM,23)=DABS(MNEU(1))
       RES(IDIM,24)=NEU(1,1)**2
       RES(IDIM,25)=NEU(1,3)**2+NEU(1,4)**2
       RES(IDIM,26)=NEU(1,5)**2

*   SUSY
       S=0d0
       DO I=1,2
        S=S+PROB(I)
       ENDDO
       DO I=20,25
        S=S+PROB(I)
       ENDDO
       RES(IDIM,28)=S

*   Higgs
       S=0d0
       DO I=3,19
        S=S+PROB(I)
       ENDDO
       S=S+PROB(41)
       RES(IDIM,29)=S

*   B Physics
       S=0d0
       DO I=32,36
        S=S+DABS(PROB(I))
       ENDDO
       DO I=38,40
        S=S+DABS(PROB(I))
       ENDDO
       RES(IDIM,30)=S

*   (g-2)_muon
       RES(IDIM,31)=PROB(37)

*   Theoretical
       S=0d0
       DO I=27,29
        S=S+PROB(I)
       ENDDO
       RES(IDIM,32)=S

*   Charged LSP
       RES(IDIM,33)=PROB(26)

*   WMAP
       RES(IDIM,34)=PROB(30)

      ELSE

       DO I=12,34
        RES(IDIM,I)=0d0
       ENDDO

      ENDIF
       
      IF(ITOT.EQ.NTOT.OR.IDIM.EQ.IMAX)THEN
       DO I=1,IDIM
        WRITE(6,10)(RES(I,J),J=1,34)
 10     FORMAT(34E16.8)
       ENDDO
       IF(IDIM.EQ.IMAX)IDIM=0
      ENDIF

      END

      SUBROUTINE ERROR(TOT,NTOT,NFAIL)

*********************************************************************      
*   Subroutine for the error file. It contains a summary of the scan:
*   Number of points that passed/failed the tests
*   and ranges for scanned parameters that passed the tests
*********************************************************************      
 
      IMPLICIT NONE

      INTEGER OMGFLAG,MAFLAG,M1FLAG,M3FLAG,I,S,TOT,NTOT,NFAIL(*)

      DOUBLE PRECISION LN,LNN,KN,KNN,TBN,TBNN,MUN,MUNN
      DOUBLE PRECISION ALN,ALNN,AKN,AKNN,MAN,MANN
      DOUBLE PRECISION M1N,M1NN,M2N,M2NN,M3N,M3NN

      COMMON/BOUNDS/LN,LNN,KN,KNN,TBN,TBNN,MUN,MUNN,
     . ALN,ALNN,AKN,AKNN,MAN,MANN,
     . M1N,M1NN,M2N,M2NN,M3N,M3NN
      COMMON/FLAGS/OMGFLAG,MAFLAG
      COMMON/SCANFLAGS/M1FLAG,M3FLAG

      WRITE(0,*)
      WRITE(0,20)"Number of points:                       "
      WRITE(0,*)
      WRITE(0,20)"  scanned                               ",NTOT
      WRITE(0,20)"  l, k, tan(beta) or mu=0               ",NFAIL(9)
      S=0
      DO I=1,7
       S=S+NFAIL(I)
      ENDDO
      WRITE(0,20)"  with mh1^2 or ma1^2 or mhc^2 < 0      ",S
      WRITE(0,20)"  with m_sfermion^2 < 0                 ",NFAIL(8)
      WRITE(0,20)"  violating phenomenological constraints",NFAIL(10)
      S=NFAIL(11)+NFAIL(12)
      WRITE(0,20)"  RGE integration problem               ",S
      WRITE(0,*)
      WRITE(0,20)"Remaining good points                   ",TOT
      IF(TOT.GT.0)THEN
       WRITE(0,*)
       WRITE(0,20)"Parameter ranges for good points:       "
       WRITE(0,*)
       WRITE(0,30)" TANB: ",TBN,TBNN
       IF(M1FLAG.EQ.1)THEN
        WRITE(0,30)" M1: ",M1N,M1NN
       ENDIF
       WRITE(0,30)" M2: ",M2N,M2NN
       IF(M3FLAG.EQ.1)THEN
        WRITE(0,30)" M3: ",M3N,M3NN
       ENDIF
       WRITE(0,30)" LAMBDA: ",LN,LNN
       WRITE(0,30)" KAPPA: ",KN,KNN
       WRITE(0,30)" MUEFF: ",MUN,MUNN
       IF(MAFLAG.EQ.0)THEN
        WRITE(0,30)" ALAMBDA: ",ALN,ALNN
        WRITE(0,30)" MA: ",MAN,MANN
        WRITE(0,20)"(MA is not an input parameter)"
       ELSE
        WRITE(0,30)" ALAMBDA: ",ALN,ALNN
        WRITE(0,20)"(ALAMBDA is not an input parameter)"
        WRITE(0,30)" MA: ",MAN,MANN
       ENDIF
       WRITE(0,30)" AKAPPA: ",AKN,AKNN
      ENDIF

 20   FORMAT(A40,I10)
 30   FORMAT(A15,2E15.4)

      END
