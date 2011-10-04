      SUBROUTINE LSIG()

*   Subroutine to calculate LHC Significances at 30 fb-1:
*    1) bbh/a -> b b tau tau
*    2) h/a -> gamma gamma
*    3) gg -> h -> ZZ -> 4 leptons
*    4) gg -> h -> WW -> 2 leptons 2 neutrinos
*    5) WW -> h -> tau tau
*    6) WW -> h -> W W
*    7) WW -> h -> gamma gamma

      IMPLICIT NONE

      INTEGER I,J

      DOUBLE PRECISION MH(5),SIG(5,7),SIGH(5,7),DM(7)
      DOUBLE PRECISION SMASS(3),SCOMP(3,3),PMASS(2),PCOMP(2,2),CMASS
      DOUBLE PRECISION BRJJ(5),BRMM(5),BRLL(5),BRSS(5),BRCC(5)
      DOUBLE PRECISION BRBB(5),BRTT(5),BRWW(3),BRZZ(3),BRGG(5)
      DOUBLE PRECISION BRZG(5),BRHHH(4),BRHAA(3,3),BRHCHC(3)
      DOUBLE PRECISION BRHAZ(3,2),BRAHA(3),BRAHZ(2,3),BRHCW(5)
      DOUBLE PRECISION BRHIGGS(5),BRNEU(5,5,5),BRCHA(5,3)
      DOUBLE PRECISION BRHSQ(3,10),BRHSL(3,7),BRASQ(2,2),BRASL(2)
      DOUBLE PRECISION BRSUSY(5),WIDTH(5)
      DOUBLE PRECISION CU(5),CD(5),CV(3),CJ(5),CG(5)
      DOUBLE PRECISION BRJJSM,BRMMSM,BRLLSM,BRSSSM,BRCCSM
      DOUBLE PRECISION BRBBSM,BRTTSM,BRWWSM,BRZZSM,BRGGSM,BRZGSM
      DOUBLE PRECISION TAUTAULOW,GAGALOW,ZZ4LLOW,WW2LLOW
      DOUBLE PRECISION WWTTLOW,WWWWLOW,WWGGLOW

      COMMON/BRN/BRJJ,BRMM,BRLL,BRSS,BRCC,BRBB,BRTT,BRWW,BRZZ,
     .      BRGG,BRZG,BRHHH,BRHAA,BRHCHC,BRHAZ,BRAHA,BRAHZ,
     .      BRHCW,BRHIGGS,BRNEU,BRCHA,BRHSQ,BRHSL,BRASQ,BRASL,
     .      BRSUSY,WIDTH
      COMMON/REDCOUP/CU,CD,CV,CJ,CG
      COMMON/HIGGSPEC/SMASS,SCOMP,PMASS,PCOMP,CMASS
      COMMON/LOWSIG/SIG

*   Loop over h1, h2, h3

      DO I=1,3

       CALL HDECAY(SMASS(I),BRJJSM,BRMMSM,BRLLSM,BRSSSM,BRCCSM,
     .     BRBBSM,BRTTSM,BRWWSM,BRZZSM,BRGGSM,BRZGSM)

*   Process 1), bbh -> bbtautau

       SIG(I,1)=TAUTAULOW(SMASS(I))*CD(I)**2*BRLL(I)/1.D-1

*   Process 2), gg -> h -> gamma gamma

       SIG(I,2)=GAGALOW(SMASS(I))*CJ(I)**2*BRGG(I)/BRGGSM

*   Process 3), gg -> h -> ZZ -> 4 leptons

       IF(BRZZSM.NE.0d0)THEN
        SIG(I,3)=ZZ4LLOW(SMASS(I))*CJ(I)**2*BRZZ(I)/BRZZSM
       ELSE
        SIG(I,3)=0d0
       ENDIF

*   Process 4), gg -> h -> WW -> 2 leptons + nu + nu

       IF(BRWWSM.NE.0d0)THEN
        SIG(I,4)=WW2LLOW(SMASS(I))*CJ(I)**2*BRWW(I)/BRWWSM
       ELSE
        SIG(I,4)=0d0
       ENDIF

*   Process 5), WW -> h -> tautau

       IF(BRLLSM.NE.0d0)THEN
        SIG(I,5)=WWTTLOW(SMASS(I))*CV(I)**2*BRLL(I)/BRLLSM
       ELSE
        SIG(I,5)=0d0
       ENDIF

*   Process 6), WW -> h -> WW

       IF(BRWWSM.NE.0d0)THEN
        SIG(I,6)=WWWWLOW(SMASS(I))*CV(I)**2*BRWW(I)/BRWWSM
       ELSE
        SIG(I,6)=0d0
       ENDIF

*   Process 7), WW -> h -> gamma gamma

       SIG(I,7)=WWGGLOW(SMASS(I))*CV(I)**2*BRGG(I)/BRGGSM

      ENDDO

*   Loop over a1, a2

      DO I=1,2

       CALL ADECAY(PMASS(I),BRJJSM,BRMMSM,BRLLSM,BRSSSM,BRCCSM,
     .     BRBBSM,BRTTSM,BRGGSM,BRZGSM)

*   Process 1), bba -> bbtautau

       SIG(I+3,1)=TAUTAULOW(PMASS(I))*CD(I+3)**2*BRLL(I+3)/1.D-1

*   Process 2), a -> gamma gamma

       SIG(I+3,2)=GAGALOW(PMASS(I))*CJ(I+3)**2*BRGG(I+3)/BRGGSM

*   Processes 5-7) = 0 (No aVV couplings)

       DO J=3,7
        SIG(I+3,J)=0d0
       ENDDO

      ENDDO

*   Modify SIG(IJ,K) if |M(I)-M(J)| < DM(K)*(M(I)+M(J))/2 :

      DM(1)=15.D-2
      DM(2)=1.D-2
      DM(3)=1.D-2
      DM(4)=1.D-2
      DM(5)=10.D-2
      DM(6)=10.D-2
      DM(7)=10.D-2

      CALL CLASSIFY(SMASS,PMASS,MH,SIG,SIGH,7)
      CALL COMB(MH,DM,SIGH,7)
      CALL DECLASSIFY(SMASS,PMASS,MH,SIG,SIGH,7)

      END


      SUBROUTINE HSIG()

*   Subroutine to calculate LHC Significances at 300 fb-1:
*    1) h/a -> gamma gamma
*    2) h/a -> gamma gamma lepton
*    3) tth/a -> bb + X
*    4) bbh/a -> bbtautau
*    5) gg -> h -> ZZ -> 4 leptons
*    6) gg -> h -> WW -> 2 leptons 2 neutrinos
*    7) WW -> h -> tautau
*    8) WW -> h -> WW
*    9) WW -> h -> invisible

      IMPLICIT NONE

      INTEGER I,J

      DOUBLE PRECISION MH(5),SIG(5,9),SIGH(5,9),DM(9),LUM
      DOUBLE PRECISION SMASS(3),SCOMP(3,3),PMASS(2),PCOMP(2,2),CMASS
      DOUBLE PRECISION BRJJ(5),BRMM(5),BRLL(5),BRSS(5),BRCC(5)
      DOUBLE PRECISION BRBB(5),BRTT(5),BRWW(3),BRZZ(3),BRGG(5)
      DOUBLE PRECISION BRZG(5),BRHHH(4),BRHAA(3,3),BRHCHC(3)
      DOUBLE PRECISION BRHAZ(3,2),BRAHA(3),BRAHZ(2,3),BRHCW(5)
      DOUBLE PRECISION BRHIGGS(5),BRNEU(5,5,5),BRCHA(5,3)
      DOUBLE PRECISION BRHSQ(3,10),BRHSL(3,7),BRASQ(2,2),BRASL(2)
      DOUBLE PRECISION BRSUSY(5),WIDTH(5),BRINV
      DOUBLE PRECISION CU(5),CD(5),CV(3),CJ(5),CG(5)
      DOUBLE PRECISION BRJJSM,BRMMSM,BRLLSM,BRSSSM,BRCCSM
      DOUBLE PRECISION BRBBSM,BRTTSM,BRWWSM,BRZZSM,BRGGSM,BRZGSM
      DOUBLE PRECISION GAGAHIG,GAGALHIG,TTBBHIG
      DOUBLE PRECISION TAUTAUHIG,ZZ4LHIG,WW2LHIG
      DOUBLE PRECISION WWTTHIG,WWWWHIG,WWINVHIG

      COMMON/BRN/BRJJ,BRMM,BRLL,BRSS,BRCC,BRBB,BRTT,BRWW,BRZZ,
     .      BRGG,BRZG,BRHHH,BRHAA,BRHCHC,BRHAZ,BRAHA,BRAHZ,
     .      BRHCW,BRHIGGS,BRNEU,BRCHA,BRHSQ,BRHSL,BRASQ,BRASL,
     .      BRSUSY,WIDTH
      COMMON/REDCOUP/CU,CD,CV,CJ,CG
      COMMON/HIGGSPEC/SMASS,SCOMP,PMASS,PCOMP,CMASS
      COMMON/HIGSIG/SIG

*   Luminosity

      LUM=300d0

*   Loop over h1, h2, h3

      DO I=1,3

       CALL HDECAY(SMASS(I),BRJJSM,BRMMSM,BRLLSM,BRSSSM,BRCCSM,
     .     BRBBSM,BRTTSM,BRWWSM,BRZZSM,BRGGSM,BRZGSM)

*   Process 1), gg -> h -> gamma gamma

       SIG(I,1)=GAGAHIG(SMASS(I))*CJ(I)**2*BRGG(I)/BRGGSM

*   Process 2), Wh + tth -> gamma gamma lepton

       SIG(I,2)=GAGALHIG(SMASS(I))*(CV(I)**2+CU(I)**2)/2d0
     .  *BRGG(I)/BRGGSM

*   Process 3), tth -> bb + X

       SIG(I,3)=TTBBHIG(SMASS(I))*CU(I)**2*BRBB(I)

*   Process 4), bbh -> bbtautau

       SIG(I,4)=TAUTAUHIG(SMASS(I))*CD(I)**2*BRLL(I)/9.D-2

*   Process 5), gg -> h -> ZZ -> 4 leptons

       IF(BRZZSM.NE.0d0)THEN
        SIG(I,5)=ZZ4LHIG(SMASS(I))*CJ(I)**2*BRZZ(I)/BRZZSM
       ELSE
        SIG(I,5)=0d0
       ENDIF

*   Process 6), gg -> h -> WW -> 2 leptons + nu + nu

       IF(BRWWSM.NE.0d0)THEN
        SIG(I,6)=WW2LHIG(SMASS(I))*CJ(I)**2*BRWW(I)/BRWWSM
       ELSE
        SIG(I,6)=0d0
       ENDIF

*   Process 7), WW -> h -> tautau

       IF(BRLLSM.NE.0d0)THEN
        SIG(I,7)=WWTTHIG(SMASS(I))*CV(I)**2*BRLL(I)/BRLLSM
       ELSE
        SIG(I,7)=0d0
       ENDIF

*   Process 8), WW -> h -> WW

       IF(BRWWSM.NE.0d0)THEN
        SIG(I,8)=WWWWHIG(SMASS(I))*CV(I)**2*BRWW(I)/BRWWSM
       ELSE
        SIG(I,8)=0d0
       ENDIF

*   Process 9), WW -> h -> invisible

       BRINV=BRNEU(I,1,1)
     .  +BRHAA(I,1)*BRNEU(4,1,1)**2
     .  +BRHAA(I,2)*BRNEU(4,1,1)*BRNEU(5,1,1)
     .  +BRHAA(I,3)*BRNEU(5,1,1)**2
       IF(I.EQ.2)
     .  BRINV=BRINV+BRHHH(1)*BRNEU(1,1,1)**2
       IF(I.EQ.3)
     .  BRINV=BRINV+BRHHH(2)*BRNEU(1,1,1)**2
     .   +BRHHH(3)*BRNEU(1,1,1)*BRNEU(2,1,1)
     .   +BRHHH(4)*BRNEU(2,1,1)**2
       SIG(I,9)=WWINVHIG(SMASS(I))*CV(I)**2*BRINV

      ENDDO

*   Loop over a1, a2

      DO I=1,2

       CALL ADECAY(PMASS(I),BRJJSM,BRMMSM,BRLLSM,BRSSSM,BRCCSM,
     .     BRBBSM,BRTTSM,BRGGSM,BRZGSM)

*   Process 1), gg -> a -> gamma gamma

       SIG(I+3,1)=GAGAHIG(PMASS(I))*CJ(I+3)**2*BRGG(I+3)/BRGGSM

*   Process 2), tta -> gamma gamma lepton

       SIG(I+3,2)=GAGALHIG(PMASS(I))*CU(I+3)**2/2d0*BRGG(I+3)/BRGGSM

*   Process 3), tta -> bb + X

       SIG(I+3,3)=TTBBHIG(PMASS(I))*CU(I+3)**2*BRBB(I+3)

*   Process 4), bba -> bbtautau

       SIG(I+3,4)=TAUTAUHIG(PMASS(I))*CD(I+3)**2*BRLL(I+3)/9.D-2

*   Processes 5-9) = 0 (No aVV couplings)

       DO J=5,9
        SIG(I+3,J)=0d0
       ENDDO

      ENDDO

*   Modify SIG(IJ,K) if |M(I)-M(J)| < DM(K)*(M(I)+M(J))/2 :

      DM(1)=1.D-2
      DM(2)=1.D-2
      DM(3)=10.D-2
      DM(4)=15.D-2
      DM(5)=1.D-2
      DM(6)=1.D-2
      DM(7)=10.D-2
      DM(8)=10.D-2
      DM(9)=10.D-2

      CALL CLASSIFY(SMASS,PMASS,MH,SIG,SIGH,9)
      CALL COMB(MH,DM,SIGH,9)
      CALL DECLASSIFY(SMASS,PMASS,MH,SIG,SIGH,9)

*   Rescaling to Luminosity LUM [fb-1]

      DO I=1,5
       SIG(I,1)=SIG(I,1)*DSQRT(LUM/1.D2)
       SIG(I,2)=SIG(I,2)*DSQRT(LUM/1.D2)
       SIG(I,3)=SIG(I,3)*DSQRT(DSQRT(LUM/1.D2))
       SIG(I,4)=SIG(I,4)*DSQRT(LUM/1.D2)
       SIG(I,5)=SIG(I,5)*DSQRT(LUM/1.D2)
       SIG(I,6)=SIG(I,6)*DSQRT(LUM/1.D2)
       SIG(I,7)=SIG(I,7)*DSQRT(DSQRT(LUM/1.D2))
       SIG(I,8)=SIG(I,8)*DSQRT(DSQRT(LUM/1.D2))
       SIG(I,9)=SIG(I,9)
      ENDDO

      END


      DOUBLE PRECISION FUNCTION TAUTAULOW(XMH)

      IMPLICIT NONE

      INTEGER NH,NHM
      PARAMETER (NHM=15)
      DOUBLE PRECISION XMH,MH(NHM),TB(NHM)

      DATA MH/107d0,110d0,120d0,125d0,130d0,135d0,140d0,200d0,
     .      240d0,250d0,285d0,295d0,330d0,500d0,800d0/
      DATA TB/55d0,50d0,30d0,26d0,26d0,17d0,16.7d0,14d0,
     .      17d0,18d0,21d0,22d0,26d0,34d0,49d0/

      TAUTAULOW=0d0

      IF(XMH.LT.MH(1) .OR. XMH.GT.MH(NHM)) RETURN

      DO 1 NH=1,NHM
       IF(MH(NH).LE.XMH) GOTO 1
       TAUTAULOW=2.5d0/TB(NH)**2+(MH(NH)-XMH)/(MH(NH)-MH(NH-1))
     .   *2.5d0*(1d0/TB(NH-1)**2-1d0/TB(NH)**2)
       GOTO 2
 1    CONTINUE
 2    RETURN

      END


      DOUBLE PRECISION FUNCTION GAGALOW(XMH)

      IMPLICIT NONE

      INTEGER NH,NHM
      PARAMETER (NHM=6)
      DOUBLE PRECISION XMH,XN(NHM),MH(NHM)

      DATA MH/100d0,110d0,120d0,130d0,140d0,150d0/
      DATA XN/2.9d0,3.9d0,4.7d0,5d0,4.5d0,3.3d0/

      GAGALOW=0d0

      IF(XMH.LT.MH(1) .OR. XMH.GT.MH(NHM)) RETURN

      DO 1 NH=1,NHM
       IF(MH(NH).LE.XMH) GOTO 1
       GAGALOW=XN(NH)+(MH(NH)-XMH)/(MH(NH)-MH(NH-1))
     .   *(XN(NH-1)-XN(NH))
       GOTO 2
 1    CONTINUE
 2    RETURN

      END


      DOUBLE PRECISION FUNCTION ZZ4LLOW(XMH)

      IMPLICIT NONE

      INTEGER NH,NHM
      PARAMETER (NHM=19)
      DOUBLE PRECISION XMH,XN(NHM),MH(NHM)

      DATA MH/115d0,120d0,130d0,140d0,150d0,160d0,165d0,170d0,
     .      180d0,190d0,200d0,250d0,300d0,350d0,400d0,450d0,
     .      500d0,550d0,600d0/
      DATA XN/2.7d0,4d0,9.1d0,13.6d0,15.5d0,9.4d0,4.7d0,4.6d0,
     .      7.9d0,15.8d0,16d0,13.3d0,13.8d0,14d0,13.3d0,11.3d0,
     .      8.8d0,6.8d0,5.7d0/

      ZZ4LLOW=0d0

      IF(XMH.LT.MH(1) .OR. XMH.GT.MH(NHM)) RETURN

      DO 1 NH=1,NHM
       IF(MH(NH).LE.XMH) GOTO 1
       ZZ4LLOW=XN(NH)+(MH(NH)-XMH)/(MH(NH)-MH(NH-1))
     .   *(XN(NH-1)-XN(NH))
       GOTO 2
 1    CONTINUE
 2    RETURN

      END


      DOUBLE PRECISION FUNCTION WW2LLOW(XMH)

      IMPLICIT NONE

      INTEGER NH,NHM
      PARAMETER (NHM=8)
      DOUBLE PRECISION XMH,XN(NHM),MH(NHM)

      DATA MH/140d0,150d0,160d0,165d0,170d0,180d0,190d0,200d0/
      DATA XN/2.7d0,7.0d0,10.3d0,11.1d0,8.4d0,7.0d0,3.1d0,1.9d0/

      WW2LLOW=0d0

      IF(XMH.LT.MH(1) .OR. XMH.GT.MH(NHM)) RETURN

      DO 1 NH=1,NHM
       IF(MH(NH).LE.XMH) GOTO 1
       WW2LLOW=XN(NH)+(MH(NH)-XMH)/(MH(NH)-MH(NH-1))
     .   *(XN(NH-1)-XN(NH))
       GOTO 2
 1    CONTINUE
 2    RETURN

      END


      DOUBLE PRECISION FUNCTION WWTTLOW(XMH)

      IMPLICIT NONE

      INTEGER NH,NHM
      PARAMETER (NHM=4)
      DOUBLE PRECISION XMH,XN(NHM),MH(NHM)

      DATA MH/115d0,125d0,135d0,145d0/
      DATA XN/4d0,3.7d0,3.9d0,2.2d0/

      WWTTLOW=0d0

      IF(XMH.LT.MH(1) .OR. XMH.GT.MH(NHM)) RETURN

      DO 1 NH=1,NHM
       IF(MH(NH).LE.XMH) GOTO 1
       WWTTLOW=XN(NH)+(MH(NH)-XMH)/(MH(NH)-MH(NH-1))
     .   *(XN(NH-1)-XN(NH))
       GOTO 2
 1    CONTINUE
 2    RETURN

      END


      DOUBLE PRECISION FUNCTION WWWWLOW(XMH)

      IMPLICIT NONE

      INTEGER NH,NHM
      PARAMETER (NHM=11)
      DOUBLE PRECISION XMH,XN(NHM),MH(NHM)

      DATA MH/130d0,140d0,150d0,160d0,170d0,180d0,190d0,
     .      200d0,210d0,220d0,250d0/
      DATA XN/2.2d0,5.3d0,7.3d0,8.1d0,8.1d0,7.3d0,6.3d0,
     .      5.3d0,4.2d0,3.4d0,2.2d0/

      WWWWLOW=0d0

      IF(XMH.LT.MH(1) .OR. XMH.GT.MH(NHM)) RETURN

      DO 1 NH=1,NHM
       IF(MH(NH).LE.XMH) GOTO 1
       WWWWLOW=XN(NH)+(MH(NH)-XMH)/(MH(NH)-MH(NH-1))
     .   *(XN(NH-1)-XN(NH))
       GOTO 2
 1    CONTINUE
 2    RETURN

      END


      DOUBLE PRECISION FUNCTION WWGGLOW(XMH)

      IMPLICIT NONE

      INTEGER NH,NHM
      PARAMETER (NHM=5)
      DOUBLE PRECISION XMH,XN(NHM),MH(NHM)

      DATA MH/115d0,120d0,130d0,140d0,150d0/
      DATA XN/2.3d0,2.3d0,2.3d0,1.7d0,1.7d0/

      WWGGLOW=0d0

      IF(XMH.LT.MH(1) .OR. XMH.GT.MH(NHM)) RETURN

      DO 1 NH=1,NHM
       IF(MH(NH).LE.XMH) GOTO 1
       WWGGLOW=XN(NH)+(MH(NH)-XMH)/(MH(NH)-MH(NH-1))
     .   *(XN(NH-1)-XN(NH))
       GOTO 2
 1    CONTINUE
 2    RETURN

      END


      DOUBLE PRECISION FUNCTION GAGAHIG(XMH)

      IMPLICIT NONE

      INTEGER NH,NHM
      PARAMETER (NHM=6)
      DOUBLE PRECISION XMH,XN(NHM),MH(NHM)

      DATA MH/100d0,110d0,120d0,130d0,140d0,150d0/
      DATA XN/4.2d0,6d0,6.8d0,8.2d0,7d0,5.2d0/

      GAGAHIG=0d0

      IF(XMH.LT.MH(1) .OR. XMH.GT.MH(NHM)) RETURN

      DO 1 NH=1,NHM
       IF(MH(NH).LE.XMH) GOTO 1
       GAGAHIG=XN(NH)+(MH(NH)-XMH)/(MH(NH)-MH(NH-1))
     .   *(XN(NH-1)-XN(NH))
       GOTO 2
 1    CONTINUE
 2    RETURN

      END


      DOUBLE PRECISION FUNCTION GAGALHIG(XMH)

      IMPLICIT NONE

      INTEGER NH,NHM
      PARAMETER (NHM=8)
      DOUBLE PRECISION XMH,XN(NHM),MH(NHM)

      DATA MH/80d0,90d0,100d0,110d0,120d0,130d0,140d0,150d0/
      DATA XN/9.4d0,10.6d0,10.9d0,14.8d0,15.7d0,13.2d0,10.4d0,8.2d0/

      GAGALHIG=0d0

      IF(XMH.LT.MH(1) .OR. XMH.GT.MH(NHM)) RETURN

      DO 1 NH=1,NHM
       IF(MH(NH).LE.XMH) GOTO 1
       GAGALHIG=XN(NH)+(MH(NH)-XMH)/(MH(NH)-MH(NH-1))
     .    *(XN(NH-1)-XN(NH))
       GOTO 2
 1    CONTINUE
 2    RETURN

      END


      DOUBLE PRECISION FUNCTION TTBBHIG(XMH)

      IMPLICIT NONE

      INTEGER NH,NHM
      PARAMETER (NHM=8)
      DOUBLE PRECISION XMH,XN(NHM),MH(NHM)

      DATA MH/80d0,90d0,100d0,110d0,120d0,130d0,140d0,150d0/
      DATA XN/17.9d0,15d0,14.1d0,12.3d0,12.7d0,13.7d0,11.3d0,10.6d0/

      TTBBHIG=0d0

      IF(XMH.LT.MH(1) .OR. XMH.GT.MH(NHM)) RETURN

      DO 1 NH=1,NHM
       IF(MH(NH).LE.XMH) GOTO 1
       TTBBHIG=XN(NH)+(MH(NH)-XMH)/(MH(NH)-MH(NH-1))
     .    *(XN(NH-1)-XN(NH))*DSQRT(30d0/100d0)
       GOTO 2
 1    CONTINUE
 2    RETURN

      END


      DOUBLE PRECISION FUNCTION TAUTAUHIG(XMH)

      IMPLICIT NONE

      INTEGER NH,NHM
      PARAMETER (NHM=13)
      DOUBLE PRECISION XMH,XN(NHM),MH(NHM)

      DATA MH/100d0,110d0,120d0,130d0,140d0,150d0,200d0,250d0,
     .      300d0,350d0,400d0,450d0,500d0/
      DATA XN/3.7D-2,4.2D-2,4.4D-2,4.5D-2,4.7D-2,4.6D-2,3.1D-2,
     .      2.1D-2,1.3D-2,1.D-2,.8D-2,.7D-2,.6D-2/

      TAUTAUHIG=0d0

      IF(XMH.LT.MH(1) .OR. XMH.GT.MH(NHM)) RETURN

      DO 1 NH=1,NHM
       IF(MH(NH).LE.XMH) GOTO 1
       TAUTAUHIG=XN(NH)+(MH(NH)-XMH)/(MH(NH)-MH(NH-1))
     .   *(XN(NH-1)-XN(NH))
       GOTO 2
 1    CONTINUE
 2    RETURN

      END


      DOUBLE PRECISION FUNCTION ZZ4LHIG(XMH)

      IMPLICIT NONE

      INTEGER NH,NHM
      PARAMETER (NHM=19)
      DOUBLE PRECISION XMH,XN(NHM),MH(NHM)

      DATA MH/100d0,120d0,130d0,140d0,150d0,160d0,170d0,180d0,
     .      190d0,200d0,250d0,275d0,350d0,400d0,500d0,600d0,
     .      700d0,800d0,1000d0/
      DATA XN/2.7d0,5.3d0,13.2d0,22.1d0,27.8d0,9.4d0,5.5d0,20.7d0,
     .      25.1d0,26.1d0,21.6d0,17.6d0,22.7d0,21.6d0,21.5d0,17.1d0,
     .      13.6d0,11.1d0,9.3d0/

      ZZ4LHIG=0d0

      IF(XMH.LT.MH(1) .OR. XMH.GT.MH(NHM)) RETURN

      DO 1 NH=1,NHM
       IF(MH(NH).LE.XMH) GOTO 1
       ZZ4LHIG=XN(NH)+(MH(NH)-XMH)/(MH(NH)-MH(NH-1))
     .   *(XN(NH-1)-XN(NH))
       GOTO 2
 1    CONTINUE
 2    RETURN

      END


      DOUBLE PRECISION FUNCTION WW2LHIG(XMH)

      IMPLICIT NONE

      INTEGER NH,NHM
      PARAMETER (NHM=13)
      DOUBLE PRECISION XMH,XN(NHM),MH(NHM)

      DATA MH/120d0,130d0,140d0,150d0,160d0,170d0,180d0,190d0,
     .      200d0,250d0,300d0,600d0,800d0/
      DATA XN/5.1d0,9.8d0,17.8d0,21.9d0,47d0,34.4d0,24.1d0,19.5d0,
     .      16.9d0,7.9d0,19.4d0,14.2d0,11.3d0/

      WW2LHIG=0d0

      IF(XMH.LT.MH(1) .OR. XMH.GT.MH(NHM)) RETURN

      DO 1 NH=1,NHM
       IF(MH(NH).LE.XMH) GOTO 1
       WW2LHIG=XN(NH)+(MH(NH)-XMH)/(MH(NH)-MH(NH-1))
     .   *(XN(NH-1)-XN(NH))
       GOTO 2
 1    CONTINUE
 2    RETURN

      END


      DOUBLE PRECISION FUNCTION WWTTHIG(XMH)

      IMPLICIT NONE

      INTEGER NH,NHM
      PARAMETER (NHM=5)
      DOUBLE PRECISION XMH,XN(NHM),MH(NHM)

      DATA MH/110d0,120d0,130d0,140d0,150d0/
      DATA XN/6.7d0,10.4d0,10.4d0,8.7d0,4.4d0/

      WWTTHIG=0d0

      IF(XMH.LT.MH(1) .OR. XMH.GT.MH(NHM)) RETURN

      DO 1 NH=1,NHM
       IF(MH(NH).LE.XMH) GOTO 1
       WWTTHIG=XN(NH)+(MH(NH)-XMH)/(MH(NH)-MH(NH-1))
     .   *(XN(NH-1)-XN(NH))
       GOTO 2
 1    CONTINUE
 2    RETURN

      END


      DOUBLE PRECISION FUNCTION WWWWHIG(XMH)

      IMPLICIT NONE

      INTEGER NH,NHM
      PARAMETER (NHM=11)
      DOUBLE PRECISION XMH,XN(NHM),MH(NHM)

      DATA MH/110d0,115d0,120d0,125d0,130d0,140d0,150d0,160d0,
     .      170d0,180d0,190d0/
      DATA XN/2.5d0,5.6d0,9.7d0,15.7d0,20.5d0,18.6d0,26.5d0,34.8d0,
     .      34.8d0,27.8d0,21.5d0/

      WWWWHIG=0d0

      IF(XMH.LT.MH(1) .OR. XMH.GT.MH(NHM)) RETURN

      DO 1 NH=1,NHM
       IF(MH(NH).LE.XMH) GOTO 1
       WWWWHIG=XN(NH)+(MH(NH)-XMH)/(MH(NH)-MH(NH-1))
     .   *(XN(NH-1)-XN(NH))
       GOTO 2
 1    CONTINUE
 2    RETURN

      END


      DOUBLE PRECISION FUNCTION WWINVHIG(XMH)

      IMPLICIT NONE

      INTEGER NH,NHM
      PARAMETER (NHM=7)
      DOUBLE PRECISION XMH,XN(NHM),MH(NHM)

      DATA MH/100d0, 150d0, 200d0, 250d0, 300d0, 350d0, 400d0/
      DATA XN/16d0, 14d0, 13d0, 11d0, 10d0, 8d0, 7d0/

      WWINVHIG=0d0

      IF(XMH.LT.MH(1) .OR. XMH.GT.MH(NHM)) RETURN

      DO 1 NH=1,NHM
       IF(MH(NH).LE.XMH) GOTO 1
       WWINVHIG=XN(NH)+(MH(NH)-XMH)/(MH(NH)-MH(NH-1))
     .   *(XN(NH-1)-XN(NH))
       GOTO 2
 1    CONTINUE
 2    RETURN

      END


      SUBROUTINE CLASSIFY(SM,PM,MH,S,SH,N)

      IMPLICIT NONE
      INTEGER N,I,J,K,L
      DOUBLE PRECISION SM(3),PM(2),MH(5),S(5,*),SH(5,*)

      I=1
      J=1
      K=1
 1    IF(I.LE.3 .AND. J.LE.2)THEN
       IF(SM(I).LE.PM(J))THEN
        MH(K)=SM(I)
        DO L=1,N
         SH(K,L)=S(I,L)
        ENDDO
        I=I+1
       ELSE
        MH(K)=PM(J)
        DO L=1,N
         SH(K,L)=S(J+3,L)
        ENDDO
        J=J+1
       ENDIF
       K=K+1
       GOTO 1
      ENDIF

 2    IF(I.LE.3)THEN
       MH(K)=SM(I)
       DO L=1,N
        SH(K,L)=S(I,L)
       ENDDO
       I=I+1
       K=K+1
       GOTO 2
      ENDIF

 3    IF(J.LE.2)THEN
       MH(K)=PM(J)
       DO L=1,N
        SH(K,L)=S(J+3,L)
       ENDDO
       J=J+1
       K=K+1
       GOTO 3
      ENDIF

      END


      SUBROUTINE COMB(MH,DM,SH,N)

      IMPLICIT NONE
      INTEGER N,I,J,K
      DOUBLE PRECISION MH(5),DM(*),SH(5,*),D(4),S

      DO I=1,N

       DO J=1,4
        D(J)=(MH(J+1)-MH(J))/DM(I)*2d0/(MH(J)+MH(J+1))
       ENDDO


       IF(D(1).LE.1d0 .AND. D(2).LE.1d0 .AND. D(3).LE.1d0
     .  .AND. D(4).LE.1d0)THEN
        S=0d0
        DO K=1,5
         S=S+SH(K,I)
        ENDDO
        DO K=1,5
         SH(K,I)=S
        ENDDO
        GOTO 1
       ENDIF

       IF(D(1).LE.1d0 .AND. D(2).LE.1d0 .AND. D(3).LE.1d0)THEN
        S=0d0
        DO K=1,4
         S=S+SH(K,I)
        ENDDO
        DO K=1,4
         SH(K,I)=S
        ENDDO
        GOTO 1
       ENDIF

       IF(D(2).LE.1d0 .AND. D(3).LE.1d0 .AND. D(4).LE.1d0)THEN
        S=0d0
        DO K=2,5
         S=S+SH(K,I)
        ENDDO
        DO K=2,5
         SH(K,I)=S
        ENDDO
        GOTO 1
       ENDIF

       IF(D(1).LE.1d0 .AND. D(2).LE.1d0)THEN
        S=0d0
        DO K=1,3
         S=S+SH(K,I)
        ENDDO
        DO K=1,3
         SH(K,I)=S
        ENDDO
        IF(D(4).LE.1d0)THEN
         S=0d0
         DO K=4,5
          S=S+SH(K,I)
         ENDDO
         DO K=4,5
          SH(K,I)=S
         ENDDO
        ENDIF
        GOTO 1
       ENDIF

       IF(D(2).LE.1d0 .AND. D(3).LE.1d0)THEN
        S=0d0
        DO K=2,4
         S=S+SH(K,I)
        ENDDO
        DO K=2,4
         SH(K,I)=S
        ENDDO
        GOTO 1
       ENDIF

       IF(D(3).LE.1d0 .AND. D(4).LE.1d0)THEN
        S=0d0
        DO K=3,5
         S=S+SH(K,I)
        ENDDO
        DO K=3,5
         SH(K,I)=S
        ENDDO
        IF(D(1).LE.1d0)THEN
         S=0d0
         DO K=1,2
          S=S+SH(K,I)
         ENDDO
         DO K=1,2
          SH(K,I)=S
         ENDDO
        ENDIF
        GOTO 1
       ENDIF

       IF(D(1).LE.1d0)THEN
        S=0d0
        DO K=1,2
         S=S+SH(K,I)
        ENDDO
        DO K=1,2
         SH(K,I)=S
        ENDDO
        IF(D(3).LE.1d0)THEN
         S=0d0
         DO K=3,4
          S=S+SH(K,I)
         ENDDO
         DO K=3,4
          SH(K,I)=S
         ENDDO
        ENDIF
        IF(D(4).LE.1d0)THEN
         S=0d0
         DO K=4,5
          S=S+SH(K,I)
         ENDDO
         DO K=4,5
          SH(K,I)=S
         ENDDO
        ENDIF
        GOTO 1
       ENDIF

       IF(D(2).LE.1d0)THEN
        S=0d0
        DO K=2,3
         S=S+SH(K,I)
        ENDDO
        DO K=2,3
         SH(K,I)=S
        ENDDO
        IF(D(4).LE.1d0)THEN
         S=0d0
         DO K=4,5
          S=S+SH(K,I)
         ENDDO
         DO K=4,5
          SH(K,I)=S
         ENDDO
        ENDIF
        GOTO 1
       ENDIF

       IF(D(3).LE.1d0)THEN
        S=0d0
        DO K=3,4
         S=S+SH(K,I)
        ENDDO
        DO K=3,4
         SH(K,I)=S
        ENDDO
        GOTO 1
       ENDIF

       IF(D(4).LE.1d0)THEN
        S=0d0
        DO K=4,5
         S=S+SH(K,I)
        ENDDO
        DO K=4,5
         SH(K,I)=S
        ENDDO
        GOTO 1
       ENDIF

 1    CONTINUE
      ENDDO

      END


      SUBROUTINE DECLASSIFY(SM,PM,MH,S,SH,N)

      IMPLICIT NONE
      INTEGER N,I,J,K,L
      DOUBLE PRECISION SM(3),PM(2),MH(5),S(5,*),SH(5,*)

      I=1
      J=1
      K=1
 1    IF(I.LE.3 .AND. J.LE.2)THEN
       IF(SM(I).EQ.MH(K))THEN
        DO L=1,N
         S(I,L)=SH(K,L)
        ENDDO
        I=I+1
       ELSE
        DO L=1,N
         S(J+3,L)=SH(K,L)
        ENDDO
        J=J+1
       ENDIF
       K=K+1
       GOTO 1
      ENDIF

 2    IF(I.LE.3)THEN
       DO L=1,N
        S(I,L)=SH(K,L)
       ENDDO
       I=I+1
       K=K+1
       GOTO 2
      ENDIF

 3    IF(J.LE.2)THEN
       DO L=1,N
        S(J+3,L)=SH(K,L)
       ENDDO
       J=J+1
       K=K+1
       GOTO 3
      ENDIF

      END
