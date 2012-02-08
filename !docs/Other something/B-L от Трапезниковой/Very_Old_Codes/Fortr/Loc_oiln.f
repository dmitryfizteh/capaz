      PROGRAM LOC_OILN
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C      INCLUDE "parixf.inc"
      INCLUDE 'parametr.inc'
      INCLUDE 'const.inc'
      INCLUDE 'common.inc'
C--------------------------------------------------------------------
C     MY PROCESSOR AND PROCESSOR COUNT
C--------------------------------------------------------------------
C      NP = nprocs()
C      MP = myprocid()
C--------------------------------------------------------------------
C     INITIALIZATION OF PROCESSOR CLOCK
C--------------------------------------------------------------------
C      call SPTIME()
C--------------------------------------------------------------------
C     INITIALIZATION OF TOPOLOGY
C--------------------------------------------------------------------
C      IF(NP.NE.4) THEN
C         WRITE(6,*) 'BAD PROCESSOR COUNT !'
C         STOP
C      ENDIF
C
C      IDTOP=makeclique(0,NP,-1,-1,-1,-1,-1,-1,IDCL,LINK)
C
C      IF(IDTOP.EQ.-1) THEN
C         WRITE(6,*) 'THE TOPOLOGY NOT INSTALLED !'
C         STOP
C      ENDIF
C      IERR=ainit(IDTOP,-1,-1)
C      IF(IERR.NE.0) THEN
C          WRITE(6,*) 'AINIT : ERR =  ', IERR
C          STOP
C      ENDIF
C
C      IDL=IDCL+1
C--------------------------------------------------------------------
C     CALCULATIONS
C--------------------------------------------------------------------
c      IF(IDL.EQ.NP) THEN
c         NJTR=4
c      ELSE
c         NJTR=35
c      ENDIF
c      write(6,*) 'idl = ',IDL,'   begin of calculation'
      NIS1=NI-1
      NIS2=NI-2
      NJS1=NJ-1
      NJS2=NJ-2
      EPS=1.D-4
C     HX=(ALEN/2)/NIH
C     HY=(DSQRT(3.D0)/2*ALEN)/NJH
      HX=5.D0
      HY=5.D0
      HT=1.D-1
      N=0
c-----------initial zero values----------
      DO 26 J=1,NJ
      DO 26 I=1,NI
         A(I,J)=0.D0
         C(I,J)=0.D0
         B(I,J)=0.D0
         ADASH(I,J)=0.D0
         BDASH(I,J)=0.D0
         F(I,J)=0.D0
 26   CONTINUE
C
C-----------INIT---------------------------------
      PA=PRESS
      DO 28 J=1,NJ
      DO 28 I=1,NI
         P(I,J)=PA
         P0(I,J)=PA
         SWOLD(I,J)=SSV
 28   CONTINUE
      DO 208 J=3,NJ-2
      DO 208 I=3,NI-2
         P0(I,J)=0.D0
 208   CONTINUE
c      P0(2,1)=PA
c      P0(2,NJ)=PA
c      P0(1,2)=PA
c      P0(NI,2)=PA
c      P0(NI-1,1)=PA
c      P0(NI-1,NJ)=PA
c      P0(1,NJ-1)=PA
c      P0(NI,NJ-1)=PA
c      P0(2,2)=PA
c      P0(2,NJ-1)=PA
c      P0(NI-1,2)=PA
c      P0(NI-1,NJ-1)=PA
C
C      IF(IDL.EQ.1)
       CALL OIL1
C 1000     ttt=APTIME()
C      WRITE(6,*) 'PROCESSOR ', IDL,':  TIME IS ',ttt
C      IRES=aexit(IDTOP)
      STOP
      END
C
      SUBROUTINE OIL1
C-----------------------------------------------------
      IMPLICIT DOUBLE PRECISION( A-H,O-Z)
C-----------------------------------------------------
C      INCLUDE "parixf.inc"
      INCLUDE 'parametr.inc'
      INCLUDE 'const.inc'
      INCLUDE 'common.inc'
c      DIMENSION BUFF1(NI),BUFF2(NI)
 102  FORMAT (33E15.8)
 101  FORMAT ('   N=',I4)
 900  FORMAT ('   SUM=',e12.5)
 104  FORMAT (5E12.5)
 105  FORMAT (5E15.8)
C 106  FORMAT (60(5E15.8\))
C
C
      DO 29 J=1,NJ
      DO 29 I=1,NI
         Q(I,J)=0.D0
         FW0(I,J)=0.D0
 29   CONTINUE
C
      QHALF=5.D-1*QN
      QDOUBLE=2.D0*QN
      Q(NIH+2,2)=-QHALF
      Q(NIH*3+2,2)=QN
      Q(NIH*5+2,2)=-QHALF
      Q(2,NJH+2)=QN
      Q(NIH*2+2,NJH+2)=-QN
      Q(NIH*4+2,NJH+2)=-QN
      Q(NIH*6+2,NJH+2)=QN
      Q(NIH+2,NJH*2+2)=-QN
      Q(NIH*3+2,NJH*2+2)=QDOUBLE
      Q(NIH*5+2,NJH*2+2)=-QN
      Q(2,NJH*3+2)=QN
      Q(NIH*2+2,NJH*3+2)=-QN
      Q(NIH*4+2,NJH*3+2)=-QN
      Q(NIH*6+2,NJH*3+2)=QN
      Q(NIH+2,NJH*4+2)=-QHALF
      Q(NIH*3+2,NJH*4+2)=QN
      Q(NIH*5+2,NJH*4+2)=-QHALF
C
      DO 30 J=1,NJ
      DO 30 I=1,NI
         IF(Q(I,J).GT.0.D0) FW0(I,J)=FW(SZ)
 30   CONTINUE
C
      DO 333 J=1,NJ
      DO 333 I=1,NI
         QTR(I,J)=Q(I,J)
         FW0TR(I,J)=FW0(I,J)
 333  CONTINUE
C
      TIME=HT
      AK=PR/(M*MUW)
 25   CONTINUE
      N=N+1
      WRITE(6,101) N
C      WRITE(6,*) 'TIME =',TIME
C
      III=1
      DO 700 J=1,NJ
      DO 700 I=1,NI
         SS(I,J)=SWOLD(I,J)
 700  CONTINUE
 702  CONTINUE
C      WRITE(6,*) 'III =',III
C
      DO 51 J=1,NJ
      DO 51 I=1,NI
         IF(QTR(I,J).LT.0.D0) FW0TR(I,J)=FW(SS(I,J))
 51   CONTINUE
C      LEN=NI*8
C         DO 500 I=1,NI
C            BUFF1(I)=SS(I,NJTR)
C            BUFF2(I)=P(I,NJTR)
C 500     CONTINUE
C      CALL send(IDTOP,LINK(IDL+1),BUFF1,LEN)
C      WRITE(6,*) IDL,': SEND SS(I,NJTR)'
C      CALL send(IDTOP,LINK(IDL+1),BUFF2,LEN)
C      WRITE(6,*) IDL,': SEND P(I,NJTR)'
C      CALL recv(IDTOP,LINK(IDL+1),BUFF1,LEN)
C      WRITE(6,*) IDL,': RECV SS(I,NJTR+1)'
C      CALL recv(IDTOP,LINK(IDL+1),BUFF2,LEN)
C      WRITE(6,*) IDL,': RECV P(I,NJTR+1)'
C         DO 501 I=1,NI
C            SS(I,NJTR+1)=BUFF1(I)
C            P(I,NJTR+1)=BUFF2(I)
C 501     CONTINUE
      DO 50 J=2,NJS1
        DO 55 I=2,NIS1
        DPDX1=(FK(SS(I+1,J))+FK(SS(I,J)))/2.D0*
     *        (P(I+1,J)-P(I,J))/HX
        DPDX2=(FK(SS(I,J))+FK(SS(I-1,J)))/2.D0*
     *        (P(I,J)-P(I-1,J))/HX
        DPDY1=(FK(SS(I,J+1))+FK(SS(I,J)))/2.D0*
     *        (P(I,J+1)-P(I,J))/HY
        DPDY2=(FK(SS(I,J))+FK(SS(I,J-1)))/2.D0*
     *        (P(I,J)-P(I,J-1))/HY
           A1=(((DPDX1+ABS(DPDX1))/2.D0-(DPDX2-ABS(DPDX2))/2.D0)*
     *          FW(SS(I,J))-(DPDX2+ABS(DPDX2))/2.D0*
     *          FW(SS(I-1,J))+(DPDX1-ABS(DPDX1))/
     *          2.D0*FW(SS(I+1,J)))/HX
           A2=(((DPDY1+ABS(DPDY1))/2.D0-(DPDY2-ABS(DPDY2))/2.D0)*
     *          FW(SS(I,J))-(DPDY2+ABS(DPDY2))/2.D0*
     *          FW(SS(I,J-1))+(DPDY1-ABS(DPDY1))/
     *          2.D0*FW(SS(I,J+1)))/HX
           IF(III.EQ.1) THEN
C           SWNEW(I,J)=HT/M*((QTR(I,J)*FW0TR(I,J))/(H*HX*HY)-(A1+A2))+
C     *               SWOLD(I,J)
           SWNEW(I,J)=SS(I,J)+1.D-2
           FN(I,J)=M/HT*(SS(I,J)-SWOLD(I,J))-QTR(I,J)/(H*HX*HY)*
     *      FW0TR(I,J)+A1+A2
           GOTO 707
           END IF
 706      CONTINUE
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          IF(ABS(SS(I,J)-SSS(I,J)).LE.1D-4)GOTO 707
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          FN(I,J)=M/HT*(SS(I,J)-SWOLD(I,J))-QTR(I,J)/(H*HX*HY)*
     *               FW0TR(I,J)+
     *                A1+A2
           SWNEW(I,J)=SS(I,J)-(SS(I,J)-SSS(I,J))*FN(I,J)/
     *                (FN(I,J)-FN1(I,J))
 707    CONTINUE
 55     CONTINUE
        SWNEW(1,J)=SWNEW(2,J)
        SWNEW(NI,J)=SWNEW(NIS1,J)
 50   CONTINUE
      DO 60 I=1,NI
         SWNEW(I,1)=SWNEW(I,2)
         SWNEW(I,NJ)=SWNEW(I,NJS1)
 60   CONTINUE
C      WRITE(6,*) IDL,':  SWNEW (III=',III,') '
      SUM=0.D0
      DO 703 J=1,NJ
      DO 703 I=1,NI
         SUM=SUM+(SWNEW(I,J)-SS(I,J))**2*HX*HY
 703  CONTINUE
C      CALL recv(IDTOP,LINK(IDL+1),SUM1,8)
C      WRITE(6,*) IDL,': RECV SUM1'
C      SUM=SUM+SUM1
      SUM=DSQRT(ABS(SUM))
C      WRITE(6,*)III
C      WRITE(6,900)SUM
      IF(SUM.LT.1.D-4)GOTO 705
      INDIC=0
C      CALL send(IDTOP,LINK(IDL+1),INDIC,4)
C      WRITE(6,*) IDL,': SEND INDIC'
      III=III+1
      DO 701 J=1,NJ
      DO 701 I=1,NI
         SSS(I,J)=SS(I,J)
         SS(I,J)=SWNEW(I,J)
         FN1(I,J)=FN(I,J)
 701  CONTINUE
      GOTO 702
 705  CONTINUE
      INDIC=1
C     CALL send(IDTOP,LINK(IDL+1),INDIC,4)
C      WRITE(6,*) IDL,': SEND INDIC'
C      write(6,*) 'swnew=', SWNEW(3*NIH+2,2*NJH+2)
C      WRITE(6,*) '"SWNEW" IS OVER'
C      OPEN(12,FILE='sw_1.res')
C       WRITE(12,*) '   sw ='
C       DO 120 J=1,NJTR
C       DO 120 I=1,NI
C           WRITE(12,*) I,J,SWNEW(I,J)
C 120   CONTINUE
C      CLOSE(12)
      CALL WRITED('oils1.bjn',NI,NJ,SWNEW)

      PI = 3.14159265368
      COSNX = COS( PI/NIS1)
      COSNX2 = COSNX*COSNX
      COSNY = COS( PI/NJS1 )
      COSNY2 = COSNY*COSNY
      COSNXY = COSNX*COSNY

C--------------------------------------------------------
C      CALL recv(IDTOP,LINK(IDL+1),BUFF1,LEN)
C      WRITE(6,*) IDL,': RECV SWNEW(I,NJTR+1)'
C         DO 503 I=1,NI
C            SWNEW(I,NJTR+1)=BUFF1(I)
C 503     CONTINUE
C         DO 502 I=1,NI
C            BUFF1(I)=SWNEW(I,NJTR)
C 502     CONTINUE
C      CALL send(IDTOP,LINK(IDL+1),BUFF1,LEN)
C      WRITE(6,*) IDL,':  SEND SWNEW(I,NJTR)'

      DO 3 J=2,NJS1
         DO 4 I=2,NIS1
            A(I,J)=(FK(SWNEW(I,J))+FK(SWNEW(I-1,J)))/(2.D0*HX**2)
            B(I,J)=(FK(SWNEW(I+1,J))+FK(SWNEW(I,J)))/(2.D0*HX**2)
            ADASH(I,J)=(FK(SWNEW(I,J))+FK(SWNEW(I,J-1)))/(2.D0*HY**2)
            BDASH(I,J)=(FK(SWNEW(I,J+1))+FK(SWNEW(I,J)))/(2.D0*HY**2)
            C(I,J)=A(I,J)+B(I,J)+ADASH(I,J)+BDASH(I,J)
            F(I,J)=-QTR(I,J)/(H*HX*HY)

            CW = A(I,J)/C(I,J)
            CE = B(I,J)/C(I,J)
            CS = ADASH(I,J)/C(I,J)
            CN = BDASH(I,J)/C(I,J)
c            CALL LRPAR
c     *              ( CS, CW, CE, CN,
c     *                COSNX2 , COSNY2 , COSNXY ,
c     *                OMG(I,J) , IFLAG                      )
            CALL LRPAR2
     *              ( CS, CW, CE, CN,
     *                OMG(I,J) , IFLAG                      )
         IF ( IFLAG .NE. 0 ) THEN
             WRITE(6,*) 'OMEGA EVALUATION ERROR!'
             GOTO 3004
         ENDIF
 4       CONTINUE
 3    CONTINUE
      CALL WRITED('omg.bjn',NI,NJ,OMG)
      OPEN(1,FILE='omg.res')
         DO 31 J=2,NJS1
         DO 31 I=2,NIS1
            WRITE(1,*) I,' ',J,'     ',OMG(I,J)
 31      CONTINUE
       CLOSE(1)
c       OPEN(2,FILE='adash.res')
c         DO 41 J=2,NJS1
c         DO 41 I=2,NIS1
c            WRITE(2,*) I,' ',J,'     ',ADASH(I,J)
c 41      CONTINUE
c      CLOSE(2)
      CALL WRITED('a.bjn',NI,NJ,A)
      CALL WRITED('adash.bjn',NI,NJ,ADASH)
      CALL WRITED('b.bjn',NI,NJ,B)
      CALL WRITED('bdash.bjn',NI,NJ,BDASH)
      CALL WRITED('c.bjn',NI,NJ,C)
C      OPEN(12,FILE='omg.res',access='append')
C         write(12,*) '  '
C         write(12,*) 'N =', N
C         write(12,*) '  '
C         write(12,*)OMG(17,20)
C      CLOSE(12)
C
C----------------------------------------------------
      CALL RELAX
C----------------------------------------------------
C
C        NPR=NPR+1
C        IF(NPR.LT.NPRINT)GOTO 850
C         OPEN(4,FILE='OIL.DAT')
C         WRITE(4,*) N
C         WRITE(4,*) ((SWNEW(I,J),I=1,NITR),J=1,NJ)
C         WRITE(4,*) ((P(I,J),I=1,NITR),J=1,NJ)
C         WRITE(4,*) ((ALF(I,J),I=1,NITR),J=1,NJ)
C         WRITE(4,*) ((BET(I,J),I=1,NITR),J=1,NJ)
C         WRITE(4,*) ((GAM(I,J),I=1,NITR),J=1,NJ)
C         WRITE(4,*) ((DEL(I,J),I=1,NITR),J=1,NJ)
C         WRITE(4,*) ((ALFDSH(I,J),I=1,NITR),J=1,NJ)
C         WRITE(4,*) ((BETDSH(I,J),I=1,NITR),J=1,NJ)
C         WRITE(4,*) ((GAMDSH(I,J),I=1,NITR),J=1,NJ)
C         WRITE(4,*) ((DELDSH(I,J),I=1,NITR),J=1,NJ)
C         DO 222 IR=2,NT
C            CALL F77_CHAN_IN_WORD(NIPR,OUPORT)
C            ALLOCATE (PR(NIPR,NJ))
C            LENPR=8*NIPR*NJ
C            DO 223 II=1,10
C               CALL F77_CHAN_IN_MESSAGE(LENPR,PR,INPORT)
C               WRITE(4,*) ((PR(I,J),I=1,NISIZE),J=1,NJ)
C 223     CONTINUE
C         DEALLOCATE (PR)
C 222     CONTINUE
C         CLOSE(4)
C         NPR=0
C 850     CONTINUE
C
C
C      DO 855 J=1,NJ
C      DO 855 I=1,NI
C         SWRES(I,J)=SWNEW(I,J)
C         PRES(I,J)=P(I,J)
C 855  CONTINUE
C      JJ=NJTR
C      DO 224 IR=2,NP
C         CALL recv(IDTOP,LINK(IDL+1),NJPR,4)
C         WRITE(6,*) IDL,': RECV NJPR=',NJPR
C         DO 225 J=1,NJPR
C            CALL recv(IDTOP,LINK(IDL+1),BUFF1,LEN)
C            DO 856 I=1,NI
C               SWRES(I,JJ+J)=BUFF1(I)
C 856        CONTINUE
C 225      CONTINUE
C          DO 226 J=1,NJPR
C            CALL recv(IDTOP,LINK(IDL+1),BUFF1,LEN)
C            DO 857 I=1,NI
C               PRES(I,JJ+J)=BUFF1(I)
C 857        CONTINUE
C 226     CONTINUE
C         JJ=JJ+NJPR
C 224  CONTINUE
      CALL WRITED('oilp1.bjn',NI,NJ,P)
      AMAXK=DFKW(SWNEW(1,1))
      DO 338 J=1,NJ
      DO 338 I=1,NI
         IF(DFKW(SWNEW(I,J)).GT.AMAXK) AMAXK=DFKW(SWNEW(I,J))
 338  CONTINUE
      AMAXPX=(P(2,1)-P(1,1))/HX
      DO 339 J=1,NJ
      DO 339 I=2,NI
         PX=(P(I,J)-P(I-1,J))/HX
         IF(PX.GT.AMAXPX) AMAXPX=PX
 339  CONTINUE
      AMAXPY=(P(1,2)-P(1,1))/HY
      DO 340 J=2,NJ
      DO 340 I=1,NI
         PY=(P(I,J)-P(I,J-1))/HY
         IF(PY.GT.AMAXPY) AMAXPY=PY
 340  CONTINUE
      AMAXP=DMAX1(AMAXPX,AMAXPY)
      HT=DMIN1(HX,HY)/(AK*AMAXK*AMAXP)
C
      DO 150 J=1,NJ
      DO 150 I=1,NI
         SWOLD(I,J)=SWNEW(I,J)
 150  CONTINUE
      WRITE(6,*) 'THE STEP IS OVER'
      IF(N.LT.1) GOTO 25

 3004 RETURN
      END
C
C
      FUNCTION FK(S)
C------------------------------------------------------
      IMPLICIT DOUBLE PRECISION( A-H,O-Z)
C------------------------------------------------------
      INCLUDE 'parametr.inc'
      INCLUDE 'const.inc'
      INCLUDE 'common.inc'
C
      FK=-PR*(FKW(S)/MUW+FKO(S)/MUO)
C
      RETURN
      END
C
      FUNCTION FKW(S)
C------------------------------------------------------
      IMPLICIT DOUBLE PRECISION( A-H,O-Z)
C------------------------------------------------------
      INCLUDE 'parametr.inc'
      INCLUDE 'const.inc'
      INCLUDE 'common.inc'
C
      IF((SSV.LE.S).AND.(S.LT.S1)) THEN
         FKW=((S-SSV)/(SZ-SSV))**N2
      ELSE IF((S1.LE.S).AND.(S.LE.1.D0)) THEN
         FKW=8.D-1*((S-SSV)/(SZ-SSV))**N3
      ELSE IF((0.D0.LE.S).AND.(S.LT.SSV)) THEN
         FKW=0.D0
      ELSE IF(S.LT.0.D0) THEN
         FKW=0.D0
      ELSE IF(S.GT.1.D0) THEN
         FKW=8.D-1*((S-SSV)/(SZ-SSV))**N3
      END IF
      RETURN
      END
C
      FUNCTION DFKW(S)
C------------------------------------------------------
      IMPLICIT DOUBLE PRECISION( A-H,O-Z)
C------------------------------------------------------
      INCLUDE 'parametr.inc'
      INCLUDE 'const.inc'
      INCLUDE 'common.inc'
C
      IF((SSV.LE.S).AND.(S.LT.S1)) THEN
         DFKW=N2*(S-SSV)**(N2-1)/(SZ-SSV)**N2
      ELSE IF((S1.LE.S).AND.(S.LE.1.D0)) THEN
         DFKW=8.D-1*N3*(S-SSV)**(N3-1)/(SZ-SSV)**N3
      ELSE IF((0.D0.LE.S).AND.(S.LT.SSV)) THEN
         DFKW=0.D0
      ELSE IF(S.LT.0.D0) THEN
         DFKW=0.D0
      ELSE IF(S.GT.1.D0) THEN
         DFKW=8.D-1*N3*(S-SSV)**(N3-1)/(SZ-SSV)**N3
      END IF
      RETURN
      END
C
      FUNCTION FKO(S)
C------------------------------------------------------
      IMPLICIT DOUBLE PRECISION( A-H,O-Z)
C------------------------------------------------------
      INCLUDE 'parametr.inc'
      INCLUDE 'const.inc'
      INCLUDE 'common.inc'
C
      IF((SSV.LE.S).AND.(S.LE.SZ)) THEN
         FKO=((SZ-S)/(SZ-SSV))**N1
      ELSE IF((0.D0.LE.S).AND.(S.LT.SSV)) THEN
         FKO=1.D0
      ELSE IF((SZ.LT.S).AND.(S.LE.1.D0)) THEN
         FKO=0.D0
      ELSE IF(S.LT.0.D0) THEN
         FKO=1.D0
      ELSE IF(S.GT.1.D0) THEN
         FKO=0.D0
      END IF
      RETURN
      END
C
      FUNCTION FW(S)
C------------------------------------------------------
      IMPLICIT DOUBLE PRECISION( A-H,O-Z)
C------------------------------------------------------
      INCLUDE 'parametr.inc'
      INCLUDE 'const.inc'
      INCLUDE 'common.inc'
C
      FW=(FKW(S)/MUW)/((FKW(S)/MUW)+(FKO(S)/MUO))
      RETURN
      END
C
      FUNCTION FO(S)
C------------------------------------------------------
      IMPLICIT DOUBLE PRECISION( A-H,O-Z)
C------------------------------------------------------
      INCLUDE 'parametr.inc'
      INCLUDE 'const.inc'
      INCLUDE 'common.inc'
C
      FO=(FKO(S)/MUO)/((FKW(S)/MUW)+(FKO(S)/MUO))
      RETURN
      END
C
      SUBROUTINE WRITED(FLNAME,N,M,Z)
C      INCLUDE "parixf.inc"
      INCLUDE 'parametr.inc'
         INTEGER*2 NN,MM
         REAL*8 Z(N,M)
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         CHARACTER*(*) FLNAME
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         REAL*4 T(NI,NJ)
         OPEN(5,FILE=FLNAME,ACCESS='DIRECT',
     *   FORM='UNFORMATTED',
     *   RECL=12+2*2+4*M*N+12)
         DO 555 I=1,N
         DO 555 J=1,M
            T(I,J)=SNGL(Z(I,J))
 555     CONTINUE
         NN=N
         MM=M
         WRITE(5,REC=1) 'JAK_BIN2_01 ',NN,MM,
     *   ((T(I,J),J=1,M),I=1,N),
     *   'JAK_BIN2_01 '
         CLOSE(5)
         RETURN
         END
C

      SUBROUTINE LRPAR

     1                  ( CS , CW , CE , CN ,
     2                    COSNX2 , COSNY2 , COSNXY ,
     3                    OMEGA , FLAG             )

      IMPLICIT DOUBLE PRECISION ( A-H , O-Z )
C======================================================================C
C     Local relaxation parameter evaluation based on:                  C
C                                                                      C
C     1. E.F.F. Botta, A.E.P. Veldman.                                 C
C        On Local Relaxation Methods and Their Application to          C
C        Convection-Diffusion Equations.                               C
C        J.Comput. Phys., 1982, v.48., N 1, pp.127-149.                C
C     2. L.W. Ehrlich.                                                 C
C        An Ad Hoc SOR Method.                                         C
C        J.Comput. Phys., 1981, v.44., N 1, pp.31-45.                  C
C     3. L.W. Ehrlich.                                                 C
C        The Ad-Hoc SOR Method: a Local Relaxation Scheme.             C
C        Elliptic Problem Solvers II                                   C
C        (Eds. G. Birkhoff and A. Schoenstadt), pp.257-269.            C
C        Academic Press, Orlando, 1984.                                C
C                                                                      C
C     for the solution of the following eqn:                           C
C                                                                      C
C         - CS*Y(i,j-1) - CW*Y(i-1,j) + Y(i,j)                         C
C         - CE*Y(i+1,j) - CN*Y(i,j+1) = F(i,j)                         C
C                                                                      C
C----------------------------------------------------------------------C

      INTEGER FLAG

C     . . . 0 < omega < 1
      OMEGAI( X2 ) = 2.D0/( 1.D0 + SQRT(1.D0+X2) )
C     . . . 1 < omega < 2
      OMEGAR( X2 ) = 2.D0/( 1.D0 + SQRT(1.D0-X2) )
C     . . . 0 < omega < 2
      OMEGAC( XR2 , XI2 ) = 2.D0/( 1.D0 + SQRT(1.D0 - XR2
     *                          + XI2/(1.D0-XR2**0.333333333333333D0)) )
C     . . . limits
      OMEGA0( X2 ) = 1.D-4
      OMEGA2( X2 ) = 2.D0 - 1.D-4
C======================================================================C

      FLAG = 0
      CWCE = CW*CE
      CSCN = CS*CN
      C = CWCE*CSCN
      IF ( C ) 10, 20, 30

C     . . . complex value . . .

  10  CONTINUE
      IF ( CWCE .GT. 0.D0 ) THEN
            RMUR2 =  4.D0*CWCE*COSNX2
            RMUI2 = -4.D0*CSCN*COSNY2
         ELSE
            RMUR2 =  4.D0*CSCN*COSNY2
            RMUI2 = -4.D0*CWCE*COSNX2
      ENDIF
      IF ( RMUR2 .GE. 1.D0 ) THEN
            OMEGA = OMEGA0( 0. )
         ELSE
            OMEGA = OMEGAC( RMUR2 , RMUI2 )
      ENDIF
      RETURN

C     . . . zero coefficients . . .

  20  CONTINUE
      IF ( CWCE ) 21, 22, 26
  21  CONTINUE
         RMUI2 =  -4.D0*CWCE*COSNX2
         OMEGA = OMEGAI( RMUI2 )
         GOTO 27
  22  CONTINUE
         IF ( CSCN ) 23, 24, 25
   23    CONTINUE
            RMUI2 = -4.D0*CSCN*COSNY2
            OMEGA = OMEGAI( RMUI2 )
            GOTO 27
   24    CONTINUE
C           . . . CWCE = CSCN = 0 !!!!!
            FLAG = 1
            GOTO 27
   25    CONTINUE
            RMUR2 =  4.D0*CSCN*COSNY2
            IF ( RMUR2 .GE. 1.D0 ) THEN
               OMEGA = OMEGA2( 2. )
            ELSE
               OMEGA = OMEGAR( RMUR2 )
            ENDIF
            GOTO 27
  26  CONTINUE
         RMUR2 =  4.D0*CWCE*COSNX2
         IF ( RMUR2 .GE. 1.D0 ) THEN
            OMEGA = OMEGA2( 2. )
         ELSE
            OMEGA = OMEGAR( RMUR2 )
         ENDIF
  27  CONTINUE
      RETURN

C     . . . singular case . . .

  30  CONTINUE
      IF ( CWCE .GT. 0.D0 ) THEN
         RMUR2 =  4.D0*
     *            ( CWCE*COSNX2 + 2.D0*SQRT(C)*COSNXY + CSCN*COSNY2 )
         IF ( RMUR2 .GE. 1.D0 ) THEN
            OMEGA = OMEGA2( 2. )
         ELSE
            OMEGA = OMEGAR( RMUR2 )
         ENDIF
      ELSE
         RMUI2 = -4.D0*
     *            ( CWCE*COSNX2 - 2.D0*SQRT(C)*COSNXY + CSCN*COSNY2 )
         OMEGA = OMEGAI( RMUI2 )
      ENDIF
      RETURN

      END
C
      SUBROUTINE LRPAR2

     1                  ( CS , CW , CE , CN ,
     2                    OMEGA , FLAG             )

      IMPLICIT DOUBLE PRECISION ( A-H , O-Z )
C======================================================================C
C     Local relaxation parameter evaluation                            C
C     (Neuman boundary conditions)                                     C
C======================================================================C

      INTEGER FLAG

C     . . . 0 < omega < 1
      OMEGAI( X2 ) = 2.D0/( 1.D0 + SQRT(1.D0+X2) )
C     . . . 1 < omega < 2
      OMEGAR( X2 ) = 2.D0/( 1.D0 + SQRT(1.D0-X2) )
C     . . . 0 < omega < 2
      OMEGAC( XR2 , XI2 ) = 2.D0/( 1.D0 + SQRT(1.D0 - XR2
     *                          + XI2/(1.D0-XR2**0.333333333333333D0)) )
C     . . . limits
      OMEGA0( X2 ) = 1.D-4
      OMEGA2( X2 ) = 2.D0 - 1.D-4
C======================================================================C

      FLAG = 0

C     . . . singular case . . .

  30  CONTINUE
         RMUR2 =  (CE + CW + CN + CS)*(CE + CW + CN + CS)
c         IF ( RMUR2 .GE. 1.D0 ) THEN
c            OMEGA = OMEGA2( 2. )
c         ELSE
            OMEGA = OMEGAR( RMUR2 )
c         ENDIF

      RETURN

      END

      SUBROUTINE RELAX
C------------------------------------------------------
      IMPLICIT DOUBLE PRECISION( A-H,O-Z)
C------------------------------------------------------
      INCLUDE 'parametr.inc'
      INCLUDE 'const.inc'
      INCLUDE 'common.inc'
C
       LLL=1
       D1=0.D0
       DO 611 I=2,NIS1
       DO 612 J=2,NJS1
          WORK1=A(I,J)*P(I-1,J)-C(I,J)*P(I,J)
     *          +B(I,J)*P(I+1,J)+ADASH(I,J)*P(I,J-1)
     *          +BDASH(I,J)*P(I,J+1)+F(I,J)
          D1=D1+WORK1*WORK1*HX*HY
 612   CONTINUE
 611   CONTINUE
       D1=DSQRT(ABS(D1))
C       IF(N.EQ.1) THEN
C          DO 5011 J=1,NJ
C          DO 5011 I=1,NI
C             P(I,J)=P0(I,J)
C 5011     CONTINUE
C       ENDIF
 602   CONTINUE
       DO 501 J=2,NJS1
       DO 551 I=2,NIS1
          IF(MOD((I+J),2).EQ.0) THEN
                OMEGA=OMG(I,J)
c                OMEGA=2.D0 - 1.D-4
                P(I,J)=(1.D0-OMEGA)*P(I,J)+
     *                OMEGA*(A(I,J)*P(I-1,J)+B(I,J)*P(I+1,J)+
     *                ADASH(I,J)*P(I,J-1)+
     *                BDASH(I,J)*P(I,J+1)+F(I,J))/C(I,J)
           ENDIF
 551   CONTINUE
 501   CONTINUE
       DO 502 J=2,NJS1
       DO 553 I=2,NIS1
          IF(MOD((I+J),2).NE.0) THEN
                OMEGA=OMG(I,J)
c                OMEGA=2.D0 - 1.D-4
                P(I,J)=(1.D0-OMEGA)*P(I,J)+
     *                OMEGA*(A(I,J)*P(I-1,J)+B(I,J)*P(I+1,J)+
     *                ADASH(I,J)*P(I,J-1)+
     *                BDASH(I,J)*P(I,J+1)+F(I,J))/C(I,J)
          ENDIF
 553   CONTINUE
 502   CONTINUE
      DO 600 J=2,NJS1
        P(1,J)=P(2,J)
        P(NI,J)=P(NIS1,J)
 600   CONTINUE
      DO 601 I=1,NI
         P(I,1)=P(I,2)
         P(I,NJ)=P(I,NJS1)
 601   CONTINUE
      D=0.D0
      DO 34 J=2,NJS1
      DO 33 I=2,NIS1
         WORK=A(I,J)*P(I-1,J)-C(I,J)*P(I,J)+B(I,J)*P(I+1,J)
     *        +ADASH(I,J)*P(I,J-1)+BDASH(I,J)*P(I,J+1)+F(I,J)
         D=D+WORK*WORK*HX*HY
 33   CONTINUE
 34   CONTINUE
      D=DSQRT(ABS(D))
      RES=D/D1
C      IF(N.EQ.NNNPR)
      WRITE(*,*) '   ITER_P  =',LLL
      WRITE(*,*) '   RESIDULE_P  =',RES
      IF(RES.GT.1.D-4) THEN
         LLL=LLL+1
         GOTO 602
      END IF
      WRITE(*,*) 'P_centr = ', P(59,72)

      RETURN
      END

