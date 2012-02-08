      PROGRAM GRIDNAPR
C----- -----------------------------------------------
      IMPLICIT DOUBLE PRECISION( A-H,O-Z)
C-----------------------------------------------------
      INCLUDE 'PARAMETR.INC'
      INCLUDE 'CONST.INC'
      INCLUDE 'COMMON.INC'
      CHARACTER*4 NUMNAME
 102  FORMAT (33E15.8)
 101  FORMAT ('   N=',I4)
 900  FORMAT ('   SUM=',e12.5)
 104  FORMAT (5E12.5)
 105  FORMAT (E23.16)
 106  FORMAT (I4)
 107  FORMAT (I4, E23.16)
 108  FORMAT (I4,2X,F6.4,2X,F10.4,2(2X,E12.5),2(2X,F10.6))
 560  FORMAT (A)
C
C--------FUNCTION-----------------------------------------------------
C
      FK(S)=-PR*(FKW(S)/MUW+FKO(S)/MUO)
C
      OPEN(2,FILE='NPRINT.DAT')
      READ(2,*)NPRINT
      WRITE(*,*)'     NPRINT='
      WRITE(*,*)NPRINT
      READ(2,*)JPR
      NNNPR=JPR
      WRITE(*,*)'     JPR=',JPR
      READ(2,*)TPR
      WRITE(*,*)'     TPR=',TPR
      CLOSE(2)
      AK=PR/(M*MUW)
      NIS1=NI-1
      NIS2=NI-2
      NJS1=NJ-1
      NJS2=NJ-2
      NPR=0
C
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
C     HX=(ALEN/2)/NIH
C     HY=(DSQRT(3.D0)/2*ALEN)/NJH
      HX=7.D0
      HY=7.D0
      WRITE(*,*) '    HX=',HX
      WRITE(*,*) '    HY=',HY
C
C-----------INIT---------------------------------
C
      PA=PRESS
      DO 29 J=1,NJ
      DO 29 I=1,NI
         Q(I,J)=0.D0
         FW0(I,J)=0.D0
 29   CONTINUE
C
      DO 30 J=1,3
            FW0(2,(J-1)*2*NJH+2)=FW(SZ)
            FW0(2*NIH+2,(J-1)*2*NJH+2)=FW(SZ)
            FW0(4*NIH+2,(J-1)*2*NJH+2)=FW(SZ)
 30   CONTINUE
C
      QHALF=5.D-1*QN
      QDOUBLE=2.D0*QN
      Q(2,2)=QHALF
      Q(NIH*2+2,2)=QN
      Q(NIH*4+2,2)=QHALF
      Q(NIH+2,NJH+2)=-QDOUBLE
      Q(NIH*3+2,NJH+2)=-QDOUBLE
      Q(2,2*NJH+2)=QN
      Q(2*NIH+2,2*NJH+2)=QDOUBLE
      Q(4*NIH+2,2*NJH+2)=QN
      Q(NIH+2,3*NJH+2)=-QDOUBLE
      Q(3*NIH+2,3*NJH+2)=-QDOUBLE
      Q(2,NJH*4+2)=QHALF
      Q(NIH*2+2,NJH*4+2)=QN
      Q(NIH*4+2,NJH*4+2)=QHALF
c
      VOIL=HX*HY*H*NIS2*NJS2*M*9.D-1
c
      WRITE(*,*) 'THE NUMBER OF N = ? (0 OR MORE THAN 0):'
      READ(*,*) N
      IF(N.GT.0) THEN
         OPEN(4,FILE='GRIDNAPR.DAT')
         READ(4,106) N
         WRITE(*,*)'   N=',N
         READ(4,106) NNNPR
         WRITE(*,*)'   NNNPR=',NNNPR
         READ(4,105) HT
         WRITE(*,*)'   HT=',HT
         READ(4,105)TIME
         WRITE(*,*)'   TIME=',TIME
         READ(4,106)NYEAR
         WRITE(*,*)'   NYEAR=',NYEAR
         READ(4,105)OILOUT
         WRITE(*,*)'   OILOUT=',OILOUT
         WRITE(*,*)'----------'
         READ(4,105) ((SWOLD(I,J),I=1,NI),J=1,NJ)
         READ(4,105) ((P(I,J),I=1,NI),J=1,NJ)
         READ(4,105) ((ALF(I,J),I=1,NI),J=1,NJ)
         READ(4,105) ((BET(I,J),I=1,NI),J=1,NJ)
         READ(4,105) ((GAM(I,J),I=1,NI),J=1,NJ)
         READ(4,105) ((DEL(I,J),I=1,NI),J=1,NJ)
         READ(4,105) ((ALFDSH(I,J),I=1,NI),J=1,NJ)
         READ(4,105) ((BETDSH(I,J),I=1,NI),J=1,NJ)
         READ(4,105) ((GAMDSH(I,J),I=1,NI),J=1,NJ)
         READ(4,105) ((DELDSH(I,J),I=1,NI),J=1,NJ)
         CLOSE(4)
         GOTO 25
      END IF
C77   OPEN(8,FILE='!S_.BJN',STATUS='NEW',IOSTAT=ISTAT)
C     WRITE(*,*) ISTAT
C     IF(ISTAT.NE.0) THEN
C       OPEN(8,FILE='!S_.BJN',STATUS='OLD',IOSTAT=ISTAT)
C       CLOSE(8,STATUS='DELETE')
C       GOTO 77
C     ENDIF
C     WRITE(8,560)'JAK_ASC0_01_'
C     WRITE(8,560)'16 0.1 0.8'
C     CLOSE(8)
C78   OPEN(7,FILE='!P_.BJN',STATUS='NEW',IOSTAT=ISTAT)
C     WRITE(*,*) ISTAT
C     IF(ISTAT.NE.0) THEN
C       OPEN(7,FILE='!P_.BJN',STATUS='OLD',IOSTAT=ISTAT)
C       CLOSE(7,STATUS='DELETE')
C       GOTO 78
C     ENDIF
C     WRITE(7,560)'JAK_ASC0_01_'
C     WRITE(7,560)'16 80 130'
C     CLOSE(7)
      HT=1.D-1
      TIME=HT
      NYEAR=0
      OILOUT=0.D0
c      OPEN(13,FILE='bal.dat',ACCESS='APPEND')
c      WRITE(13,*) '  N    HT        TIME       SW   SW_O    BAL    OILK'
c      CLOSE(13)
      DO 27 J=1,NJ
      DO 27 I=1,NI
         ALF(I,J)=0.D0
         BET(I,J)=0.D0
         GAM(I,J)=0.D0
         DEL(I,J)=0.D0
         ALFDSH(I,J)=0.D0
         BETDSH(I,J)=0.D0
         GAMDSH(I,J)=0.D0
         DELDSH(I,J)=0.D0
 27   CONTINUE
      DO 28 J=1,NJ
      DO 28 I=1,NI
         P(I,J)=PA
         SWOLD(I,J)=SSV
 28   CONTINUE
      DO 100 J=2,NJS1
         ALF(2,J)=1.D0
         BET(2,J)=0.D0
         GAM(NIS1,J)=1.D0
         DEL(NIS1,J)=0.D0
 100  CONTINUE
      DO 200 I=2,NIS1
         ALFDSH(I,2)=1.D0
         BETDSH(I,2)=0.D0
         GAMDSH(I,NJS1)=1.D0
         DELDSH(I,NJS1)=0.D0
 200  CONTINUE
      ALF(2,2)=0.D0
      BET(2,2)=PA
      ALF(2,NJS1)=0.D0
      BET(2,NJS1)=PA
      ALFDSH(2,2)=0.D0
      ALFDSH(NIS1,2)=0.D0
      BETDSH(2,2)=PA
      BETDSH(NIS1,2)=PA
      GAM(NI-1,2)=0.D0
      DEL(NI-1,2)=PA
      GAM(NI-1,NJS1)=0.D0
      DEL(NI-1,NJS1)=PA
      GAMDSH(2,NJ-1)=0.D0
      GAMDSH(NIS1,NJ-1)=0.D0
      DELDSH(2,NJ-1)=PA
      DELDSH(NIS1,NJ-1)=PA
C
 25   CONTINUE
      N=N+1
c*************************************
      IF(N.EQ.NNNPR) WRITE(*,101) N
      IF(N.EQ.NNNPR)WRITE(*,*) 'HT=',HT
      IF(N.EQ.NNNPR) WRITE(*,*) '   TIME=',TIME
c*************************************
C
      DO 51 J=1,2
            FW0(NIH+2,(2*J-1)*NJH+2)=FW(SWOLD(NIH+2,(2*J-1)*NJH+2))
            FW0(NIH*3+2,(2*J-1)*NJH+2)=FW(SWOLD(NIH*3+2,(2*J-1)*NJH+2))
 51   CONTINUE
      DO 50 J=2,NJS1
        DO 55 I=2,NIS1
          DPDX1=(FK(SWOLD(I+1,J))+FK(SWOLD(I,J)))/2.D0*
     *          (P(I+1,J)-P(I,J))/HX
          DPDX2=(FK(SWOLD(I,J))+FK(SWOLD(I-1,J)))/2.D0*
     *          (P(I,J)-P(I-1,J))/HX
          DPDY1=(FK(SWOLD(I,J+1))+FK(SWOLD(I,J)))/2.D0*
     *          (P(I,J+1)-P(I,J))/HY
          DPDY2=(FK(SWOLD(I,J))+FK(SWOLD(I,J-1)))/2.D0*
     *          (P(I,J)-P(I,J-1))/HY
          A1=(((DPDX1+ABS(DPDX1))/2.D0-(DPDX2-ABS(DPDX2))/2.D0)*
     *          FW(SWOLD(I,J))-(DPDX2+ABS(DPDX2))/2.D0*
     *          FW(SWOLD(I-1,J))+(DPDX1-ABS(DPDX1))/
     *          2.D0*FW(SWOLD(I+1,J)))/HX
           A2=(((DPDY1+ABS(DPDY1))/2.D0-(DPDY2-ABS(DPDY2))/2.D0)*
     *          FW(SWOLD(I,J))-(DPDY2+ABS(DPDY2))/2.D0*
     *          FW(SWOLD(I,J-1))+(DPDY1-ABS(DPDY1))/
     *          2.D0*FW(SWOLD(I,J+1)))/HX
           SWNEW(I,J)=HT/M*((Q(I,J)*FW0(I,J))/(H*HX*HY)-(A1+A2))+
     *               SWOLD(I,J)
 55     CONTINUE
        SWNEW(1,J)=SWNEW(2,J)
        SWNEW(NI,J)=SWNEW(NIS1,J)
 50   CONTINUE
      DO 60 I=1,NI
         SWNEW(I,1)=SWNEW(I,2)
         SWNEW(I,NJ)=SWNEW(I,NJS1)
 60   CONTINUE
c**************************
      IF(N.EQ.NNNPR)WRITE(*,*) 'SWNEW=',SWNEW(2*NIH+2,2*NJH+2)
c**************************
      DO 3 J=2,NJS1
      DO 4 I=2,NIS1
            A(I,J)=(FK(SWNEW(I,J))+FK(SWNEW(I-1,J)))/(2.D0*HX**2)
            B(I,J)=(FK(SWNEW(I+1,J))+FK(SWNEW(I,J)))/(2.D0*HX**2)
            ADASH(I,J)=(FK(SWNEW(I,J))+FK(SWNEW(I,J-1)))/(2.D0*HY**2)
            BDASH(I,J)=(FK(SWNEW(I,J+1))+FK(SWNEW(I,J)))/(2.D0*HY**2)
         C(I,J)=A(I,J)+B(I,J)+ADASH(I,J)+BDASH(I,J)
         F(I,J)=-Q(I,J)/(H*HX*HY)
 4    CONTINUE
 3    CONTINUE
C----------------------------------------------------
      CALL ALBET
C----------------------------------------------------
c      IF(N.EQ.NNNPR)CALL WRITED('S_GRID.BJN',NI,NJ,SWNEW)
c      IF(N.EQ.NNNPR)CALL WRITED('P_GRID.BJN',NI,NJ,P)
c      IF(N.NE.NNNPR) GOTO 800
c      OPEN(10,FILE='SJ22GRID.bjn')
c      DO 777 I=1,NI
c         WRITE(10,107) I,SWNEW(I,2*NJH+2)
c 777  CONTINUE
c      CLOSE(10)
c      OPEN(11,FILE='SI22GRID.bjn')
c      DO 778 J=1,NJ
c         WRITE(11,107) J,SWNEW(2*NIH+2,J)
c 778  CONTINUE
c      CLOSE(11)
 800  continue
      OILIN=0.D0
      DO 151 J=2,NJS1
      DO 151 I=2,NIS1
         OILIN=OILIN+HX*HY*(1.D0-SWNEW(I,J))*H*M
 151  CONTINUE
      DO 152 J=2,NJS1
      DO 153 I=2,NIS1
      IF(Q(I,J).GE.0) GOTO 153
         OILOUT=OILOUT-(1.D0-FW(SWNEW(I,J)))*Q(I,J)*HT
 153  CONTINUE
 152  CONTINUE
      OILCONTR=VOIL-OILIN-OILOUT
      BAL=OILCONTR/VOIL*100
      OILK=OILOUT/VOIL
c      IF(N.EQ.NNNPR) THEN
c        WRITE(*,*) '  VOIL =',VOIL
c        WRITE(*,*) '  OILIN =',OILIN
c        WRITE(*,*) '  OILOUT =',OILOUT
c        WRITE(*,*) '  OILCONTR =',OILCONTR
c        WRITE(*,*) '  BAL =',BAL
c        WRITE(*,*) '  OILK =',OILK
c        OPEN(13,FILE='bal.dat',ACCESS='APPEND')
c      WRITE(13,108) N,HT,TIME,SWNEW(2*NIH+2,2*NJH+2),
c     *              SWNEW(3*NIH+2,3*NJH+2),BAL,OILK
c        CLOSE(13)
c      ENDIF
c      IF((TIME/365).GE.NYEAR) THEN
c         NEX=NYEAR
c         NUMNAME=' '
c         NCOUNT=0
c 114      IF(NEX.NE.0) THEN
c            NFL=MOD(NEX,10)+48
c            NUMNAME=CHAR(NFL)//NUMNAME
c            NCOUNT=NCOUNT+1
c            NEX=NEX/10
c            GOTO 114
c         ENDIF
c        CALL WRITED('S_'//NUMNAME(1:NCOUNT)//'.BJN',NI,NJ,SWNEW)
c        CALL WRITED('P_'//NUMNAME(1:NCOUNT)//'.BJN',NI,NJ,P)
c        OPEN(17,FILE='SJ22_'//NUMNAME(1:NCOUNT)//'.bjn')
c        DO 877 I=1,NI
c           WRITE(17,107) I,SWNEW(I,2*NJH+2)
c 877    CONTINUE
c        CLOSE(17)
c        OPEN(18,FILE='SI22_'//NUMNAME(1:NCOUNT)//'.bjn')
c        DO 878 J=1,NJ
c           WRITE(18,107) J,SWNEW(2*NIH+2,J)
c 878    CONTINUE
c        CLOSE(18)
c        NYEAR=NYEAR+1
c      ENDIF

      IF(MOD(N,200).EQ.0) THEN
         NEX=N
         NUMNAME=' '
         NCOUNT=0
 114      IF(NEX.NE.0) THEN
            NFL=MOD(NEX,10)+48
            NUMNAME=CHAR(NFL)//NUMNAME
            NCOUNT=NCOUNT+1
            NEX=NEX/10
            GOTO 114
         ENDIF
C
        OPEN(17,FILE='SG'//NUMNAME(1:NCOUNT)//'.dat')
        REWIND (17)
        WRITE(17, * ) 'TITLE = "WATER SATURATION, N = ', N ,
     *                ' TIME = ', TIME, '"'
        WRITE  (17, * ) 'VARIABLES = "X", "Y", "SWNEW"'
        WRITE  (17, * )
     *   'ZONE T ="BIG ZONE", I =', NI, ' J =', NJ,  ' F = POINT'
      DO 630 J = 1, NJ
      DO 631 I = 1, NI
         WRITE (17,*)
     *         REAL((I-1)*HX), REAL((J-1)*HY), REAL(SWNEW(I,J))
  631 CONTINUE
  630 CONTINUE
      CLOSE  (17)
      ENDIF

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
      HT=TPR*DMIN1(HX,HY)/(AK*AMAXK*AMAXP)
      TIME=TIME+HT
      IF(N.EQ.NNNPR)NNNPR=NNNPR+JPR
      NPR=NPR+1
      IF(NPR.LT.NPRINT)GOTO 850
        OPEN(4,FILE='GRIDNAPR.DAT')
          WRITE(4,106) N
          WRITE(4,106) NNNPR
          WRITE(4,105) HT
          WRITE(4,105)TIME
          WRITE(4,106)NYEAR
          WRITE(4,105)OILOUT
          WRITE(4,105) ((SWNEW(I,J),I=1,NI),J=1,NJ)
          WRITE(4,105) ((P(I,J),I=1,NI),J=1,NJ)
          WRITE(4,105) ((ALF(I,J),I=1,NI),J=1,NJ)
          WRITE(4,105) ((BET(I,J),I=1,NI),J=1,NJ)
          WRITE(4,105) ((GAM(I,J),I=1,NI),J=1,NJ)
          WRITE(4,105) ((DEL(I,J),I=1,NI),J=1,NJ)
          WRITE(4,105) ((ALFDSH(I,J),I=1,NI),J=1,NJ)
          WRITE(4,105) ((BETDSH(I,J),I=1,NI),J=1,NJ)
          WRITE(4,105) ((GAMDSH(I,J),I=1,NI),J=1,NJ)
          WRITE(4,105) ((DELDSH(I,J),I=1,NI),J=1,NJ)
        CLOSE(4)
        NPR=0
 850  CONTINUE
      DO 150 J=1,NJ
      DO 150 I=1,NI
         SWOLD(I,J)=SWNEW(I,J)
 150  CONTINUE
      WRITE(*,*) 'THE STEP IS OVER'
      IF(N.EQ.NNNPR)WRITE(*,*) '--------------------------------------'
      IF(N.LE.10200)  GOTO 25
c      IF(NYEAR.LT.7) GOTO 25
      STOP
      END
C
      FUNCTION FKW(S)
C------------------------------------------------------
      IMPLICIT DOUBLE PRECISION( A-H,O-Z)
C------------------------------------------------------
      INCLUDE 'PARAMETR.INC'
      INCLUDE 'CONST.INC'
      INCLUDE 'COMMON.INC'
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
      INCLUDE 'PARAMETR.INC'
      INCLUDE 'CONST.INC'
      INCLUDE 'COMMON.INC'
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
      INCLUDE 'PARAMETR.INC'
      INCLUDE 'CONST.INC'
      INCLUDE 'COMMON.INC'
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
      INCLUDE 'PARAMETR.INC'
      INCLUDE 'CONST.INC'
      INCLUDE 'COMMON.INC'
C
      FW=(FKW(S)/MUW)/((FKW(S)/MUW)+(FKO(S)/MUO))
      RETURN
      END
C
      FUNCTION FO(S)
C------------------------------------------------------
      IMPLICIT DOUBLE PRECISION( A-H,O-Z)
C------------------------------------------------------
      INCLUDE 'PARAMETR.INC'
      INCLUDE 'CONST.INC'
      INCLUDE 'COMMON.INC'
C
      FO=(FKO(S)/MUO)/((FKW(S)/MUW)+(FKO(S)/MUO))
      RETURN
      END
C
      SUBROUTINE WRITED(FLNAME,N,M,Z)
      INCLUDE 'PARAMETR.INC'
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
       SUBROUTINE ALBET
C-----------------------------------------------------
       IMPLICIT DOUBLE PRECISION( A-H,O-Z)
C-----------------------------------------------------
      INCLUDE 'PARAMETR.INC'
      INCLUDE 'COMMON.INC'
C
 100  FORMAT(5X,' ITRGAM=',I3,' DGAM=',E17.11)
 110  FORMAT(5X,' ITRDEL=',I3,' RES=',E17.11)
 111  FORMAT(5X,'ITERATIONS ALFA IN',I3,'ITERATIONS ARE NOT CONVERGED',
     *          'ITRGAM=',I3,' DGAM=',E17.11,' EPSGAM=',E17.11)
 112  FORMAT(5X,'ITERATIONS BETA IN',I3,'ITERATIONS ARE NOT CONVERGED',
     *          'ITRDEL=',I3,' D=',E17.11)
C
C
      open(3,FILE='omega.dat')
      READ(3,500)OMEGA
 500  FORMAT(e12.5)
C      WRITE(*,*) 'OMEGA=',OMEGA
      CLOSE(3)
      NIS1=NI-1
      NIS2=NI-2
      NJS1=NJ-1
      NJS2=NJ-2
      MAXGAM=20
      EPSGAM=1.D-10
      MAXDEL=1000
      EPS=1.D-4
C
      DO 600 I=2,NIS1
      DO 601 J=2,NJS1
         WORK1=A(I,J)*P(I-1,J)-C(I,J)*P(I,J)+B(I,J)*P(I+1,J)
     *         +ADASH(I,J)*P(I,J-1)+BDASH(I,J)*P(I,J+1)+F(I,J)
         D1=D1+WORK1*WORK1*HX*HY
 601  CONTINUE
 600  CONTINUE
      D1=DSQRT(ABS(D1))
C--------ALFA-ITERATIONS-----------------------------------------------
      ITRGAM=0
 10   CONTINUE
C
C-------ITERATIONS UP--------------------------------------------------
      DO 14 J=2,NJS1
         CALL ALFJ(J,2,NIS1)
         CALL GAMJ(J,2,NIS1)
      DO 13 I=2,NIS1
         CALL ALFD(I,J)
 13   CONTINUE
 14   CONTINUE
C
C-------ITERATIONS DOWN------------------------------------------------
      DGAM=0.D0
      DO 18 K=1,NJS2
         J=NJ-K
         CALL ALFJ(J,2,NIS1)
         CALL GAMJ(J,2,NIS1)
C
      DO 17 I=2,NIS1
         GAMOLD=GAMDSH(I,J-1)
         CALL GAMD(I,J)
         WORK=GAMDSH(I,J-1)-GAMOLD
         DGAM=DGAM+HX*HY*WORK*WORK
 17   CONTINUE
 18   CONTINUE
      ITRGAM=ITRGAM+1
      DGAM=DSQRT(ABS(DGAM))
      IF(ITRGAM.GT.MAXGAM) GO TO 21
      IF(DGAM.LT.EPSGAM) GO TO 19
      INDIC=0
      GO TO 10
 21   INDIC=2
      GO TO 41
 19   INDIC=1
C      WRITE(*,*) 'END OF ALFA'
C      WRITE(*,*) 'BETA'
C---------BETA-ITERATIONS----------------------------------------------
C
      ITRDEL=0
 20   CONTINUE
C--------ITERATIONS UP-------------------------------------------------
      DO 24 J=2,NJS1
         CALL BETJ(J,2,NIS1)
         CALL DELJ(J,2,NIS1)
C
      DO 23 I=2,NIS1
         CALL BETD(I,J)
 23   CONTINUE
 24   CONTINUE
C--------ITERATIONS DOWN-----------------------------------------------
      DO 28 K=1,NJS2
         J=NJ-K
         CALL BETJ(J,2,NIS1)
         CALL DELJ(J,2,NIS1)
      DO 27 I=2,NIS1
         CALL DELD(I,J)
 27   CONTINUE
 28   CONTINUE
      ITRDEL=ITRDEL+1
      DO 30 I=2,NIS1
         P(I,1)=(ALFDSH(I,2)*DELDSH(I,1)+BETDSH(I,2))/
     *              (1.D0-ALFDSH(I,2)*GAMDSH(I,1))
 30   CONTINUE
      DO 31 J=2,NJS1
      DO 32 I=2,NIS1
         P(I,J)=GAMDSH(I,J-1)*P(I,J-1)+DELDSH(I,J-1)
 32   CONTINUE
 31   CONTINUE
      DO 827 I=2,NIS1
c        P(I,1)=P(I,2)
         P(I,NJ)=P(I,NJS1)
 827  CONTINUE
      DO 821 J=1,NJ
         P(1,J)=P(2,J)
         P(NI,J)=P(NIS1,J)
 821  CONTINUE
c     P(NIH*3+2,1)=P(NIH*3+2,NJ)
      P(2,1)=PA
      P(2,NJ)=PA
      P(1,2)=PA
      P(NI,2)=PA
      P(NI-1,1)=PA
      P(NI-1,NJ)=PA
      P(1,NJ-1)=PA
      P(NI,NJ-1)=PA
C      P(NIH*3+2,1)=PA
C      P(1,NJH+2)=PA
C      P(1,3*NJH+2)=PA
C      P(NI,NJH+2)=PA
C      P(NI,3*NJH+2)=PA
C      P(3*NIH+2,NJ)=PA
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
      IF(N.EQ.NNNPR) WRITE(*,110)ITRDEL,RES
      IF(ITRDEL.GT.MAXDEL)GO TO 35
      IF(RES.LT.EPS)GO TO 36
      INDIC=0
      GO TO 20
 35   INDIC=2
      GO TO 42
 36   INDIC=1
      GO TO 43
 41   WRITE(*,111)MAXGAM,ITRGAM,DGAM,EPSGAM
c     WRITE(2,111)MAXGAM,ITRGAM,DGAM,EPSGAM
      STOP
 42   WRITE(*,112)MAXDEL,ITRDEL,D
c     WRITE(2,112)MAXDEL,ITRDEL,D
 43   CONTINUE
c     CALL APRINT('P ',P,NI,1,NI,1,1,NJ,1)

      RETURN
      END



      SUBROUTINE ALFA(I,J)
C-----------------------------------------------------
      IMPLICIT DOUBLE PRECISION( A-H,O-Z)
C-----------------------------------------------------
      INCLUDE 'PARAMETR.INC'
      INCLUDE 'COMMON.INC'
C
         ALF(I+1,J)=B(I,J)/(C(I,J)-ALF(I,J)*A(I,J)
     *              -ALFDSH(I,J)*ADASH(I,J)-GAMDSH(I,J)*BDASH(I,J))
      RETURN
      END


      SUBROUTINE BETA(I,J)
C-----------------------------------------------------
      IMPLICIT DOUBLE PRECISION( A-H,O-Z)
C-----------------------------------------------------
      INCLUDE 'PARAMETR.INC'
      INCLUDE 'COMMON.INC'
C
         BET(I+1,J)=(F(I,J)+BET(I,J)*A(I,J)+BETDSH(I,J)*ADASH(I,J)
     *              +DELDSH(I,J)*BDASH(I,J))*ALF(I+1,J)/B(I,J)
      RETURN
      END


      SUBROUTINE GAMMA(I,J)
C-----------------------------------------------------
      IMPLICIT DOUBLE PRECISION( A-H,O-Z)
C-----------------------------------------------------
      INCLUDE 'PARAMETR.INC'
      INCLUDE 'COMMON.INC'
C
         GAM(I-1,J)=A(I,J)/(C(I,J)-GAM(I,J)*B(I,J)
     *              -ALFDSH(I,J)*ADASH(I,J)-GAMDSH(I,J)*BDASH(I,J))
      RETURN
      END


      SUBROUTINE DELTA(I,J)
C-----------------------------------------------------
      IMPLICIT DOUBLE PRECISION( A-H,O-Z)
C-----------------------------------------------------
      INCLUDE 'PARAMETR.INC'
      INCLUDE 'COMMON.INC'
C
         DEL(I-1,J)=(F(I,J)+DEL(I,J)*B(I,J)+BETDSH(I,J)*ADASH(I,J)
     *              +DELDSH(I,J)*BDASH(I,J))*GAM(I-1,J)/A(I,J)
      RETURN
      END


      SUBROUTINE ALFD(I,J)
C-----------------------------------------------------
      IMPLICIT DOUBLE PRECISION( A-H,O-Z)
C-----------------------------------------------------
      INCLUDE 'PARAMETR.INC'
      INCLUDE 'COMMON.INC'
C
         ALFDSH(I,J+1)=BDASH(I,J)/(C(I,J)-ALFDSH(I,J)*ADASH(I,J)
     *              -ALF(I,J)*A(I,J)-GAM(I,J)*B(I,J))
      RETURN
      END


      SUBROUTINE BETD(I,J)
C-----------------------------------------------------
      IMPLICIT DOUBLE PRECISION( A-H,O-Z)
C-----------------------------------------------------
      INCLUDE 'PARAMETR.INC'
      INCLUDE 'COMMON.INC'
C
C        BETDSH(I,J+1)=(F(I,J)+BET(I,J)*A(I,J)+BETDSH(I,J)*ADASH(I,J)
C    *              +DEL(I,J)*B(I,J))*ALFDSH(I,J+1)/BDASH(I,J)
C
      BETDSH(I,J+1)=OMEGA*(F(I,J)+BET(I,J)*A(I,J)+BETDSH(I,J)*ADASH(I,J)
     *              +DEL(I,J)*B(I,J))*ALFDSH(I,J+1)/BDASH(I,J)
     *              +(1.D0-OMEGA)*BETDSH(I,J+1)
      RETURN
      END


      SUBROUTINE GAMD(I,J)
C-----------------------------------------------------
      IMPLICIT DOUBLE PRECISION( A-H,O-Z)
C-----------------------------------------------------
      INCLUDE 'PARAMETR.INC'
      INCLUDE 'COMMON.INC'
C
         GAMDSH(I,J-1)=ADASH(I,J)/(C(I,J)-GAMDSH(I,J)*BDASH(I,J)
     *              -ALF(I,J)*A(I,J)-GAM(I,J)*B(I,J))
      RETURN
      END


      SUBROUTINE DELD(I,J)
C-----------------------------------------------------
      IMPLICIT DOUBLE PRECISION( A-H,O-Z)
C-----------------------------------------------------
      INCLUDE 'PARAMETR.INC'
      INCLUDE 'COMMON.INC'
C
C        DELDSH(I,J-1)=(F(I,J)+DELDSH(I,J)*BDASH(I,J)+BET(I,J)*A(I,J)
C    *              +DEL(I,J)*B(I,J))*GAMDSH(I,J-1)/ADASH(I,J)
C
      DELDSH(I,J-1)=OMEGA*(F(I,J)+DELDSH(I,J)*BDASH(I,J)+BET(I,J)*A(I,J)
     *              +DEL(I,J)*B(I,J))*GAMDSH(I,J-1)/ADASH(I,J)
     *              +(1.D0-OMEGA)*DELDSH(I,J-1)
      RETURN
      END

      SUBROUTINE ALFJ(J,I1,I2)
C-----------------------------------------------------------------
            DO 1 I=I1,I2
               CALL ALFA(I,J)
 1    CONTINUE
      RETURN
      END


      SUBROUTINE GAMJ(J,I1,I2)
C---------------------------------------------------------------
         DO 1 L=I1,I2
            I=I2+I1-L
            CALL GAMMA(I,J)
 1    CONTINUE
      RETURN
      END


      SUBROUTINE BETJ(J,I1,I2)
C----------------------------------------------------------------
         DO 1 I=I1,I2
            CALL BETA(I,J)
 1    CONTINUE
      RETURN
      END


      SUBROUTINE DELJ(J,I1,I2)
C---------------------------------------------------------------
         DO 1 L=I1,I2
            I=I2+I1-L
            CALL DELTA(I,J)
 1    CONTINUE
      RETURN
      END

      SUBROUTINE APRINT
     1                  ( NAME ,    A ,    NI ,
     2                    IBEG , IEND , ISTEP ,
     3                    JBEG , JEND , JSTEP )
      IMPLICIT REAL*8 ( A-H , O-Z )
C________________________________________________
C
C_________INPUT: PARAMETERS FOR PRINT
C
      CHARACTER*(*) NAME
      DIMENSION A(NI,1)
C
  1   FORMAT (//,57X,'*** ',A,' ***',/)
  2   FORMAT (1X,(11I11))
  3   FORMAT (1X,I4, 1P11D11.3 )
C
C________________________________________________
C
      JE=JBEG
C___________________________________
  10  JB=JE
      IF ( JB .GE. JEND ) RETURN
C
      JE=JB + 10*JSTEP
      IF ( JE .GT. JEND ) JE=JEND
C
      WRITE(2,1)NAME
      WRITE(2,2)(J,J=JB,JE,JSTEP)
      DO 20 I=IBEG,IEND,ISTEP
      IDASH=IEND - I + IBEG
      WRITE(2,3)IDASH,(A(IDASH,J),J=JB,JE,JSTEP)
  20  CONTINUE
C
      GOTO 10
C
      END
