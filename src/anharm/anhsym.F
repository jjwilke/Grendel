C     PROGRAM ANHSYM                                                    ANH00010
C**********************************************************************
C***LAST UPDATED ON January 10, 2006 BY DR. YUKIO YAMAGUCHI         ***
C   Modification for a stand-alone version                          ***
C**********************************************************************
C***LAST UPDATED ON April 20, 2004 BY YUKIO YAMAGUCHI               ***
C   Modification                                                    ***
C   for anharmonic zero-point vibrational energy (ZPVE)             ***
C   and vibrational G_0 energy                                      ***
C********************************************************************** ANH00020
C***ERRORS WERE CORRECED FOR X(TT') AND G(TT')                        * ANH00030
C***NOTIFIED BY WALTER THIEL                                          * ANH00040
C***ON OCTOBER 20, 1989 BY YUKIO YAMAGUCHI                            * ANH00050
C********************************************************************** ANH00060
C***LAST UPDATED ON JULY 03, 1986 BY YUKIO YAMAGUCHI                  * ANH00070
C********************************************************************** ANH00080
      SUBROUTINE LTYPSY(ZETA,AXY,F3Q,EE,CLTYP)                          ANH00090
C   THE CALCULATION OF L-TYPE DOUBLING CONSTANTS                        ANH00100
C   FOR A SYMMETRIC TOP MOLECULE                                        ANH00110
      IMPLICIT REAL*8 (A-H,O-Z)                                         ANH00120
      DIMENSION ZETA(N3N,N3N,3),AXY(3,3,N3N),F3Q(N3N,N3N,N3N)           ANH00130
      DIMENSION EE(N3N,N3N),CLTYP(N3N,N3N)                              ANH00140
      COMMON/VIB101/NATOM,N3N,NATRI,ILIN,NVIB                           ANH00150
      COMMON/VIB106/ROTAA(3),ROTGC(3),ROTCM(3),ROTMH(3)                 ANH00160
      COMMON/VIB108/IAXIS(3),NDEG(150),NDAB(150,5),IMAG(150)            ANH00170
      COMMON/VIB109/CLIMIT,FLIM1,FLIM2                                  ANH00180
      COMMON/VIB203/IOFF(150),IPRNT                                     ANH00190
      COMMON/VIB204/SQM(150),ROOT(150),FREQ(150)                        ANH00200
      COMMON/VIB205/IFREQ,NFRQ(150)                                     ANH00210
      DATA PH,AVN / 6.626176D+00 , 6.022045D+00 /                       ANH00220
      DATA PI,CL / 3.1415926536D+00 , 2.99792458D+00 /                  ANH00230
      DATA A00,HALF,ONE,TWO / 0.0D+00 , 0.5D+00 , 1.0D+00 , 2.0D+00 /   ANH00240
      DATA ALIMIT / 1.0D-04 /                                           ANH00250
      DATA WLIMIT / 1.0D-02 /                                           ANH00260
    1 FORMAT(//,2X,' ***********************************'/              ANH00270
     1          2X,' ***THE L-TYPE DOUBLING CONSTANTS***'/              ANH00280
     2          2X,' ***********************************')              ANH00290
    2 FORMAT(//,2X,' ZETA MATRIX, IABC = ',I5/)                         ANH00300
    3 FORMAT(//,2X,' F3(QI*QJ*QK) MATRIX (IN CM-1), IVIB = ',I5/)       ANH00310
    4 FORMAT(//,2X,' AXY MATRIX, IVIB = ',I5/)                          ANH00320
    5 FORMAT(2X,I5,4F10.5)                                              ANH00330
    6 FORMAT(//,2X,' QE(T) VALU'/)                                      ANH00340
    7 FORMAT(2X,' QVIB',5X,' IVIB',7X,' FREQENCY',9X,' HARMONIC',       ANH00290
     1      11X,' ANHARMONIC',9X,' TOTAL',11X,' TOTAL'/                 ANH00300
     2      24X,' IN CM-1',10X,' IN CM-1',12X,' IN CM-1',               ANH00310
     3      12X,' IN CM-1',9X,' IN MHZ'/)                               ANH00320
    8 FORMAT(2X,I5,5X,I5,8X,F10.5,F20.10,F20.10,F20.10,F20.10)          ANH00330
C                                                                       ANH00400
      WRITE(6,1)                                                        ANH00410
      RCONST=1.0D-01*PI*DSQRT(CL/(PH*AVN))                              ANH00420
C                                                                       ANH00430
C   NOTE : THE SYMMETRY AXIS IS CHOSEN AS THE Z AXIS                    ANH00440
      IXX=1                                                             ANH00450
      IYY=2                                                             ANH00460
      IZZ=3                                                             ANH00470
C                                                                       ANH00480
      IF(IPRNT.LE.2) GO TO 201                                          ANH00490
      DO 101 IABC=1,3                                                   ANH00500
      WRITE(6,2) IABC                                                   ANH00510
      CALL MATOUT(ZETA(1,1,IABC),N3N,N3N,N3N,N3N,6)                     ANH00520
  101 CONTINUE                                                          ANH00530
      DO 103 III=1,NVIB                                                 ANH00540
      IVIB=NFRQ(III)                                                    ANH00550
      WRITE(6,3) IVIB                                                   ANH00560
      DO 102 I=1,NVIB                                                   ANH00570
      II=NFRQ(I)                                                        ANH00580
      DO 102 J=1,NVIB                                                   ANH00590
      JJ=NFRQ(J)                                                        ANH00600
      EE(I,J)=F3Q(II,JJ,IVIB)                                           ANH00610
  102 CONTINUE                                                          ANH00620
      CALL MATOUT(EE,N3N,N3N,NVIB,NVIB,6)                               ANH00630
  103 CONTINUE                                                          ANH00640
      DO 104 IVIB=1,NVIB                                                ANH00650
      WRITE(6,4) IVIB                                                   ANH00660
      CALL MATOUT(AXY(1,1,IVIB),3,3,3,3,6)                              ANH00670
  104 CONTINUE                                                          ANH00680
C                                                                       ANH00690
C:::::::::::::::::::::::::::::::::::::::::::::::::::                    ANH00700
C:::ONLY FOR DEGENERATE MODES (FOR C3V MOLECULES):::                    ANH00710
C:::::::::::::::::::::::::::::::::::::::::::::::::::                    ANH00720
  201 CONTINUE                                                          ANH00730
      DO 120 IVIB=1,NVIB                                                ANH00740
      IF(NDEG(IVIB).EQ.1) GO TO 120                                     ANH00750
      FRQI=FREQ(IVIB)                                                   ANH00760
      FRQI2=FRQI*FRQI                                                   ANH00770
      AXZVAL=AXY(IXX,IZZ,IVIB)                                          ANH00780
      AYYVAL=AXY(IYY,IYY,IVIB)                                          ANH00790
      IF(DABS(AXZVAL).LE.ALIMIT.AND.DABS(AYYVAL).LE.ALIMIT) GO TO 202   ANH00800
C                                                                       ANH00810
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~                                    ANH00820
C~~~FOR AN A-DEGENERATE (T1) MODE~~~                                    ANH00830
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~                                    ANH00840
      VALU1=-0.75D+00*(AXZVAL**2/ROTAA(IZZ))                            ANH00850
      VALU2=A00                                                         ANH00860
      VALU3=A00                                                         ANH00870
C   SUM OVER NON-DEGENERATE MODES                                       ANH00880
      DO 112 JVIB=1,NVIB                                                ANH00890
      IF(NDEG(JVIB).NE.1) GO TO 112                                     ANH00900
      FRQJ=FREQ(JVIB)                                                   ANH00910
      FRQJ2=FRQJ*FRQJ                                                   ANH00920
      FACT2=(FRQI2*3.0D+00+FRQJ2)/(FRQI2-FRQJ2)                         ANH00930
      VALU2=VALU2+(ZETA(IVIB,JVIB,IYY)**2-ZETA(IVIB,JVIB,IXX)**2)*FACT2 ANH00940
  112 CONTINUE                                                          ANH00950
      DO 113 JVIB=1,NVIB                                                ANH00960
      IF(NDEG(JVIB).EQ.1) GO TO 113                                     ANH00970
      FRQJ=FREQ(JVIB)                                                   ANH00980
      FRQJ2=FRQJ*FRQJ                                                   ANH00990
      AXZVAL=AXY(IXX,IZZ,JVIB)                                          ANH01000
      AYYVAL=AXY(IYY,IYY,JVIB)                                          ANH01010
      IF(DABS(AXZVAL).LE.ALIMIT.AND.DABS(AYYVAL).LE.ALIMIT) GO TO 113   ANH01020
      FACT3=(FRQI*F3Q(IVIB,IVIB,JVIB)*RCONST*TWO)/(FRQJ**(1.5D+00))     ANH01030
      VALU3=VALU3+AYYVAL*FACT3                                          ANH01040
  113 CONTINUE                                                          ANH01050
      GO TO 203                                                         ANH01060
C                                                                       ANH01070
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~                                     ANH01080
C~~~FOR A B-DEGENERATE (T2) MODE~~~                                     ANH01090
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~                                     ANH01100
C   NOTE : CUBIC FORCE CONSTANTS HAVE THE FOLLOWING RELATIONSHIP        ANH01110
C          (TA,TA,TA')=-(TB,TB,TA')                                     ANH01120
C          (TB,TB,TB)=0                                                 ANH01130
C          (TB,TB,TB')=0                                                ANH01140
  202 CONTINUE                                                          ANH01150
      AYZVAL=AXY(IYY,IZZ,IVIB)                                          ANH01160
      AXYVAL=AXY(IXX,IYY,IVIB)                                          ANH01170
      IF(DABS(AYZVAL).LE.ALIMIT.AND.DABS(AXYVAL).LE.ALIMIT) GO TO 120   ANH01180
      VALU1=-0.75D+00*(AYZVAL**2/ROTAA(IZZ))                            ANH01190
      VALU2=A00                                                         ANH01200
      VALU3=A00                                                         ANH01210
      DO 116 JVIB=1,NVIB                                                ANH01220
      IF(NDEG(JVIB).NE.1) GO TO 116                                     ANH01230
      FRQJ=FREQ(JVIB)                                                   ANH01240
      FRQJ2=FRQJ*FRQJ                                                   ANH01250
      FACT2=(FRQI2*3.0D+00+FRQJ2)/(FRQI2-FRQJ2)                         ANH01260
      VALU2=VALU2-(ZETA(IVIB,JVIB,IYY)**2-ZETA(IVIB,JVIB,IXX)**2)*FACT2 ANH01270
  116 CONTINUE                                                          ANH01280
      DO 117 JVIB=1,NVIB                                                ANH01290
      IF(NDEG(JVIB).EQ.1) GO TO 117                                     ANH01300
      FRQJ=FREQ(JVIB)                                                   ANH01310
      FRQJ2=FRQJ*FRQJ                                                   ANH01320
      AXZVAL=AXY(IXX,IZZ,JVIB)                                          ANH01330
      AYYVAL=AXY(IYY,IYY,JVIB)                                          ANH01340
      IF(DABS(AXZVAL).LE.ALIMIT.AND.DABS(AYYVAL).LE.ALIMIT) GO TO 117   ANH01350
      FACT3=(FRQI*F3Q(IVIB,IVIB,JVIB)*RCONST*TWO)/(FRQJ**(1.5D+00))     ANH01360
      VALU3=VALU3-AYYVAL*FACT3                                          ANH01370
  117 CONTINUE                                                          ANH01380
C                                                                       ANH01390
  203 CONTINUE                                                          ANH01400
      FACTB=(ROTCM(IXX)*ROTCM(IXX)*TWO)/FRQI                            ANH01410
      VALUH=(VALU1+VALU2)*FACTB                                         ANH01420
      VALUA=VALU3*FACTB                                                 ANH01430
      VALUT=VALUH+VALUA                                                 ANH01440
      CLTYP(IVIB,1)=VALUH                                               ANH01450
      CLTYP(IVIB,2)=VALUA                                               ANH01460
      CLTYP(IVIB,3)=VALUT                                               ANH01470
      IF(IPRNT.LE.2) GO TO 120                                          ANH01480
      WRITE(6,5) IVIB,FRQI,VALUH,VALUA,VALUT                            ANH01490
  120 CONTINUE                                                          ANH01500
C                                                                       ANH01510
      WRITE(6,6)                                                        ANH01520
      WRITE(6,7)                                                        ANH01530
      DO 125 II=1,NVIB                                                  ANH01540
      IVIB=NFRQ(II)                                                     ANH01550
      IF(NDEG(IVIB).LE.1) GO TO 125                                     ANH01560
      FRQI=FREQ(IVIB)                                                   ANH01570
      VALUH=CLTYP(IVIB,1)                                               ANH01580
      VALUA=CLTYP(IVIB,2)                                               ANH01590
      VALUT=CLTYP(IVIB,3)                                               ANH01600
      VALUM=VALUT*CL*1.0D+04                                            ANH01610
      WRITE(6,8) II,IVIB,FRQI,VALUH,VALUA,VALUT,VALUM                   ANH01620
  125 CONTINUE                                                          ANH01630
C                                                                       ANH01640
      RETURN                                                            ANH01650
      END                                                               ANH01660
      SUBROUTINE ROVSYM(ZETA,AXY,F3Q,EE,RVCNST)                         ANH01670
C   THE CALCULATION OF VIBRATION-ROTATION INTERACTION CONSTANTS         ANH01680
C   FOR A SYMMETRIC TOP MOLECULE                                        ANH01690
      IMPLICIT REAL*8 (A-H,O-Z)                                         ANH01700
      CHARACTER*1 RAXIS                                                 ANH01710
      CHARACTER*25 RTYPE                                                ANH01720
      DIMENSION ZETA(N3N,N3N,3),AXY(3,3,N3N),F3Q(N3N,N3N,N3N)           ANH01730
      DIMENSION EE(N3N,N3N),RVCNST(NVIB,3)                              ANH01740
      COMMON/VIB101/NATOM,N3N,NATRI,ILIN,NVIB                           ANH01750
      COMMON/VIB103/PARA,WAVE,CONST,CYCL,CONV                           ANH01760
      COMMON/VIB106/ROTAA(3),ROTGC(3),ROTCM(3),ROTMH(3)                 ANH01770
      COMMON/VIB108/IAXIS(3),NDEG(150),NDAB(150,5),IMAG(150)            ANH01780
      COMMON/VIB109/CLIMIT,FLIM1,FLIM2                                  ANH01790
      COMMON/VIB203/IOFF(150),IPRNT                                     ANH01800
      COMMON/VIB204/SQM(150),ROOT(150),FREQ(150)                        ANH01810
      COMMON/VIB205/IFREQ,NFRQ(150)                                     ANH01820
      COMMON/VIB206/RAXIS(3),RTYPE(6)                                   ANH01830
      COMMON/VIB209/SALP(3)                                             ANH01840
      DATA PH,AVN / 6.626176D+00 , 6.022045D+00 /                       ANH01850
      DATA PI,CL / 3.1415926536D+00 , 2.99792458D+00 /                  ANH01860
      DATA A00,HALF,ONE,TWO / 0.0D+00 , 0.5D+00 , 1.0D+00 , 2.0D+00 /   ANH01870
      DATA ALIMIT / 1.0D-04 /                                           ANH01880
      DATA WLIMIT / 1.0D-02 /                                           ANH01890
    1 FORMAT(//,2X,' **************************************************'ANH01900
     1        /,2X,' ***THE VIBRATION-ROTATION INTERACTION CONSTANTS***'ANH01910
     2        /,2X,' **************************************************'ANH01920
     3        )                                                         ANH01930
    2 FORMAT(//,2X,' ZETA MATRIX, IABC = ',I5/)                         ANH01940
    3 FORMAT(//,2X,' F3(QI*QJ*QK) MATRIX (IN CM-1), IVIB = ',I5/)       ANH01950
    4 FORMAT(//,2X,' AXY MATRIX, IVIB = ',I5/)                          ANH01960
    5 FORMAT(//,2X,' QVIB',5X,' IVIB',10X,' A',21X,' B',21X,' C'/       ANH01970
     1          27X,' IN CM-1',15X,' IN CM-1',15X,' IN CM-1'/)          ANH01980
    6 FORMAT(2X,I5,5X,I5,3X,F20.10,3X,F20.10,3X,F20.10)                 ANH01990
    7 FORMAT(//,2X,' PRINCIPAL MOMENTS OF INERTIA WITH ZERO-POINT VIBRATANH02000
     1IONAL CORRECTION'/                                                ANH02010
     2 4X,5H AXIS,9X,11H IN AMU.A+2,9X,10H IN G.CM+2,10X,               ANH02020
     3 8H IN CM-1,15X,7H IN MHZ/38X,8H (*D+39)/)                        ANH02030
    8 FORMAT(8X,A1,3X,3F20.10,F20.5)                                    ANH02040
    9 FORMAT(//,2X,' ::::::::::::::::::::::::::::::::::::'/             ANH02050
     1          2X,' :::CORIOLIS RESONANCE IS OBSERVED:::'/             ANH02060
     2          2X,' ::::::::::::::::::::::::::::::::::::'/)            ANH02070
   10 FORMAT(2X,' RESONANCE LIMIT IS SET TO ',F10.3,' CM-1'/)           ANH02080
C                                                                       ANH02090
      WRITE(6,1)                                                        ANH02100
      RCONST=1.0D-01*PI*DSQRT(CL/(PH*AVN))                              ANH02110
C                                                                       ANH02120
C   NOTE : THE SYMMETRY AXIS IS CHOSEN AS THE Z AXIS                    ANH02130
      IXX=1                                                             ANH02140
      IYY=2                                                             ANH02150
      IZZ=3                                                             ANH02160
      NCORIO=0                                                          ANH02170
C                                                                       ANH02180
      IF(IPRNT.LE.2) GO TO 201                                          ANH02190
      DO 101 IABC=1,3                                                   ANH02200
      WRITE(6,2) IABC                                                   ANH02210
      CALL MATOUT(ZETA(1,1,IABC),N3N,N3N,N3N,N3N,6)                     ANH02220
  101 CONTINUE                                                          ANH02230
      DO 103 III=1,NVIB                                                 ANH02240
      IVIB=NFRQ(III)                                                    ANH02250
      WRITE(6,3) IVIB                                                   ANH02260
      DO 102 I=1,NVIB                                                   ANH02270
      II=NFRQ(I)                                                        ANH02280
      DO 102 J=1,NVIB                                                   ANH02290
      JJ=NFRQ(J)                                                        ANH02300
      EE(I,J)=F3Q(II,JJ,IVIB)                                           ANH02310
  102 CONTINUE                                                          ANH02320
      CALL MATOUT(EE,N3N,N3N,NVIB,NVIB,6)                               ANH02330
  103 CONTINUE                                                          ANH02340
      DO 104 IVIB=1,NVIB                                                ANH02350
      WRITE(6,4) IVIB                                                   ANH02360
      CALL MATOUT(AXY(1,1,IVIB),3,3,3,3,6)                              ANH02370
  104 CONTINUE                                                          ANH02380
C                                                                       ANH02390
C::::::::::::::::::::::::::::::                                         ANH02400
C:::FOR NON-DEGENERATE MODES:::                                         ANH02410
C::::::::::::::::::::::::::::::                                         ANH02420
  201 CONTINUE                                                          ANH02430
      CALL ZERO(RVCNST,NVIB*3)                                          ANH02440
      DO 110 IVIB=1,NVIB                                                ANH02450
      IF(NDEG(IVIB).NE.1) GO TO 110                                     ANH02460
      FRQI=FREQ(IVIB)                                                   ANH02470
      FRQI2=FRQI*FRQI                                                   ANH02480
      VALXS1=((AXY(IXX,IXX,IVIB)**2+AXY(IXX,IYY,IVIB)**2)/ROTAA(IXX))   ANH02490
     1       *0.75D+00                                                  ANH02500
      VALYS1=((AXY(IYY,IYY,IVIB)**2+AXY(IYY,IXX,IVIB)**2)/ROTAA(IYY))   ANH02510
     1       *0.75D+00                                                  ANH02520
      VALZS1=(AXY(IZZ,IZZ,IVIB)**2/ROTAA(IZZ))*0.75D+00                 ANH02530
      VALXS2=A00                                                        ANH02540
      VALYS2=A00                                                        ANH02550
      VALZS2=A00                                                        ANH02560
      VALXS3=A00                                                        ANH02570
      VALYS3=A00                                                        ANH02580
      VALZS3=A00                                                        ANH02590
C   SUM OVER NON-DEGENERATE MODES                                       ANH02600
      DO 107 JVIB=1,NVIB                                                ANH02610
      IF(NDEG(JVIB).NE.1) GO TO 107                                     ANH02620
      FRQJ=FREQ(JVIB)                                                   ANH02630
      FRQJ2=FRQJ*FRQJ                                                   ANH02640
      IF(IVIB.EQ.JVIB) GO TO 202                                        ANH02650
      FACT2=(FRQI2*3.0D+00+FRQJ2)/(FRQI2-FRQJ2)                         ANH02660
      VALZS2=VALZS2+(ZETA(IVIB,JVIB,IZZ)**2)*FACT2                      ANH02670
  202 CONTINUE                                                          ANH02680
      FACT3=(FRQI*F3Q(IVIB,IVIB,JVIB)*RCONST)/(FRQJ**(1.5D+00))         ANH02690
      VALXS3=VALXS3+AXY(IXX,IXX,JVIB)*FACT3                             ANH02700
      VALYS3=VALYS3+AXY(IYY,IYY,JVIB)*FACT3                             ANH02710
      VALZS3=VALZS3+AXY(IZZ,IZZ,JVIB)*FACT3                             ANH02720
  107 CONTINUE                                                          ANH02730
C   SUM OVER DEGENERATE MODES                                           ANH02740
      DO 108 JVIB=1,NVIB                                                ANH02750
      IF(NDEG(JVIB).EQ.1) GO TO 108                                     ANH02760
      FRQJ=FREQ(JVIB)                                                   ANH02770
      FRQJ2=FRQJ*FRQJ                                                   ANH02780
      FACT2A=0.5D+00*(FRQI2*3.0D+00+FRQJ2)/(FRQI2-FRQJ2)                ANH02790
      FACT2B=ZETA(IVIB,JVIB,IXX)**2+ZETA(IVIB,JVIB,IYY)**2              ANH02800
      VALXS2=VALXS2+FACT2A*FACT2B                                       ANH02810
      VALYS2=VALYS2+FACT2A*FACT2B                                       ANH02820
  108 CONTINUE                                                          ANH02830
      ROTA=-(VALXS1+VALXS2+VALXS3)*(ROTCM(IXX)**2*TWO)/FRQI             ANH02840
      ROTB=-(VALYS1+VALYS2+VALYS3)*(ROTCM(IYY)**2*TWO)/FRQI             ANH02850
      ROTC=-(VALZS1+VALZS2+VALZS3)*(ROTCM(IZZ)**2*TWO)/FRQI             ANH02860
      RVCNST(IVIB,1)=ROTA                                               ANH02870
      RVCNST(IVIB,2)=ROTB                                               ANH02880
      RVCNST(IVIB,3)=ROTC                                               ANH02890
  110 CONTINUE                                                          ANH02900
C                                                                       ANH02910
C::::::::::::::::::::::::::                                             ANH02920
C:::FOR DEGENERATE MODES:::                                             ANH02930
C::::::::::::::::::::::::::                                             ANH02940
C   NOTE : RELATIONSHIPS AMONG CUBIC FORCE CONSTANTS                    ANH02950
C     (S,TA,TA)=(S,TB,TB)                                               ANH02960
C     (S,TA,TA')=(S,TB,TB')                                             ANH02970
C   NOTE : RELATIONSHIPS AMONG ZETA CONSTANTS                           ANH02980
C     (TA,TA,IZZ)=(TB,TB,IZZ)=0                                         ANH02990
C     (TA,TA',IZZ)=(TB,TB',IZZ)=0                                       ANH03000
      DO 120 IVIB=1,NVIB                                                ANH03010
      IF(NDEG(IVIB).EQ.1) GO TO 120                                     ANH03020
      FRQI=FREQ(IVIB)                                                   ANH03030
      FRQI2=FRQI*FRQI                                                   ANH03040
      VALXT1=A00                                                        ANH03050
      VALYT1=A00                                                        ANH03060
      VALZT1=A00                                                        ANH03070
      VALXT2=A00                                                        ANH03080
      VALYT2=A00                                                        ANH03090
      VALZT2=A00                                                        ANH03100
      VALXT3=A00                                                        ANH03110
      VALYT3=A00                                                        ANH03120
      VALZT3=A00                                                        ANH03130
      VALXT4=A00                                                        ANH03140
      VALYT4=A00                                                        ANH03150
      VALZT4=A00                                                        ANH03160
      AXZVAL=AXY(IXX,IZZ,IVIB)                                          ANH03170
      AYYVAL=AXY(IYY,IYY,IVIB)                                          ANH03180
      IF(DABS(AXZVAL).LE.ALIMIT.AND.DABS(AYYVAL).LE.ALIMIT) GO TO 205   ANH03190
C                                                                       ANH03200
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~                                    ANH03210
C~~~FOR AN A-DEGENERATE (T1) MODE~~~                                    ANH03220
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~                                    ANH03230
      AXZSQ=AXZVAL**2                                                   ANH03240
      AXXSQ=AXY(IXX,IXX,IVIB)**2                                        ANH03250
      AYYSQ=AXY(IYY,IYY,IVIB)**2                                        ANH03260
      VALXT1=(0.375D+00*AXZSQ)/ROTAA(IZZ)                               ANH03270
     1      +(0.750D+00*AXXSQ)/ROTAA(IXX)                               ANH03280
      VALYT1=(0.375D+00*AXZSQ)/ROTAA(IZZ)                               ANH03290
     1      +(0.750D+00*AYYSQ)/ROTAA(IYY)                               ANH03300
      VALZT1=(0.750D+00*AXZSQ)/ROTAA(IXX)                               ANH03310
C   SUM OVER NON-DEGENERATE MODES                                       ANH03320
      DO 115 JVIB=1,NVIB                                                ANH03330
      IF(NDEG(JVIB).NE.1) GO TO 115                                     ANH03340
      FRQJ=FREQ(JVIB)                                                   ANH03350
      FRQJ2=FRQJ*FRQJ                                                   ANH03360
      FACT2A=0.50D+00*(FRQI2*3.0D+00+FRQJ2)/(FRQI2-FRQJ2)               ANH03370
      FACT2B=(ZETA(IVIB,JVIB,IXX)**2+ZETA(IVIB,JVIB,IYY)**2)            ANH03380
      VALXT2=VALXT2+FACT2A*FACT2B                                       ANH03390
      VALYT2=VALYT2+FACT2A*FACT2B                                       ANH03400
      FACT4=(FRQI*F3Q(IVIB,IVIB,JVIB)*RCONST)/(FRQJ**(1.5D+00))         ANH03410
      VALXT4=VALXT4+AXY(IXX,IXX,JVIB)*FACT4                             ANH03420
      VALYT4=VALYT4+AXY(IYY,IYY,JVIB)*FACT4                             ANH03430
      VALZT4=VALZT4+AXY(IZZ,IZZ,JVIB)*FACT4                             ANH03440
  115 CONTINUE                                                          ANH03450
C   SUM OVER DEGENERATE MODES                                           ANH03460
      DO 116 JVIB=1,NVIB                                                ANH03470
      IF(NDEG(JVIB).EQ.1) GO TO 116                                     ANH03480
      FRQJ=FREQ(JVIB)                                                   ANH03490
      FRQJ2=FRQJ*FRQJ                                                   ANH03500
      IF(IVIB.EQ.JVIB) GO TO 116                                        ANH03510
      IF(DABS(FRQI-FRQJ).LE.WLIMIT) GO TO 116                           ANH03520
      FACT3A=(FRQI2*3.0D+00+FRQJ2)/(FRQI2-FRQJ2)                        ANH03530
      FACT3B=(ZETA(IVIB,JVIB,IXX)**2+ZETA(IVIB,JVIB,IYY)**2)*HALF       ANH03540
      VALXT3=VALXT3+FACT3A*FACT3B                                       ANH03550
      VALYT3=VALYT3+FACT3A*FACT3B                                       ANH03560
      VALZT3=VALZT3+(ZETA(IVIB,JVIB,IZZ)**2)*FACT3A                     ANH03570
  116 CONTINUE                                                          ANH03580
      GO TO 210                                                         ANH03590
C                                                                       ANH03600
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~                                     ANH03610
C~~~FOR A B-DEGENERATE (T2) MODE~~~                                     ANH03620
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~                                     ANH03630
  205 CONTINUE                                                          ANH03640
      AYZVAL=AXY(IYY,IZZ,IVIB)                                          ANH03650
      AXYVAL=AXY(IXX,IYY,IVIB)                                          ANH03660
      IF(DABS(AYZVAL).LE.ALIMIT.AND.DABS(AXYVAL).LE.ALIMIT) GO TO 210   ANH03670
      AYZSQ=AYZVAL**2                                                   ANH03680
      AXYSQ=AXYVAL**2                                                   ANH03690
      VALXT1=(0.375D+00*AYZSQ)/ROTAA(IZZ)                               ANH03700
     1      +(0.750D+00*AXYSQ)/ROTAA(IXX)                               ANH03710
      VALYT1=(0.375D+00*AYZSQ)/ROTAA(IZZ)                               ANH03720
     1      +(0.750D+00*AXYSQ)/ROTAA(IYY)                               ANH03730
      VALZT1=(0.750D+00*AYZSQ)/ROTAA(IXX)                               ANH03740
      DO 117 JVIB=1,NVIB                                                ANH03750
      IF(NDEG(JVIB).NE.1) GO TO 117                                     ANH03760
      FRQJ=FREQ(JVIB)                                                   ANH03770
      FRQJ2=FRQJ*FRQJ                                                   ANH03780
      FACT2A=0.50D+00*(FRQI2*3.0D+00+FRQJ2)/(FRQI2-FRQJ2)               ANH03790
      FACT2B=(ZETA(IVIB,JVIB,IXX)**2+ZETA(IVIB,JVIB,IYY)**2)            ANH03800
      VALXT2=VALXT2+FACT2A*FACT2B                                       ANH03810
      VALYT2=VALYT2+FACT2A*FACT2B                                       ANH03820
      FACT4=(FRQI*F3Q(IVIB,IVIB,JVIB)*RCONST)/(FRQJ**(1.5D+00))         ANH03830
      VALXT4=VALXT4+AXY(IXX,IXX,JVIB)*FACT4                             ANH03840
      VALYT4=VALYT4+AXY(IYY,IYY,JVIB)*FACT4                             ANH03850
      VALZT4=VALZT4+AXY(IZZ,IZZ,JVIB)*FACT4                             ANH03860
  117 CONTINUE                                                          ANH03870
C   SUM OVER DEGENERATE MODES                                           ANH03880
      DO 118 JVIB=1,NVIB                                                ANH03890
      IF(NDEG(JVIB).EQ.1) GO TO 118                                     ANH03900
      FRQJ=FREQ(JVIB)                                                   ANH03910
      FRQJ2=FRQJ*FRQJ                                                   ANH03920
      IF(IVIB.EQ.JVIB) GO TO 118                                        ANH03930
      IF(DABS(FRQI-FRQJ).LE.WLIMIT) GO TO 118                           ANH03940
      FACT3A=(FRQI2*3.0D+00+FRQJ2)/(FRQI2-FRQJ2)                        ANH03950
      FACT3B=(ZETA(IVIB,JVIB,IXX)**2+ZETA(IVIB,JVIB,IYY)**2)*HALF       ANH03960
      VALXT3=VALXT3+FACT3A*FACT3B                                       ANH03970
      VALYT3=VALYT3+FACT3A*FACT3B                                       ANH03980
      VALZT3=VALZT3+(ZETA(IVIB,JVIB,IZZ)**2)*FACT3A                     ANH03990
  118 CONTINUE                                                          ANH04000
C                                                                       ANH04010
  210 CONTINUE                                                          ANH04020
      ROTA=-(VALXT1+VALXT2+VALXT3+VALXT4)*(ROTCM(IXX)**2*TWO)/FRQI      ANH04030
      ROTB=-(VALYT1+VALYT2+VALYT3+VALYT4)*(ROTCM(IYY)**2*TWO)/FRQI      ANH04040
      ROTC=-(VALZT1+VALZT2+VALZT3+VALZT4)*(ROTCM(IZZ)**2*TWO)/FRQI      ANH04050
      RVCNST(IVIB,1)=ROTA                                               ANH04060
      RVCNST(IVIB,2)=ROTB                                               ANH04070
      RVCNST(IVIB,3)=ROTC                                               ANH04080
C                                                                       ANH04090
  120 CONTINUE                                                          ANH04100
C                                                                       ANH04110
C   CALCULATE THE ROTATIONAL CONSTANTS                                  ANH04120
C   WITH ZERO-POINT VIBRATIONAL CORRECTION                              ANH04130
      CALL ZERO(SALP,3)                                                 ANH04140
      WRITE(6,5)                                                        ANH04150
      DO 121 III=1,NVIB                                                 ANH04160
      IVIB=NFRQ(III)                                                    ANH04170
      ROTA=RVCNST(IVIB,1)                                               ANH04180
      ROTB=RVCNST(IVIB,2)                                               ANH04190
      ROTC=RVCNST(IVIB,3)                                               ANH04200
      SALP(1)=SALP(1)+ROTA                                              ANH04210
      SALP(2)=SALP(2)+ROTB                                              ANH04220
      SALP(3)=SALP(3)+ROTC                                              ANH04230
      WRITE(6,6) III,IVIB,ROTA,ROTB,ROTC                                ANH04240
  121 CONTINUE                                                          ANH04250
      WRITE(6,7)                                                        ANH04260
      DO 122 I=1,3                                                      ANH04270
      PB=ROTCM(I)-SALP(I)*HALF                                          ANH04280
      PT=CONST/PB                                                       ANH04290
      PA=PT*PARA                                                        ANH04300
      PC=CYCL/PT                                                        ANH04310
      WRITE(6,8) RAXIS(I),PT,PA,PB,PC                                   ANH04320
  122 CONTINUE                                                          ANH04330
C                                                                       ANH04340
      IF(NCORIO.NE.0) THEN                                              ANH04350
        WRITE(6,9)                                                      ANH04360
        WRITE(6,10)                                                     ANH04370
      END IF                                                            ANH04380
C                                                                       ANH04390
      RETURN                                                            ANH04400
      END                                                               ANH04410
      SUBROUTINE CENTSY(F3Q,EE)                                         ANH04420
C   THE CALCULATION OF CENTRIFUGAL DISTORTION CONSTANTS                 ANH04430
C   FOR A SYMMETRIC MOLECULE                                            ANH04440
      IMPLICIT REAL*8 (A-H,O-Z)                                         ANH04450
      DIMENSION F3Q(N3N,N3N,N3N),EE(N3N,N3N)                            ANH04460
      COMMON/VIB101/NATOM,N3N,NATRI,ILIN,NVIB                           ANH04470
      COMMON/VIB102/ITHREE,IFOUR,N3TOT,N4TOT                            ANH04480
      COMMON/VIB106/ROTAA(3),ROTGC(3),ROTCM(3),ROTMH(3)                 ANH04490
      COMMON/VIB108/IAXIS(3),NDEG(150),NDAB(150,5),IMAG(150)            ANH04500
      COMMON/VIB203/IOFF(150),IPRNT                                     ANH04510
      COMMON/VIB204/SQM(150),ROOT(150),FREQ(150)                        ANH04520
      COMMON/VIB205/IFREQ,NFRQ(150)                                     ANH04530
      COMMON/VIB207/TABCD(3,3,3,3),PABC(3,3,3)                          ANH04540
      COMMON/VIB208/TAABB(3,3),TABAB(3,3)                               ANH04550
      DATA WLIMIT / 1.0D-2 /                                            ANH04560
      DATA PH,AVN / 6.626176D+00 , 6.022045D+00 /                       ANH04570
      DATA PI,CL / 3.1415926536D+00 , 2.99792458D+00 /                  ANH04580
      DATA A00,HALF,ONE,TWO / 0.0D+00 , 0.5D+00 , 1.0D+00 , 2.0D+00 /   ANH04590
    1 FORMAT(//,2X,' ******************************************'/       ANH04600
     1          2X,' ***THE CENTRIFUGAL DISTORTION CONSTANTS***'/       ANH04610
     2          2X,' ******************************************')       ANH04620
    2 FORMAT(//,2X,' F3(QI*QJ*QK) MATRIX (IN CM-1), IVIB = ',I5/)       ANH04630
    3 FORMAT(//,2X,' TAABB MATRIX FROM BXY MATRICES (IN CM-1)'/)        ANH04640
    4 FORMAT(//,2X,' TABAB MATRIX FROM BXY MATRICES (IN CM-1)'/)        ANH04650
    5 FORMAT(//,2X,' FOR A SYMMETRIC TOP MOLECULE'/                     ANH04660
     1          2X,' DJ   = ',F15.5,' D-06 CM-1',                       ANH04670
     2          8X,' DJ   = ',F15.5,' MHZ'/                             ANH04680
     3          2X,' DK   = ',F15.5,' D-06 CM-1',                       ANH04690
     4          8X,' DK   = ',F15.5,' MHZ'/                             ANH04700
     5          2X,' DJK  = ',F15.5,' D-06 CM-1',                       ANH04710
     6          8X,' DJK  = ',F15.5,' MHZ'/)                            ANH04720
    6 FORMAT(//,2X,' FOR A SYMMETRIC TOP MOLECULE'/                     ANH04730
     1          2X,' HJ   = ',F15.5,' D-06 CM-1',                       ANH04740
     2          8X,' HJ   = ',F15.5,' MHZ'/                             ANH04750
     3          2X,' HJK  = ',F15.5,' D-06 CM-1',                       ANH04760
     4          8X,' HJK  = ',F15.5,' MHZ'/                             ANH04770
     5          2X,' HKJ  = ',F15.5,' D-06 CM-1',                       ANH04780
     6          8X,' HKJ  = ',F15.5,' MHZ'/                             ANH04790
     7          2X,' HK   = ',F15.5,' D-06 CM-1',                       ANH04800
     8          8X,' HK   = ',F15.5,' MHZ'/                             ANH04810
     9          2X,' DEL3 = ',F15.5,' D-06 CM-1',                       ANH04820
     A          8X,' DEL3 = ',F15.5,' MHZ'/)                            ANH04830
C                                                                       ANH04840
      WRITE(6,1)                                                        ANH04850
C                                                                       ANH04860
C   ELEMENTS OF F3Q ARE IN CM-1                                         ANH04870
      IF(IPRNT.LE.2) GO TO 201                                          ANH04880
      DO 102 III=1,NVIB                                                 ANH04890
      IVIB=NFRQ(III)                                                    ANH04900
      WRITE(6,2) IVIB                                                   ANH04910
      DO 101 I=1,NVIB                                                   ANH04920
      II=NFRQ(I)                                                        ANH04930
      DO 101 J=1,NVIB                                                   ANH04940
      JJ=NFRQ(J)                                                        ANH04950
      EE(I,J)=F3Q(II,JJ,IVIB)                                           ANH04960
  101 CONTINUE                                                          ANH04970
      CALL MATOUT(EE,N3N,N3N,N3N,N3N,6)                                 ANH04980
  102 CONTINUE                                                          ANH04990
C                                                                       ANH05000
C   ELEMENTS OF TAABB AND TABAB MATRICES ARE IN CM-1                    ANH05010
  201 CONTINUE                                                          ANH05020
      DO 103 IA=1,3                                                     ANH05030
      DO 103 IB=1,3                                                     ANH05040
      TAABB(IA,IB)=TABCD(IA,IA,IB,IB)                                   ANH05050
      TABAB(IA,IB)=TABCD(IA,IB,IA,IB)                                   ANH05060
  103 CONTINUE                                                          ANH05070
      IF(IPRNT.LE.2) GO TO 202                                          ANH05080
      WRITE(6,3)                                                        ANH05090
      CALL MATOUT(TAABB,3,3,3,3,6)                                      ANH05100
      WRITE(6,4)                                                        ANH05110
      CALL MATOUT(TABAB,3,3,3,3,6)                                      ANH05120
C                                                                       ANH05130
C   FOR A SYMMETRIC TOP MOLECULE                                        ANH05140
  202 CONTINUE                                                          ANH05150
      IXX=1                                                             ANH05160
      IYY=2                                                             ANH05170
      IZZ=3                                                             ANH05180
      DJ=-(TAABB(IXX,IXX)*3.0D+00+TAABB(IYY,IYY)*3.0D+00                ANH05190
     1    +TAABB(IXX,IYY)*2.0D+00+TABAB(IXX,IYY)*4.0D+00)/32.0D+00      ANH05200
      DK=DJ-(TAABB(IZZ,IZZ)-TAABB(IXX,IZZ)-TAABB(IYY,IZZ)               ANH05210
     1      -TABAB(IXX,IZZ)*2.0D+00-TABAB(IYY,IZZ)*2.0D+00)/4.0D+00     ANH05220
      DJK=-DJ-DK-TAABB(IZZ,IZZ)/4.0D+00                                 ANH05230
      DJCM=DJ*1.0D+06                                                   ANH05240
      DKCM=DK*1.0D+06                                                   ANH05250
      DJKCM=DJK*1.0D+06                                                 ANH05260
      DJMH=DJ*CL*1.0D+04                                                ANH05270
      DKMH=DK*CL*1.0D+04                                                ANH05280
      DJKMH=DJK*CL*1.0D+04                                              ANH05290
      WRITE(6,5) DJCM,DJMH,DKCM,DKMH,DJKCM,DJKMH                        ANH05300
C                                                                       ANH05310
C   HIGHER ORDER CENTRIFUGAL DISTORTION CONSTANTS                       ANH05320
      HJ=(PABC(IXX,IXX,IXX)+PABC(IYY,IYY,IYY))*HALF                     ANH05330
      HJK=PABC(IZZ,IXX,IXX)*TWO-HJ*3.0D+00                              ANH05340
      HKJ=PABC(IZZ,IZZ,IXX)*TWO-HJK*TWO-HJ*3.0D+00                      ANH05350
      HK=PABC(IZZ,IZZ,IZZ)-HKJ-HJK-HJ                                   ANH05360
      DEL3=(PABC(IXX,IXX,IXX)-PABC(IYY,IYY,IYY))*HALF                   ANH05370
      HJCM=HJ*1.0D+06                                                   ANH05380
      HJKCM=HJK*1.0D+06                                                 ANH05390
      HKJCM=HKJ*1.0D+06                                                 ANH05400
      HKCM=HK*1.0D+06                                                   ANH05410
      DEL3CM=DEL3*1.0D+06                                               ANH05420
      HJMH=HJ*CL*1.0D+04                                                ANH05430
      HJKMH=HJK*CL*1.0D+04                                              ANH05440
      HKJMH=HKJ*CL*1.0D+04                                              ANH05450
      HKMH=HK*CL*1.0D+04                                                ANH05460
      DEL3MH=DEL3*CL*1.0D+04                                            ANH05470
      WRITE(6,6) HJCM,HJMH,HJKCM,HJKMH,HKJCM,HKJMH,HKCM,HKMH,           ANH05480
     1           DEL3CM,DEL3MH                                          ANH05490
C                                                                       ANH05500
      RETURN                                                            ANH05510
      END                                                               ANH05520
      SUBROUTINE VLTYSY(ZETA,AXY,F3Q,F4Q,CLTYP)                         ANH05530
C   THE CALCULATION OF VIBRATIONAL L-TYPE DOUBLING CONSTANTS            ANH05540
C   FOR A SYMMETRIC TOP MOLECULE                                        ANH05550
      IMPLICIT REAL*8 (A-H,O-Z)                                         ANH05560
      DIMENSION CLTYP(N3N,N3N,4)                                        ANH05570
      DIMENSION ZETA(N3N,N3N,3),AXY(3,3,N3N)                            ANH05580
      DIMENSION F3Q(N3N,N3N,N3N),F4Q(N3N,N3N,N3N,N3N)                   ANH05590
      COMMON/VIB101/NATOM,N3N,NATRI,ILIN,NVIB                           ANH05600
      COMMON/VIB106/ROTAA(3),ROTGC(3),ROTCM(3),ROTMH(3)                 ANH05610
      COMMON/VIB108/IAXIS(3),NDEG(150),NDAB(150,5),IMAG(150)            ANH05620
      COMMON/VIB109/CLIMIT,FLIM1,FLIM2                                  ANH05630
      COMMON/VIB203/IOFF(150),IPRNT                                     ANH05640
      COMMON/VIB204/SQM(150),ROOT(150),FREQ(150)                        ANH05650
      COMMON/VIB205/IFREQ,NFRQ(150)                                     ANH05660
      DATA A00,HALF,ONE,TWO / 0.0D+00 , 0.5D+00 , 1.0D+00 , 2.0D+00 /   ANH05670
      DATA ALIMIT / 1.0D-04 /                                           ANH05680
      DATA WLIMIT / 1.0D-02 /                                           ANH05690
      DATA DLIMIT / 1.0D-02 /                                           ANH05700
      DATA PI,CL / 3.1415926536D+00 , 2.99792458D+00 /                  ANH05710
    1 FORMAT(//,2X,' ***********************************************'/  ANH05720
     1          2X,' ***THE VIBRATIONAL L-TYPE DOUBLING CONSTANTS***'/  ANH05730
     2          2X,' ***********************************************')  ANH05740
    2 FORMAT(2X,2I5,5X,7F12.5)                                          ANH05750
    3 FORMAT(//,2X,' DOUBLING CONSTANTS IN CM-1'/                       ANH05760
     1          4X,' QI',5X,' QJ',8X,' FREQ(I)',5X,' FREQ(J)',          ANH05770
     2          9X,' QUADR.',8X,' CUBIC',9X,' QUARTIC',7X,' TOTAL'/)    ANH05780
    4 FORMAT(2X,I5,3X,I5,3X,F13.3,F13.3,7X,F10.5,5X,F10.5,5X,F10.5,     ANH05790
     1       5X,F10.5)                                                  ANH05800
C                                                                       ANH05810
      WRITE(6,1)                                                        ANH05820
C                                                                       ANH05830
C   NOTE : THE MOLECULAR AXIS IS CHOSEN AS THE Z AXIS                   ANH05840
      IXX=1                                                             ANH05850
      IYY=2                                                             ANH05860
      IZZ=3                                                             ANH05870
C                                                                       ANH05880
C:::::::::::::::::::::::::::::::                                        ANH05890
C:::ONLY FOR DEGENERATE MODES:::                                        ANH05900
C:::::::::::::::::::::::::::::::                                        ANH05910
      CALL ZERO(CLTYP,N3N*N3N*4)                                        ANH05920
      DO 110 IVIB=1,NVIB                                                ANH05930
      IF(NDEG(IVIB).EQ.1) GO TO 110                                     ANH05940
      FRQI=FREQ(IVIB)                                                   ANH05950
      FRQI2=FRQI*FRQI                                                   ANH05960
      IA=NDAB(IVIB,1)                                                   ANH05970
      IB=NDAB(IVIB,2)                                                   ANH05980
      DO 109 JVIB=1,NVIB                                                ANH05990
      IF(NDEG(JVIB).EQ.1) GO TO 109                                     ANH06000
      FRQJ=FREQ(JVIB)                                                   ANH06010
      FRQJ2=FRQJ*FRQJ                                                   ANH06020
      JA=NDAB(JVIB,1)                                                   ANH06030
      JB=NDAB(JVIB,2)                                                   ANH06040
      IF(IVIB.EQ.JVIB) GO TO 109                                        ANH06050
      IF(DABS(FRQI-FRQJ).LE.DLIMIT) GO TO 109                           ANH06060
C                                                                       ANH06070
      VALU1=A00                                                         ANH06080
      VALU2=A00                                                         ANH06090
      VALU3=A00                                                         ANH06100
      VALU4=A00                                                         ANH06110
C                                                                       ANH06120
C   SUM OVER NON-DEGENERATE MODES                                       ANH06130
C   NOTE : RELATIONSHIP AMONG CUBIC FORCE CONSTANTS                     ANH06140
C     (S,TA,TA')=(S,TB,TB')                                             ANH06150
C     (S,TA,TB')=(S,TB,TA')=0                                           ANH06160
      DO 101 KVIB=1,NVIB                                                ANH06170
      IF(NDEG(KVIB).NE.1) GO TO 101                                     ANH06180
      FRQK=FREQ(KVIB)                                                   ANH06190
      FRQK2=FRQK*FRQK                                                   ANH06200
      DNUM=(FRQK2-FRQI2-FRQJ2)*FRQK                                     ANH06210
      FACT1=F3Q(IA,JA,KVIB)**2+F3Q(IA,JB,KVIB)**2                       ANH06220
     1     +F3Q(IB,JA,KVIB)**2+F3Q(IB,JB,KVIB)**2                       ANH06230
      WPPP=FRQI+FRQJ+FRQK                                               ANH06240
      WPPM=FRQI+FRQJ-FRQK                                               ANH06250
      WPMP=FRQI-FRQJ+FRQK                                               ANH06260
      WPMM=FRQI-FRQJ-FRQK                                               ANH06270
      DIJK=WPPP*WPPM*WPMP*WPMM                                          ANH06280
      VALU1=VALU1-(0.5D+00*FACT1*DNUM)/(DIJK*TWO)                       ANH06290
  101 CONTINUE                                                          ANH06300
C                                                                       ANH06310
C  SUM OVER DEGENERATE MODES                                            ANH06320
      DO 102 KVIB=1,NVIB                                                ANH06330
      IF(NDEG(KVIB).LE.1) GO TO 102                                     ANH06340
      FRQK=FREQ(KVIB)                                                   ANH06350
      VALU2=VALU2-F3Q(IVIB,IVIB,KVIB)*F3Q(JVIB,JVIB,KVIB)/FRQK          ANH06360
  102 CONTINUE                                                          ANH06370
      VALU2=VALU2*HALF                                                  ANH06380
      VALU3=(F4Q(IA,IA,JA,JA)-F4Q(IA,IA,JB,JB))*0.5D+00                 ANH06390
      FACT4=ZETA(IA,JA,IZZ)**2+ZETA(IA,JB,IZZ)**2                       ANH06400
     1     +ZETA(IB,JA,IZZ)**2+ZETA(IB,JB,IZZ)**2                       ANH06410
      VALU4=-(HALF*FACT4*ROTCM(IZZ)*(FRQI2+FRQJ2))/(FRQI*FRQJ)          ANH06420
      VALUT=VALU1+VALU2+VALU3+VALU4                                     ANH06430
      CLTYP(IVIB,JVIB,1)=VALU4                                          ANH06440
      CLTYP(IVIB,JVIB,2)=VALU1+VALU2                                    ANH06450
      CLTYP(IVIB,JVIB,3)=VALU3                                          ANH06460
      CLTYP(IVIB,JVIB,4)=VALUT                                          ANH06470
      IF(IPRNT.LE.2) GO TO 109                                          ANH06480
      WRITE(6,2) IVIB,JVIB,FRQI,FRQJ,VALU1,VALU2,VALU3,VALU4,VALUT      ANH06490
  109 CONTINUE                                                          ANH06500
  110 CONTINUE                                                          ANH06510
C                                                                       ANH06520
      WRITE(6,3)                                                        ANH06530
      DO 112 II=1,NVIB                                                  ANH06540
      IVIB=NFRQ(II)                                                     ANH06550
      IF(NDEG(IVIB).LE.1) GO TO 112                                     ANH06560
      FRQI=FREQ(IVIB)                                                   ANH06570
      DO 111 JJ=1,NVIB                                                  ANH06580
      JVIB=NFRQ(JJ)                                                     ANH06590
      IF(NDEG(JVIB).LE.1) GO TO 111                                     ANH06600
      FRQJ=FREQ(JVIB)                                                   ANH06610
      IF(DABS(FRQI-FRQJ).LE.DLIMIT) GO TO 111                           ANH06620
      VALU2=CLTYP(IVIB,JVIB,1)                                          ANH06630
      VALU3=CLTYP(IVIB,JVIB,2)                                          ANH06640
      VALU4=CLTYP(IVIB,JVIB,3)                                          ANH06650
      VALUT=CLTYP(IVIB,JVIB,4)                                          ANH06660
      WRITE(6,4) II,JJ,FRQI,FRQJ,VALU2,VALU3,VALU4,VALUT                ANH06670
  111 CONTINUE                                                          ANH06680
  112 continue
C                                                                       ANH06690
      RETURN                                                            ANH06700
      END                                                               ANH06710
      SUBROUTINE SYMTOP(ZETA,F3Q,F4Q,XIJ,GIJ)                           ANH06720
C   THE ANHARMONICITY CONTSTANTS FOR A SYMMETRIC TOP                    ANH06730
      IMPLICIT REAL*8 (A-H,O-Z)                                         ANH06740
      DIMENSION ZETA(N3N,N3N,3),F3Q(N3N,N3N,N3N),F4Q(N3N,N3N,N3N,N3N)   ANH06750
      DIMENSION XIJ(N3N,N3N),GIJ(N3N,N3N)                               ANH06760
      COMMON/VIB101/NATOM,N3N,NATRI,ILIN,NVIB                           ANH06770
      COMMON/VIB106/ROTAA(3),ROTGC(3),ROTCM(3),ROTMH(3)                 ANH06780
      COMMON/VIB108/IAXIS(3),NDEG(150),NDAB(150,5),IMAG(150)            ANH06790
      COMMON/VIB109/CLIMIT,FLIM1,FLIM2                                  ANH06800
      COMMON/VIB111/ANHTOT,HZPVE,AZPVE
      COMMON/VIB203/IOFF(150),IPRNT                                     ANH06810
      COMMON/VIB204/SQM(150),ROOT(150),FREQ(150)                        ANH06820
      COMMON/VIB205/IFREQ,NFRQ(150)                                     ANH06830
      DATA A00,ONE / 0.0D+00 , 1.0D+00 /                                ANH06840
      DATA DLIMIT / 1.0D-02 /                                           ANH06850
      DATA WLIMIT / 3.0D+00 /                                           ANH06860
      DATA PLIMIT / 1.0D-05 /                                           ANH06870
    1 FORMAT(//,2X,' *********************************************'/    ANH06880
     1          2X,' ***THE VIBRATIONAL ANHARMONICITY CONSTANTS***'/    ANH06890
     2          2X,' *********************************************'/)   ANH06900
    2 FORMAT(//,2X,' ANHARMONIC CONSTANTS IN CM-1'/                     ANH06910
     1          2X,' IVIB',3X,' JVIB',8X,' FREQ(I)',5X,' FREQ(J)',      ANH06920
     2          7X,' VALU1',7X,' VALU2',7X,' VALU3',7X,' VALU4',        ANH06930
     3          7X,' XIJ',9X,' GIJ'/)                                   ANH06940
    3 FORMAT(2X,I5,3X,I5,3X,8(F13.3))                                   ANH06950
    4 FORMAT(//,2X,' ANHARMONIC CONSTANTS IN CM-1'/                     ANH06960
     1          4X,' QI',5X,' QJ',8X,' FREQ(I)',5X,' FREQ(J)',          ANH06970
     2          7X,' XIJ',9X,' GIJ'/)                                   ANH06980
    5 FORMAT(2X,I5,3X,I5,3X,4(F13.3))                                   ANH06990
    6 FORMAT(//,2X,' THEORETICAL FREQUENCIES IN CM-1'/                  ANH07000
     1          2X,' NO.',9X,' HARMONIC',14X,' DELTA',                  ANH07010
     2         13X,' FUNDAMENTAL'/)                                     ANH07020
    7 FORMAT(2X,I5,4X,F13.3,7X,F13.3,12X,F13.3)                         ANH07030
    8 FORMAT(//,2X,' :::::::::::::::::::::::::::::::::::::::::::::::'/  ANH07040
     1          2X,' :::ANHARMONIC RESONANCE (TYPE I) IS OBSERVED:::'/  ANH07050
     2          2X,' :::::::::::::::::::::::::::::::::::::::::::::::'/) ANH07060
    9 FORMAT(//,2X,' ::::::::::::::::::::::::::::::::::::::::::::::::'/ ANH07070
     1          2X,' :::ANHARMONIC RESONANCE (TYPE II) IS OBSERVED:::'/ ANH07080
     2          2X,' ::::::::::::::::::::::::::::::::::::::::::::::::'/)ANH07090
   10 FORMAT(2X,' RESONANCE LIMIT IS SET TO ',F10.3,' CM-1'/)           ANH07100
C                                                                       ANH07110
C:::::::::::::::::::::::::::::::::::::::                                ANH07120
C:::CALCULATE ANHARMONICITY CONSTANTS:::                                ANH07130
C:::::::::::::::::::::::::::::::::::::::                                ANH07140
      IANH1=0                                                           ANH07150
      IANH2=0                                                           ANH07160
C                                                                       ANH07170
      WRITE(6,1)                                                        ANH07180
      IF(IPRNT.GT.2)                                                    ANH07190
     *WRITE(6,2)                                                        ANH07200
      IXX=1                                                             ANH07210
      IYY=2                                                             ANH07220
      IZZ=3                                                             ANH07230
C                                                                       ANH07240
      DO 120 IVIB=1,N3N                                                 ANH07250
      FRQI=FREQ(IVIB)                                                   ANH07260
      FRQI2=FRQI*FRQI                                                   ANH07270
      IVA=NDAB(IVIB,1)                                                  ANH07280
      IVB=NDAB(IVIB,2)                                                  ANH07290
      DO 120 JVIB=1,N3N                                                 ANH07300
      FRQJ=FREQ(JVIB)                                                   ANH07310
      FRQJ2=FRQJ*FRQJ                                                   ANH07320
      JVA=NDAB(JVIB,1)                                                  ANH07330
      JVB=NDAB(JVIB,2)                                                  ANH07340
C                                                                       ANH07350
      VALU1=A00                                                         ANH07360
      VALU2=A00                                                         ANH07370
      VALU3=A00                                                         ANH07380
      VALU4=A00                                                         ANH07390
      VALU5=A00                                                         ANH07400
      GVAL1=A00                                                         ANH07410
      GVAL2=A00                                                         ANH07420
      GVAL3=A00                                                         ANH07430
      GVAL4=A00                                                         ANH07440
      XRS=A00                                                           ANH07450
      GRS=A00                                                           ANH07460
      IF(DABS(FRQI).LE.WLIMIT) GO TO 250                                ANH07470
      IF(DABS(FRQJ).LE.WLIMIT) GO TO 250                                ANH07480
      IF(IVIB.NE.JVIB) GO TO 210                                        ANH07490
C                                                                       ANH07500
C+++++++++++++++++++++++++++++++++++++++++                              ANH07510
C+++DIAGONAL CONSTANTS (NON-DEGENERATE)+++                              ANH07520
C+++++++++++++++++++++++++++++++++++++++++                              ANH07530
      IF(NDEG(IVIB).NE.1) GO TO 205                                     ANH07540
      VALU1=F4Q(IVIB,IVIB,IVIB,IVIB)                                    ANH07550
      DO 101 KVIB=1,N3N                                                 ANH07560
C   SUM OVER ONLY NON-DEGENERATE MODES                                  ANH07570
      FRQK=FREQ(KVIB)                                                   ANH07580
      IF(DABS(FRQK).LE.WLIMIT) GO TO 101                                ANH07590
      IF(NDEG(KVIB).NE.1) GO TO 101                                     ANH07600
      FRQK2=FRQK*FRQK                                                   ANH07610
      FAC=F3Q(IVIB,IVIB,KVIB)                                           ANH07620
      FACT=FAC*FAC                                                      ANH07630
      IF(DABS(FRQI+FRQI-FRQK).LE.FLIM1) GO TO 201                       ANH07640
      DNUM=FRQI2*8.0D+00-FRQK2*3.0D+00                                  ANH07650
      DENM=FRQI2*4.0D+00-FRQK2                                          ANH07660
      VALU2=VALU2+(DNUM*FACT)/(FRQK*DENM)                               ANH07670
      GO TO 101                                                         ANH07680
C/////////////////////////////////                                      ANH07690
C///FERMI RESONANCE CORRECTIONS///                                      ANH07700
C/////////////////////////////////                                      ANH07710
  201 CONTINUE                                                          ANH07720
      DNUM=ONE/((FRQI+FRQI+FRQK)*2.0D+00)+2.0D+00/FRQK                  ANH07730
      IF(DABS(FAC).GT.PLIMIT) THEN                                      ANH07740
        IANH1=IANH1+1                                                   ANH07750
      END IF                                                            ANH07760
      VALU2=VALU2+FACT*DNUM                                             ANH07770
C////////////////////////////                                           ANH07780
C///END OF FERMI RESONANCE///                                           ANH07790
C////////////////////////////                                           ANH07800
  101 CONTINUE                                                          ANH07810
      XRS=(VALU1-VALU2)/16.0D+00                                        ANH07820
      GO TO 250                                                         ANH07830
C                                                                       ANH07840
C+++++++++++++++++++++++++++++++++++++                                  ANH07850
C+++DIAGONAL CONSTANTS (DEGENERATE)+++                                  ANH07860
C+++++++++++++++++++++++++++++++++++++                                  ANH07870
C   NOTE : RELATIONSHIPS AMONG CUBIC FORCE CONSTANTS                    ANH07880
C     (S,TA,TA)=(S,TB,TB)                                               ANH07890
C     (S,TA,TA')=(S,TB,TB')                                             ANH07900
C     (S,TA,TB)=0                                                       ANH07910
C     (S,TA,TB')=0                                                      ANH07920
C     (S,TB,TA')=0                                                      ANH07930
C   NOTE : RELATIONSHIPS AMONG QUARTIC FORCE CONSTANTS                  ANH07940
C     (TA,TA,TA,TA)=(TB,TB,TB,TB)                                       ANH07950
C                                                                       ANH07960
  205 CONTINUE                                                          ANH07970
      VALU1=F4Q(IVIB,IVIB,IVIB,IVIB)                                    ANH07980
      GVAL1=-VALU1                                                      ANH07990
C   SUM OVER NON-DEGENERATE MODES                                       ANH08000
      DO 102 KVIB=1,N3N                                                 ANH08010
      FRQK=FREQ(KVIB)                                                   ANH08020
      IF(DABS(FRQK).LE.WLIMIT) GO TO 102                                ANH08030
      FRQK2=FRQK*FRQK                                                   ANH08040
      FAC=F3Q(IVIB,IVIB,KVIB)                                           ANH08050
      FACT=FAC*FAC                                                      ANH08060
      IF(NDEG(KVIB).NE.1) GO TO 102                                     ANH08070
      IF(DABS(FRQI+FRQI-FRQK).LE.FLIM1) GO TO 206                       ANH08080
      DNUM=FRQI2*8.0D+00-FRQK2*3.0D+00                                  ANH08090
      DENM=FRQI2*4.0D+00-FRQK2                                          ANH08100
      VALU2=VALU2+(DNUM*FACT)/(FRQK*DENM)                               ANH08110
      GVAL2=GVAL2+(FACT*FRQK)/DENM                                      ANH08120
      GO TO 102                                                         ANH08130
C/////////////////////////////////                                      ANH08140
C///FERMI RESONANCE CORRECTIONS///                                      ANH08150
C/////////////////////////////////                                      ANH08160
  206 CONTINUE                                                          ANH08170
      DNUM=ONE/((FRQI+FRQI+FRQK)*2.0D+00)+2.0D+00/FRQK                  ANH08180
      IF(DABS(FAC).LE.PLIMIT) THEN                                      ANH08190
        IANH1=IANH1+1                                                   ANH08200
      END IF                                                            ANH08210
      VALU2=VALU2+FACT*DNUM                                             ANH08220
      GNUM=ONE/((FRQI+FRQI+FRQK)*2.0D+00)                               ANH08230
      GVAL2=GVAL2-FACT*GNUM                                             ANH08240
C////////////////////////////                                           ANH08250
C///END OF FERMI RESONANCE///                                           ANH08260
C////////////////////////////                                           ANH08270
  102 CONTINUE                                                          ANH08280
C   SUM OVER DEGENERATE MODES                                           ANH08290
      DO 103 KVIB=1,N3N                                                 ANH08300
      FRQK=FREQ(KVIB)                                                   ANH08310
      IF(DABS(FRQK).LE.WLIMIT) GO TO 103                                ANH08320
      FRQK2=FRQK*FRQK                                                   ANH08330
      FAC=F3Q(IVIB,IVIB,KVIB)                                           ANH08340
      FACT=FAC*FAC                                                      ANH08350
      IF(NDEG(KVIB).EQ.1) GO TO 103                                     ANH08360
      IF(DABS(FRQI+FRQI-FRQK).LE.FLIM1) GO TO 208                       ANH08370
      DNUM=FRQI2*8.0D+00-FRQK2*3.0D+00                                  ANH08380
      DNUG=FRQI2*8.0D+00-FRQK2                                          ANH08390
      DENM=FRQI2*4.0D+00-FRQK2                                          ANH08400
      VALU3=VALU3+(DNUM*FACT)/(FRQK*DENM)                               ANH08410
      GVAL3=GVAL3+(DNUG*FACT)/(FRQK*DENM)                               ANH08420
      GO TO 103                                                         ANH08430
C/////////////////////////////////                                      ANH08440
C///FERMI RESONANCE CORRECTIONS///                                      ANH08450
C/////////////////////////////////                                      ANH08460
  208 CONTINUE                                                          ANH08470
      DNUM=ONE/((FRQI+FRQI+FRQK)*2.0D+00)+2.0D+00/FRQK                  ANH08480
      IF(DABS(FAC).LE.PLIMIT) THEN                                      ANH08490
        IANH1=IANH1+1                                                   ANH08500
      END IF                                                            ANH08510
      VALU2=VALU2+FACT*DNUM                                             ANH08520
      GNUM=-ONE/((FRQI+FRQI+FRQK)*2.0D+00)+2.0D+00/FRQK                 ANH08530
      GVAL2=GVAL2+FACT*GNUM                                             ANH08540
C////////////////////////////                                           ANH08550
C///END OF FERMI RESONANCE///                                           ANH08560
C////////////////////////////                                           ANH08570
  103 CONTINUE                                                          ANH08580
      FACT4=ZETA(IVA,IVB,IZZ)**2                                        ANH08590
      GVAL4=FACT4*ROTCM(IZZ)                                            ANH08600
      XRS=(VALU1-VALU2-VALU3)/16.0D+00                                  ANH08610
      GRS=(GVAL1/3.0D+00-GVAL2+GVAL3)/16.0D+00+GVAL4                    ANH08620
      GO TO 250                                                         ANH08630
C                                                                       ANH08640
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++         ANH08650
C+++OFF-DIAGONAL CONSTANTS (NON-DEGENERATE & NON-DEGENERATE)+++         ANH08660
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++         ANH08670
  210 CONTINUE                                                          ANH08680
      IF(NDEG(IVIB).NE.1) GO TO 215                                     ANH08690
      IF(NDEG(JVIB).NE.1) GO TO 215                                     ANH08700
      VALU1=F4Q(IVIB,IVIB,JVIB,JVIB)                                    ANH08710
      DO 104 KVIB=1,NVIB                                                ANH08720
C   SUM OVER ONLY NON-DEGENERATE MODES                                  ANH08730
      FRQK=FREQ(KVIB)                                                   ANH08740
      IF(DABS(FRQK).LE.WLIMIT) GO TO 104                                ANH08750
      IF(NDEG(KVIB).NE.1) GO TO 104                                     ANH08760
      FRQK2=FRQK*FRQK                                                   ANH08770
      WPPP=FRQI+FRQJ+FRQK                                               ANH08780
      WPPM=FRQI+FRQJ-FRQK                                               ANH08790
      WPMP=FRQI-FRQJ+FRQK                                               ANH08800
      WPMM=FRQI-FRQJ-FRQK                                               ANH08810
      DIJK=WPPP*WPPM*WPMP*WPMM                                          ANH08820
      DNUM2=F3Q(IVIB,IVIB,KVIB)*F3Q(KVIB,JVIB,JVIB)                     ANH08830
      DNUM3=(FRQK2-FRQI2-FRQJ2)*FRQK                                    ANH08840
      FAC=F3Q(IVIB,JVIB,KVIB)                                           ANH08850
      FACT=FAC*FAC                                                      ANH08860
      VALU2=VALU2+DNUM2/FRQK                                            ANH08870
      IF(DABS(WPPM).LE.FLIM2) GO TO 211                                 ANH08880
      IF(DABS(WPMP).LE.FLIM2) GO TO 212                                 ANH08890
      IF(DABS(WPMM).LE.FLIM2) GO TO 213                                 ANH08900
      VALU3=VALU3+(DNUM3*FACT)/DIJK                                     ANH08910
      GO TO 104                                                         ANH08920
C/////////////////////////////////                                      ANH08930
C///FERMI RESONANCE CORRECTIONS///                                      ANH08940
C/////////////////////////////////                                      ANH08950
C   WI+WJ=WK                                                            ANH08960
  211 CONTINUE                                                          ANH08970
      DNUM3=ONE/WPPP+ONE/WPMP-ONE/WPMM                                  ANH08980
      IF(DABS(FAC).GT.PLIMIT) THEN                                      ANH08990
        IANH2=IANH2+1                                                   ANH09000
      END IF                                                            ANH09010
      VALU3=VALU3+(DNUM3*FACT)/4.0D+00                                  ANH09020
      GO TO 104                                                         ANH09030
C   WI+WK=WJ                                                            ANH09040
  212 CONTINUE                                                          ANH09050
      DNUM3=ONE/WPPP-ONE/WPPM-ONE/WPMM                                  ANH09060
      IF(DABS(FAC).GT.PLIMIT) THEN                                      ANH09070
        IANH2=IANH2+1                                                   ANH09080
      END IF                                                            ANH09090
      VALU3=VALU3+(DNUM3*FACT)/4.0D+00                                  ANH09100
      GO TO 104                                                         ANH09110
C   WJ+WK=WI                                                            ANH09120
  213 CONTINUE                                                          ANH09130
      DNUM3=ONE/WPPP-ONE/WPPM+ONE/WPMP                                  ANH09140
      IF(DABS(FAC).GT.PLIMIT) THEN                                      ANH09150
        IANH2=IANH2+1                                                   ANH09160
      END IF                                                            ANH09170
      VALU3=VALU3+(DNUM3*FACT)/4.0D+00                                  ANH09180
C////////////////////////////                                           ANH09190
C///END OF FERMI RESONANCE///                                           ANH09200
C////////////////////////////                                           ANH09210
  104 CONTINUE                                                          ANH09220
      ZFACT=ROTCM(IZZ)*ZETA(IVIB,JVIB,IZZ)**2                           ANH09230
      VALU4=(FRQI/FRQJ+FRQJ/FRQI)*ZFACT                                 ANH09240
      XRS=(VALU1-VALU2-VALU3*2.0D+00)/4.0D+00+VALU4                     ANH09250
      GO TO 250                                                         ANH09260
C                                                                       ANH09270
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++             ANH09280
C+++OFF-DIAGONAL CONSTANTS (DEGENERATE & NON-DEGENERATE)+++             ANH09290
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++             ANH09300
  215 CONTINUE                                                          ANH09310
      IF(NDEG(IVIB).NE.1.AND.NDEG(JVIB).NE.1) GO TO 225                 ANH09320
      VALU1=F4Q(IVIB,IVIB,JVIB,JVIB)                                    ANH09330
C   SUM OVER NON-DEGENERATE MODES                                       ANH09340
      DO 105 KVIB=1,NVIB                                                ANH09350
      FRQK=FREQ(KVIB)                                                   ANH09360
      IF(DABS(FRQK).LE.WLIMIT) GO TO 105                                ANH09370
      IF(NDEG(KVIB).NE.1) GO TO 105                                     ANH09380
      DNUM2=F3Q(IVIB,IVIB,KVIB)*F3Q(KVIB,JVIB,JVIB)                     ANH09390
      VALU2=VALU2+DNUM2/FRQK                                            ANH09400
  105 CONTINUE                                                          ANH09410
C  SUM OVER DEGENERATE MODES                                            ANH09420
      DO 106 KVIB=1,NVIB                                                ANH09430
      FRQK=FREQ(KVIB)                                                   ANH09440
      IF(DABS(FRQK).LE.WLIMIT) GO TO 106                                ANH09450
      FRQK2=FRQK*FRQK                                                   ANH09460
      IF(NDEG(KVIB).EQ.1) GO TO 106                                     ANH09470
      WPPP=FRQI+FRQJ+FRQK                                               ANH09480
      WPPM=FRQI+FRQJ-FRQK                                               ANH09490
      WPMP=FRQI-FRQJ+FRQK                                               ANH09500
      WPMM=FRQI-FRQJ-FRQK                                               ANH09510
      DIJK=WPPP*WPPM*WPMP*WPMM                                          ANH09520
      DNUM3=(FRQK2-FRQI2-FRQJ2)*FRQK                                    ANH09530
      FAC=F3Q(IVIB,JVIB,KVIB)                                           ANH09540
      FACT=FAC*FAC                                                      ANH09550
      IF(DABS(WPPM).LE.FLIM2) GO TO 217                                 ANH09560
      IF(DABS(WPMP).LE.FLIM2) GO TO 218                                 ANH09570
      IF(DABS(WPMM).LE.FLIM2) GO TO 219                                 ANH09580
      VALU3=VALU3+(DNUM3*FACT)/DIJK                                     ANH09590
      GO TO 106                                                         ANH09600
C/////////////////////////////////                                      ANH09610
C///FERMI RESONANCE CORRECTIONS///                                      ANH09620
C/////////////////////////////////                                      ANH09630
C   WI+WJ=WK                                                            ANH09640
  217 CONTINUE                                                          ANH09650
      DNUM3=ONE/WPPP+ONE/WPMP-ONE/WPMM                                  ANH09660
      IF(DABS(FAC).GT.PLIMIT) THEN                                      ANH09670
        IANH2=IANH2+1                                                   ANH09680
      END IF                                                            ANH09690
      VALU3=VALU3+(DNUM3*FACT)/4.0D+00                                  ANH09700
      GO TO 106                                                         ANH09710
C   WI+WK=WJ                                                            ANH09720
  218 CONTINUE                                                          ANH09730
      DNUM3=ONE/WPPP-ONE/WPPM-ONE/WPMM                                  ANH09740
      IF(DABS(FAC).GT.PLIMIT) THEN                                      ANH09750
        IANH2=IANH2+1                                                   ANH09760
      END IF                                                            ANH09770
      VALU3=VALU3+(DNUM3*FACT)/4.0D+00                                  ANH09780
      GO TO 106                                                         ANH09790
C   WJ+WK=WI                                                            ANH09800
  219 CONTINUE                                                          ANH09810
      DNUM3=ONE/WPPP-ONE/WPPM+ONE/WPMP                                  ANH09820
      IF(DABS(FAC).GT.PLIMIT) THEN                                      ANH09830
        IANH2=IANH2+1                                                   ANH09840
      END IF                                                            ANH09850
      VALU3=VALU3+(DNUM3*FACT)/4.0D+00                                  ANH09860
C//////////////////////                                                 ANH09870
C///END OF RESONANCE///                                                 ANH09880
C//////////////////////                                                 ANH09890
  106 CONTINUE                                                          ANH09900
      ZFACT=(ZETA(IVIB,JVIB,IXX)**2+ZETA(IVIB,JVIB,IYY)**2)*ROTCM(IXX)  ANH09910
      VALU4=(FRQI/FRQJ+FRQJ/FRQI)*ZFACT                                 ANH09920
      XRS=(VALU1-VALU2-VALU3*2.0D+00)/4.0D+00+VALU4                     ANH09930
      GO TO 250                                                         ANH09940
C                                                                       ANH09950
C++++++++++++++++++++++++++++++++++++++++++++++++++++++                 ANH09960
C+++OFF-DIAGONAL CONSTANTS (DEGENERATE & DEGENERATE)+++                 ANH09970
C++++++++++++++++++++++++++++++++++++++++++++++++++++++                 ANH09980
  225 CONTINUE                                                          ANH09990
      IF(DABS(FREQ(IVIB)-FREQ(JVIB)).LE.DLIMIT) GO TO 250               ANH10000
      VALU1=(F4Q(IVIB,IVIB,JVA,JVA)+F4Q(IVIB,IVIB,JVB,JVB))/2.0D+00     ANH10010
C  SUM OVER NON-DEGENERATE MODES                                        ANH10020
      DO 109 KVIB=1,NVIB                                                ANH10030
      FRQK=FREQ(KVIB)                                                   ANH10040
      IF(DABS(FRQK).LE.WLIMIT) GO TO 109                                ANH10050
      IF(NDEG(KVIB).NE.1) GO TO 109                                     ANH10060
      FRQK2=FRQK*FRQK                                                   ANH10070
      WPPP=FRQI+FRQJ+FRQK                                               ANH10080
      WPPM=FRQI+FRQJ-FRQK                                               ANH10090
      WPMP=FRQI-FRQJ+FRQK                                               ANH10100
      WPMM=FRQI-FRQJ-FRQK                                               ANH10110
      DIJK=WPPP*WPPM*WPMP*WPMM                                          ANH10120
      DNUM2=F3Q(IVIB,IVIB,KVIB)*F3Q(KVIB,JVIB,JVIB)                     ANH10130
      DNUM3=(FRQK2-FRQI2-FRQJ2)*FRQK                                    ANH10140
      FAC=F3Q(IVIB,JVA,KVIB)+F3Q(IVIB,JVB,KVIB)                         ANH10150
      FACT=FAC*FAC                                                      ANH10160
      VALU2=VALU2+DNUM2/FRQK                                            ANH10170
      IF(DABS(WPPM).LE.FLIM2) GO TO 226                                 ANH10180
      IF(DABS(WPMP).LE.FLIM2) GO TO 227                                 ANH10190
      IF(DABS(WPMM).LE.FLIM2) GO TO 228                                 ANH10200
      VALU3=VALU3+(DNUM3*FACT)/DIJK                                     ANH10210
      GVAL1=GVAL1+(FRQI*FRQJ*FRQK*FACT)/DIJK                            ANH10220
      GO TO 109                                                         ANH10230
C/////////////////////////////////                                      ANH10240
C///FERMI RESONANCE CORRECTIONS///                                      ANH10250
C/////////////////////////////////                                      ANH10260
C   WI+WJ=WK                                                            ANH10270
  226 CONTINUE                                                          ANH10280
      DNUM3=ONE/WPPP+ONE/WPMP-ONE/WPMM                                  ANH10290
      IF(DABS(FAC).GT.PLIMIT) THEN                                      ANH10300
        IANH2=IANH2+1                                                   ANH10310
      END IF                                                            ANH10320
      VALU3=VALU3+(DNUM3*FACT)/4.0D+00                                  ANH10330
      GFACT=ONE/WPPP-ONE/WPMP+ONE/WPMM                                  ANH10340
      GVAL1=GVAL1+(GFACT*FACT)/8.0D+00                                  ANH10350
      GO TO 109                                                         ANH10360
C   WI+WK=WJ                                                            ANH10370
  227 CONTINUE                                                          ANH10380
      DNUM3=ONE/WPPP-ONE/WPPM-ONE/WPMM                                  ANH10390
      IF(DABS(FAC).GT.PLIMIT) THEN                                      ANH10400
        IANH2=IANH2+1                                                   ANH10410
      END IF                                                            ANH10420
      VALU3=VALU3+(DNUM3*FACT)/4.0D+00                                  ANH10430
      GFACT=ONE/WPPP-ONE/WPPM+ONE/WPMM                                  ANH10440
      GVAL1=GVAL1+(GFACT*FACT)/8.0D+00                                  ANH10450
      GO TO 109                                                         ANH10460
C   WJ+WK=WI                                                            ANH10470
  228 CONTINUE                                                          ANH10480
      DNUM3=ONE/WPPP-ONE/WPPM+ONE/WPMP                                  ANH10490
      IF(DABS(FAC).GT.PLIMIT) THEN                                      ANH10500
        IANH2=IANH2+1                                                   ANH10510
      END IF                                                            ANH10520
      VALU3=VALU3+(DNUM3*FACT)/4.0D+00                                  ANH10530
      GFACT=ONE/WPPP-ONE/WPPM-ONE/WPMP                                  ANH10540
      GVAL1=GVAL1+(GFACT*FACT)/8.0D+00                                  ANH10550
C//////////////////////                                                 ANH10560
C///END OF RESONANCE///                                                 ANH10570
C//////////////////////                                                 ANH10580
  109 CONTINUE                                                          ANH10590
C                                                                       ANH10600
C   SUM OVER DEGENERATE MODES                                           ANH10610
      DO 110 KVIB=1,NVIB                                                ANH10620
      FRQK=FREQ(KVIB)                                                   ANH10630
      IF(DABS(FRQK).LE.WLIMIT) GO TO 110                                ANH10640
      IF(NDEG(KVIB).EQ.1) GO TO 110                                     ANH10650
      FRQK2=FRQK*FRQK                                                   ANH10660
      KVA=NDAB(KVIB,1)                                                  ANH10670
      KVB=NDAB(KVIB,2)                                                  ANH10680
      WPPP=FRQI+FRQJ+FRQK                                               ANH10690
      WPPM=FRQI+FRQJ-FRQK                                               ANH10700
      WPMP=FRQI-FRQJ+FRQK                                               ANH10710
      WPMM=FRQI-FRQJ-FRQK                                               ANH10720
      DIJK=WPPP*WPPM*WPMP*WPMM                                          ANH10730
      DNUM4=(FRQK2-FRQI2-FRQJ2)*FRQK                                    ANH10740
      FAC=(DABS(F3Q(IVA,JVA,KVIB))+DABS(F3Q(IVA,JVB,KVIB))              ANH10750
     1    +DABS(F3Q(IVB,JVA,KVIB))+DABS(F3Q(IVB,JVB,KVIB)))/2.0D+00     ANH10760
      FACT=FAC*FAC                                                      ANH10770
      IF(DABS(WPPM).LE.FLIM2) GO TO 230                                 ANH10780
      IF(DABS(WPMP).LE.FLIM2) GO TO 231                                 ANH10790
      IF(DABS(WPMM).LE.FLIM2) GO TO 232                                 ANH10800
      VALU4=VALU4+(DNUM4*FACT)/DIJK                                     ANH10810
      GVAL2=GVAL2-(FRQI*FRQJ*FRQK*FACT)/DIJK                            ANH10820
      GO TO 110                                                         ANH10830
C/////////////////////////////////                                      ANH10840
C///FERMI RESONANCE CORRECTIONS///                                      ANH10850
C/////////////////////////////////                                      ANH10860
C   WI+WJ=WK                                                            ANH10870
  230 CONTINUE                                                          ANH10880
      DNUM4=ONE/WPPP+ONE/WPMP-ONE/WPMM                                  ANH10890
      IF(DABS(FAC).GT.PLIMIT) THEN                                      ANH10900
        IANH2=IANH2+1                                                   ANH10910
      END IF                                                            ANH10920
      VALU4=VALU4+(DNUM4*FACT)/4.0D+00                                  ANH10930
      GFACT=ONE/WPPP-ONE/WPMP+ONE/WPMM                                  ANH10940
      GVAL2=GVAL2+(GFACT*FACT)/8.0D+00                                  ANH10950
      GO TO 110                                                         ANH10960
C   WI+WK=WJ                                                            ANH10970
  231 CONTINUE                                                          ANH10980
      DNUM4=ONE/WPPP-ONE/WPPM-ONE/WPMM                                  ANH10990
      IF(DABS(FAC).GT.PLIMIT) THEN                                      ANH11000
        IANH2=IANH2+1                                                   ANH11010
      END IF                                                            ANH11020
      VALU4=VALU4+(DNUM4*FACT)/4.0D+00                                  ANH11030
      GFACT=ONE/WPPP-ONE/WPPM+ONE/WPMM                                  ANH11040
      GVAL2=GVAL2+(GFACT*FACT)/8.0D+00                                  ANH11050
      GO TO 110                                                         ANH11060
C   WJ+WK=WI                                                            ANH11070
  232 CONTINUE                                                          ANH11080
      DNUM4=ONE/WPPP-ONE/WPPM+ONE/WPMP                                  ANH11090
      IF(DABS(FAC).GT.PLIMIT) THEN                                      ANH11100
        IANH2=IANH2+1                                                   ANH11110
      END IF                                                            ANH11120
      VALU4=VALU4+(DNUM4*FACT)/4.0D+00                                  ANH11130
      GFACT=ONE/WPPP-ONE/WPPM-ONE/WPMP                                  ANH11140
      GVAL2=GVAL2+(GFACT*FACT)/8.0D+00                                  ANH11150
C//////////////////////                                                 ANH11160
C///END OF RESONANCE///                                                 ANH11170
C//////////////////////                                                 ANH11180
  110 CONTINUE                                                          ANH11190
      ZETAZ=ZETA(IVIB,JVA,IZZ)+ZETA(IVIB,JVB,IZZ)                       ANH11200
      ZETAY=ZETA(IVIB,JVA,IYY)+ZETA(IVIB,JVB,IYY)                       ANH11210
      ZFACT=ZETAZ**2*0.5D+00*ROTCM(IZZ)                                 ANH11220
     1     +ZETAY**2*ROTCM(IYY)                                         ANH11230
      VALU5=(FRQI/FRQJ+FRQJ/FRQI)*ZFACT                                 ANH11240
      GVAL3=-ZETAY**2*ROTCM(IYY)*2.0D+00                                ANH11250
     1      +ZETAZ**2*ROTCM(IZZ)                                        ANH11260
     2      +ZETA(IVA,IVB,IZZ)*ZETA(JVA,JVB,IZZ)*ROTCM(IZZ)*2.0D+00     ANH11270
C*10/20/89  XRS=(VALU1-VALU2-VALU3-VALU4*2.0D+00)/4.0D+00+VALU5         ANH11280
      XRS=(VALU1-VALU2-VALU3-VALU4)/4.0D+00+VALU5                       ANH11290
C*10/20/89  GRS=GVAL1*0.5D+00+GVAL2+GVAL3                               ANH11300
      GRS=GVAL1*0.5D+00+GVAL2*0.5D+00+GVAL3                             ANH11310
C                                                                       ANH11320
  250 CONTINUE                                                          ANH11330
      XIJ(IVIB,JVIB)=XRS                                                ANH11340
      GIJ(IVIB,JVIB)=GRS                                                ANH11350
      IF(IPRNT.LE.2) GO TO 120                                          ANH11360
      WRITE(6,3) IVIB,JVIB,FRQI,FRQJ,VALU1,VALU2,VALU3,VALU4,XRS,GRS    ANH11370
  120 CONTINUE                                                          ANH11380
      IF(IANH1.NE.0) THEN                                               ANH11390
        WRITE(6,8)                                                      ANH11400
        WRITE(6,10) FLIM1                                               ANH11410
      END IF                                                            ANH11420
      IF(IANH2.NE.0) THEN                                               ANH11430
        WRITE(6,9)                                                      ANH11440
        WRITE(6,10) FLIM2                                               ANH11450
      END IF                                                            ANH11460
C                                                                       ANH11470
      IF(IPRNT.GT.3) GO TO 251                                          ANH11480
      WRITE(6,4)                                                        ANH11490
      DO 122 II=1,NVIB                                                  ANH11500
      IVIB=NFRQ(II)                                                     ANH11510
      FRQI=FREQ(IVIB)                                                   ANH11520
      DO 122 JJ=1,NVIB                                                  ANH11530
      JVIB=NFRQ(JJ)                                                     ANH11540
      FRQJ=FREQ(JVIB)                                                   ANH11550
      XRS=XIJ(IVIB,JVIB)                                                ANH11560
      GRS=GIJ(IVIB,JVIB)                                                ANH11570
      WRITE(6,5) II,JJ,FRQI,FRQJ,XRS,GRS                                ANH11580
  122 CONTINUE                                                          ANH11590
C                                                                       ANH11600
C:::::::::::::::::::::::::::::::::::::::                                ANH11610
C:::CALCULATE FUNDAMENTAL FREQUENCIES:::                                ANH11620
C:::::::::::::::::::::::::::::::::::::::                                ANH11630
  251 CONTINUE                                                          ANH11640
      WRITE(6,6)                                                        ANH11650
      DO 125 II=1,NVIB                                                  ANH11660
      IVIB=NFRQ(II)                                                     ANH11670
      DFAC=DBLE(NDEG(IVIB))                                             ANH11680
      VALU1=FREQ(IVIB)                                                  ANH11690
      VALU2=XIJ(IVIB,IVIB)*(ONE+DFAC)                                   ANH11700
      VALU3=A00                                                         ANH11710
      VALU4=GIJ(IVIB,IVIB)                                              ANH11720
      DO 124 JJ=1,NVIB                                                  ANH11730
      JVIB=NFRQ(JJ)                                                     ANH11740
      IF(IVIB.EQ.JVIB) GO TO 124                                        ANH11750
      VALU3=VALU3+XIJ(IVIB,JVIB)                                        ANH11760
  124 CONTINUE                                                          ANH11770
      DELTA=VALU2+VALU3*0.5D+00+VALU4                                   ANH11780
      VALUT=VALU1+DELTA                                                 ANH11790
      WRITE(6,7) II,VALU1,DELTA,VALUT                                   ANH11800
  125 CONTINUE                                                          ANH11810
C
C::::::::::::::::::::::::::::::::::::::::::::::::::
C:::CALCULATE ANHARMONIC VIBRATIONAL ENERGY :::::::
C::::::::::::::::::::::::::::::::::::::::::::::::::
C   ZPVE = 1/2 Sum_(i)^{unique} d_i W_i
C        + 1/4 Sum_(i>=j)^{unique} d_i d_j X_(ij)
C::::::::::::::::::::::::::::::::::::::::::::::::::
      TERM1=A00
      TERM2=A00
      TERM3=A00
      TERM4=A00
      TERM5=A00
      TERM6=A00
      TERM23=A00
      TERM45=A00
      TERM456=A00
      DO 127 IVIB=1,NVIB
      FRQI=FREQ(IVIB)
      NDII=NDEG(IVIB)
      FACI=DFLOAT(NDII)
      TERM1=TERM1+FRQI/(2.0D+00)
      XRR=XIJ(IVIB,IVIB)
      IF(NDII.EQ.1) THEN
        TERM2=TERM2+XRR/(4.0D+00)
      END IF
      IF(NDII.EQ.2) THEN
        IVA=NDAB(IVIB,1)
        IVB=NDAB(IVIB,2)
        TERM3=TERM3+XRR/(2.0D+00)
      END IF
      TERM23=TERM23+(FACI*XRR)/(4.0D+00)
      DO 126 JVIB=1,NVIB
      IF(IVIB.EQ.JVIB) GO TO 126
      NDJJ=NDEG(JVIB)
      FACJ=DFLOAT(NDJJ)
      IF(NDJJ.EQ.2) THEN
        JVA=NDAB(JVIB,1)
        JVB=NDAB(JVIB,2)
      END IF
      XRS=XIJ(IVIB,JVIB)
      IF(NDII.EQ.1 .AND. NDJJ.EQ.1) THEN
        TERM4=TERM4+XRS/(8.0D+00)
      END IF
      IF(NDII.EQ.1 .AND. NDJJ.EQ.2) THEN
        TERM5=TERM5+XRS/(8.0D+00)
      END IF
      IF(NDII.EQ.2 .AND. NDJJ.EQ.1) THEN
        TERM5=TERM5+XRS/(8.0D+00)
      END IF
      IF(NDII.EQ.2 .AND. NDJJ.EQ.2) THEN
        IF(IVA.NE.JVA .AND. IVB.NE.JVB) THEN
          TERM6=TERM6+XRS/(8.0D+00)
        END IF
      ELSE
        TERM45=TERM45+XRS/(8.0D+00)
      END IF
      TERM456=TERM456+XRS/(8.0D+00)
  126 CONTINUE
  127 CONTINUE
C
      ANHTOT=TERM1+TERM2+TERM3
     *      +TERM4+TERM5+TERM6
      IF(IPRNT.GT.10) THEN
        WRITE(6,*) ' TERM1 (Q2)        = ',TERM1
        WRITE(6,*) ' TERM2 (XRR-ND)    = ',TERM2
        WRITE(6,*) ' TERM3 (XRR-DG)    = ',TERM3
        WRITE(6,*) ' TERM4 (XRS-ND_ND) = ',TERM4
        WRITE(6,*) ' TERM5 (XRS-ND_DG) = ',TERM5
        WRITE(6,*) ' TERM6 (XRS-DG_DG) = ',TERM6
        WRITE(6,*) ' ANHTOT            = ',ANHTOT
        WRITE(6,*) ' TERM2+TERM3       = ',TERM2+TERM3
        WRITE(6,*) ' TERM4+TERM5       = ',TERM4+TERM5
        WRITE(6,*) ' TERM4+TERM5+TERM6 = ',TERM4+TERM5+TERM6
        WRITE(6,*) ' TERM23  (2+3)     = ',TERM23
        WRITE(6,*) ' TERM45  (4+5)     = ',TERM45
        WRITE(6,*) ' TERM456 (4+5+6)   = ',TERM456
      END IF
C
      RETURN
      END
      SUBROUTINE SYMZPE(ZETA,F3Q,F4Q) 
C   THE ANHARMONIC ZPVE FOR A SYMMETRIC TOP     
      IMPLICIT REAL*8 (A-H,O-Z)                                      
      DIMENSION ZETA(N3N,N3N,3)                                     
      DIMENSION F3Q(N3N,N3N,N3N),F4Q(N3N,N3N,N3N,N3N)              
      COMMON/VIB101/NATOM,N3N,NATRI,ILIN,NVIB                     
      COMMON/VIB106/ROTAA(3),ROTGC(3),ROTCM(3),ROTMH(3)               
      COMMON/VIB108/IAXIS(3),NDEG(150),NDAB(150,5),IMAG(150)         
      COMMON/VIB111/ANHTOT,HZPVE,AZPVE
      COMMON/VIB203/IOFF(150),IPRNT                                
      COMMON/VIB204/SQM(150),ROOT(150),FREQ(150)                  
      DATA A00,ONE / 0.0D+00 , 1.0D+00 /                            
    1 FORMAT(//,2X,' **************************************************' 
     1        /,2X,' ***THE ANHARMONIC ZERO-POINT VIBRATIONAL ENERGY***'   
     2        /,2X,' **************************************************'  
     3      )
    2 FORMAT(//,2X,' FOR A SYMMETRIC TOP MOLECULE'/)
    3 FORMAT(   2X,' VALU1      = ',F20.10,' CM-1',/
     1          2X,' VALU2      = ',F20.10,' CM-1',/
     2          2X,' VALU3      = ',F20.10,' CM-1',/
     3          2X,' VALU4      = ',F20.10,' CM-1',/
     4          2X,' VALU5      = ',F20.10,' CM-1',/)
    4 FORMAT(   2X,' HARM. (Q2) = ',F20.10,' CM-1',/
     1          2X,' CUBIC (Q3) = ',F20.10,' CM-1',/
     2          2X,' QUART.(Q4) = ',F20.10,' CM-1',/
     3          2X,' Q3+Q4      = ',F20.10,' CM-1',/
     4          2X,' KINET.(K)  = ',F20.10,' CM-1',/
     5          2X,' TOTAL      = ',F20.10,' CM-1',/
     6          2X,' (Q2+Q3+Q4+K)'/)
C
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C:::CALCULATE ANHARMONIC VIBRATIONAL ENERGY ::::::::::::::::::
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C   ZPVE = 1/2 Sum_(i)^{all} W_i
C        - 1/32 Sum_(ijk)^{all} (phi_{iik})(phi_{kjj})/(W_k)
C        - 1/48 Sum_(ijk)^{all} (phi_{ijk})^{2}/(W_i+W_j+W_k)
C        + 1/32 Sum_(ij)^{all} phi_{iijj}
C        + Z_{kinetic}
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      WRITE(6,1)
      VALU1=A00
      VALU2=A00
      VALU3=A00
      VALU4=A00
      VALU5=A00
      DO 103 IVIB=1,NVIB
      FRQI=FREQ(IVIB)
      VALU1=VALU1+FRQI
      DO 102 JVIB=1,NVIB
      FRQJ=FREQ(JVIB)
      FRIIJJ=F4Q(IVIB,IVIB,JVIB,JVIB)
      VALU4=VALU4+FRIIJJ
      DO 101 KVIB=1,NVIB
      FRQK=FREQ(KVIB)
      FRIIK=F3Q(IVIB,IVIB,KVIB)
      FRKJJ=F3Q(KVIB,JVIB,JVIB)
      FRIJK=F3Q(IVIB,JVIB,KVIB)
      VALU2=VALU2+(FRIIK*FRKJJ)/(FRQK)
      VALU3=VALU3+(FRIJK*FRIJK)/(FRQI+FRQJ+FRQK)
  101 CONTINUE
  102 CONTINUE
  103 CONTINUE
C
C==========================================
C   Following codes are equivalent to 
C   other codes described below
C==========================================
c     DO 107 IABC=1,3
c     II=IAXIS(IABC)
c     BE=ROTCM(IABC)
c     TERM=ONE
c     DO 106 IVIB=1,NVIB
c     FRQI=FREQ(IVIB)
c     DO 105 JVIB=1,NVIB
c     IF(IVIB.LE.JVIB) GO TO 105
c     FRQJ=FREQ(JVIB)
c     ZET=ZETA(IVIB,JVIB,II)
c     ZETSQ=ZET*ZET
c     WIPJ=FRQI+FRQJ
c     WIMJ=FRQI-FRQJ
c     WIXJ=FRQI*FRQJ
c     WIMJSQ=WIMJ*WIMJ
c     TERM=TERM+ZETSQ*(BE*WIPJ-WIMJSQ)/(WIXJ)
c 105 CONTINUE
c 106 CONTINUE
c     VALU5=VALU5-TERM*BE
c 107 CONTINUE
C==========================================
C==========================================
      DO 107 IABC=1,3
      II=IAXIS(IABC)
      BE=ROTCM(IABC)
      TERM=-ONE
      DO 106 IVIB=1,NVIB
      FRQI=FREQ(IVIB)
      DO 105 JVIB=1,NVIB
      FRQJ=FREQ(JVIB)
      ZET=ZETA(IVIB,JVIB,II)
      ZETSQ=ZET*ZET
      WIMJ=FRQI-FRQJ
      TERM=TERM+ZETSQ*(WIMJ-BE)/(FRQJ)
  105 CONTINUE
  106 CONTINUE
      VALU5=VALU5+TERM*BE
  107 CONTINUE
C==========================================
C
      VALU1 =+VALU1/(2.0D+00)
      VALU2 =-VALU2/(32.0D+00)
      VALU3 =-VALU3/(48.0D+00)
      VALU4 =+VALU4/(32.0D+00)
      VALU5 =+VALU5/(4.0D+00)
      HZPVE = VALU1
      ZCUBIC= VALU2+VALU3
      ZQUART= VALU4
      ZANHR = ZCUBIC+ZQUART
      ZKINET= VALU5
      AZPVE = HZPVE+ZANHR+ZKINET
      WRITE(6,2)
      IF(IPRNT.GT.3)
     *WRITE(6,3) VALU1,VALU2,VALU3,VALU4,VALU5
      WRITE(6,4) HZPVE,ZCUBIC,ZQUART,ZANHR,ZKINET,AZPVE
C
      RETURN
      END
      SUBROUTINE SYMG00(ZETA,F3Q,F4Q)
C   THE VIBRATIONAL G0 VALUE FOR A SYMMETRIC TOP
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION ZETA(N3N,N3N,3)
      DIMENSION F3Q(N3N,N3N,N3N),F4Q(N3N,N3N,N3N,N3N)
      COMMON/VIB101/NATOM,N3N,NATRI,ILIN,NVIB
      COMMON/VIB106/ROTAA(3),ROTGC(3),ROTCM(3),ROTMH(3)
      COMMON/VIB108/IAXIS(3),NDEG(150),NDAB(150,5),IMAG(150)
      COMMON/VIB109/CLIMIT,FLIM1,FLIM2
      COMMON/VIB111/ANHTOT,HZPVE,AZPVE
      COMMON/VIB203/IOFF(150),IPRNT
      COMMON/VIB204/SQM(150),ROOT(150),FREQ(150)
      DATA A00,ONE / 0.0D+00 , 1.0D+00 /
      DATA WLIMIT / 1.0D+00 /
      DATA PLIMIT / 1.0D-05 /
    1 FORMAT(//,2X,' ********************************'
     1        /,2X,' ***THE G_0 VIBRATIONAL ENERGY***'
     2        /,2X,' ********************************'
     3      )
    2 FORMAT(//,2X,' :::::::::::::::::::::::::::::::::::::::::::::::'/
     1          2X,' :::ANHARMONIC RESONANCE (TYPE I) IS OBSERVED:::'/
     2          2X,' :::::::::::::::::::::::::::::::::::::::::::::::'/)
    3 FORMAT(//,2X,' ::::::::::::::::::::::::::::::::::::::::::::::::'/
     1          2X,' :::ANHARMONIC RESONANCE (TYPE II) IS OBSERVED:::'/
     2          2X,' ::::::::::::::::::::::::::::::::::::::::::::::::'/)
    4 FORMAT(2X,' NUMBER OF OCCURENCES IS ',I5/
     1       2X,' RESONANCE LIMIT IS SET TO ',F10.3,' CM-1'/)
    5 FORMAT(//,2X,' FOR A SYMMETRIC TOP MOLECULE')
    6 FORMAT(   2X,' VALU1      = ',F20.10,' CM-1',/
     1          2X,' VALU2      = ',F20.10,' CM-1',/
     2          2X,' VALU3      = ',F20.10,' CM-1',/
     3          2X,' VALU4      = ',F20.10,' CM-1',/
     4          2X,' VALU5      = ',F20.10,' CM-1',/ 
     5          2X,' VALU6      = ',F20.10,' CM-1',/ 
     6          2X,' VALU7      = ',F20.10,' CM-1',/ 
     7          2X,' VALU8      = ',F20.10,' CM-1',/ 
     8          2X,' VALU9      = ',F20.10,' CM-1',/
     9          2X,' VALU10     = ',F20.10,' CM-1',/
     A          2X,' VALU11     = ',F20.10,' CM-1',/)
    7 FORMAT(   2X,' CUBIC (Q3) = ',F20.10,' CM-1',/
     1          2X,' QUART.(Q4) = ',F20.10,' CM-1',/
     2          2X,' KINET.(K)  = ',F20.10,' CM-1',/
     3          2X,' G_0        = ',F20.10,' CM-1',/
     4          2X,' (Q3+Q4+K)'/)
    8 FORMAT(   2X,' INTERNAL CONSISTENCY TEST',/
     1          2X,' ZPVE(Q)    = ',F20.10,' CM-1',/
     2          2X,' ZPVE(TOTAL)= ',F20.10,' CM-1',//
     3          2X,' ZPVE(XRS)  = ',F20.10,' CM-1',/
     4          2X,' ZPVE(ALL)  = ',F20.10,' CM-1',/
     5          2X,' (XRS+G_0)'/)
C
      WRITE(6,1)
C=============================================
C  For non-degenerate * non-degenerate modes
C=============================================
      IANH1=0
      IANH2=0
      VALU1=A00
      VALU2=A00
      VALU3=A00
      VALU4=A00
      DO 103 IVIB=1,NVIB
      IF(NDEG(IVIB).NE.1) GO TO 103
      FRQI=FREQ(IVIB)
      FRIII=F3Q(IVIB,IVIB,IVIB)
      FRIIII=F4Q(IVIB,IVIB,IVIB,IVIB)
      VALU1=VALU1+FRIIII
      VALU2=VALU2+(FRIII*FRIII)/(FRQI)
      FRISQ4=FRQI*FRQI*4.0D+00
      DO 102 JVIB=1,NVIB
      IF(NDEG(JVIB).NE.1) GO TO 102
      FRQJ=FREQ(JVIB)
      FRIIJ=F3Q(IVIB,IVIB,JVIB)
      FIIJSQ=FRIIJ*FRIIJ
      FACT3=FRIIJ*FRIIJ*FRQJ
      IF(DABS(FRQI+FRQI-FRQJ).LE.FLIM1) GO TO 201
      FRQJSQ=FRQJ*FRQJ
      DENOM=FRISQ4-FRQJSQ
      VALU3=VALU3+(FACT3)/(DENOM)
      GO TO 202
C/////////////////////////////////
C///FERMI RESONANCE CORRECTIONS///
C/////////////////////////////////
C   WI+WI=WJ                          
  201 CONTINUE
      IF(DABS(FRIIJ).GT.PLIMIT) THEN
        IANH1=IANH1+1
      END IF
      DNUM=ONE/(FRQI+FRQI+FRQJ)
      VALU3=VALU3-(FIIJSQ*DNUM)/(2.0D+00)
  202 CONTINUE
      DO 101 KVIB=1,NVIB
      IF(NDEG(KVIB).NE.1) GO TO 101
      IF(IVIB.LE.JVIB) GO TO 101
      IF(IVIB.LE.KVIB) GO TO 101
      IF(JVIB.LE.KVIB) GO TO 101
      FRQK=FREQ(KVIB)
      FRIJK=F3Q(IVIB,JVIB,KVIB)
      FACT4=FRIJK*FRIJK*FRQI*FRQJ*FRQK
      FACTSQ=FRIJK*FRIJK
      WPPP=FRQI+FRQJ+FRQK
      WPPM=FRQI+FRQJ-FRQK
      WPMP=FRQI-FRQJ+FRQK
      WPMM=FRQI-FRQJ-FRQK
      IF(DABS(WPPM).LE.FLIM2) GO TO 203
      IF(DABS(WPMP).LE.FLIM2) GO TO 204
      IF(DABS(WPMM).LE.FLIM2) GO TO 205
      DIJK=WPPP*WPPM*WPMP*WPMM
      VALU4=VALU4+(FACT4)/(DIJK)
      GO TO 101
C/////////////////////////////////       
C///FERMI RESONANCE CORRECTIONS///      
C/////////////////////////////////     
C   WI+WJ=WK                          
  203 CONTINUE
      DNUM4=ONE/WPPP-ONE/WPMP+ONE/WPMM
      IF(DABS(FRIJK).GT.PLIMIT) THEN
        IANH2=IANH2+1
      END IF
      VALU4=VALU4-(DNUM4*FACTSQ)/(8.0D+00)
      GO TO 101
C   WI+WK=WJ                          
  204 CONTINUE
      DNUM4=ONE/WPPP-ONE/WPPM+ONE/WPMM
      IF(DABS(FRIJK).GT.PLIMIT) THEN
        IANH2=IANH2+1
      END IF
      VALU4=VALU4-(DNUM4*FACTSQ)/(8.0D+00)
      GO TO 101
C   WJ+WK=WI
  205 CONTINUE
      DNUM4=ONE/WPPP-ONE/WPPM-ONE/WPMP
      IF(DABS(FRIJK).GT.PLIMIT) THEN
        IANH2=IANH2+1
      END IF
      VALU4=VALU4-(DNUM4*FACTSQ)/(8.0D+00)
  101 CONTINUE
  102 CONTINUE
  103 CONTINUE
C
      IF(IANH1.NE.0) THEN
        WRITE(6,2) 
        WRITE(6,4) IANH1,FLIM1
      END IF
      IF(IANH2.NE.0) THEN
        WRITE(6,3)
        WRITE(6,4) IANH2,FLIM2
      END IF
C
C=============================================
C  For degenerate * non-degenerate modes
C=============================================
      IANH1=0
      IANH2=0
      VALU5=A00
      VALU6=A00
      VALU7=A00
      VALU8=A00
      DO 106 IVIB=1,NVIB
      IF(NDEG(IVIB).NE.2) GO TO 106
      FRQI=FREQ(IVIB)
      FRIII=F3Q(IVIB,IVIB,IVIB)
      FRIIII=F4Q(IVIB,IVIB,IVIB,IVIB)
      VALU5=VALU5+FRIIII
      VALU6=VALU6+(FRIII*FRIII)/(FRQI)
      FRISQ4=FRQI*FRQI*4.0D+00
      DO 105 JVIB=1,NVIB
      IF(NDEG(JVIB).NE.1) GO TO 105
      FRQJ=FREQ(JVIB)
      FRIIJ=F3Q(IVIB,IVIB,JVIB)
      FIIJSQ=FRIIJ*FRIIJ
      FACT7=FRIIJ*FRIIJ*FRQJ
      IF(DABS(FRQI+FRQI-FRQJ).LE.FLIM1) GO TO 206
      FRQJSQ=FRQJ*FRQJ
      DENOM=FRISQ4-FRQJSQ
      VALU7=VALU7+(FACT7)/(DENOM)
      GO TO 207
C/////////////////////////////////
C///FERMI RESONANCE CORRECTIONS///
C/////////////////////////////////
C   WI+WI=WJ
  206 CONTINUE
      IF(DABS(FRIIJ).GT.PLIMIT) THEN
        IANH1=IANH1+1
      END IF
      DNUM=ONE/(FRQI+FRQI+FRQJ)
      VALU7=VALU7-(FIIJSQ*DNUM)/(2.0D+00)
  207 CONTINUE
      DO 104 KVIB=1,NVIB
      IF(NDEG(KVIB).NE.2) GO TO 104
      IF(IVIB.LE.KVIB) GO TO 104
      FRQK=FREQ(KVIB)
      FRIJK=F3Q(IVIB,JVIB,KVIB)
      FACT8=FRIJK*FRIJK*FRQI*FRQJ*FRQK
      FACTSQ=FRIJK*FRIJK
      WPPP=FRQI+FRQJ+FRQK
      WPPM=FRQI+FRQJ-FRQK
      WPMP=FRQI-FRQJ+FRQK
      WPMM=FRQI-FRQJ-FRQK
      IF(DABS(WPPM).LE.FLIM2) GO TO 208
      IF(DABS(WPMP).LE.FLIM2) GO TO 209
      IF(DABS(WPMM).LE.FLIM2) GO TO 210
      DIJK=WPPP*WPPM*WPMP*WPMM
      VALU8=VALU8+(FACT8)/(DIJK)
      GO TO 104
C/////////////////////////////////
C///FERMI RESONANCE CORRECTIONS///
C/////////////////////////////////
C   WI+WJ=WK
  208 CONTINUE
      DNUM8=ONE/WPPP-ONE/WPMP+ONE/WPMM
      IF(DABS(FRIJK).GT.PLIMIT) THEN
        IANH2=IANH2+1
      END IF
      VALU8=VALU8-(DNUM8*FACTSQ)/(8.0D+00)
      GO TO 104
C   WI+WK=WJ
  209 CONTINUE
      DNUM8=ONE/WPPP-ONE/WPPM+ONE/WPMM
      IF(DABS(FRIJK).GT.PLIMIT) THEN
        IANH2=IANH2+1
      END IF
      VALU8=VALU8-(DNUM8*FACTSQ)/(8.0D+00)
      GO TO 104
C   WJ+WK=WI
  210 CONTINUE
      DNUM8=ONE/WPPP-ONE/WPPM-ONE/WPMP
      IF(DABS(FRIJK).GT.PLIMIT) THEN
        IANH2=IANH2+1
      END IF
      VALU8=VALU8-(DNUM8*FACTSQ)/(8.0D+00)
  104 CONTINUE
  105 CONTINUE
  106 CONTINUE
C
      IF(IANH1.NE.0) THEN
        WRITE(6,2)
        WRITE(6,4) IANH1,FLIM1
      END IF
      IF(IANH2.NE.0) THEN
        WRITE(6,3)
        WRITE(6,4) IANH2,FLIM2
      END IF
C
C=============================================
C  For degenerate * degenerate modes
C=============================================
      IANH1=0
      IANH2=0
      VALU9=A00
      VALU10=A00
      DO 109 IVIB=1,NVIB
      IF(NDEG(IVIB).NE.2) GO TO 109
      IVA=NDAB(IVIB,1)
      IVB=NDAB(IVIB,2)
      FRQI=FREQ(IVIB)
      FRQISQ=FRQI*FRQI
      FRISQ8=FRQISQ*8.0D+00
      FRISQ4=FRQISQ*4.0D+00
      DO 108 JVIB=1,NVIB
      IF(NDEG(JVIB).NE.2) GO TO 108
      JVA=NDAB(JVIB,1)
      JVB=NDAB(JVIB,2)
      FRQJ=FREQ(JVIB)
      FRIIJ=F3Q(IVIB,IVIB,JVIB)
      FIIJSQ=FRIIJ*FRIIJ
      FRQJSQ=FRQJ*FRQJ
      FACT9=FRIIJ*FRIIJ*(FRISQ8+FRQJSQ)
      IF(DABS(FRQI+FRQI-FRQJ).LE.FLIM1) GO TO 211
      DENOM=(FRISQ4-FRQJSQ)*FRQJ
      VALU9=VALU9+(FACT9)/(DENOM)
      GO TO 212
C/////////////////////////////////
C///FERMI RESONANCE CORRECTIONS///
C/////////////////////////////////
C   WI+WI=WJ
  211 CONTINUE
      IF(DABS(FRIIJ).GT.PLIMIT) THEN
        IANH1=IANH1+1
      END IF
      DNUM=(2.0D+00)/(FRQJ)
     *    -(3.0D+00)/(4.0D+00*FRQI+2.0D+00*FRQJ)
      VALU9=VALU9-(FIIJSQ*DNUM)/(2.0D+00)
  212 CONTINUE
      DO 107 KVIB=1,NVIB
      IF(NDEG(KVIB).NE.2) GO TO 107
      KVA=NDAB(KVIB,1)
      KVB=NDAB(KVIB,2)
      IF(IVIB.LE.JVIB) GO TO 107
      IF(IVA.EQ.JVA .OR. IVA.EQ.JVB) GO TO 107
      IF(IVB.EQ.JVA .OR. IVB.EQ.JVB) GO TO 107
      IF(IVIB.LE.KVIB) GO TO 107
      IF(IVA.EQ.KVA .OR. IVA.EQ.KVB) GO TO 107
      IF(IVB.EQ.KVA .OR. IVB.EQ.KVB) GO TO 107
      IF(JVIB.LE.KVIB) GO TO 107
      IF(JVA.EQ.KVA .OR. JVA.EQ.KVB) GO TO 107
      IF(JVB.EQ.KVA .OR. JVB.EQ.KVB) GO TO 107
      FRQK=FREQ(KVIB)
      FRIJK=F3Q(IVIB,JVIB,KVIB)
      FACT10=FRIJK*FRIJK*FRQI*FRQJ*FRQK
      FACTSQ=FRIJK*FRIJK
      WPPP=FRQI+FRQJ+FRQK
      WPPM=FRQI+FRQJ-FRQK
      WPMP=FRQI-FRQJ+FRQK
      WPMM=FRQI-FRQJ-FRQK
      IF(DABS(WPPM).LE.FLIM2) GO TO 213
      IF(DABS(WPMP).LE.FLIM2) GO TO 214
      IF(DABS(WPMM).LE.FLIM2) GO TO 215
      DIJK=WPPP*WPPM*WPMP*WPMM
      VALU10=VALU10+(FACT10)/(DIJK)
      GO TO 107
C/////////////////////////////////
C///FERMI RESONANCE CORRECTIONS///
C/////////////////////////////////
C   WI+WJ=WK
  213 CONTINUE
      DNUM10=ONE/WPPP-ONE/WPMP+ONE/WPMM
      IF(DABS(FRIJK).GT.PLIMIT) THEN
        IANH2=IANH2+1
      END IF
      VALU10=VALU10-(DNUM10*FACTSQ)/(8.0D+00)
      GO TO 107
C   WI+WK=WJ
  214 CONTINUE
      DNUM10=ONE/WPPP-ONE/WPPM+ONE/WPMM
      IF(DABS(FRIJK).GT.PLIMIT) THEN
        IANH2=IANH2+1
      END IF
      VALU10=VALU10-(DNUM10*FACTSQ)/(8.0D+00)
      GO TO 107
C   WJ+WK=WI
  215 CONTINUE
      DNUM10=ONE/WPPP-ONE/WPPM-ONE/WPMP
      IF(DABS(FRIJK).GT.PLIMIT) THEN
        IANH2=IANH2+1
      END IF
      VALU10=VALU10-(DNUM10*FACTSQ)/(8.0D+00)
  107 CONTINUE
  108 CONTINUE
  109 CONTINUE
C
      IF(IANH1.NE.0) THEN
        WRITE(6,2)
        WRITE(6,4) IANH1,FLIM1
      END IF
      IF(IANH2.NE.0) THEN
        WRITE(6,3)
        WRITE(6,4) IANH2,FLIM2
      END IF
C
      VALU11=A00
      DO 113 IABC=1,3
      II=IAXIS(IABC)
      BE=ROTCM(IABC)
      TERM=ONE
      DO 112 IVIB=1,NVIB
      DO 111 JVIB=1,NVIB
      ZET=ZETA(IVIB,JVIB,II)
      ZETSQ=ZET*ZET
      TERM=TERM+ZETSQ
  111 CONTINUE
  112 CONTINUE
      VALU11=VALU11+TERM*BE
  113 CONTINUE
C
      VALU1 =+VALU1/(64.0D+00)
      VALU2 =-VALU2/(36.0D+00)
      VALU3 =+(VALU3*3.0D+00)/(64.0D+00)
      VALU4 =-VALU4/(4.0D+00)
      VALU5 =+VALU5/(96.0D+00)
      VALU6 =-VALU6/(9.0D+00)
      VALU7 =+VALU7/(32.0D+00)
      VALU8 =-VALU8/(4.0D+00)
      VALU9 =+VALU9/(32.0D+00)
      VALU10=-VALU10/(2.0D+00)
      VALU11=-VALU11/(4.0D+00)
C
      GCUBIC= VALU2+VALU3+VALU4
     *       +VALU6+VALU7+VALU8 
     *       +VALU9+VALU10
      GQUART= VALU1+VALU5
      GKINET= VALU11
      G00   = GCUBIC+GQUART+GKINET
      ANHG00=ANHTOT+G00
      WRITE(6,5)
      IF(IPRNT.GT.3)
     *WRITE(6,6) VALU1,VALU2,VALU3,VALU4,
     *           VALU5,VALU6,VALU7,VALU8,
     *           VALU9,VALU10,VALU11
      WRITE(6,7) GCUBIC,GQUART,GKINET,G00
      WRITE(6,8) HZPVE,AZPVE,ANHTOT,ANHG00
C
      RETURN
      END
