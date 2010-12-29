C     PROGRAM ANHASY                                                
C**********************************************************************
C***LAST UPDATED ON January 10, 2006 BY DR. YUKIO YAMAGUCHI         ***
C   Modification for a stand-alone version                          ***
C**********************************************************************
C***LAST UPDATED ON April 20, 2004 BY DR. YUKIO YAMAGUCHI           ***
C   Modification                                                    ***
C   for anharmonic zero-point vibrational energy (ZPVE)             ***
C   and vibrational G_0 energy                                      ***
C**********************************************************************
C***LAST UPDATED ON SEPTEMBER 07, 1987 BY DR. YUKIO YAMAGUCHI       ***  
C**********************************************************************
      SUBROUTINE ROVASY(ZETA,AXY,F3Q,EE,RVCNST)                         ANH00050
C   THE CALCULATION OF VIBRATION-ROTATION INTERACTION CONSTANTS         ANH00060
C   FOR AN ASYMMETRIC TOP MOLECULE                                      ANH00070
      IMPLICIT REAL*8 (A-H,O-Z)                                         ANH00080
      CHARACTER*1 RAXIS                                                 ANH00090
      CHARACTER*25 RTYPE                                                ANH00100
      DIMENSION ZETA(N3N,N3N,3),AXY(3,3,N3N)                            ANH00110
      DIMENSION F3Q(N3N,N3N,N3N),EE(N3N,N3N),RVCNST(NVIB,3)             ANH00120
      COMMON/VIB101/NATOM,N3N,NATRI,ILIN,NVIB                           ANH00130
      COMMON/VIB102/ITHREE,IFOUR,N3TOT,N4TOT                            ANH00140
      COMMON/VIB103/PARA,WAVE,CONST,CYCL,CONV                           ANH00150
      COMMON/VIB106/ROTAA(3),ROTGC(3),ROTCM(3),ROTMH(3)                 ANH00160
      COMMON/VIB108/IAXIS(3),NDEG(150),NDAB(150,5),IMAG(150)            ANH00170
      COMMON/VIB109/CLIMIT,FLIM1,FLIM2                                  ANH00180
      COMMON/VIB203/IOFF(150),IPRNT                                     ANH00190
      COMMON/VIB204/SQM(150),ROOT(150),FREQ(150)                        ANH00200
      COMMON/VIB205/IFREQ,NFRQ(150)                                     ANH00210
      COMMON/VIB206/RAXIS(3),RTYPE(6)                                   ANH00220
      COMMON/VIB209/SALP(3)                                             ANH00230
      DATA PH,AVN / 6.626176D+00 , 6.022045D+00 /                       ANH00240
      DATA PI,CL / 3.1415926536D+00 , 2.99792458D+00 /                  ANH00250
      DATA A00,HALF,ONE,TWO / 0.0D+00 , 0.5D+00 , 1.0D+00 , 2.0D+00 /   ANH00260
      DATA WLIMIT / 1.0D+00 /                                           ANH00270
      DATA ZLIMIT / 1.0D-05 /                                           ANH00280
    1 FORMAT(//,2X,' **************************************************'ANH00290
     1        /,2X,' ***THE VIBRATION-ROTATION INTERACTION CONSTANTS***'ANH00300
     2        /,2X,' **************************************************'ANH00310
     3        )                                                         ANH00320
    2 FORMAT(//,2X,' ZETA MATRIX, IABC = ',I5/)                         ANH00330
    3 FORMAT(//,2X,' F3(QI*QJ*QK) MATRIX (IN CM-1), IVIB = ',I5/)       ANH00340
    4 FORMAT(//,2X,' AXY MATRIX, IVIB = ',I5/)                          ANH00350
    5 FORMAT(//,2X,' QVIB',5X,' IVIB',10X,' A',21X,' B',21X,' C'/       ANH00360
     1          27X,' IN CM-1',15X,' IN CM-1',15X,' IN CM-1'/)          ANH00370
    6 FORMAT(2X,I5,5X,I5,3X,F20.10,3X,F20.10,3X,F20.10)                 ANH00380
    7 FORMAT(//,2X,' PRINCIPAL MOMENTS OF INERTIA WITH ZERO-POINT VIBRATANH00390
     1IONAL CORRECTION'/                                                ANH00400
     2 4X,5H AXIS,9X,11H IN AMU.A+2,9X,10H IN G.CM+2,10X,               ANH00410
     3 8H IN CM-1,15X,7H IN MHZ/38X,8H (*D+39)/)                        ANH00420
    8 FORMAT(8X,A1,3X,3F20.10,F20.5)                                    ANH00430
    9 FORMAT(//,2X,' ::::::::::::::::::::::::::::::::::::'/             ANH00440
     1          2X,' :::CORIOLIS RESONANCE IS OBSERVED:::'/             ANH00450
     2          2X,' ::::::::::::::::::::::::::::::::::::'/)            ANH00460
   10 FORMAT(2X,' NUMBER OF OCCURENCES IS ',I5/                         ANH00470
     1       2X,' RESONANCE LIMIT IS SET TO ',F10.3,' CM-1'/)           ANH00480
C                                                                       ANH00490
      WRITE(6,1)                                                        ANH00500
      RCONST=1.0D-01*PI*DSQRT(CL/(PH*AVN))                              ANH00510
      IAA=IAXIS(1)                                                      ANH00520
      IBB=IAXIS(2)                                                      ANH00530
      ICC=IAXIS(3)                                                      ANH00540
      NCORIO=0                                                          ANH00550
C                                                                       ANH00560
      IF(IPRNT.LE.2) GO TO 201                                          ANH00570
      DO 101 IABC=1,3                                                   ANH00580
      WRITE(6,2) IABC                                                   ANH00590
      CALL MATOUT(ZETA(1,1,IABC),N3N,N3N,N3N,N3N,6)                     ANH00600
  101 CONTINUE                                                          ANH00610
      DO 103 III=1,NVIB                                                 ANH00620
      IVIB=NFRQ(III)                                                    ANH00630
      WRITE(6,3) IVIB                                                   ANH00640
      DO 102 I=1,NVIB                                                   ANH00650
      II=NFRQ(I)                                                        ANH00660
      DO 102 J=1,NVIB                                                   ANH00670
      JJ=NFRQ(J)                                                        ANH00680
      EE(I,J)=F3Q(II,JJ,IVIB)                                           ANH00690
  102 CONTINUE                                                          ANH00700
      CALL MATOUT(EE,N3N,N3N,NVIB,NVIB,6)                               ANH00710
  103 CONTINUE                                                          ANH00720
      DO 104 IVIB=1,NVIB                                                ANH00730
      WRITE(6,4) IVIB                                                   ANH00740
      CALL MATOUT(AXY(1,1,IVIB),3,3,3,3,6)                              ANH00750
  104 CONTINUE                                                          ANH00760
C                                                                       ANH00770
  201 CONTINUE                                                          ANH00780
      DO 110 IVIB=1,NVIB                                                ANH00790
      FRQI=FREQ(IVIB)                                                   ANH00800
      FRQI2=FRQI*FRQI                                                   ANH00810
      VALUA1=A00                                                        ANH00820
      VALUB1=A00                                                        ANH00830
      VALUC1=A00                                                        ANH00840
      DO 106 I=1,3                                                      ANH00850
      VALUA1=VALUA1+(AXY(IAA,IAXIS(I),IVIB)**2/ROTAA(I))*0.75D+00       ANH00860
      VALUB1=VALUB1+(AXY(IBB,IAXIS(I),IVIB)**2/ROTAA(I))*0.75D+00       ANH00870
      VALUC1=VALUC1+(AXY(ICC,IAXIS(I),IVIB)**2/ROTAA(I))*0.75D+00       ANH00880
  106 CONTINUE                                                          ANH00890
      VALUA2=A00                                                        ANH00900
      VALUB2=A00                                                        ANH00910
      VALUC2=A00                                                        ANH00920
      VALUA3=A00                                                        ANH00930
      VALUB3=A00                                                        ANH00940
      VALUC3=A00                                                        ANH00950
      DO 107 JVIB=1,NVIB                                                ANH00960
      FRQJ=FREQ(JVIB)                                                   ANH00970
      FRQJ2=FRQJ*FRQJ                                                   ANH00980
      IF(IVIB.EQ.JVIB) GO TO 205                                        ANH00990
      IF(DABS(FRQI-FRQJ).LE.CLIMIT) GO TO 203                           ANH01000
      FACT2=(FRQI2*3.0D+00+FRQJ2)/(FRQI2-FRQJ2)                         ANH01010
      GO TO 204                                                         ANH01020
C////////////////////////////////////                                   ANH01030
C///CORIOLIS RESONANCE CORRECTIONS///                                   ANH01040
C////////////////////////////////////                                   ANH01050
  203 CONTINUE                                                          ANH01060
      FACT2=-(0.5D+00*(FRQI-FRQJ)**2)/(FRQJ*(FRQI+FRQJ))                ANH01070
      IF(DABS(ZETA(IVIB,JVIB,IAA)).GT.ZLIMIT) THEN                      ANH01080
        NCORIO=NCORIO+1                                                 ANH01090
      END IF                                                            ANH01100
      IF(DABS(ZETA(IVIB,JVIB,IBB)).GT.ZLIMIT) THEN                      ANH01110
        NCORIO=NCORIO+1                                                 ANH01120
      END IF                                                            ANH01130
      IF(DABS(ZETA(IVIB,JVIB,ICC)).GT.ZLIMIT) THEN                      ANH01140
        NCORIO=NCORIO+1                                                 ANH01150
      END IF                                                            ANH01160
C///////////////////////////////                                        ANH01170
C///END OF CORIOLIS RESONANCE///                                        ANH01180
C///////////////////////////////                                        ANH01190
  204 CONTINUE                                                          ANH01200
      VALUA2=VALUA2+(ZETA(IVIB,JVIB,IAA)**2)*FACT2                      ANH01210
      VALUB2=VALUB2+(ZETA(IVIB,JVIB,IBB)**2)*FACT2                      ANH01220
      VALUC2=VALUC2+(ZETA(IVIB,JVIB,ICC)**2)*FACT2                      ANH01230
  205 CONTINUE                                                          ANH01240
      FACT3=(FRQI*F3Q(IVIB,IVIB,JVIB)*RCONST)/(FRQJ**(1.5D+00))         ANH01250
      VALUA3=VALUA3+AXY(IAA,IAA,JVIB)*FACT3                             ANH01260
      VALUB3=VALUB3+AXY(IBB,IBB,JVIB)*FACT3                             ANH01270
      VALUC3=VALUC3+AXY(ICC,ICC,JVIB)*FACT3                             ANH01280
  107 CONTINUE                                                          ANH01290
      ROTA=-(VALUA1+VALUA2+VALUA3)*(ROTCM(1)*ROTCM(1)*TWO)/FRQI         ANH01300
      ROTB=-(VALUB1+VALUB2+VALUB3)*(ROTCM(2)*ROTCM(2)*TWO)/FRQI         ANH01310
      ROTC=-(VALUC1+VALUC2+VALUC3)*(ROTCM(3)*ROTCM(3)*TWO)/FRQI         ANH01320
      RVCNST(IVIB,1)=ROTA                                               ANH01330
      RVCNST(IVIB,2)=ROTB                                               ANH01340
      RVCNST(IVIB,3)=ROTC                                               ANH01350
  110 CONTINUE                                                          ANH01360
C                                                                       ANH01370
C   CALCULATE THE ROTATIONAL CONSTANTS                                  ANH01380
C   WITH ZERO-POINT VIBRATIONAL CORRECTION                              ANH01390
      CALL ZERO(SALP,3)                                                 ANH01400
      WRITE(6,5)                                                        ANH01410
      DO 115 III=1,NVIB                                                 ANH01420
      IVIB=NFRQ(III)                                                    ANH01430
      ROTA=RVCNST(IVIB,1)                                               ANH01440
      ROTB=RVCNST(IVIB,2)                                               ANH01450
      ROTC=RVCNST(IVIB,3)                                               ANH01460
      SALP(1)=SALP(1)+ROTA                                              ANH01470
      SALP(2)=SALP(2)+ROTB                                              ANH01480
      SALP(3)=SALP(3)+ROTC                                              ANH01490
      WRITE(6,6) III,IVIB,ROTA,ROTB,ROTC                                ANH01500
  115 CONTINUE                                                          ANH01510
      WRITE(6,7)                                                        ANH01520
      DO 116 I=1,3                                                      ANH01530
      PB=ROTCM(I)-SALP(I)*HALF                                          ANH01540
      PT=CONST/PB                                                       ANH01550
      PA=PT*PARA                                                        ANH01560
      PC=CYCL/PT                                                        ANH01570
      WRITE(6,8) RAXIS(I),PT,PA,PB,PC                                   ANH01580
  116 CONTINUE                                                          ANH01590
C                                                                       ANH01600
      IF(NCORIO.NE.0) THEN                                              ANH01610
        WRITE(6,9)                                                      ANH01620
        WRITE(6,10) NCORIO,CLIMIT                                       ANH01630
      END IF                                                            ANH01640
C                                                                       ANH01650
      RETURN                                                            ANH01660
      END                                                               ANH01670
      SUBROUTINE CENTAS(F3Q,EE)                                         ANH01680
C   THE CALCULATION OF CENTRIFUGAL DISTORTION CONSTANTS                 ANH01690
C   FOR AN ASYMMETRIC TOP MOLECULE                                      ANH01700
      IMPLICIT REAL*8 (A-H,O-Z)                                         ANH01710
      CHARACTER*1 RAXIS                                                 ANH01720
      CHARACTER*25 RTYPE                                                ANH01730
      DIMENSION RCENT(3),RTOTC(3),RTOTM(3)                              ANH01740
      DIMENSION F3Q(N3N,N3N,N3N),EE(N3N,N3N)                            ANH01750
      COMMON/VIB101/NATOM,N3N,NATRI,ILIN,NVIB                           ANH01760
      COMMON/VIB102/ITHREE,IFOUR,N3TOT,N4TOT                            ANH01770
      COMMON/VIB103/PARA,WAVE,CONST,CYCL,CONV                           ANH01780
      COMMON/VIB106/ROTAA(3),ROTGC(3),ROTCM(3),ROTMH(3)                 ANH01790
      COMMON/VIB108/IAXIS(3),NDEG(150),NDAB(150,5),IMAG(150)            ANH01800
      COMMON/VIB110/ISIGMA                                              ANH01810
      COMMON/VIB203/IOFF(150),IPRNT                                     ANH01820
      COMMON/VIB205/IFREQ,NFRQ(150)                                     ANH01830
      COMMON/VIB206/RAXIS(3),RTYPE(6)                                   ANH01840
      COMMON/VIB207/TABCD(3,3,3,3),PABC(3,3,3)                          ANH01850
      COMMON/VIB208/TAABB(3,3),TABAB(3,3)                               ANH01860
      COMMON/VIB209/SALP(3)                                             ANH01870
      DATA CL / 2.99792458D+00 /                                        ANH01880
      DATA A00,HALF,ONE,TWO / 0.0D+00 , 0.5D+00 , 1.0D+00 , 2.0D+00 /   ANH01890
    1 FORMAT(//,2X,' ******************************************'/       ANH01900
     1          2X,' ***THE CENTRIFUGAL DISTORTION CONSTANTS***'/       ANH01910
     2          2X,' ******************************************')       ANH01920
    2 FORMAT(//,2X,' F3(QI*QJ*QK) MATRIX (IN CM-1), IVIB = ',I5/)       ANH01930
    3 FORMAT(//,2X,' TAABB MATRIX FROM BXY MATRICES (IN CM-1)'/)        ANH01940
    4 FORMAT(//,2X,' TABAB MATRIX FROM BXY MATRICES (IN CM-1)'/)        ANH01950
    5 FORMAT(//,2X,' FOR AN ASYMMETRIC TOP MOLECULE'/                   ANH01960
     1          2X,' DJ       = ',F20.10,' D-06 CM-1',                  ANH01970
     2          8X,' DJ       = ',F20.10,' MHZ'/                        ANH01980
     3          2X,' DK       = ',F20.10,' D-06 CM-1',                  ANH01990
     4          8X,' DK       = ',F20.10,' MHZ'/                        ANH02000
     5          2X,' DJK      = ',F20.10,' D-06 CM-1',                  ANH02010
     6          8X,' DJK      = ',F20.10,' MHZ'/                        ANH02020
     7          2X,' R5       = ',F20.10,' D-06 CM-1',                  ANH02030
     8          8X,' R5       = ',F20.10,' MHZ'/                        ANH02040
     9          2X,' R6       = ',F20.10,' D-06 CM-1',                  ANH02050
     A          8X,' R6       = ',F20.10,' MHZ'/                        ANH02060
     B          2X,' DELJ     = ',F20.10,' D-06 CM-1',                  ANH02070
     C          8X,' DELJ     = ',F20.10,' MHZ'/)                       ANH02080
    6 FORMAT(//,2X,' PRINCIPAL MOMENTS OF INERTIA AT EQUILIBRIUM GEOMETRANH02090
     1Y'/                                                               ANH02100
     2 4X,5H AXIS,9X,11H IN AMU.A+2,9X,10H IN G.CM+2,10X,               ANH02110
     3 8H IN CM-1,15X,7H IN MHZ/38X,8H (*D+39)/)                        ANH02120
    7 FORMAT(8X,A1,3X,3F20.10,F20.5)                                    ANH02130
    8 FORMAT(//,2X,' PRINCIPAL MOMENTS OF INERTIA WITH ZERO-POINT VIBRATANH02140
     1IONAL CORRECTION'/                                                ANH02150
     2 4X,5H AXIS,9X,11H IN AMU.A+2,9X,10H IN G.CM+2,10X,               ANH02160
     3 8H IN CM-1,15X,7H IN MHZ/38X,8H (*D+39)/)                        ANH02170
    9 FORMAT(//,2X,' PRINCIPAL MOMENTS OF INERTIA WITH ZERO-POINT VIBRATANH02180
     1IONAL AND CENTRIFUGAL DISTORTION CORRECTIONS'/                    ANH02190
     2 4X,5H AXIS,9X,11H IN AMU.A+2,9X,10H IN G.CM+2,10X,               ANH02200
     3 8H IN CM-1,15X,7H IN MHZ/38X,8H (*D+39)/)                        ANH02210
   10 FORMAT(/,2X,' SIGMA (0'') = ',F20.10/                             ANH02220
     1         2X,' SIGMA (E)  = ',F20.10/)                             ANH02230
   11 FORMAT(/,2X,' SIGMA (0'') WILL BE USED TO CALCULATE CENTRIFUGAL DIANH02240
     1STORTION CONSTANTS'/)                                             ANH02250
   12 FORMAT(/,2X,' SIGMA (E) WILL BE USED TO CALCULATE CENTRIFUGAL DISTANH02260
     1ORTION CONSTANTS'/)                                               ANH02270
   13 FORMAT(/,2X,' A (TILDE)  = ',F20.10,' CM-1',                      ANH02280
     1        13X,' A (TILDE)  = ',F20.10,' MHZ'/                       ANH02290
     2         2X,' B (TILDE)  = ',F20.10,' CM-1',                      ANH02300
     3        13X,' B (TILDE)  = ',F20.10,' MHZ'/                       ANH02310
     4         2X,' C (TILDE)  = ',F20.10,' CM-1',                      ANH02320
     5        13X,' C (TILDE)  = ',F20.10,' MHZ'/)                      ANH02330
   14 FORMAT(//,2X,' FOR AN ASYMMETRIC TOP MOLECULE'/                   ANH02340
     1          2X,' DELTA J    = ',F20.10,' D-06 CM-1',                ANH02350
     2          8X,' DELTA J    = ',F20.10,' MHZ'/                      ANH02360
     3          2X,' DELTA K    = ',F20.10,' D-06 CM-1',                ANH02370
     4          8X,' DELTA K    = ',F20.10,' MHZ'/                      ANH02380
     5          2X,' DELTA JK   = ',F20.10,' D-06 CM-1',                ANH02390
     6          8X,' DELTA JK   = ',F20.10,' MHZ'/                      ANH02400
     7          2X,' DELK       = ',F20.10,' D-06 CM-1',                ANH02410
     8          8X,' DELK       = ',F20.10,' MHZ'/                      ANH02420
     9          2X,' DELJ       = ',F20.10,' D-06 CM-1',                ANH02430
     A          8X,' DELJ       = ',F20.10,' MHZ'/)                     ANH02440
   15 FORMAT(//,2X,' FOR AN ASYMMETRIC TOP MOLECULE'/                   ANH02450
     1          2X,' DDJ        = ',F20.10,' D-06 CM-1',                ANH02460
     2          8X,' DDJ        = ',F20.10,' MHZ'/                      ANH02470
     3          2X,' DDK        = ',F20.10,' D-06 CM-1',                ANH02480
     4          8X,' DDK        = ',F20.10,' MHZ'/                      ANH02490
     5          2X,' DDJK       = ',F20.10,' D-06 CM-1',                ANH02500
     6          8X,' DDJK       = ',F20.10,' MHZ'/                      ANH02510
     7          2X,' DDWJ       = ',F20.10,' D-06'/                     ANH02520
     8          2X,' DDWK       = ',F20.10,' D-06'/)                    ANH02530
   16 FORMAT(//,2X,' FOR AN ASYMMETRIC TOP MOLECULE (S111=0)'/          ANH02540
     1          2X,' HJ         = ',F20.10,' D-06 CM-1',                ANH02550
     2          8X,' HJ         = ',F20.10,' MHZ'/                      ANH02560
     3          2X,' HJK        = ',F20.10,' D-06 CM-1',                ANH02570
     4          8X,' HJK        = ',F20.10,' MHZ'/                      ANH02580
     5          2X,' HKJ        = ',F20.10,' D-06 CM-1',                ANH02590
     6          8X,' HKJ        = ',F20.10,' MHZ'/                      ANH02600
     7          2X,' HK         = ',F20.10,' D-06 CM-1',                ANH02610
     8          8X,' HK         = ',F20.10,' MHZ'/                      ANH02620
     9          2X,' HHJ        = ',F20.10,' D-06 CM-1',                ANH02630
     A          8X,' HHJ        = ',F20.10,' MHZ'/                      ANH02640
     B          2X,' HHJK       = ',F20.10,' D-06 CM-1',                ANH02650
     C          8X,' HHJK       = ',F20.10,' MHZ'/                      ANH02660
     D          2X,' HHK        = ',F20.10,' D-06 CM-1',                ANH02670
     E          8X,' HHK        = ',F20.10,' MHZ'/)                     ANH02680
   17 FORMAT(//,2X,' FOR AN ASYMMETRIC TOP MOLECULE (S111=-4*R6/(B-C))'/ANH02690
     1          2X,' HJ         = ',F20.10,' D-06 CM-1',                ANH02700
     2          8X,' HJ         = ',F20.10,' MHZ'/                      ANH02710
     3          2X,' HJK        = ',F20.10,' D-06 CM-1',                ANH02720
     4          8X,' HJK        = ',F20.10,' MHZ'/                      ANH02730
     5          2X,' HKJ        = ',F20.10,' D-06 CM-1',                ANH02740
     6          8X,' HKJ        = ',F20.10,' MHZ'/                      ANH02750
     7          2X,' HK         = ',F20.10,' D-06 CM-1',                ANH02760
     8          8X,' HK         = ',F20.10,' MHZ'/                      ANH02770
     9          2X,' HHJ        = ',F20.10,' D-06 CM-1',                ANH02780
     A          8X,' HHJ        = ',F20.10,' MHZ'/                      ANH02790
     B          2X,' HHJK       = ',F20.10,' D-06 CM-1',                ANH02800
     C          8X,' HHJK       = ',F20.10,' MHZ'/                      ANH02810
     D          2X,' HHK        = ',F20.10,' D-06 CM-1',                ANH02820
     E          8X,' HHK        = ',F20.10,' MHZ'/)                     ANH02830
C                                                                       ANH02840
      WRITE(6,1)                                                        ANH02850
C   ELEMENTS OF F3Q ARE IN CM-1                                         ANH02860
      IF(IPRNT.LE.2) GO TO 201                                          ANH02870
      DO 102 III=1,NVIB                                                 ANH02880
      IVIB=NFRQ(III)                                                    ANH02890
      WRITE(6,2) IVIB                                                   ANH02900
      DO 101 I=1,NVIB                                                   ANH02910
      II=NFRQ(I)                                                        ANH02920
      DO 101 J=1,NVIB                                                   ANH02930
      JJ=NFRQ(J)                                                        ANH02940
      EE(I,J)=F3Q(II,JJ,IVIB)                                           ANH02950
  101 CONTINUE                                                          ANH02960
      CALL MATOUT(EE,N3N,N3N,NVIB,NVIB,6)                               ANH02970
  102 CONTINUE                                                          ANH02980
C                                                                       ANH02990
  201 CONTINUE                                                          ANH03000
      DO 103 I=1,3                                                      ANH03010
      DO 103 J=1,3                                                      ANH03020
      TAABB(I,J)=TABCD(I,I,J,J)                                         ANH03030
      TABAB(I,J)=TABCD(I,J,I,J)                                         ANH03040
  103 CONTINUE                                                          ANH03050
      IF(IPRNT.LE.2) GO TO 202                                          ANH03060
      WRITE(6,3)                                                        ANH03070
      CALL MATOUT(TAABB,3,3,3,3,6)                                      ANH03080
      WRITE(6,4)                                                        ANH03090
      CALL MATOUT(TABAB,3,3,3,3,6)                                      ANH03100
C                                                                       ANH03110
C   FOR AN ASYMMETRIC TOP MOLECULE                                      ANH03120
  202 CONTINUE                                                          ANH03130
      IAA=IAXIS(1)                                                      ANH03140
      IBB=IAXIS(2)                                                      ANH03150
      ICC=IAXIS(3)                                                      ANH03160
C                                                                       ANH03170
      DJ=-(TAABB(IBB,IBB)*3.0D+00+TAABB(ICC,ICC)*3.0D+00                ANH03180
     1    +TAABB(IBB,ICC)*2.0D+00+TABAB(IBB,ICC)*4.0D+00)/32.0D+00      ANH03190
      DK=DJ-(TAABB(IAA,IAA)-TAABB(IAA,IBB)-TAABB(ICC,IAA)               ANH03200
     1      -TABAB(IAA,IBB)*2.0D+00-TABAB(ICC,IAA)*2.0D+00)/4.0D+00     ANH03210
      DJK=-DJ-DK-TAABB(IAA,IAA)/4.0D+00                                 ANH03220
      R5=-(TAABB(IBB,IBB)-TAABB(ICC,ICC)                                ANH03230
     1    -TAABB(IAA,IBB)*2.0D+00-TABAB(IAA,IBB)*4.0D+00                ANH03240
     2    +TAABB(ICC,IAA)*2.0D+00+TABAB(ICC,IAA)*4.0D+00)/32.0D+00      ANH03250
      R6=(TAABB(IBB,IBB)+TAABB(ICC,ICC)                                 ANH03260
     1   -TAABB(IBB,ICC)*2.0D+00-TABAB(IBB,ICC)*4.0D+00)/64.0D+00       ANH03270
      DELJ=-(TAABB(IBB,IBB)-TAABB(ICC,ICC))/16.0D+00                    ANH03280
      DJCM=DJ*1.0D+06                                                   ANH03290
      DKCM=DK*1.0D+06                                                   ANH03300
      DJKCM=DJK*1.0D+06                                                 ANH03310
      R5CM=R5*1.0D+06                                                   ANH03320
      R6CM=R6*1.0D+06                                                   ANH03330
      DELJCM=DELJ*1.0D+06                                               ANH03340
      DJMH=DJ*CL*1.0D+04                                                ANH03350
      DKMH=DK*CL*1.0D+04                                                ANH03360
      DJKMH=DJK*CL*1.0D+04                                              ANH03370
      R5MH=R5*CL*1.0D+04                                                ANH03380
      R6MH=R6*CL*1.0D+04                                                ANH03390
      DELJMH=DELJ*CL*1.0D+04                                            ANH03400
      WRITE(6,5) DJCM,DJMH,DKCM,DKMH,DJKCM,DJKMH,R5CM,R5MH,R6CM,R6MH,   ANH03410
     1           DELJCM,DELJMH                                          ANH03420
C                                                                       ANH03430
C   ROTATIONAL CONSTANTS WITH CENTRIFUGAL DISTORTION CORRECTIONS        ANH03440
      WRITE(6,6)                                                        ANH03450
      DO 105 I=1,3                                                      ANH03460
      PB=ROTCM(I)                                                       ANH03470
      PT=CONST/PB                                                       ANH03480
      PA=PT*PARA                                                        ANH03490
      PC=CYCL/PT                                                        ANH03500
      WRITE(6,7) RAXIS(I),PT,PA,PB,PC                                   ANH03510
  105 CONTINUE                                                          ANH03520
      WRITE(6,8)                                                        ANH03530
      DO 106 I=1,3                                                      ANH03540
      PB=ROTCM(I)-SALP(I)*HALF                                          ANH03550
      PT=CONST/PB                                                       ANH03560
      PA=PT*PARA                                                        ANH03570
      PC=CYCL/PT                                                        ANH03580
      WRITE(6,7) RAXIS(I),PT,PA,PB,PC                                   ANH03590
  106 CONTINUE                                                          ANH03600
      WRITE(6,9)                                                        ANH03610
      RCENT(1)=(TABAB(IBB,ICC)*3.0D+00-TABAB(ICC,IAA)*TWO               ANH03620
     1                                -TABAB(IAA,IBB)*TWO)/4.0D+00      ANH03630
      RCENT(2)=(TABAB(ICC,IAA)*3.0D+00-TABAB(IAA,IBB)*TWO               ANH03640
     1                                -TABAB(IBB,ICC)*TWO)/4.0D+00      ANH03650
      RCENT(3)=(TABAB(IAA,IBB)*3.0D+00-TABAB(IBB,ICC)*TWO               ANH03660
     1                                -TABAB(ICC,IAA)*TWO)/4.0D+00      ANH03670
C                                                                       ANH03680
      DO 107 I=1,3                                                      ANH03690
      PB=ROTCM(I)-SALP(I)*HALF+RCENT(I)                                 ANH03700
      RTOTC(I)=PB                                                       ANH03710
      PT=CONST/PB                                                       ANH03720
      PA=PT*PARA                                                        ANH03730
      PC=CYCL/PT                                                        ANH03740
      RTOTM(I)=PC                                                       ANH03750
      WRITE(6,7) RAXIS(I),PT,PA,PB,PC                                   ANH03760
  107 CONTINUE                                                          ANH03770
C                                                                       ANH03780
C   ASYMMETRY PARAMETER                                                 ANH03790
C####################################################                   ANH03800
C#  NOTE :                                         ##                   ANH03810
C#       IF(ISIGMA.EQ.0) THEN SIGMA IS CALCULATED  ##                   ANH03820
C#       USING A(0') , B(0') AND C(0') CONSTANTS   ##                   ANH03830
C#                                                 ##                   ANH03840
C#       IF(ISIGMA.NE.0) THEN SIGMA IS CALCULATED  ##                   ANH03850
C#       USING A(E) , B(E) AND C(E) CONSTANTS      ##                   ANH03860
C#                                                 ##                   ANH03870
C####################################################                   ANH03880
      SIGMA=(RTOTC(IAA)*TWO-RTOTC(IBB)-RTOTC(ICC))                      ANH03890
     1     /(RTOTC(IBB)-RTOTC(ICC))                                     ANH03900
      SIGMB=(ROTCM(IAA)*TWO-ROTCM(IBB)-ROTCM(ICC))                      ANH03910
     1     /(ROTCM(IBB)-ROTCM(ICC))                                     ANH03920
      WRITE(6,10) SIGMA,SIGMB                                           ANH03930
      IF(ISIGMA.EQ.0) THEN                                              ANH03940
        WRITE(6,11)                                                     ANH03950
      ELSE                                                              ANH03960
        SIGMA=SIGMB                                                     ANH03970
        WRITE(6,12)                                                     ANH03980
      END IF                                                            ANH03990
C                                                                       ANH04000
C   DEFINE ATILDE , BTILDE AND CTILDE                                   ANH04010
      ATILDE=RTOTC(IAA)+16.0D+00*R6                                     ANH04020
      BTILDE=RTOTC(IBB)-8.0D+00*R6*(SIGMA+1.0D+00)                      ANH04030
      CTILDE=RTOTC(ICC)+8.0D+00*R6*(SIGMA-1.0D+00)                      ANH04040
C*************************                                              ANH04050
C*  NOTE :  ISIGMA.NE.0 **                                              ANH04060
C*************************                                              ANH04070
      IF(ISIGMA.NE.0) THEN                                              ANH04080
        ATILDE=ROTCM(IAA)+16.0D+00*R6                                   ANH04090
        BTILDE=ROTCM(IBB)-8.0D+00*R6*(SIGMA+1.0D+00)                    ANH04100
        CTILDE=ROTCM(ICC)+8.0D+00*R6*(SIGMA-1.0D+00)                    ANH04110
      END IF                                                            ANH04120
      ATILDM=ATILDE*CL*1.0D+04                                          ANH04130
      BTILDM=BTILDE*CL*1.0D+04                                          ANH04140
      CTILDM=CTILDE*CL*1.0D+04                                          ANH04150
      WRITE(6,13) ATILDE,ATILDM,BTILDE,BTILDM,CTILDE,CTILDM             ANH04160
C                                                                       ANH04170
C   DEFINE CETRIFUGAL CONSTANTS (DELTA ETC.)                            ANH04180
      DLTJCM=DJCM-R6CM*2.0D+00                                          ANH04190
      DLTKCM=DKCM-R6CM*10.0D+00                                         ANH04200
      DLTJKC=DJKCM+R6CM*12.0D+00                                        ANH04210
      DELK=-R5*2.0D+00-4.0D+00*R6*SIGMA                                 ANH04220
      DELKCM=DELK*1.0D+06                                               ANH04230
      DLTJMH=DJMH-R6MH*2.0D+00                                          ANH04240
      DLTKMH=DKMH-R6MH*10.0D+00                                         ANH04250
      DLTJKM=DJKMH+R6MH*12.0D+00                                        ANH04260
      DELKMH=DELK*CL*1.0D+04                                            ANH04270
      WRITE(6,14) DLTJCM,DLTJMH,DLTKCM,DLTKMH,DLTJKC,DLTJKM,            ANH04280
     1            DELKCM,DELKMH,DELJCM,DELJMH                           ANH04290
C                                                                       ANH04300
C   CALCULATE D(J) , D(JK) , D(K) , D(WJ) , AND D(WK)                   ANH04310
      TILMM=BTILDE-CTILDE                                               ANH04320
      TILBC=(BTILDE+CTILDE)/TILMM                                       ANH04330
      TILAA=(2.0D+00*ATILDE-BTILDE-CTILDE)/TILMM                        ANH04340
      DDJCM=DLTJCM-DELJCM*2.0D+00*TILBC                                 ANH04350
      DDKCM=DLTKCM-DELKCM*2.0D+00*TILAA                                 ANH04360
      DDJKCM=DLTJKC-DELKCM*2.0D+00*TILBC-DELJCM*2.0D+00*TILAA           ANH04370
      DDWJ=(DELJ*4.0D+00*1.0D+06)/TILMM                                 ANH04380
      DDWK=(DELK*4.0D+00*1.0D+06)/TILMM                                 ANH04390
      DDJMH=DLTJMH-DELJMH*2.0D+00*TILBC                                 ANH04400
      DDKMH=DLTKMH-DELKMH*2.0D+00*TILAA                                 ANH04410
      DDJKMH=DLTJKM-DELKMH*2.0D+00*TILBC-DELJMH*2.0D+00*TILAA           ANH04420
      WRITE(6,15) DDJCM,DDJMH,DDKCM,DDKMH,DDJKCM,DDJKMH,                ANH04430
     1            DDWJ,DDWK                                             ANH04440
C                                                                       ANH04450
C   CALCULATE THE SEXTIC DISTIRTION CONSTANTS                           ANH04460
      P600=(PABC(IBB,IBB,IBB)+PABC(ICC,ICC,ICC))*(5.0D+00/16.0D+00)     ANH04470
     1    +(PABC(IBB,IBB,ICC)+PABC(ICC,ICC,IBB))/8.0D+00                ANH04480
      P420=(PABC(IAA,IBB,IBB)+PABC(IAA,ICC,ICC))*(3.0D+00/4.0D+00)      ANH04490
     1    +PABC(IAA,IBB,ICC)/4.0D+00-P600*3.0D+00                       ANH04500
      P240=PABC(IAA,IAA,IBB)+PABC(IAA,IAA,ICC)-P420*TWO-P600*3.0D+00    ANH04510
      P060=PABC(IAA,IAA,IAA)-P240-P420-P600                             ANH04520
      P402=(PABC(IBB,IBB,IBB)-PABC(ICC,ICC,ICC))*(15.0D+00/64.0D+00)    ANH04530
     1    +(PABC(IBB,IBB,ICC)-PABC(ICC,ICC,IBB))/32.0D+00               ANH04540
      P222=(PABC(IAA,IBB,IBB)-PABC(IAA,ICC,ICC))*HALF-P402*TWO          ANH04550
      P042=(PABC(IAA,IAA,IBB)-PABC(IAA,IAA,ICC))*HALF-P222-P402         ANH04560
      P204=(PABC(IBB,IBB,IBB)+PABC(ICC,ICC,ICC))*(3.0D+00/32.0D+00)     ANH04570
     1    -(PABC(IBB,IBB,ICC)+PABC(ICC,ICC,IBB))/16.0D+00               ANH04580
      P024=(PABC(IAA,IBB,IBB)+PABC(IAA,ICC,ICC)-PABC(IAA,IBB,ICC))      ANH04590
     1     /8.0D+00-P204                                                ANH04600
      P006=(PABC(IBB,IBB,IBB)-PABC(ICC,ICC,ICC))/64.0D+00               ANH04610
     1    -(PABC(IBB,IBB,ICC)-PABC(ICC,ICC,IBB))/32.0D+00               ANH04620
C                                                                       ANH04630
C   DEFINE DIFFERENCES IN ROTATIONAL CONSTANS                           ANH04640
      D2ABC=RTOTC(IAA)*TWO-RTOTC(IBB)-RTOTC(ICC)                        ANH04650
      DBC=RTOTC(IBB)-RTOTC(ICC)                                         ANH04660
C*************************                                              ANH04670
C*  NOTE :  ISIGMA.NE.0 **                                              ANH04680
C*************************                                              ANH04690
      IF(ISIGMA.NE.0) THEN                                              ANH04700
        D2ABC=ROTCM(IAA)*TWO-ROTCM(IBB)-ROTCM(ICC)                      ANH04710
        DBC=ROTCM(IBB)-ROTCM(ICC)                                       ANH04720
      END IF                                                            ANH04730
C                                                                       ANH04740
C   S111=0                                                              ANH04750
      SIGMA2=SIGMA*SIGMA                                                ANH04760
      S111=A00                                                          ANH04770
      S111SQ=S111*S111                                                  ANH04780
      HJ=P600+P204*TWO                                                  ANH04790
      HJK=P420-P204*12.0D+00+P024*TWO+P006*SIGMA*16.0D+00               ANH04800
     1   -(R5-SIGMA*R6*TWO)*S111*16.0D+00                               ANH04810
     2   +D2ABC*S111SQ*8.0D+00                                          ANH04820
      HKJ=P240+P420*(10.0D+00/3.0D+00)-P204*30.0D+00                    ANH04830
     1        -HJK*(10.0D+00/3.0D+00)                                   ANH04840
      HK=P060-P420*(7.0D+00/3.0D+00)+P204*28.0D+00+HJK*(7.0D+00/3.0D+00)ANH04850
      HHJ=P402+P006                                                     ANH04860
      HHJK=P222-P006*10.0D+00+P204*SIGMA*4.0D+00                        ANH04870
     1    +(DJK-SIGMA*DELJ*TWO-R6*4.0D+00)*S111*TWO                     ANH04880
     2    -DBC*S111SQ*4.0D+00                                           ANH04890
      HHK=P042+P024*(SIGMA*4.0D+00/3.0D+00)                             ANH04900
     1   +(9.0D+00+SIGMA2*32.0D+00/3.0D+00)*P006                        ANH04910
     2   +(DK-(SIGMA*R5*TWO)/3.0D+00+TWO*(ONE+8.0D+00*SIGMA2/3.0D+00)*R6ANH04920
     3    )*S111*4.0D+00+(6.0D+00+10.0D+00*SIGMA2/3.0D+00)              ANH04930
     4   *DBC*S111SQ                                                    ANH04940
C                                                                       ANH04950
      HJCM=HJ*1.0D+06                                                   ANH04960
      HJKCM=HJK*1.0D+06                                                 ANH04970
      HKJCM=HKJ*1.0D+06                                                 ANH04980
      HKCM=HK*1.0D+06                                                   ANH04990
      HHJCM=HHJ*1.0D+06                                                 ANH05000
      HHJKCM=HHJK*1.0D+06                                               ANH05010
      HHKCM=HHK*1.0D+06                                                 ANH05020
      HJMH=HJ*CL*1.0D+04                                                ANH05030
      HJKMH=HJK*CL*1.0D+04                                              ANH05040
      HKJMH=HKJ*CL*1.0D+04                                              ANH05050
      HKMH=HK*CL*1.0D+04                                                ANH05060
      HHJMH=HHJ*CL*1.0D+04                                              ANH05070
      HHJKMH=HHJK*CL*1.0D+04                                            ANH05080
      HHKMH=HHK*CL*1.0D+04                                              ANH05090
      WRITE(6,16) HJCM,HJMH,HJKCM,HJKMH,HKJCM,HKJMH,HKCM,HKMH,          ANH05100
     1            HHJCM,HHJMH,HHJKCM,HHJKMH,HHKCM,HHKMH                 ANH05110
C                                                                       ANH05120
C   S111=-(4*R6)/(B-C)                                                  ANH05130
      SIGMA2=SIGMA*SIGMA                                                ANH05140
      S111=-(4.0D+00*R6)/DBC                                            ANH05150
      S111SQ=S111*S111                                                  ANH05160
      HJ=P600+P204*TWO                                                  ANH05170
      HJK=P420-P204*12.0D+00+P024*TWO+P006*SIGMA*16.0D+00               ANH05180
     1   -(R5-SIGMA*R6*TWO)*S111*16.0D+00                               ANH05190
     2   +D2ABC*S111SQ*8.0D+00                                          ANH05200
      HKJ=P240+P420*(10.0D+00/3.0D+00)-P204*30.0D+00                    ANH05210
     1        -HJK*(10.0D+00/3.0D+00)                                   ANH05220
      HK=P060-P420*(7.0D+00/3.0D+00)+P204*28.0D+00+HJK*(7.0D+00/3.0D+00)ANH05230
      HHJ=P402+P006                                                     ANH05240
      HHJK=P222-P006*10.0D+00+P204*SIGMA*4.0D+00                        ANH05250
     1    +(DJK-SIGMA*DELJ*TWO-R6*4.0D+00)*S111*TWO                     ANH05260
     2    -DBC*S111SQ*4.0D+00                                           ANH05270
      HHK=P042+P024*(SIGMA*4.0D+00/3.0D+00)                             ANH05280
     1   +(9.0D+00+SIGMA2*32.0D+00/3.0D+00)*P006                        ANH05290
     2   +(DK-(SIGMA*R5*TWO)/3.0D+00+TWO*(ONE+8.0D+00*SIGMA2/3.0D+00)*R6ANH05300
     3    )*S111*4.0D+00+(6.0D+00+10.0D+00*SIGMA2/3.0D+00)              ANH05310
     4   *DBC*S111SQ                                                    ANH05320
C                                                                       ANH05330
      HJCM=HJ*1.0D+06                                                   ANH05340
      HJKCM=HJK*1.0D+06                                                 ANH05350
      HKJCM=HKJ*1.0D+06                                                 ANH05360
      HKCM=HK*1.0D+06                                                   ANH05370
      HHJCM=HHJ*1.0D+06                                                 ANH05380
      HHJKCM=HHJK*1.0D+06                                               ANH05390
      HHKCM=HHK*1.0D+06                                                 ANH05400
      HJMH=HJ*CL*1.0D+04                                                ANH05410
      HJKMH=HJK*CL*1.0D+04                                              ANH05420
      HKJMH=HKJ*CL*1.0D+04                                              ANH05430
      HKMH=HK*CL*1.0D+04                                                ANH05440
      HHJMH=HHJ*CL*1.0D+04                                              ANH05450
      HHJKMH=HHJK*CL*1.0D+04                                            ANH05460
      HHKMH=HHK*CL*1.0D+04                                              ANH05470
      WRITE(6,17) HJCM,HJMH,HJKCM,HJKMH,HKJCM,HKJMH,HKCM,HKMH,          ANH05480
     1            HHJCM,HHJMH,HHJKCM,HHJKMH,HHKCM,HHKMH                 ANH05490
C                                                                       ANH05500
      RETURN                                                            ANH05510
      END                                                               ANH05520
      SUBROUTINE ASYTOP(ZETA,F3Q,F4Q,XIJ)                               ANH05530
C   THE ANHARMONICITY CONTSTANTS FOR AN ASYMMETRIC TOP                  ANH05540
      IMPLICIT REAL*8 (A-H,O-Z)                                         ANH05550
      DIMENSION ZETA(N3N,N3N,3)                                         ANH05560
      DIMENSION F3Q(N3N,N3N,N3N),F4Q(N3N,N3N,N3N,N3N)                   ANH05570
      DIMENSION XIJ(N3N,N3N)                                            ANH05580
      COMMON/VIB101/NATOM,N3N,NATRI,ILIN,NVIB                           ANH05590
      COMMON/VIB106/ROTAA(3),ROTGC(3),ROTCM(3),ROTMH(3)                 ANH05600
      COMMON/VIB108/IAXIS(3),NDEG(150),NDAB(150,5),IMAG(150)            ANH05610
      COMMON/VIB109/CLIMIT,FLIM1,FLIM2                                  ANH05620
      COMMON/VIB111/ANHTOT,HZPVE,AZPVE
      COMMON/VIB203/IOFF(150),IPRNT                                     ANH05630
      COMMON/VIB204/SQM(150),ROOT(150),FREQ(150)                        ANH05640
      COMMON/VIB205/IFREQ,NFRQ(150)                                     ANH05650
      DATA A00,ONE / 0.0D+00 , 1.0D+00 /                                ANH05660
      DATA WLIMIT / 1.0D+00 /                                           ANH05670
      DATA PLIMIT / 1.0D-05 /                                           ANH05680
    1 FORMAT(//,2X,' ********************************************'      ANH05690
     1        /,2X,' ***THE VIBRATIONAL ANHARMONICTY CONSTANTS***'      ANH05700
     2        /,2X,' ********************************************')     ANH05710
    2 FORMAT(//,2X,' ANHARMONIC CONSTANTS IN CM-1'/                     ANH05720
     1          2X,' IVIB',3X,' JVIB',8X,' FREQ(I)',5X,' FREQ(J)',      ANH05730
     2          7X,' VALU1',7X,' VALU2',7X,' VALU3',7X,' VALU4',        ANH05740
     3          7X,' XIJ'/)                                             ANH05750
    3 FORMAT(2X,I5,3X,I5,3X,7(F13.3))                                   ANH05760
    4 FORMAT(//,2X,' ANHARMONIC CONSTANTS IN CM-1'/                     ANH05770
     1          4X,' QI',5X,' QJ',8X,' FREQ(I)',5X,' FREQ(J)',          ANH05780
     2          7X,' XIJ'/)                                             ANH05790
    5 FORMAT(2X,I5,3X,I5,3X,3(F13.3))                                   ANH05800
    6 FORMAT(//,2X,' THEORETICAL FREQUENCIES IN CM-1'/                  ANH05810
     1          2X,' NO.',9X,' HARMONIC',14X,' DELTA',                  ANH05820
     2         13X,' FUNDAMENTAL'/)                                     ANH05830
    7 FORMAT(2X,I5,4X,F13.3,7X,F13.3,12X,F13.3)                         ANH05840
    8 FORMAT(//,2X,' :::::::::::::::::::::::::::::::::::::::::::::::'/  ANH05850
     1          2X,' :::ANHARMONIC RESONANCE (TYPE I) IS OBSERVED:::'/  ANH05860
     2          2X,' :::::::::::::::::::::::::::::::::::::::::::::::'/) ANH05870
    9 FORMAT(//,2X,' ::::::::::::::::::::::::::::::::::::::::::::::::'/ ANH05880
     1          2X,' :::ANHARMONIC RESONANCE (TYPE II) IS OBSERVED:::'/ ANH05890
     2          2X,' ::::::::::::::::::::::::::::::::::::::::::::::::'/)ANH05900
   10 FORMAT(2X,' NUMBER OF OCCURENCES IS ',I5/                         ANH05910
     1       2X,' RESONANCE LIMIT IS SET TO ',F10.3,' CM-1'/)           ANH05920
C                                                                       ANH05930
C:::::::::::::::::::::::::::::::::::::::                                ANH05940
C:::CALCULATE ANHARMONICITY CONSTANTS:::                                ANH05950
C:::::::::::::::::::::::::::::::::::::::                                ANH05960
      IANH1=0                                                           ANH05970
      IANH2=0                                                           ANH05980
C                                                                       ANH05990
      WRITE(6,1)                                                        ANH06000
      IF(IPRNT.GT.2)                                                    ANH06010
     *WRITE(6,2)                                                        ANH06020
      DO 110 IVIB=1,NVIB                                                ANH06030
      FRQI=FREQ(IVIB)                                                   ANH06040
      FRQI2=FRQI*FRQI                                                   ANH06050
      DO 110 JVIB=1,NVIB                                                ANH06060
      FRQJ=FREQ(JVIB)                                                   ANH06070
      FRQJ2=FRQJ*FRQJ                                                   ANH06080
C                                                                       ANH06090
      VALU1=A00                                                         ANH06100
      VALU2=A00                                                         ANH06110
      VALU3=A00                                                         ANH06120
      VALU4=A00                                                         ANH06130
      XRS=A00                                                           ANH06140
      IF(DABS(FRQI).LE.WLIMIT) GO TO 210                                ANH06150
      IF(DABS(FRQJ).LE.WLIMIT) GO TO 210                                ANH06160
      IF(IVIB.NE.JVIB) GO TO 205                                        ANH06170
C                                                                       ANH06180
C++++++++++++++++++++++++                                               ANH06190
C+++DIAGONAL CONSTANTS+++                                               ANH06200
C++++++++++++++++++++++++                                               ANH06210
      VALU1=F4Q(IVIB,IVIB,IVIB,IVIB)                                    ANH06220
      DO 101 KVIB=1,NVIB                                                ANH06230
      FRQK=FREQ(KVIB)                                                   ANH06240
      IF(DABS(FRQK).LE.WLIMIT) GO TO 101                                ANH06250
      FRQK2=FRQK*FRQK                                                   ANH06260
      FAC=F3Q(IVIB,IVIB,KVIB)                                           ANH06270
      FACT=FAC*FAC                                                      ANH06280
      IF(DABS(FRQI+FRQI-FRQK).LE.FLIM1) GO TO 201                       ANH06290
      DNUM=FRQI2*8.0D+00-FRQK2*3.0D+00                                  ANH06300
      DENM=FRQI2*4.0D+00-FRQK2                                          ANH06310
      VALU2=VALU2+(DNUM*FACT)/(FRQK*DENM)                               ANH06320
C***********************************************************************ANH06330
C***FOLLWOING CODES ARE FOR AN ALTERNATE EXPRESSION                    *ANH06340
C***THEY ARE TESTED TO GIVE THE SAME RESUTLS                           *ANH06350
C***********************************************************************ANH06360
C*    DNUM=ONE/((FRQI+FRQI+FRQK)*2.0D+00)-ONE/((FRQI+FRQI-FRQK)*2.0D+00)ANH06370
C*   1    +2.0D+00/FRQK                                                *ANH06380
C*    VALU2=VALU2+(DNUM*FACT)                                          *ANH06390
C***********************************************************************ANH06400
      GO TO 101                                                         ANH06410
C/////////////////////////////////                                      ANH06420
C///FERMI RESONANCE CORRECTIONS///                                      ANH06430
C/////////////////////////////////                                      ANH06440
  201 CONTINUE                                                          ANH06450
      DNUM=ONE/((FRQI+FRQI+FRQK)*2.0D+00)+2.0D+00/FRQK                  ANH06460
      IF(DABS(FAC).GT.PLIMIT) THEN                                      ANH06470
        write(6,*) 'Fermi resonance'
        IANH1=IANH1+1                                                   ANH06480
      END IF                                                            ANH06490
      VALU2=VALU2+FACT*DNUM                                             ANH06500
C////////////////////////////                                           ANH06510
C///END OF FERMI RESONANCE///                                           ANH06520
C////////////////////////////                                           ANH06530
  101 CONTINUE                                                          ANH06540
      XRS=(VALU1-VALU2)/16.0D+00                                        ANH06550
      GO TO 210                                                         ANH06560
C                                                                       ANH06570
C++++++++++++++++++++++++++++                                           ANH06580
C+++OFF-DIAGONAL CONSTANTS+++                                           ANH06590
C++++++++++++++++++++++++++++                                           ANH06600
  205 CONTINUE                                                          ANH06610
      VALU1=F4Q(IVIB,IVIB,JVIB,JVIB)                                    ANH06620
      DO 102 KVIB=1,NVIB                                                ANH06630
      FRQK=FREQ(KVIB)                                                   ANH06640
      IF(DABS(FRQK).LE.WLIMIT) GO TO 102                                ANH06650
      FRQK2=FRQK*FRQK                                                   ANH06660
      WPPP=FRQI+FRQJ+FRQK                                               ANH06670
      WPPM=FRQI+FRQJ-FRQK                                               ANH06680
      WPMP=FRQI-FRQJ+FRQK                                               ANH06690
      WPMM=FRQI-FRQJ-FRQK                                               ANH06700
      DIJK=WPPP*WPPM*WPMP*WPMM                                          ANH06710
      DNUM2=F3Q(IVIB,IVIB,KVIB)*F3Q(KVIB,JVIB,JVIB)                     ANH06720
      VALU2=VALU2+DNUM2/FRQK                                            ANH06730
      FAC=F3Q(IVIB,JVIB,KVIB)                                           ANH06740
      FACT=FAC*FAC                                                      ANH06750
      IF(DABS(WPPM).LE.FLIM2) GO TO 206                                 ANH06760
      IF(DABS(WPMP).LE.FLIM2) GO TO 207                                 ANH06770
      IF(DABS(WPMM).LE.FLIM2) GO TO 208                                 ANH06780
      DNUM3=(FRQK2-FRQI2-FRQJ2)*FRQK                                    ANH06790
      VALU3=VALU3+(DNUM3*FACT)/DIJK                                     ANH06800
C***********************************************************************ANH06810
C***FOLLOWING CODES ARE FOR AN ALTERNATE EXPRESSION                    *ANH06820
C***THEY ARE TESTED TO GIVE THE SAME RESULTS                           *ANH06830
C***********************************************************************ANH06840
C*    DNUM3=ONE/WPPP-ONE/WPPM+ONE/WPMP-ONE/WPMM                        *ANH06850
C*    VALU3=VALU3+(DNUM3*FACT)/4.0D+00                                 *ANH06860
C***********************************************************************ANH06870
      GO TO 102                                                         ANH06880
C/////////////////////////////////                                      ANH06890
C///FERMI RESONANCE CORRECTIONS///                                      ANH06900
C/////////////////////////////////                                      ANH06910
C   WI+WJ=WK                                                            ANH06920
  206 CONTINUE                                                          ANH06930
      DNUM3=ONE/WPPP+ONE/WPMP-ONE/WPMM                                  ANH06940
      IF(DABS(FAC).GT.PLIMIT) THEN                                      ANH06950
        IANH2=IANH2+1                                                   ANH06960
      END IF                                                            ANH06970
      VALU3=VALU3+(DNUM3*FACT)/4.0D+00                                  ANH06980
      GO TO 102                                                         ANH06990
C   WI+WK=WJ                                                            ANH07000
  207 CONTINUE                                                          ANH07010
      DNUM3=ONE/WPPP-ONE/WPPM-ONE/WPMM                                  ANH07020
      IF(DABS(FAC).GT.PLIMIT) THEN                                      ANH07030
        IANH2=IANH2+1                                                   ANH07040
      END IF                                                            ANH07050
      VALU3=VALU3+(DNUM3*FACT)/4.0D+00                                  ANH07060
      GO TO 102                                                         ANH07070
C   WJ+WK=WI                                                            ANH07080
  208 CONTINUE                                                          ANH07090
      DNUM3=ONE/WPPP-ONE/WPPM+ONE/WPMP                                  ANH07100
      IF(DABS(FAC).GT.PLIMIT) THEN                                      ANH07110
        IANH2=IANH2+1                                                   ANH07120
      END IF                                                            ANH07130
      VALU3=VALU3+(DNUM3*FACT)/4.0D+00                                  ANH07140
C////////////////////////////                                           ANH07150
C///END OF FERMI RESONANCE///                                           ANH07160
C////////////////////////////                                           ANH07170
  102 CONTINUE                                                          ANH07180
      ZFACT=A00                                                         ANH07190
      DO 103 IABC=1,3                                                   ANH07200
      II=IAXIS(IABC)                                                    ANH07210
      ZET=ZETA(IVIB,JVIB,II)                                            ANH07220
      ZFACT=ZFACT+ZET*ZET*ROTCM(IABC)                                   ANH07230
  103 CONTINUE                                                          ANH07240
      VALU4=(FRQI/FRQJ+FRQJ/FRQI)*ZFACT                                 ANH07250
      XRS=(VALU1-VALU2)/4.0D+00-VALU3/2.0D+00+VALU4                     ANH07260
C                                                                       ANH07270
  210 CONTINUE                                                          ANH07280
      XIJ(IVIB,JVIB)=XRS                                                ANH07290
      IF(IPRNT.LE.2) GO TO 110                                          ANH07300
      WRITE(6,3) IVIB,JVIB,FRQI,FRQJ,VALU1,VALU2,VALU3,VALU4,XRS        ANH07310
  110 CONTINUE                                                          ANH07320
      IF(IANH1.NE.0) THEN                                               ANH07330
        write(6,*) 'wtf'
        WRITE(6,8)                                                      ANH07340
        WRITE(6,10) IANH1,FLIM1                                         ANH07350
      END IF                                                            ANH07360
      IF(IANH2.NE.0) THEN                                               ANH07370
        WRITE(6,8)                                                      ANH07380
        WRITE(6,10) IANH2,FLIM2                                         ANH07390
      END IF                                                            ANH07400
C                                                                       ANH07410
      IF(IPRNT.GT.3) GO TO 211                                          ANH07420
      WRITE(6,4)                                                        ANH07430
      DO 112 II=1,NVIB                                                  ANH07440
      IVIB=NFRQ(II)                                                     ANH07450
      FRQI=FREQ(IVIB)                                                   ANH07460
      DO 112 JJ=1,II                                                    ANH07470
      JVIB=NFRQ(JJ)                                                     ANH07480
      FRQJ=FREQ(JVIB)                                                   ANH07490
      XRS=XIJ(IVIB,JVIB)                                                ANH07500
      WRITE(6,5) II,JJ,FRQI,FRQJ,XRS                                    ANH07510
  112 CONTINUE                                                          ANH07520
C                                                                       ANH07530
C:::::::::::::::::::::::::::::::::::::::                                ANH07540
C:::CALCULATE FUNDAMENTAL FREQUENCIES:::                                ANH07550
C:::::::::::::::::::::::::::::::::::::::                                ANH07560
  211 CONTINUE                                                          ANH07570
      WRITE(6,6)                                                        ANH07580
      DO 115 II=1,NVIB                                                  ANH07590
      IVIB=NFRQ(II)                                                     ANH07600
      VALU1=FREQ(IVIB)                                                  ANH07610
      VALU2=XIJ(IVIB,IVIB)*2.0D+00                                      ANH07620
      VALU3=A00                                                         ANH07630
      DO 114 JJ=1,NVIB                                                  ANH07640
      JVIB=NFRQ(JJ)                                                     ANH07650
      IF(IVIB.EQ.JVIB) GO TO 114                                        ANH07660
      VALU3=VALU3+XIJ(IVIB,JVIB)                                        ANH07670
  114 CONTINUE                                                          ANH07680
      DELTA=VALU2+VALU3*0.5D+00                                         ANH07690
      VALUT=VALU1+DELTA                                                 ANH07700
      WRITE(6,7) II,VALU1,DELTA,VALUT                                   ANH07710
  115 CONTINUE                                                          ANH07720
C                                                                       ANH07730
C::::::::::::::::::::::::::::::::::::::::::::::::::
C:::CALCULATE ANHARMONIC VIBRATIONAL ENERGY :::::::
C::::::::::::::::::::::::::::::::::::::::::::::::::
C   ZPVE = 1/2 Sum_(i)^{unique} d_i W_i
C        + 1/4 Sum_(i>=j)^{unique} d_i d_j X_(ij)
C::::::::::::::::::::::::::::::::::::::::::::::::::
      TERM1=A00
      TERM2=A00
      TERM3=A00
      DO 117 IVIB=1,NVIB
      FRQI=FREQ(IVIB)
      XRR=XIJ(IVIB,IVIB)
      TERM1=TERM1+FRQI/(2.0D+00)
      TERM2=TERM2+XRR/(4.0D+00)
      DO 116 JVIB=1,NVIB
      IF(IVIB.EQ.JVIB) GO TO 116
      XRS=XIJ(IVIB,JVIB)
      TERM3=TERM3+XRS/(8.0D+00)
  116 CONTINUE
  117 CONTINUE
      ANHTOT=TERM1+TERM2+TERM3
      IF(IPRNT.GT.10) THEN
        WRITE(6,*) ' TERM1 (Q2)        = ',TERM1
        WRITE(6,*) ' TERM2 (XRR-ND)    = ',TERM2
        WRITE(6,*) ' TERM3 (XRS-ND_ND) = ',TERM3
        WRITE(6,*) ' ANHTOT            = ',ANHTOT
      END IF
C
      RETURN                                              
      END                                                
      SUBROUTINE ASYZPE(ZETA,F3Q,F4Q) 
C   THE ANHARMONIC ZPVE FOR AN ASYMMETRIC TOP     
      IMPLICIT REAL*8 (A-H,O-Z)                                      
      DIMENSION ZETA(N3N,N3N,3)                                     
      DIMENSION F3Q(N3N,N3N,N3N),F4Q(N3N,N3N,N3N,N3N)              
      COMMON/VIB101/NATOM,N3N,NATRI,ILIN,NVIB                     
      COMMON/VIB106/ROTAA(3),ROTGC(3),ROTCM(3),ROTMH(3)               
      COMMON/VIB108/IAXIS(3),NDEG(150),NDAB(150,5),IMAG(150)         
      COMMON/VIB203/IOFF(150),IPRNT                                
      COMMON/VIB204/SQM(150),ROOT(150),FREQ(150)                  
      COMMON/VIB111/ANHTOT,HZPVE,AZPVE
      DATA A00,ONE / 0.0D+00 , 1.0D+00 /                            
    1 FORMAT(//,2X,' **************************************************' 
     1        /,2X,' ***THE ANHARMONIC ZERO-POINT VIBRATIONAL ENERGY***'   
     2        /,2X,' **************************************************'  
     3      )
    2 FORMAT(//,2X,' FOR AN ASYMMETRIC TOP MOLECULE')
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
C   Following codes for a kinetic term
C   are equivalent to the codes 
C   described below
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
      IF(IPRNT.GT.2)
     *WRITE(6,3) VALU1,VALU2,VALU3,VALU4,VALU5 
      WRITE(6,4) HZPVE,ZCUBIC,ZQUART,ZANHR,ZKINET,AZPVE
C
      RETURN
      END
      SUBROUTINE ASYG00(ZETA,F3Q,F4Q)
C   THE VIBRATIONAL G0 VALUE FOR AN ASYMMETRIC TOP
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
    5 FORMAT(//,2X,' FOR AN ASYMMETRIC TOP MOLECULE')
    6 FORMAT(   2X,' VALU1      = ',F20.10,' CM-1',/
     1          2X,' VALU2      = ',F20.10,' CM-1',/
     2          2X,' VALU3      = ',F20.10,' CM-1',/
     3          2X,' VALU4      = ',F20.10,' CM-1',/
     4          2X,' VALU5      = ',F20.10,' CM-1',/)
    7 FORMAT(   2X,' CUBIC (Q3) = ',F20.10,' CM-1',/
     1          2X,' QUART.(Q4) = ',F20.10,' CM-1',/
     2          2X,' KINET.(K)  = ',F20.10,' CM-1',/
     3          2X,' G_0        = ',F20.10,' CM-1',/
     4          2X,' (Q3+Q4+K)'/)
    8 FORMAT(   2X,' INTERNAL CONSISTENCY TEST',/
     1          2X,' ZPVE(Q2)   = ',F20.10,' CM-1',/
     2          2X,' ZPVE(TOTAL)= ',F20.10,' CM-1',//
     3          2X,' ZPVE(XRS)  = ',F20.10,' CM-1',/
     4          2X,' ZPVE(ALL)  = ',F20.10,' CM-1',/
     5          2X,' (XRS+G_0)'/)
C
      WRITE(6,1)
      IANH1=0
      IANH2=0
      VALU1=A00
      VALU2=A00
      VALU3=A00
      VALU4=A00
      VALU5=A00
      DO 103 IVIB=1,NVIB
      FRQI=FREQ(IVIB)
      FRIII=F3Q(IVIB,IVIB,IVIB)
      FRIIII=F4Q(IVIB,IVIB,IVIB,IVIB)
      VALU1=VALU1+FRIIII
      VALU2=VALU2+(FRIII*FRIII)/(FRQI)
      FRISQ4=FRQI*FRQI*4.0D+00
      DO 102 JVIB=1,NVIB
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
        write(6, *) 'Fermi resonance'
        IANH1=IANH1+1
      END IF
      DNUM=ONE/(FRQI+FRQI+FRQJ)
      VALU3=VALU3-(FIIJSQ*DNUM)/(2.0D+00)
  202 CONTINUE
      DO 101 KVIB=1,NVIB
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
      DO 107 IABC=1,3
      II=IAXIS(IABC)
      BE=ROTCM(IABC)
      TERM=ONE
      DO 106 IVIB=1,NVIB
      DO 105 JVIB=1,NVIB
      ZET=ZETA(IVIB,JVIB,II)
      ZETSQ=ZET*ZET
      TERM=TERM+ZETSQ
  105 CONTINUE
  106 CONTINUE
      VALU5=VALU5+TERM*BE
  107 CONTINUE
C
      VALU1 =+VALU1/(64.0D+00)
      VALU2 =-VALU2/(36.0D+00)
      VALU3 =+(VALU3*3.0D+00)/(64.0D+00)
      VALU4 =-VALU4/(4.0D+00)
      VALU5 =-VALU5/(4.0D+00)
      GCUBIC=VALU2+VALU3+VALU4
      GQUART=VALU1
      GKINET=VALU5
      G00   = GCUBIC+GQUART+GKINET
      ANHG00=ANHTOT+G00
      WRITE(6,5)
      IF(IPRNT.GT.2)
     *WRITE(6,6) VALU1,VALU2,VALU3,VALU4,VALU5
      WRITE(6,7) GCUBIC,GQUART,GKINET,G00
      WRITE(6,8) HZPVE,AZPVE,ANHTOT,ANHG00
C
      RETURN
      END
