C     PROGRAM ANHSPH                                                    ANH00010
C******************************************************                 ANH00020
C***LAST UPDATED ON JUNE 05, 1986 BY YUKIO YAMAGUCHI***                 ANH00030
C******************************************************                 ANH00040
      SUBROUTINE LTYPSP(ZETA,AXY,F3Q,EE,CLTYP)                          ANH00050
C   THE CALCULATION OF L-TYPE DOUBLING CONSTANTS                        ANH00060
C   FOR A SYMMETRIC TOP MOLECULE                                        ANH00070
      IMPLICIT REAL*8 (A-H,O-Z)                                         ANH00080
      DIMENSION ZETA(N3N,N3N,3),AXY(3,3,N3N),F3Q(N3N,N3N,N3N)           ANH00090
      DIMENSION EE(N3N,N3N),CLTYP(N3N,N3N)                              ANH00100
      COMMON/VIB101/NATOM,N3N,NATRI,ILIN,NVIB                           ANH00110
      COMMON/VIB106/ROTAA(3),ROTGC(3),ROTCM(3),ROTMH(3)                 ANH00120
      COMMON/VIB108/IAXIS(3),NDEG(150),NDAB(150,5),IMAG(150)            ANH00130
      COMMON/VIB109/CLIMIT,FLIM1,FLIM2                                  ANH00140
      COMMON/VIB203/IOFF(150),IPRNT                                     ANH00150
      COMMON/VIB204/SQM(150),ROOT(150),FREQ(150)                        ANH00160
      COMMON/VIB205/IFREQ,NFRQ(150)                                     ANH00170
      DATA PH,AVN / 6.626176D+00 , 6.022045D+00 /                       ANH00180
      DATA PI,CL / 3.1415926536D+00 , 2.99792458D+00 /                  ANH00190
      DATA A00,HALF,ONE,TWO / 0.0D+00 , 0.5D+00 , 1.0D+00 , 2.0D+00 /   ANH00200
      DATA ALIMIT / 1.0D-04 /                                           ANH00210
      DATA WLIMIT / 1.0D-02 /                                           ANH00220
    1 FORMAT(//,2X,' ***********************************'/              ANH00230
     1          2X,' ***THE L-TYPE DOUBLING CONSTANTS***'/              ANH00240
     2          2X,' ***********************************')              ANH00250
    2 FORMAT(//,2X,' ZETA MATRIX, IABC = ',I5/)                         ANH00260
    3 FORMAT(//,2X,' F3(QI*QJ*QK) MATRIX (IN CM-1), IVIB = ',I5/)       ANH00270
    4 FORMAT(//,2X,' AXY MATRIX, IVIB = ',I5/)                          ANH00280
    5 FORMAT(2X,I5,4F10.5)                                              ANH00290
    6 FORMAT(//,2X,' QVIB',5X,' IVIB',7X,' FREQENCY',9X,' HARMONIC',    ANH00300
     1          9X,' ANHARMONIC',7X,' TOTAL',9X,' TOTAL'/               ANH00310
     2         24X,' IN CM-1',10X,' IN CM-1',10X,' IN CM-1',            ANH00320
     3         10X,' IN CM-1',7X,' IN MHZ'/)                            ANH00330
    7 FORMAT(2X,I5,5X,I5,8X,F10.5,F15.5,3X,F15.5,3X,F15.5,3X,F15.5)     ANH00340
C                                                                       ANH00350
      WRITE(6,1)                                                        ANH00360
      RCONST=1.0D-01*PI*DSQRT(CL/(PH*AVN))                              ANH00370
C                                                                       ANH00380
C   NOTE : THE SYMMETRY AXIS IS CHOSEN AS THE Z AXIS                    ANH00390
      IXX=1                                                             ANH00400
      IYY=2                                                             ANH00410
      IZZ=3                                                             ANH00420
C                                                                       ANH00430
      IF(IPRNT.LE.2) GO TO 201                                          ANH00440
      DO 101 IABC=1,3                                                   ANH00450
      WRITE(6,2) IABC                                                   ANH00460
      CALL MATOUT(ZETA(1,1,IABC),N3N,N3N,N3N,N3N,6)                     ANH00470
  101 CONTINUE                                                          ANH00480
      DO 103 III=1,NVIB                                                 ANH00490
      IVIB=NFRQ(III)                                                    ANH00500
      WRITE(6,3) IVIB                                                   ANH00510
      DO 102 I=1,NVIB                                                   ANH00520
      II=NFRQ(I)                                                        ANH00530
      DO 102 J=1,NVIB                                                   ANH00540
      JJ=NFRQ(J)                                                        ANH00550
      EE(I,J)=F3Q(II,JJ,IVIB)                                           ANH00560
  102 CONTINUE                                                          ANH00570
      CALL MATOUT(EE,N3N,N3N,NVIB,NVIB,6)                               ANH00580
  103 CONTINUE                                                          ANH00590
      DO 104 IVIB=1,NVIB                                                ANH00600
      WRITE(6,4) IVIB                                                   ANH00610
      CALL MATOUT(AXY(1,1,IVIB),3,3,3,3,6)                              ANH00620
  104 CONTINUE                                                          ANH00630
C                                                                       ANH00640
C:::::::::::::::::::::::::::::::                                        ANH00650
C:::ONLY FOR DEGENERATE MODES:::                                        ANH00660
C:::::::::::::::::::::::::::::::                                        ANH00670
  201 CONTINUE                                                          ANH00680
      DO 120 IVIB=1,NVIB                                                ANH00690
      IF(NDEG(IVIB).EQ.1) GO TO 120                                     ANH00700
      FRQI=FREQ(IVIB)                                                   ANH00710
      FRQI2=FRQI*FRQI                                                   ANH00720
      AXZVAL=AXY(IXX,IZZ,IVIB)                                          ANH00730
      AYYVAL=AXY(IYY,IYY,IVIB)                                          ANH00740
      IF(DABS(AXZVAL).LE.ALIMIT.AND.DABS(AYYVAL).LE.ALIMIT) GO TO 202   ANH00750
C                                                                       ANH00760
C   FOR AN A-DEGENERATE (T1) MODE                                       ANH00770
      VALU1=-0.75D+00*(AXZVAL**2/ROTAA(IZZ))                            ANH00780
      VALU2=A00                                                         ANH00790
      VALU3=A00                                                         ANH00800
C   SUM OVER NON-DEGENERATE MODES                                       ANH00810
      DO 112 JVIB=1,NVIB                                                ANH00820
      IF(NDEG(JVIB).NE.1) GO TO 112                                     ANH00830
      FRQJ=FREQ(JVIB)                                                   ANH00840
      FRQJ2=FRQJ*FRQJ                                                   ANH00850
      FACT2=(FRQI2*3.0D+00+FRQJ2)/(FRQI2-FRQJ2)                         ANH00860
      VALU2=VALU2+(ZETA(IVIB,JVIB,IYY)**2-ZETA(IVIB,JVIB,IXX)**2)*FACT2 ANH00870
  112 CONTINUE                                                          ANH00880
      DO 113 JVIB=1,NVIB                                                ANH00890
      IF(NDEG(JVIB).EQ.1) GO TO 113                                     ANH00900
      FRQJ=FREQ(JVIB)                                                   ANH00910
      FRQJ2=FRQJ*FRQJ                                                   ANH00920
      AXZVAL=AXY(IXX,IZZ,JVIB)                                          ANH00930
      AYYVAL=AXY(IYY,IYY,JVIB)                                          ANH00940
      IF(DABS(AXZVAL).LE.ALIMIT.AND.DABS(AYYVAL).LE.ALIMIT) GO TO 113   ANH00950
      FACT3=(FRQI*F3Q(IVIB,IVIB,JVIB)*RCONST*TWO)/(FRQJ**(1.5D+00))     ANH00960
      VALU3=VALU3+AYYVAL*FACT3                                          ANH00970
  113 CONTINUE                                                          ANH00980
      GO TO 203                                                         ANH00990
C                                                                       ANH01000
C   FOR A B-DEGENERATE (T2) MODE                                        ANH01010
  202 CONTINUE                                                          ANH01020
      AYZVAL=AXY(IYY,IZZ,IVIB)                                          ANH01030
      AXYVAL=AXY(IXX,IYY,IVIB)                                          ANH01040
      IF(DABS(AYZVAL).LE.ALIMIT.AND.DABS(AXYVAL).LE.ALIMIT) GO TO 120   ANH01050
      VALU1=-0.75D+00*(AYZVAL**2/ROTAA(IZZ))                            ANH01060
      VALU2=A00                                                         ANH01070
      VALU3=A00                                                         ANH01080
      DO 116 JVIB=1,NVIB                                                ANH01090
      IF(NDEG(JVIB).NE.1) GO TO 116                                     ANH01100
      FRQJ=FREQ(JVIB)                                                   ANH01110
      FRQJ2=FRQJ*FRQJ                                                   ANH01120
      FACT2=(FRQI2*3.0D+00+FRQJ2)/(FRQI2-FRQJ2)                         ANH01130
      VALU2=VALU2-(ZETA(IVIB,JVIB,IYY)**2-ZETA(IVIB,JVIB,IXX)**2)*FACT2 ANH01140
  116 CONTINUE                                                          ANH01150
      DO 117 JVIB=1,NVIB                                                ANH01160
      IF(NDEG(JVIB).EQ.1) GO TO 117                                     ANH01170
      FRQJ=FREQ(JVIB)                                                   ANH01180
      FRQJ2=FRQJ*FRQJ                                                   ANH01190
      AXZVAL=AXY(IXX,IZZ,JVIB)                                          ANH01200
      AYYVAL=AXY(IYY,IYY,JVIB)                                          ANH01210
      IF(DABS(AXZVAL).LE.ALIMIT.AND.DABS(AYYVAL).LE.ALIMIT) GO TO 117   ANH01220
      FACT3=(FRQI*F3Q(IVIB,IVIB,JVIB)*RCONST*TWO)/(FRQJ**(1.5D+00))     ANH01230
      VALU3=VALU3-AYYVAL*FACT3                                          ANH01240
  117 CONTINUE                                                          ANH01250
C                                                                       ANH01260
  203 CONTINUE                                                          ANH01270
      VALUH=(VALU1+VALU2)*(ROTCM(IXX)*ROTCM(IXX)*TWO)/FRQI              ANH01280
      VALUA=(VALU3)*(ROTCM(IXX)*ROTCM(IXX)*TWO)/FRQI                    ANH01290
      VALUT=VALUH+VALUA                                                 ANH01300
      CLTYP(IVIB,1)=VALUH                                               ANH01310
      CLTYP(IVIB,2)=VALUA                                               ANH01320
      CLTYP(IVIB,3)=VALUT                                               ANH01330
      IF(IPRNT.LE.2) GO TO 120                                          ANH01340
      WRITE(6,5) IVIB,FRQI,VALUH,VALUA,VALUT                            ANH01350
  120 CONTINUE                                                          ANH01360
C                                                                       ANH01370
      WRITE(6,6)                                                        ANH01380
      DO 125 II=1,NVIB                                                  ANH01390
      IVIB=NFRQ(II)                                                     ANH01400
      IF(NDEG(IVIB).LE.1) GO TO 125                                     ANH01410
      FRQI=FREQ(IVIB)                                                   ANH01420
      VALUH=CLTYP(IVIB,1)                                               ANH01430
      VALUA=CLTYP(IVIB,2)                                               ANH01440
      VALUT=CLTYP(IVIB,3)                                               ANH01450
      VALUM=VALUT*CL*1.0D+04                                            ANH01460
      WRITE(6,7) II,IVIB,FRQI,VALUH,VALUA,VALUT,VALUM                   ANH01470
  125 CONTINUE                                                          ANH01480
C                                                                       ANH01490
      RETURN                                                            ANH01500
      END                                                               ANH01510
      SUBROUTINE ROVSPH(ZETA,AXY,F3Q,EE,RVCNST)                         ANH01520
C   THE CALCULATION OF VIBRATION-ROTATION INTERACTION CONSTANTS         ANH01530
C   FOR THE SPHERICAL TOP MOLECULES WITH ONLY REAL FREQUENCIES          ANH01540
      IMPLICIT REAL*8 (A-H,O-Z)                                         ANH01550
      CHARACTER*1 RAXIS                                                 ANH01560
      CHARACTER*25 RTYPE                                                ANH01570
      DIMENSION ZETA(N3N,N3N,3),AXY(3,3,N3N),F3Q(N3N,N3N,N3N)           ANH01580
      DIMENSION EE(N3N,N3N),RVCNST(NVIB,3)                              ANH01590
      COMMON/VIB101/NATOM,N3N,NATRI,ILIN,NVIB                           ANH01600
      COMMON/VIB103/PARA,WAVE,CONST,CYCL,CONV                           ANH01610
      COMMON/VIB106/ROTAA(3),ROTGC(3),ROTCM(3),ROTMH(3)                 ANH01620
      COMMON/VIB108/IAXIS(3),NDEG(150),NDAB(150,5),IMAG(150)            ANH01630
      COMMON/VIB109/CLIMIT,FLIM1,FLIM2                                  ANH01640
      COMMON/VIB203/IOFF(150),IPRNT                                     ANH01650
      COMMON/VIB204/SQM(150),ROOT(150),FREQ(150)                        ANH01660
      COMMON/VIB205/IFREQ,NFRQ(150)                                     ANH01670
      COMMON/VIB206/RAXIS(3),RTYPE(6)                                   ANH01680
      COMMON/VIB209/SALP(3)                                             ANH01690
      DATA PH,AVN / 6.626176D+00 , 6.022045D+00 /                       ANH01700
      DATA PI,CL / 3.1415926536D+00 , 2.99792458D+00 /                  ANH01710
      DATA A00,ONE,TWO / 0.0D+00 , 1.0D+00 , 2.0D+00 /                  ANH01720
      DATA WLIMIT / 1.0D-02 /                                           ANH01730
    1 FORMAT(//,2X,' **************************************************'ANH01740
     1        /,2X,' ***THE VIBRATION-ROTATION INTERACTION CONSTANTS***'ANH01750
     2        /,2X,' **************************************************'ANH01760
     3        )                                                         ANH01770
    2 FORMAT(//,2X,' ZETA MATRIX, IABC = ',I5/)                         ANH01780
    3 FORMAT(//,2X,' F3(QI*QJ*QK) MATRIX (IN CM-1), IVIB = ',I5/)       ANH01790
    4 FORMAT(//,2X,' AXY MATRIX, IVIB = ',I5/)                          ANH01800
    5 FORMAT(//,2X,' QVIB',5X,' IVIB',7X,' A',16X,' B',16X,' C'/        ANH01810
     1          24X,' IN CM-1',10X,' IN CM-1',10X,' IN CM-1'/)          ANH01820
    6 FORMAT(2X,I5,5X,I5,5X,F10.5,8X,F10.5,8X,F10.5)                    ANH01830
    7 FORMAT(//,2X,' PRINCIPAL MOMENTS OF INERTIA WITH ZERO-POINT VIBRATANH01840
     1IONAL CORRECTION'/                                                ANH01850
     2 4X,5H AXIS,9X,11H IN AMU.A+2,9X,10H IN G.CM+2,10X,               ANH01860
     3 8H IN CM-1,15X,7H IN MHZ/38X,8H (*D+39)/)                        ANH01870
    8 FORMAT(8X,A1,3X,3F20.10,F20.5)                                    ANH01880
    9 FORMAT(//,2X,' ::::::::::::::::::::::::::::::::::::'/             ANH01890
     1          2X,' :::CORIOLIS RESONANCE IS OBSERVED:::'/             ANH01900
     2          2X,' ::::::::::::::::::::::::::::::::::::'/)            ANH01910
   10 FORMAT(2X,' RESONANCE LIMIT IS SET TO ',F10.3,' CM-1'/)           ANH01920
C                                                                       ANH01930
      WRITE(6,1)                                                        ANH01940
      RCONST=1.0D-01*PI*DSQRT(CL/(PH*AVN))                              ANH01950
      IAA=IAXIS(1)                                                      ANH01960
      IBB=IAXIS(2)                                                      ANH01970
      ICC=IAXIS(3)                                                      ANH01980
      NCORIO=0                                                          ANH01990
C                                                                       ANH02000
      IF(IPRNT.LE.2) GO TO 201                                          ANH02010
      DO 101 IABC=1,3                                                   ANH02020
      WRITE(6,2) IABC                                                   ANH02030
      CALL MATOUT(ZETA(1,1,IABC),N3N,N3N,N3N,N3N,6)                     ANH02040
  101 CONTINUE                                                          ANH02050
      DO 103 III=1,NVIB                                                 ANH02060
      IVIB=NFRQ(III)                                                    ANH02070
      WRITE(6,3) IVIB                                                   ANH02080
      DO 102 I=1,NVIB                                                   ANH02090
      II=NFRQ(I)                                                        ANH02100
      DO 102 J=1,NVIB                                                   ANH02110
      JJ=NFRQ(J)                                                        ANH02120
      EE(I,J)=F3Q(II,JJ,IVIB)                                           ANH02130
  102 CONTINUE                                                          ANH02140
      CALL MATOUT(EE,N3N,N3N,NVIB,NVIB,6)                               ANH02150
  103 CONTINUE                                                          ANH02160
      DO 104 IVIB=1,NVIB                                                ANH02170
      WRITE(6,4) IVIB                                                   ANH02180
      CALL MATOUT(AXY(1,1,IVIB),3,3,3,3,6)                              ANH02190
  104 CONTINUE                                                          ANH02200
C                                                                       ANH02210
C   NOTE : THE SYMMETRY AXIS IS CHOSEN AS THE Z AXIS                    ANH02220
  201 CONTINUE                                                          ANH02230
      IXX=1                                                             ANH02240
      IYY=2                                                             ANH02250
      IZZ=3                                                             ANH02260
C:::FOR NON-DEGENERATE MODES:::                                         ANH02270
      DO 113 IVIB=1,NVIB                                                ANH02280
      IF(NDEG(IVIB).NE.1) GO TO 113                                     ANH02290
      FRQI=FREQ(IVIB)                                                   ANH02300
      FRQI2=FRQI*FRQI                                                   ANH02310
      VALXS1=((AXY(IXX,IXX,IVIB)**2+AXY(IXX,IYY,IVIB)**2)/ROTAA(IXX))   ANH02320
     1       *0.75D+00                                                  ANH02330
      VALYS1=((AXY(IYY,IYY,IVIB)**2+AXY(IYY,IXX,IVIB)**2)/ROTAA(IYY))   ANH02340
     1       *0.75D+00                                                  ANH02350
      VALZS1=(AXY(IZZ,IZZ,IVIB)**2/ROTAA(IZZ))*0.75D+00                 ANH02360
      VALXS2=A00                                                        ANH02370
      VALYS2=A00                                                        ANH02380
      VALZS2=A00                                                        ANH02390
      VALUS2=A00                                                        ANH02400
      VALXS3=A00                                                        ANH02410
      VALYS3=A00                                                        ANH02420
      VALZS3=A00                                                        ANH02430
C   SUM OVER NON-DEGENERATE MODES                                       ANH02440
      DO 111 JVIB=1,NVIB                                                ANH02450
      IF(NDEG(JVIB).NE.1) GO TO 111                                     ANH02460
      FRQJ=FREQ(JVIB)                                                   ANH02470
      FRQJ2=FRQJ*FRQJ                                                   ANH02480
      IF(IVIB.EQ.JVIB) GO TO 202                                        ANH02490
      FACT2=(FRQI2*3.0D+00+FRQJ2)/(FRQI2-FRQJ2)                         ANH02500
      VALZS2=VALZS2+(ZETA(IVIB,JVIB,IZZ)**2)*FACT2                      ANH02510
  202 CONTINUE                                                          ANH02520
      FACT3=(FRQI*F3Q(IVIB,IVIB,JVIB)*RCONST)/(FRQJ**(1.5D+00))         ANH02530
      VALXS3=VALXS3+AXY(IXX,IXX,JVIB)*FACT3                             ANH02540
      VALYS3=VALYS3+AXY(IYY,IYY,JVIB)*FACT3                             ANH02550
      VALZS3=VALZS3+AXY(IZZ,IZZ,JVIB)*FACT3                             ANH02560
  111 CONTINUE                                                          ANH02570
C   SUM OVER DEGENERATE MODES                                           ANH02580
      DO 112 JVIB=1,NVIB                                                ANH02590
      IF(NDEG(JVIB).EQ.1) GO TO 112                                     ANH02600
      FRQJ=FREQ(JVIB)                                                   ANH02610
      FRQJ2=FRQJ*FRQJ                                                   ANH02620
      FACT2=0.5D+00*(FRQI2*3.0D+00+FRQJ2)/(FRQI2-FRQJ2)                 ANH02630
      VALUS2=VALUS2+(ZETA(IVIB,JVIB,IXX)**2+ZETA(IVIB,JVIB,IYY)**2)     ANH02640
     1              *FACT2                                              ANH02650
  112 CONTINUE                                                          ANH02660
      VALXS2=VALUS2                                                     ANH02670
      VALYS2=VALUS2                                                     ANH02680
      ROTA=-(VALXS1+VALXS2+VALXS3)*(ROTCM(IXX)**2*TWO)/FRQI             ANH02690
      ROTB=-(VALYS1+VALYS2+VALYS3)*(ROTCM(IYY)**2*TWO)/FRQI             ANH02700
      ROTC=-(VALZS1+VALZS2+VALZS3)*(ROTCM(IZZ)**2*TWO)/FRQI             ANH02710
      RVCNST(IVIB,1)=ROTA                                               ANH02720
      RVCNST(IVIB,2)=ROTB                                               ANH02730
      RVCNST(IVIB,3)=ROTC                                               ANH02740
  113 CONTINUE                                                          ANH02750
C                                                                       ANH02760
C:::FOR DEGENERATE MODES:::                                             ANH02770
      DO 118 IVIB=1,NVIB                                                ANH02780
      IF(NDEG(IVIB).EQ.1) GO TO 118                                     ANH02790
      FRQI=FREQ(IVIB)                                                   ANH02800
      FRQI2=FRQI*FRQI                                                   ANH02810
      NII=NDEG(IVIB)                                                    ANH02820
      FTXYZZ=(AXY(IXX,IZZ,IVIB)**2+AXY(IYY,IZZ,IVIB)**2)                ANH02830
      FTXXYY=(AXY(IXX,IXX,IVIB)**2+AXY(IYY,IYY,IVIB)**2                 ANH02840
     1       +AXY(IXX,IYY,IVIB)**2+AXY(IYY,IXX,IVIB)**2)                ANH02850
      VALXT1=(0.375D+00*FTXYZZ)/ROTAA(IZZ)                              ANH02860
     1      +(0.750D+00*FTXXYY)/(ROTAA(IXX)+ROTAA(IYY))                 ANH02870
      VALYT1=(0.375D+00*FTXYZZ)/ROTAA(IZZ)                              ANH02880
     1      +(0.750D+00*FTXXYY)/(ROTAA(IXX)+ROTAA(IYY))                 ANH02890
      VALZT1=(1.500D+00*FTXYZZ)/(ROTAA(IXX)+ROTAA(IYY))                 ANH02900
      VALXT2=A00                                                        ANH02910
      VALYT2=A00                                                        ANH02920
      VALZT2=A00                                                        ANH02930
      VALUT2=A00                                                        ANH02940
      VALXT3=A00                                                        ANH02950
      VALYT3=A00                                                        ANH02960
      VALZT3=A00                                                        ANH02970
      VALUT3=A00                                                        ANH02980
      VALXT4=A00                                                        ANH02990
      VALYT4=A00                                                        ANH03000
      VALZT4=A00                                                        ANH03010
C   SUM OVER NON-DEGENERATE MODES                                       ANH03020
      DO 115 JVIB=1,NVIB                                                ANH03030
      IF(NDEG(JVIB).NE.1) GO TO 115                                     ANH03040
      FRQJ=FREQ(JVIB)                                                   ANH03050
      FRQJ2=FRQJ*FRQJ                                                   ANH03060
      FACT2=0.250D+00*(FRQI2*3.0D+00+FRQJ2)/(FRQI2-FRQJ2)               ANH03070
      DO 114 KVIB=1,NII                                                 ANH03080
      LVIB=NDAB(IVIB,KVIB)                                              ANH03090
      VALUT2=VALUT2+(ZETA(LVIB,JVIB,IXX)**2+ZETA(LVIB,JVIB,IYY)**2)     ANH03100
     1             *FACT2                                               ANH03110
  114 CONTINUE                                                          ANH03120
      FACT4=(FRQI*F3Q(IVIB,IVIB,JVIB)*RCONST)/(FRQJ**(1.5D+00))         ANH03130
      VALXT4=VALXT4+AXY(IXX,IXX,JVIB)*FACT4                             ANH03140
      VALYT4=VALYT4+AXY(IYY,IYY,JVIB)*FACT4                             ANH03150
      VALZT4=VALZT4+AXY(IZZ,IZZ,JVIB)*FACT4                             ANH03160
  115 CONTINUE                                                          ANH03170
      VALXT2=VALUT2                                                     ANH03180
      VALYT2=VALUT2                                                     ANH03190
C   SUM OVER DEGENERATE MODES                                           ANH03200
      DO 116 JVIB=1,NVIB                                                ANH03210
      IF(NDEG(JVIB).EQ.1) GO TO 116                                     ANH03220
      FRQJ=FREQ(JVIB)                                                   ANH03230
      FRQJ2=FRQJ*FRQJ                                                   ANH03240
      IF(IVIB.EQ.JVIB) GO TO 116                                        ANH03250
      IF(DABS(FRQI-FRQJ).LE.WLIMIT) GO TO 116                           ANH03260
      FACT3=(FRQI2*3.0D+00+FRQJ2)/(FRQI2-FRQJ2)                         ANH03270
      VALUT3=VALUT3+(ZETA(IVIB,JVIB,IXX)**2+ZETA(IVIB,JVIB,IYY)**2)     ANH03280
     1       *FACT3                                                     ANH03290
      VALZT3=VALZT3+(ZETA(IVIB,JVIB,IZZ)**2)*FACT3                      ANH03300
  116 CONTINUE                                                          ANH03310
      VALXT3=VALUT3*0.5D+00                                             ANH03320
      VALYT3=VALUT3*0.5D+00                                             ANH03330
CCC   VALXT3=VALUT3*0.25D+00                                            ANH03340
CCC   VALYT3=VALUT3*0.25D+00                                            ANH03350
      ROTA=-(VALXT1+VALXT2+VALXT3+VALXT4)*(ROTCM(IXX)**2*TWO)/FRQI      ANH03360
      ROTB=-(VALYT1+VALYT2+VALYT3+VALYT4)*(ROTCM(IYY)**2*TWO)/FRQI      ANH03370
      ROTC=-(VALZT1+VALZT2+VALZT3+VALZT4)*(ROTCM(IZZ)**2*TWO)/FRQI      ANH03380
      RVCNST(IVIB,1)=ROTA                                               ANH03390
      RVCNST(IVIB,2)=ROTB                                               ANH03400
      RVCNST(IVIB,3)=ROTC                                               ANH03410
C                                                                       ANH03420
CCC   FACTX=(ROT(IXX)**2*TWO)/FRQI                                      ANH03430
CCC   FACTY=(ROT(IYY)**2*TWO)/FRQI                                      ANH03440
CCC   FACTZ=(ROT(IZZ)**2*TWO)/FRQI                                      ANH03450
CCC   VALXT1=VALXT1*FACTX                                               ANH03460
CCC   VALXT2=VALXT2*FACTX                                               ANH03470
CCC   VALXT3=VALXT3*FACTX                                               ANH03480
CCC   VALXT4=VALXT4*FACTX                                               ANH03490
CCC   VALYT1=VALYT1*FACTY                                               ANH03500
CCC   VALYT2=VALYT2*FACTY                                               ANH03510
CCC   VALYT3=VALYT3*FACTY                                               ANH03520
CCC   VALYT4=VALYT4*FACTY                                               ANH03530
CCC   VALZT1=VALZT1*FACTZ                                               ANH03540
CCC   VALZT2=VALZT2*FACTZ                                               ANH03550
CCC   VALZT3=VALZT3*FACTZ                                               ANH03560
CCC   VALZT4=VALZT4*FACTZ                                               ANH03570
CCC   WRITE(6,*) ' IVIB ',IVIB                                          ANH03580
CCC   WRITE(6,*) ' VALXT1,VALXT2,VALXT3,VALXT4,FACTX,ROTA'              ANH03590
CCC   WRITE(6,*)   VALXT1,VALXT2,VALXT3,VALXT4,FACTX,ROTA               ANH03600
CCC   WRITE(6,*) ' VALYT1,VALYT2,VALYT3,VALYT4,FACTY,ROTB'              ANH03610
CCC   WRITE(6,*)   VALYT1,VALYT2,VALYT3,VALYT4,FACTY,ROTB               ANH03620
CCC   WRITE(6,*) ' VALZT1,VALZT2,VALZT3,VALZT4,FACTZ,ROTC'              ANH03630
CCC   WRITE(6,*)   VALZT1,VALZT2,VALZT3,VALZT4,FACTZ,ROTC               ANH03640
  118 CONTINUE                                                          ANH03650
C                                                                       ANH03660
C   CALCULATE THE ROTATIONAL CONSTANTS                                  ANH03670
C   WITH ZERO-POINT VIBRATIONAL CORRECTION                              ANH03680
      CALL ZERO(SALP,3)                                                 ANH03690
      WRITE(6,5)                                                        ANH03700
      DO 120 III=1,NVIB                                                 ANH03710
      IVIB=NFRQ(III)                                                    ANH03720
      ROTA=RVCNST(IVIB,1)                                               ANH03730
      ROTB=RVCNST(IVIB,2)                                               ANH03740
      ROTC=RVCNST(IVIB,3)                                               ANH03750
      SALP(1)=SALP(1)+ROTA                                              ANH03760
      SALP(2)=SALP(2)+ROTB                                              ANH03770
      SALP(3)=SALP(3)+ROTC                                              ANH03780
      WRITE(6,6) III,IVIB,ROTA,ROTB,ROTC                                ANH03790
  120 CONTINUE                                                          ANH03800
      WRITE(6,7)                                                        ANH03810
      DO 122 I=1,3                                                      ANH03820
      PB=ROTCM(I)-SALP(I)*HALF                                          ANH03830
      PT=CONST/PB                                                       ANH03840
      PA=PT*PARA                                                        ANH03850
      PC=CYCL/PT                                                        ANH03860
      WRITE(6,8) RAXIS(I),PT,PA,PB,PC                                   ANH03870
  122 CONTINUE                                                          ANH03880
C                                                                       ANH03890
      IF(NCORIO.NE.0) THEN                                              ANH03900
        WRITE(6,9)                                                      ANH03910
        WRITE(6,10)                                                     ANH03920
      END IF                                                            ANH03930
C                                                                       ANH03940
      RETURN                                                            ANH03950
      END                                                               ANH03960
      SUBROUTINE CENTSP(F3Q,EE)                                         ANH03970
C   THE CALCULATION OF CENTRIFUGAL DISTORTION CONSTANTS                 ANH03980
      IMPLICIT REAL*8 (A-H,O-Z)                                         ANH03990
      DIMENSION F3Q(N3N,N3N,N3N),EE(N3N,N3N)                            ANH04000
      COMMON/VIB101/NATOM,N3N,NATRI,ILIN,NVIB                           ANH04010
      COMMON/VIB102/ITHREE,IFOUR,N3TOT,N4TOT                            ANH04020
      COMMON/VIB106/ROTAA(3),ROTGC(3),ROTCM(3),ROTMH(3)                 ANH04030
      COMMON/VIB108/IAXIS(3),NDEG(150),NDAB(150,5),IMAG(150)            ANH04040
      COMMON/VIB203/IOFF(150),IPRNT                                     ANH04050
      COMMON/VIB204/SQM(150),ROOT(150),FREQ(150)                        ANH04060
      COMMON/VIB205/IFREQ,NFRQ(150)                                     ANH04070
      COMMON/VIB207/TABCD(3,3,3,3),PABC(3,3,3)                          ANH04080
      COMMON/VIB208/TAABB(3,3),TABAB(3,3)                               ANH04090
      COMMON/VIB209/SALP(3)                                             ANH04100
      DATA WLIMIT / 1.0D-2 /                                            ANH04110
      DATA PH,AVN / 6.626176D+00 , 6.022045D+00 /                       ANH04120
      DATA PI,CL / 3.1415926536D+00 , 2.99792458D+00 /                  ANH04130
      DATA A00,HALF,ONE,TWO / 0.0D+00 , 0.5D+00 , 1.0D+00 , 2.0D+00 /   ANH04140
    1 FORMAT(//,2X,' ******************************************'/       ANH04150
     1          2X,' ***THE CENTRIFUGAL DISTORTION CONSTANTS***'/       ANH04160
     2          2X,' ******************************************')       ANH04170
    2 FORMAT(//,2X,' F3(QI*QJ*QK) MATRIX (IN CM-1), IVIB = ',I5/)       ANH04180
    3 FORMAT(//,2X,' TAABB MATRIX FROM BXY MATRICES (IN CM-1)'/)        ANH04190
    4 FORMAT(//,2X,' TABAB MATRIX FROM BXY MATRICES (IN CM-1)'/)        ANH04200
    5 FORMAT(//,2X,' FOR A SPHERICAL TOP MOLECULE'/                     ANH04210
     1          2X,' DJ   = ',F15.5,' D-06 CM-1',                       ANH04220
     2         10X,' DJ   = ',F15.5,' MHZ'/                             ANH04230
     3          2X,' DK   = ',F15.5,' D-06 CM-1',                       ANH04240
     4         10X,' DK   = ',F15.5,' MHZ'/                             ANH04250
     5          2X,' DJK  = ',F15.5,' D-06 CM-1',                       ANH04260
     6         10X,' DJK  = ',F15.5,' MHZ'/)                            ANH04270
    6 FORMAT(//,2X,' FOR A SPHERICAL TOP MOLECULE'/                     ANH04280
     1          2X,' H    = ',F15.5,' D-06 CM-1',                       ANH04290
     2         10X,' H    = ',F15.5,' MHZ'/                             ANH04300
     3          2X,' H4T  = ',F15.5,' D-06 CM-1',                       ANH04310
     4         10X,' H4T  = ',F15.5,' MHZ'/                             ANH04320
     5          2X,' H6T  = ',F15.5,' D-06 CM-1',                       ANH04330
     6         10X,' H6T  = ',F15.5,' MHZ'/)                            ANH04340
C                                                                       ANH04350
      WRITE(6,1)                                                        ANH04360
C                                                                       ANH04370
C   ELEMENTS OF F3Q ARE IN CM-1                                         ANH04380
      IF(IPRNT.LE.2) GO TO 201                                          ANH04390
      DO 102 III=1,NVIB                                                 ANH04400
      IVIB=NFRQ(III)                                                    ANH04410
      WRITE(6,2) IVIB                                                   ANH04420
      DO 101 I=1,NVIB                                                   ANH04430
      II=NFRQ(I)                                                        ANH04440
      DO 101 J=1,NVIB                                                   ANH04450
      JJ=NFRQ(J)                                                        ANH04460
      EE(I,J)=F3Q(II,JJ,IVIB)                                           ANH04470
  101 CONTINUE                                                          ANH04480
      CALL MATOUT(EE,N3N,N3N,NVIB,NVIB,6)                               ANH04490
  102 CONTINUE                                                          ANH04500
C                                                                       ANH04510
  201 CONTINUE                                                          ANH04520
      DO 103 IA=1,3                                                     ANH04530
      DO 103 IB=1,3                                                     ANH04540
      TAABB(IA,IB)=TABCD(IA,IA,IB,IB)                                   ANH04550
      TABAB(IA,IB)=TABCD(IA,IB,IA,IB)                                   ANH04560
  103 CONTINUE                                                          ANH04570
      IF(IPRNT.LE.2) GO TO 202                                          ANH04580
      WRITE(6,3)                                                        ANH04590
      CALL MATOUT(TAABB,3,3,3,3,6)                                      ANH04600
      WRITE(6,4)                                                        ANH04610
      CALL MATOUT(TABAB,3,3,3,3,6)                                      ANH04620
C                                                                       ANH04630
C   FOR A SYMMETRIC TOP MOLECULE                                        ANH04640
  202 CONTINUE                                                          ANH04650
      IXX=1                                                             ANH04660
      IYY=2                                                             ANH04670
      IZZ=3                                                             ANH04680
      DJ=-(TAABB(IXX,IXX)*3.0D+00+TAABB(IYY,IYY)*3.0D+00                ANH04690
     1    +TAABB(IXX,IYY)*2.0D+00+TABAB(IXX,IYY)*4.0D+00)/32.0D+00      ANH04700
      DK=DJ-(TAABB(IZZ,IZZ)-TAABB(IXX,IZZ)-TAABB(IYY,IZZ)               ANH04710
     1      -TABAB(IXX,IZZ)*2.0D+00-TABAB(IYY,IZZ)*2.0D+00)/4.0D+00     ANH04720
      DJK=-DJ-DK-TAABB(IZZ,IZZ)/4.0D+00                                 ANH04730
      DJCM=DJ*1.0D+06                                                   ANH04740
      DKCM=DK*1.0D+06                                                   ANH04750
      DJKCM=DJK*1.0D+06                                                 ANH04760
      DJMH=DJ*CL*1.0D+04                                                ANH04770
      DKMH=DK*CL*1.0D+04                                                ANH04780
      DJKMH=DJK*CL*1.0D+04                                              ANH04790
      WRITE(6,5) DJCM,DJMH,DKCM,DKMH,DJKCM,DJKMH                        ANH04800
C                                                                       ANH04810
C   HIGHER ORDER CENTRIFUGAL DISTORTION CONSTANTS                       ANH04820
      H=(PABC(IXX,IXX,IXX)*45.0D+00+PABC(IXX,IXX,IYY)*36.0D+00          ANH04830
     1  +PABC(IXX,IYY,IZZ)*2.0D+00)/105.0D+00                           ANH04840
      H4T=(-PABC(IXX,IXX,IXX)*15.0D+00+PABC(IXX,IXX,IYY)*8.0D+00        ANH04850
     1     +PABC(IXX,IYY,IZZ))/110.0D+00                                ANH04860
      H6T=(PABC(IXX,IXX,IXX)**6.0D+00-PABC(IXX,IXX,IYY)*12.0D+00        ANH04870
     1    +PABC(IXX,IYY,IZZ))/231.0D+00                                 ANH04880
      HCM=H*1.0D+06                                                     ANH04890
      H4TCM=H4T*1.0D+06                                                 ANH04900
      H6TCM=H6T*1.0D+06                                                 ANH04910
      HMH=H*CL*1.0D+04                                                  ANH04920
      H4TMH=H46*CL*1.0D+04                                              ANH04930
      H6TMH=H6T*CL*1.0D+04                                              ANH04940
      WRITE(6,6) HCM,HMH,H4TCM,H4TMH,H6TCM,H6TMH                        ANH04950
C                                                                       ANH04960
      RETURN                                                            ANH04970
      END                                                               ANH04980
      SUBROUTINE VLTYSP(F3Q,F4Q,CLTYP)                                  ANH04990
C   THE CALCULATION OF VIBRATIONAL L-TYPE DOUBLING CONSTANTS            ANH05000
C   FOR A SPHERICAL TOP MOLECULE                                        ANH05010
      IMPLICIT REAL*8 (A-H,O-Z)                                         ANH05020
      DIMENSION CLTYP(N3N,N3N)                                          ANH05030
      DIMENSION F3Q(N3N,N3N,N3N),F4Q(N3N,N3N,N3N,N3N)                   ANH05040
      COMMON/VIB101/NATOM,N3N,NATRI,ILIN,NVIB                           ANH05050
      COMMON/VIB108/IAXIS(3),NDEG(150),NDAB(150,5),IMAG(150)            ANH05060
      COMMON/VIB109/CLIMIT,FLIM1,FLIM2                                  ANH05070
      COMMON/VIB203/IOFF(150),IPRNT                                     ANH05080
      COMMON/VIB204/SQM(150),ROOT(150),FREQ(150)                        ANH05090
      COMMON/VIB205/IFREQ,NFRQ(150)                                     ANH05100
      DATA A00,HALF,ONE,TWO / 0.0D+00 , 0.5D+00 , 1.0D+00 , 2.0D+00 /   ANH05110
      DATA WLIMIT / 1.0D-02 /                                           ANH05120
      DATA PI,CL / 3.1415926536D+00 , 2.99792458D+00 /                  ANH05130
    1 FORMAT(//,2X,' ***********************************************'/  ANH05140
     1          2X,' ***THE VIBRATIONAL L-TYPE DOUBLING CONSTANTS***'/  ANH05150
     2          2X,' ***********************************************')  ANH05160
    2 FORMAT(2X,I5,4F10.5)                                              ANH05170
    3 FORMAT(//,2X,' QVIB',5X,' IVIB',7X,' FREQENCY',9X,' CUBIC',       ANH05180
     1         12X,' QUARTIC',10X,' TOTAL',9X,' TOTAL'/                 ANH05190
     2         24X,' IN CM-1',10X,' IN CM-1',10X,' IN CM-1',            ANH05200
     3         10X,' IN CM-1',7X,' IN MHZ'/)                            ANH05210
    4 FORMAT(2X,I5,5X,I5,8X,F10.5,F15.5,3X,F15.5,3X,F15.5,3X,F15.5)     ANH05220
C                                                                       ANH05230
      WRITE(6,1)                                                        ANH05240
C                                                                       ANH05250
C   NOTE : THE MOLECULAR AXIS IS CHOSEN AS THE Z AXIS                   ANH05260
      IXX=1                                                             ANH05270
      IYY=2                                                             ANH05280
      IZZ=3                                                             ANH05290
C                                                                       ANH05300
C:::::::::::::::::::::::::::::::                                        ANH05310
C:::ONLY FOR DEGENERATE MODES:::                                        ANH05320
C:::::::::::::::::::::::::::::::                                        ANH05330
      DO 105 IVIB=1,NVIB                                                ANH05340
      IF(NDEG(IVIB).EQ.1) GO TO 105                                     ANH05350
      FRQI=FREQ(IVIB)                                                   ANH05360
      FRQI2=FRQI*FRQI                                                   ANH05370
      DO 104 JVIB=1,NVIB                                                ANH05380
      IF(NDEG(JVIB).EQ.1) GO TO 104                                     ANH05390
      FRQJ=FREQ(JVIB)                                                   ANH05400
      FRQJ2=FRQJ*FRQJ                                                   ANH05410
      IF(IVIB.EQ.JVIB) GO TO 104                                        ANH05420
C                                                                       ANH05430
C   SUM OVER NON-DEGENERATE MODES                                       ANH05440
      VALU3=A00                                                         ANH05450
      DO 103 KVIB=1,NVIB                                                ANH05460
      IF(NDEG(KVIB).NE.1) GO TO 103                                     ANH05470
      FRQK=FREQ(JVIB)                                                   ANH05480
      FRQK2=FRQK*FRQK                                                   ANH05490
      DNUM=(FRQK2-FRQI2-FRQJ2)*FRQK                                     ANH05500
      FAC=F3Q(IVIB,JVIB,KVIB)                                           ANH05510
      FACT=FAC*FAC                                                      ANH05520
      WPPP=FRQI+FRQJ+FRQK                                               ANH05530
      WPPM=FRQI+FRQJ-FRQK                                               ANH05540
      WPMP=FRQI-FRQJ+FRQK                                               ANH05550
      WPMM=FRQI-FRQJ-FRQK                                               ANH05560
      DIJK=WPPP*WPPM*WPMP*WPMM                                          ANH05570
      VALU3=VALU3-(0.5D+00*FACT*DNUM)/(DIJK)                            ANH05580
  103 CONTINUE                                                          ANH05590
      VALU4=(F4Q(IVIB,IVIB,JVIB,JVIB)-F4Q(IVIB,IVIB,JVIB,JVIB))*0.5D+00 ANH05600
      VALUT=VALU3+VALU4                                                 ANH05610
      CLTYP(IVIB,1)=VALU3                                               ANH05620
      CLTYP(IVIB,2)=VALU4                                               ANH05630
      CLTYP(IVIB,3)=VALUT                                               ANH05640
      IF(IPRNT.LE.2) GO TO 105                                          ANH05650
      WRITE(6,2) IVIB,FRQI,VALU3,VALU4,VALUT                            ANH05660
  104 CONTINUE                                                          ANH05670
  105 CONTINUE                                                          ANH05680
C                                                                       ANH05690
      WRITE(6,3)                                                        ANH05700
      DO 106 II=1,NVIB                                                  ANH05710
      IVIB=NFRQ(II)                                                     ANH05720
      IF(NDEG(IVIB).LE.1) GO TO 106                                     ANH05730
      FRQI=FREQ(IVIB)                                                   ANH05740
      VALU3=CLTYP(IVIB,1)                                               ANH05750
      VALU4=CLTYP(IVIB,2)                                               ANH05760
      VALUT=CLTYP(IVIB,3)                                               ANH05770
      VALUM=VALUT*CL*1.0D+04                                            ANH05780
      WRITE(6,4) II,IVIB,FRQI,VALU3,VALU4,VALUT,VALUM                   ANH05790
  106 CONTINUE                                                          ANH05800
C                                                                       ANH05810
      RETURN                                                            ANH05820
      END                                                               ANH05830
      SUBROUTINE SPHTOP(ZETA,F3Q,F4Q,XIJ,GIJ)                           ANH05840
C   THE ANHARMONICITY CONTSTANTS FOR A SYMMETRIC TOP                    ANH05850
      IMPLICIT REAL*8 (A-H,O-Z)                                         ANH05860
      DIMENSION ZETA(N3N,N3N,3),F3Q(N3N,N3N,N3N),F4Q(N3N,N3N,N3N,N3N)   ANH05870
      DIMENSION XIJ(N3N,N3N),GIJ(N3N,N3N)                               ANH05880
      COMMON/VIB101/NATOM,N3N,NATRI,ILIN,NVIB                           ANH05890
      COMMON/VIB106/ROTAA(3),ROTGC(3),ROTCM(3),ROTMH(3)                 ANH05900
      COMMON/VIB108/IAXIS(3),NDEG(150),NDAB(150,5),IMAG(150)            ANH05910
      COMMON/VIB109/CLIMIT,FLIM1,FLIM2                                  ANH05920
      COMMON/VIB203/IOFF(150),IPRNT                                     ANH05930
      COMMON/VIB204/SQM(150),ROOT(150),FREQ(150)                        ANH05940
      COMMON/VIB205/IFREQ,NFRQ(150)                                     ANH05950
      DATA A00,ONE / 0.0D+00 , 1.0D+00 /                                ANH05960
      DATA WLIMIT / 3.0D+00 /                                           ANH05970
    1 FORMAT(//,2X,' ANHARMONIC CONSTANTS IN CM-1'/                     ANH05980
     1          2X,' IVIB',3X,' JVIB',8X,' FREQ(I)',5X,' FREQ(J)',      ANH05990
     2          7X,' VALU1',7X,' VALU2',7X,' VALU3',7X,' VALU4',        ANH06000
     3          7X,' XIJ',9X,' GIJ'/)                                   ANH06010
    2 FORMAT(2X,I5,3X,I5,3X,8(F13.3))                                   ANH06020
    3 FORMAT(//,2X,' ANHARMONIC CONSTANTS IN CM-1'/                     ANH06030
     1          4X,' QI',5X,' QJ',8X,' FREQ(I)',5X,' FREQ(J)',          ANH06040
     2          7X,' XIJ',9X,' GIJ'/)                                   ANH06050
    4 FORMAT(2X,I5,3X,I5,3X,4(F13.3))                                   ANH06060
    5 FORMAT(//,2X,' THEORETICAL FREQUENCIES IN CM-1'/                  ANH06070
     1          2X,' NO.',9X,' HARMONIC',14X,' DELTA',                  ANH06080
     2         13X,' FUNDAMENTAL'/)                                     ANH06090
    6 FORMAT(2X,I5,4X,F13.3,7X,F13.3,12X,F13.3)                         ANH06100
    7 FORMAT(//,2X,' :::::::::::::::::::::::::::::::::::::::::::::::'/  ANH06110
     1          2X,' :::ANHARMONIC RESONANCE (TYPE I) IS OBSERVED:::'/  ANH06120
     2          2X,' :::::::::::::::::::::::::::::::::::::::::::::::'/) ANH06130
    8 FORMAT(//,2X,' ::::::::::::::::::::::::::::::::::::::::::::::::'/ ANH06140
     1          2X,' :::ANHARMONIC RESONANCE (TYPE II) IS OBSERVED:::'/ ANH06150
     2          2X,' ::::::::::::::::::::::::::::::::::::::::::::::::'/)ANH06160
    9 FORMAT(2X,' RESONANCE LIMIT IS SET TO ',F10.3,' CM-1'/)           ANH06170
C                                                                       ANH06180
C:::::::::::::::::::::::::::::::::::::::                                ANH06190
C:::CALCULATE ANHARMONICITY CONSTANTS:::                                ANH06200
C:::::::::::::::::::::::::::::::::::::::                                ANH06210
      IANH1=0                                                           ANH06220
      IANH2=0                                                           ANH06230
C                                                                       ANH06240
      IPRNT=5                                                           ANH06250
      IF(IPRNT.GT.2)                                                    ANH06260
     *WRITE(6,1)                                                        ANH06270
      IAA=IAXIS(1)                                                      ANH06280
      IBB=IAXIS(2)                                                      ANH06290
      ICC=IAXIS(3)                                                      ANH06300
C                                                                       ANH06310
      DO 120 IVIB=1,N3N                                                 ANH06320
      FRQI=FREQ(IVIB)                                                   ANH06330
      FRQI2=FRQI*FRQI                                                   ANH06340
      DO 120 JVIB=1,N3N                                                 ANH06350
      FRQJ=FREQ(JVIB)                                                   ANH06360
      FRQJ2=FRQJ*FRQJ                                                   ANH06370
C                                                                       ANH06380
      VALU1=A00                                                         ANH06390
      VALU2=A00                                                         ANH06400
      VALU3=A00                                                         ANH06410
      VALU4=A00                                                         ANH06420
      VALU5=A00                                                         ANH06430
      GVAL1=A00                                                         ANH06440
      GVAL2=A00                                                         ANH06450
      GVAL3=A00                                                         ANH06460
      GVAL4=A00                                                         ANH06470
      XRS=A00                                                           ANH06480
      GRS=A00                                                           ANH06490
      IF(DABS(FRQI).LE.WLIMIT) GO TO 230                                ANH06500
      IF(DABS(FRQJ).LE.WLIMIT) GO TO 230                                ANH06510
      IF(IVIB.NE.JVIB) GO TO 210                                        ANH06520
C                                                                       ANH06530
C+++++++++++++++++++++++++++++++++++++++++                              ANH06540
C+++DIAGONAL CONSTANTS (NON-DEGENERATE)+++                              ANH06550
C+++++++++++++++++++++++++++++++++++++++++                              ANH06560
      IF(NDEG(IVIB).NE.1) GO TO 205                                     ANH06570
      VALU1=F4Q(IVIB,IVIB,IVIB,IVIB)                                    ANH06580
      DO 101 KVIB=1,N3N                                                 ANH06590
C   SUM OVER ONLY NON-DEGENERATE MODES                                  ANH06600
      FRQK=FREQ(KVIB)                                                   ANH06610
      IF(DABS(FRQK).LE.WLIMIT) GO TO 101                                ANH06620
      IF(NDEG(KVIB).NE.1) GO TO 101                                     ANH06630
      FRQK2=FRQK*FRQK                                                   ANH06640
      FAC=F3Q(IVIB,IVIB,KVIB)                                           ANH06650
      FACT=FAC*FAC                                                      ANH06660
      IF(DABS(FRQI+FRQI-FRQK).LE.FLIM1) GO TO 201                       ANH06670
      DNUM=FRQI2*8.0D+00-FRQK2*3.0D+00                                  ANH06680
      DENM=FRQI2*4.0D+00-FRQK2                                          ANH06690
      VALU2=VALU2+(DNUM*FACT)/(FRQK*DENM)                               ANH06700
      GO TO 101                                                         ANH06710
C   FERMI RESONANCE CORRECTIONS                                         ANH06720
  201 CONTINUE                                                          ANH06730
      IANH1=IANH1+1                                                     ANH06740
      DNUM=ONE/(FRQK+FRQK)+ONE/((FRQI+FRQI+FRQK)*8.0D+00)               ANH06750
      VALU2=VALU2+4.0D+00*FACT*DNUM                                     ANH06760
  101 CONTINUE                                                          ANH06770
      XRS=(VALU1-VALU2)/16.0D+00                                        ANH06780
      GO TO 230                                                         ANH06790
C                                                                       ANH06800
C+++++++++++++++++++++++++++++++++++++                                  ANH06810
C+++DIAGONAL CONSTANTS (DEGENERATE)+++                                  ANH06820
C+++++++++++++++++++++++++++++++++++++                                  ANH06830
  205 CONTINUE                                                          ANH06840
      VALU1=F4Q(IVIB,IVIB,IVIB,IVIB)                                    ANH06850
      GVAL1=-VALU1                                                      ANH06860
      DO 102 KVIB=1,N3N                                                 ANH06870
C   SUM OVER NON-DEGENERATE MODES                                       ANH06880
      FRQK=FREQ(KVIB)                                                   ANH06890
      IF(DABS(FRQK).LE.WLIMIT) GO TO 102                                ANH06900
      FRQK2=FRQK*FRQK                                                   ANH06910
      FAC=F3Q(IVIB,IVIB,KVIB)                                           ANH06920
      FACT=FAC*FAC                                                      ANH06930
      IF(NDEG(KVIB).NE.1) GO TO 207                                     ANH06940
      IF(DABS(FRQI+FRQI-FRQK).LE.FLIM1) GO TO 206                       ANH06950
      DNUM=FRQI2*8.0D+00-FRQK2*3.0D+00                                  ANH06960
      DENM=FRQI2*4.0D+00-FRQK2                                          ANH06970
      VALU2=VALU2+(DNUM*FACT)/(FRQK*DENM)                               ANH06980
      GVAL2=GVAL2+(FACT*FRQK)/DENM                                      ANH06990
C     GO TO 102                                                         ANH07000
      GO TO 2222                                                        ANH07010
C   FERMI RESONANCE CORRECTIONS                                         ANH07020
  206 CONTINUE                                                          ANH07030
      IANH1=IANH1+1                                                     ANH07040
      DNUM=ONE/(FRQK+FRQK)+ONE/((FRQI+FRQI+FRQK)*8.0D+00)               ANH07050
      VALU2=VALU2+4.0D+00*FACT*DNUM                                     ANH07060
C     GO TO 102                                                         ANH07070
      GO TO 2222                                                        ANH07080
C   SUM OVER DEGENERATE MODES                                           ANH07090
  207 CONTINUE                                                          ANH07100
      DNUM=FRQI2*8.0D+00-FRQK2*3.0D+00                                  ANH07110
      DENM=FRQI2*4.0D+00-FRQK2                                          ANH07120
      DFAC=DBLE(NDEG(KVIB))                                             ANH07130
      VALU3=VALU3+(DNUM*FACT)/(FRQK*DENM*DFAC)                          ANH07140
      GVAL3=GVAL3+(DNUM*FACT)/(FRQK*DENM*DFAC)                          ANH07150
      GVAL4=ZETA(IVIB,IVIB,ICC)**2*ROTCM(IAA)                           ANH07160
 2222 CONTINUE                                                          ANH07170
C     WRITE(6,*) ' IVIB,JVIB,KVIB'                                      ANH07180
C     WRITE(6,*)   IVIB,JVIB,KVIB                                       ANH07190
C     WRITE(6,*) ' DNUM,DENM,DFAC'                                      ANH07200
C     WRITE(6,*)   DNUM,DENM,DFAC                                       ANH07210
C     WRITE(6,*) ' VALU1,VALU2,VALU3'                                   ANH07220
C     WRITE(6,*)   VALU1,VALU2,VALU3                                    ANH07230
  102 CONTINUE                                                          ANH07240
      XRS=(VALU1-VALU2-VALU3)/16.0D+00                                  ANH07250
      GRS=(GVAL1/3.0D+00-GVAL2+GVAL3)/16.0D+00+GVAL4                    ANH07260
      GO TO 230                                                         ANH07270
C                                                                       ANH07280
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++         ANH07290
C+++OFF-DIAGONAL CONSTANTS (NON-DEGENERATE & NON-DEGENERATE)+++         ANH07300
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++         ANH07310
  210 CONTINUE                                                          ANH07320
      IF(NDEG(IVIB).NE.1) GO TO 215                                     ANH07330
      IF(NDEG(JVIB).NE.1) GO TO 215                                     ANH07340
      VALU1=F4Q(IVIB,IVIB,JVIB,JVIB)                                    ANH07350
      DO 103 KVIB=1,NVIB                                                ANH07360
C   SUM OVER ONLY NON-DEGENERATE MODES                                  ANH07370
      FRQK=FREQ(KVIB)                                                   ANH07380
      IF(DABS(FRQK).LE.WLIMIT) GO TO 103                                ANH07390
      IF(NDEG(KVIB).NE.1) GO TO 103                                     ANH07400
      FRQK2=FRQK*FRQK                                                   ANH07410
      WPPP= FRQI+FRQJ+FRQK                                              ANH07420
      WPMM= FRQI-FRQJ-FRQK                                              ANH07430
      WMPM=-FRQI+FRQJ-FRQK                                              ANH07440
      WMMP=-FRQI-FRQJ+FRQK                                              ANH07450
      DIJK=WPPP*WPMM*WMPM*WMMP                                          ANH07460
      DNUM2=F3Q(IVIB,IVIB,KVIB)*F3Q(KVIB,JVIB,JVIB)                     ANH07470
      DNUM3=(FRQK2-FRQI2-FRQJ2)*FRQK                                    ANH07480
      FAC=F3Q(IVIB,JVIB,KVIB)                                           ANH07490
      FACT=FAC*FAC                                                      ANH07500
      VALU2=VALU2+DNUM2/FRQK                                            ANH07510
      IF(DABS(WPMM).LE.FLIM2) GO TO 211                                 ANH07520
      IF(DABS(WMPM).LE.FLIM2) GO TO 212                                 ANH07530
      IF(DABS(WMMP).LE.FLIM2) GO TO 213                                 ANH07540
      VALU3=VALU3+(DNUM3*FACT)/DIJK                                     ANH07550
      GO TO 103                                                         ANH07560
C   FERMI RESONANCE CORRECTIONS                                         ANH07570
  211 CONTINUE                                                          ANH07580
      IANH2=IANH2+1                                                     ANH07590
      DNUM3=ONE/WPPP-ONE/WMMP-ONE/WMPM                                  ANH07600
      VALU3=VALU3-(DNUM3*FACT)/4.0D+00                                  ANH07610
      GO TO 103                                                         ANH07620
  212 CONTINUE                                                          ANH07630
      IANH2=IANH2+1                                                     ANH07640
      DNUM3=ONE/WPPP-ONE/WMMP-ONE/WPMM                                  ANH07650
      VALU3=VALU3-(DNUM3*FACT)/4.0D+00                                  ANH07660
      GO TO 103                                                         ANH07670
  213 CONTINUE                                                          ANH07680
      IANH2=IANH2+1                                                     ANH07690
      DNUM3=ONE/WPPP-ONE/WMPM-ONE/WPMM                                  ANH07700
      VALU3=VALU3-(DNUM3*FACT)/4.0D+00                                  ANH07710
  103 CONTINUE                                                          ANH07720
      ZFACT=ZETA(IVIB,JVIB,ICC)**2*ROTCM(IAA)                           ANH07730
      VALU4=(FRQI/FRQJ+FRQJ/FRQI)*ZFACT                                 ANH07740
      XRS=(VALU1-VALU2-VALU3*2.0D+00)/4.0D+00+VALU4                     ANH07750
      GO TO 230                                                         ANH07760
C                                                                       ANH07770
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++             ANH07780
C+++OFF-DIAGONAL CONSTANTS (DEGENERATE & NON-DEGENERATE)+++             ANH07790
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++             ANH07800
  215 CONTINUE                                                          ANH07810
      IF(NDEG(IVIB).NE.1.AND.NDEG(JVIB).NE.1) GO TO 225                 ANH07820
      VALU1=F4Q(IVIB,IVIB,JVIB,JVIB)                                    ANH07830
      DO 105 KVIB=1,NVIB                                                ANH07840
C   SUM OVER NON-DEGENERATE MODES                                       ANH07850
      FRQK=FREQ(KVIB)                                                   ANH07860
      IF(DABS(FRQK).LE.WLIMIT) GO TO 105                                ANH07870
      FRQK2=FRQK*FRQK                                                   ANH07880
      IF(NDEG(KVIB).NE.1) GO TO 216                                     ANH07890
      DNUM2=F3Q(IVIB,IVIB,KVIB)*F3Q(KVIB,JVIB,JVIB)                     ANH07900
      VALU2=VALU2+DNUM2/FRQK                                            ANH07910
      GO TO 105                                                         ANH07920
C  SUM OVER DEGENERATE MODES                                            ANH07930
  216 CONTINUE                                                          ANH07940
      WPPP= FRQI+FRQJ+FRQK                                              ANH07950
      WPMM= FRQI-FRQJ-FRQK                                              ANH07960
      WMPM=-FRQI+FRQJ-FRQK                                              ANH07970
      WMMP=-FRQI-FRQJ+FRQK                                              ANH07980
      DIJK=WPPP*WPMM*WMPM*WMMP                                          ANH07990
      DNUM3=(FRQK2-FRQI2-FRQJ2)*FRQK                                    ANH08000
      FAC=F3Q(IVIB,JVIB,KVIB)                                           ANH08010
      FACT=FAC*FAC                                                      ANH08020
      IF(DABS(WPMM).LE.FLIM2) GO TO 217                                 ANH08030
      IF(DABS(WMPM).LE.FLIM2) GO TO 218                                 ANH08040
      IF(DABS(WMMP).LE.FLIM2) GO TO 219                                 ANH08050
      VALU3=VALU3+(DNUM3*FACT)/DIJK                                     ANH08060
      GO TO 105                                                         ANH08070
C   FERMI RESONANCE CORRECTIONS                                         ANH08080
  217 CONTINUE                                                          ANH08090
      IANH2=IANH2+1                                                     ANH08100
      DNUM3=ONE/WPPP-ONE/WMMP-ONE/WMPM                                  ANH08110
      VALU3=VALU3-(DNUM3*FACT)/4.0D+00                                  ANH08120
      GO TO 105                                                         ANH08130
  218 CONTINUE                                                          ANH08140
      IANH2=IANH2+1                                                     ANH08150
      DNUM3=ONE/WPPP-ONE/WMMP-ONE/WPMM                                  ANH08160
      VALU3=VALU3-(DNUM3*FACT)/4.0D+00                                  ANH08170
      GO TO 105                                                         ANH08180
  219 CONTINUE                                                          ANH08190
      IANH2=IANH2+1                                                     ANH08200
      DNUM3=ONE/WPPP-ONE/WMPM-ONE/WPMM                                  ANH08210
      VALU3=VALU3-(DNUM3*FACT)/4.0D+00                                  ANH08220
  105 CONTINUE                                                          ANH08230
      ZFACT=(ZETA(IVIB,JVIB,IAA)**2+ZETA(IVIB,JVIB,IBB)**2)*ROTCM(IBB)  ANH08240
      VALU4=(FRQI/FRQJ+FRQJ/FRQI)*ZFACT                                 ANH08250
      XRS=(VALU1-VALU2-VALU3*2.0D+00)/4.0D+00+VALU4                     ANH08260
      GO TO 230                                                         ANH08270
C                                                                       ANH08280
C++++++++++++++++++++++++++++++++++++++++++++++++++++++                 ANH08290
C+++OFF-DIAGONAL CONSTANTS (DEGENERATE & DEGENERATE)+++                 ANH08300
C++++++++++++++++++++++++++++++++++++++++++++++++++++++                 ANH08310
  225 CONTINUE                                                          ANH08320
      VALU1=(F4Q(IVIB,IVIB,JVIB,JVIB)+F4Q(IVIB,IVIB,JVIB,JVIB))*0.5D+00 ANH08330
      DO 109 KVIB=1,NVIB                                                ANH08340
C   SUM OVER NON-DEGENERATE MODES                                       ANH08350
      FRQK=FREQ(KVIB)                                                   ANH08360
      IF(DABS(FRQK).LE.WLIMIT) GO TO 109                                ANH08370
      FRQK2=FRQK*FRQK                                                   ANH08380
      WPPP= FRQI+FRQJ+FRQK                                              ANH08390
      WPMM= FRQI-FRQJ-FRQK                                              ANH08400
      WMPM=-FRQI+FRQJ-FRQK                                              ANH08410
      WMMP=-FRQI-FRQJ+FRQK                                              ANH08420
      DIJK=WPPP*WPMM*WMPM*WMMP                                          ANH08430
      DNUM2=F3Q(IVIB,IVIB,KVIB)*F3Q(KVIB,JVIB,JVIB)                     ANH08440
      DNUM3=(FRQK2-FRQI2-FRQJ2)*FRQK                                    ANH08450
      FAC=F3Q(IVIB,JVIB,KVIB)                                           ANH08460
      FACT=FAC*FAC                                                      ANH08470
      IF(NDEG(KVIB).NE.1) GO TO 229                                     ANH08480
      VALU2=VALU2+DNUM2/FRQK                                            ANH08490
      IF(DABS(WPMM).LE.FLIM2) GO TO 226                                 ANH08500
      IF(DABS(WMPM).LE.FLIM2) GO TO 227                                 ANH08510
      IF(DABS(WMMP).LE.FLIM2) GO TO 228                                 ANH08520
      VALU3=VALU3+(DNUM3*FACT)/DIJK                                     ANH08530
      GVAL1=GVAL1+(FRQI*FRQJ*FRQK*FACT*0.5D+00)/DIJK                    ANH08540
      GO TO 109                                                         ANH08550
C   FERMI RESONANCE CORRECTIONS                                         ANH08560
  226 CONTINUE                                                          ANH08570
      IANH2=IANH2+1                                                     ANH08580
      DNUM3=ONE/WPPP-ONE/WMMP-ONE/WMPM                                  ANH08590
      VALU3=VALU3-(DNUM3*FACT)/4.0D+00                                  ANH08600
      GO TO 109                                                         ANH08610
  227 CONTINUE                                                          ANH08620
      IANH2=IANH2+1                                                     ANH08630
      DNUM3=ONE/WPPP-ONE/WMMP-ONE/WPMM                                  ANH08640
      VALU3=VALU3-(DNUM3*FACT)/4.0D+00                                  ANH08650
      GO TO 109                                                         ANH08660
  228 CONTINUE                                                          ANH08670
      IANH2=IANH2+1                                                     ANH08680
      DNUM3=ONE/WPPP-ONE/WMPM-ONE/WPMM                                  ANH08690
      VALU3=VALU3-(DNUM3*FACT)/4.0D+00                                  ANH08700
      GO TO 109                                                         ANH08710
  229 CONTINUE                                                          ANH08720
C   SUM OVER DEGENERATE MODES                                           ANH08730
      DFAC=DBLE(NDEG(KVIB))                                             ANH08740
      VALU4=VALU4+(DNUM3*FACT)/(DIJK*DFAC)                              ANH08750
      GVAL2=GVAL2-(FRQI*FRQJ*FRQK*FACT)/(DIJK*DFAC)                     ANH08760
  109 CONTINUE                                                          ANH08770
      ZFACT=ZETA(IVIB,JVIB,ICC)**2*0.5D+00*ROTCM(IAA)                   ANH08780
     1     +ZETA(IVIB,JVIB,IBB)**2*ROTCM(IBB)                           ANH08790
      VALU5=(FRQI/FRQJ+FRQJ/FRQI)*ZFACT                                 ANH08800
      GVAL3=-ZETA(IVIB,JVIB,IBB)**2*ROTCM(IBB)*2.0D+00                  ANH08810
     1      +ZETA(IVIB,JVIB,ICC)**2*ROTCM(IAA)                          ANH08820
     2      +ZETA(IVIB,JVIB,ICC)*ZETA(IVIB,JVIB,ICC)*ROTCM(IAA)*2.0D+00 ANH08830
      XRS=(VALU1-VALU2-VALU3-VALU4*2.0D+00)/4.0D+00+VALU5               ANH08840
      GRS=GVAL1+GVAL2+GVAL3                                             ANH08850
C                                                                       ANH08860
  230 CONTINUE                                                          ANH08870
      XIJ(IVIB,JVIB)=XRS                                                ANH08880
      GIJ(IVIB,JVIB)=GRS                                                ANH08890
      IF(IPRNT.LE.2) GO TO 120                                          ANH08900
      WRITE(6,2) IVIB,JVIB,FRQI,FRQJ,VALU1,VALU2,VALU3,VALU4,XRS,GRS    ANH08910
  120 CONTINUE                                                          ANH08920
      IF(IANH1.NE.0) THEN                                               ANH08930
        WRITE(6,7)                                                      ANH08940
        WRITE(6,9) FLIM1                                                ANH08950
      END IF                                                            ANH08960
      IF(IANH2.NE.0) THEN                                               ANH08970
        WRITE(6,8)                                                      ANH08980
        WRITE(6,9) FLIM2                                                ANH08990
      END IF                                                            ANH09000
C                                                                       ANH09010
      IF(IPRNT.GT.2) GO TO 231                                          ANH09020
      WRITE(6,3)                                                        ANH09030
      DO 122 II=1,NVIB                                                  ANH09040
      IVIB=NFRQ(II)                                                     ANH09050
      FRQI=FREQ(IVIB)                                                   ANH09060
      DO 122 JJ=1,II                                                    ANH09070
      JVIB=NFRQ(JJ)                                                     ANH09080
      FRQJ=FREQ(JVIB)                                                   ANH09090
      XRS=XIJ(IVIB,JVIB)                                                ANH09100
      GRS=GIJ(IVIB,JVIB)                                                ANH09110
      WRITE(6,4) II,JJ,FRQI,FRQJ,XRS,GRS                                ANH09120
  122 CONTINUE                                                          ANH09130
C                                                                       ANH09140
C:::::::::::::::::::::::::::::::::::::::                                ANH09150
C:::CALCULATE FUNDAMENTAL FREQUENCIES:::                                ANH09160
C:::::::::::::::::::::::::::::::::::::::                                ANH09170
  231 CONTINUE                                                          ANH09180
      WRITE(6,5)                                                        ANH09190
      DO 125 II=1,NVIB                                                  ANH09200
      IVIB=NFRQ(II)                                                     ANH09210
      VALU1=FREQ(IVIB)                                                  ANH09220
      VALU2=XIJ(IVIB,IVIB)*2.0D+00                                      ANH09230
      VALU3=A00                                                         ANH09240
      DO 124 JJ=1,NVIB                                                  ANH09250
      JVIB=NFRQ(JJ)                                                     ANH09260
      IF(IVIB.EQ.JVIB) GO TO 124                                        ANH09270
      VALU3=VALU3+XIJ(IVIB,JVIB)                                        ANH09280
  124 CONTINUE                                                          ANH09290
      DELTA=VALU2+VALU3*0.5D+00                                         ANH09300
      VALUT=VALU1+DELTA                                                 ANH09310
      WRITE(6,6) II,VALU1,DELTA,VALUT                                   ANH09320
  125 CONTINUE                                                          ANH09330
C                                                                       ANH09340
      RETURN                                                            ANH09350
      END                                                               ANH09360
