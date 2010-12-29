C     PROGRAM ANHLIN                                           
C**********************************************************************
C***LAST UPDATED ON January 10, 2006 BY DR. YUKIO YAMAGUCHI         ***
C   Modification for a stand-alone version                          ***
C**********************************************************************
C***LAST UPDATED ON April 20, 2004 BY YUKIO YAMAGUCHI               ***
C   Modification                                                    ***
C   for anharmonic zero-point vibrational energy (ZPVE)             ***
C   and vibrational G_0 energy                                      ***
C**********************************************************************     
C***LAST UPDATED ON JANUARY 25, 1989 BY YUKIO YAMAGUCHI             ***             
C**********************************************************************  
      SUBROUTINE LTYPLI(ZETA,CXY,F3Q,CLTYP)                             ANH00050
C   THE CALCULATION OF L-TYPE DOUBLING CONSTANTS                        ANH00060
C   FOR A LINEAR POLYATOMIC MOLECULE                                    ANH00070
      IMPLICIT REAL*8 (A-H,O-Z)                                         ANH00080
      DIMENSION ZETA(N3N,N3N,3),CXY(3,3,N3N),F3Q(N3N,N3N,N3N)           ANH00090
      DIMENSION CLTYP(N3N,3,3)                                          ANH00100
      DIMENSION CCC(3,3)                                                ANH00110
      COMMON/VIB101/NATOM,N3N,NATRI,ILIN,NVIB                           ANH00120
      COMMON/VIB106/ROTAA(3),ROTGC(3),ROTCM(3),ROTMH(3)                 ANH00130
      COMMON/VIB108/IAXIS(3),NDEG(150),NDAB(150,5),IMAG(150)            ANH00140
      COMMON/VIB109/CLIMIT,FLIM1,FLIM2                                  ANH00150
      COMMON/VIB203/IOFF(150),IPRNT                                     ANH00160
      COMMON/VIB204/SQM(150),ROOT(150),FREQ(150)                        ANH00170
      COMMON/VIB205/IFREQ,NFRQ(150)                                     ANH00180
      DATA A00,HALF,ONE,TWO / 0.0D+00 , 0.5D+00 , 1.0D+00 , 2.0D+00 /   ANH00190
      DATA WLIMIT / 1.0D-02 /                                           ANH00200
      DATA DLIMIT / 1.0D-02 /                                           ANH00210
      DATA PI,CL / 3.1415926536D+00 , 2.99792458D+00 /                  ANH00220
    1 FORMAT(//,2X,' ***********************************'/              ANH00230
     1          2X,' ***THE L-TYPE DOUBLING CONSTANTS***'/              ANH00240
     2          2X,' ***********************************')              ANH00250
    2 FORMAT(//,2X,' ZETA MATRIX, IABC = ',I5/)                         ANH00260
    3 FORMAT(2X,I5,4F10.5)                                              ANH00270
    4 FORMAT(//,2X,' QE(T) VALU'/)                                      ANH00280
    5 FORMAT(2X,' QVIB',5X,' IVIB',7X,' FREQENCY',9X,' HARMONIC',       ANH00290
     1      11X,' ANHARMONIC',9X,' TOTAL',11X,' TOTAL'/                 ANH00300
     2      24X,' IN CM-1',10X,' IN CM-1',12X,' IN CM-1',               ANH00310
     3      12X,' IN CM-1',9X,' IN MHZ'/)                               ANH00320
    6 FORMAT(2X,I5,5X,I5,8X,F10.5,F20.10,F20.10,F20.10,F20.10)          ANH00330
    7 FORMAT(//,2X,' DE = ',F20.10/                                     ANH00340
     1          2X,' BE = ',F20.10/)                                    ANH00350
    8 FORMAT(//,2X,' QET = ',F20.10/                                    ANH00360
     1          2X,' AET = ',F20.10/)                                   ANH00370
    9 FORMAT(//,2X,' VALU4H = ',F20.10/                                 ANH00380
     1          2X,' VALU4A = ',F20.10/                                 ANH00390
     2          2X,' VALU4  = ',F20.10/                                 ANH00400
     3          2X,' VALU4X = ',F20.10/)                                ANH00410
   10 FORMAT(//,2X,' VALU8T = ',F20.10/                                 ANH00420
     1          2X,' VALU8S = ',F20.10/                                 ANH00430
     2          2X,' VALU8  = ',F20.10/)                                ANH00440
   11 FORMAT(//,2X,' VALU9T = ',F20.10/                                 ANH00450
     1          2X,' VALU9S = ',F20.10/                                 ANH00460
     2          2X,' VALU9SP= ',F20.10/                                 ANH00470
     3          2X,' VALU9  = ',F20.10/)                                ANH00480
   12 FORMAT(//,2X,' IVIB   =  ',I5/                                    ANH00490
     1          2X,' VALU1  = ',F20.10/                                 ANH00500
     2          2X,' VALU2  = ',F20.10/                                 ANH00510
     3          2X,' VALU3  = ',F20.10/                                 ANH00520
     4          2X,' VALU4  = ',F20.10/                                 ANH00530
     5          2X,' VALU5  = ',F20.10/                                 ANH00540
     6          2X,' VALU6  = ',F20.10/                                 ANH00550
     7          2X,' VALU7  = ',F20.10/                                 ANH00560
     8          2X,' VALU8  = ',F20.10/                                 ANH00570
     9          2X,' VALU9  = ',F20.10/                                 ANH00580
     A          2X,' VALU10 = ',F20.10/                                 ANH00590
     B          2X,' QJT    = ',F20.10/                                 ANH00600
     C          2X,' VALUH  = ',F20.10/                                 ANH00610
     D          2X,' VALUA  = ',F20.10/                                 ANH00620
     E          2X,' VALUT  = ',F20.10/)                                ANH00630
   13 FORMAT(//,2X,' QJ(T) VALU'/)                                      ANH00640
   14 FORMAT(2X,' QVIB',5X,' IVIB',7X,' FREQENCY',9X,' HARMONIC',       ANH00650
     1      11X,' ANHARMONIC',9X,' TOTAL',11X,' TOTAL'/                 ANH00660
     2      24X,' IN CM-1',10X,' IN D-07 CM-1',7X,' IN D-07 CM-1',      ANH00670
     3       7X,' IN D-07 CM-1',4X,' IN KHZ'/)                          ANH00680
   15 FORMAT(//,2X,' QK(T) VALU'/)                                      ANH00690
C                                                                       ANH00700
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ANH00710
C%  HIGHER-ORDER L-TYPE DOUBLING CONSTANTS FOR LINEAR MOLEUCLES        %ANH00720
C%  REFERENCE : J.K.G. WATSON, J. MOL. SPECTROSC. 101,83 (1983).       %ANH00730
C%                                                                     %ANH00740
C%  NOTE :                                                             %ANH00750
C%         ZETA(S,T)=-ZETA(T,S)                                        %ANH00760
C%         IN THIS SUBROUTINE ZETA(S,T) HAS BEEN USED.                 %ANH00770
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ANH00780
C                                                                       ANH00790
      WRITE(6,1)                                                        ANH00800
C                                                                       ANH00810
C   NOTE : THE MOLECULAR AXIS IS CHOSEN AS THE Z AXIS                   ANH00820
      IXX=1                                                             ANH00830
      IYY=2                                                             ANH00840
      IZZ=3                                                             ANH00850
      IF(IPRNT.LE.2) GO TO 201                                          ANH00860
      DO 101 IABC=1,3                                                   ANH00870
      WRITE(6,2) IABC                                                   ANH00880
      CALL MATOUT(ZETA(1,1,IABC),N3N,N3N,NVIB,NVIB,6)                   ANH00890
  101 CONTINUE                                                          ANH00900
      IF(IPRNT.GT.2) THEN                                               ANH00910
        DO 150 IABC=1,NVIB                                              ANH00920
        WRITE(6,*) ' IABC = ',IABC                                      ANH00930
        DO 149 I=1,3                                                    ANH00940
        DO 149 J=1,3                                                    ANH00950
        CCC(I,J)=CXY(I,J,IABC)*1.0D+07                                  ANH00960
  149   CONTINUE                                                        ANH00970
        CALL MATOUT(CCC,3,3,3,3,6)                                      ANH00980
  150   CONTINUE                                                        ANH00990
      END IF                                                            ANH01000
C                                                                       ANH01010
C:::::::::::::::::::::::::::::::                                        ANH01020
C:::ONLY FOR DEGENERATE MODES:::                                        ANH01030
C:::::::::::::::::::::::::::::::                                        ANH01040
  201 CONTINUE                                                          ANH01050
C   CALCULATION OF Q(E)                                                 ANH01060
      CALL ZERO(CLTYP,N3N*9)                                            ANH01070
C***********************************************************************ANH01080
C***FOLLOWINGS ARE FOR AN ALTERNATE EXPRESSION                         *ANH01090
C***THEY ARE TESTED TO GIVE THE SAME RESULTS                           *ANH01100
C***********************************************************************ANH01110
C*    DO 103 IVIB=1,NVIB                                               *ANH01120
C*    IF(NDEG(IVIB).EQ.1) GO TO 103                                    *ANH01130
C*    FRQI=FREQ(IVIB)                                                  *ANH01140
C*    FRQI2=FRQI*FRQI                                                  *ANH01150
C                                                                      *ANH01160
C   SUM OVER NON-DEGENERATE MODES                                      *ANH01170
C*    VALU2=A00                                                        *ANH01180
C*    DO 102 JVIB=1,NVIB                                               *ANH01190
C*    IF(NDEG(JVIB).NE.1) GO TO 102                                    *ANH01200
C*    FRQJ=FREQ(JVIB)                                                  *ANH01210
C*    FRQJ2=FRQJ*FRQJ                                                  *ANH01220
C*    FACT2=(FRQI2*3.0D+00+FRQJ2)/(FRQI2-FRQJ2)                        *ANH01230
C*    VALU2=VALU2+(ZETA(IVIB,JVIB,IXX)**2+ZETA(IVIB,JVIB,IYY)**2)*FACT2*ANH01240
C*102 CONTINUE                                                         *ANH01250
C*    VALUH=-(VALU2*ROTCM(IXX)*ROTCM(IXX)*TWO)/FRQI                    *ANH01260
C*    VALUA=A00                                                        *ANH01270
C*    VALUT=VALUH+VALUA                                                *ANH01280
C*    CLTYP(IVIB,1,1)=VALUH                                            *ANH01290
C*    CLTYP(IVIB,2,1)=VALUA                                            *ANH01300
C*    CLTYP(IVIB,3,1)=VALUT                                            *ANH01310
C*    IF(IPRNT.LE.2) GO TO 103                                         *ANH01320
C*    WRITE(6,3) IVIB,FRQI,VALUH,VALUA,VALUT                           *ANH01330
C*103 CONTINUE                                                         *ANH01340
C***********************************************************************ANH01350
C                                                                       ANH01360
      BE=ROTCM(IXX)                                                     ANH01370
      BE2=BE*BE                                                         ANH01380
      DO 105 IVIB=1,NVIB                                                ANH01390
      IF(NDEG(IVIB).EQ.1) GO TO 105                                     ANH01400
      FRQI=FREQ(IVIB)                                                   ANH01410
      FRQI2=FRQI*FRQI                                                   ANH01420
C                                                                       ANH01430
C   SUM OVER NON-DEGENERATE MODES                                       ANH01440
      VALU2=ONE                                                         ANH01450
      DO 104 JVIB=1,NVIB                                                ANH01460
      IF(NDEG(JVIB).NE.1) GO TO 104                                     ANH01470
      FRQJ=FREQ(JVIB)                                                   ANH01480
      FRQJ2=FRQJ*FRQJ                                                   ANH01490
      FACT2=(FRQI2*4.0D+00)/(FRQJ2-FRQI2)                               ANH01500
      VALU2=VALU2+(ZETA(IVIB,JVIB,IXX)**2+ZETA(IVIB,JVIB,IYY)**2)*FACT2 ANH01510
  104 CONTINUE                                                          ANH01520
      VALUH=(VALU2*BE2*TWO)/FRQI                                        ANH01530
      VALUA=A00                                                         ANH01540
      VALUT=VALUH+VALUA                                                 ANH01550
      CLTYP(IVIB,1,1)=VALUH                                             ANH01560
      CLTYP(IVIB,2,1)=VALUA                                             ANH01570
      CLTYP(IVIB,3,1)=VALUT                                             ANH01580
      IF(IPRNT.LE.2) GO TO 105                                          ANH01590
      WRITE(6,3) IVIB,FRQI,VALUH,VALUA,VALUT                            ANH01600
  105 CONTINUE                                                          ANH01610
C                                                                       ANH01620
      WRITE(6,4)                                                        ANH01630
      WRITE(6,5)                                                        ANH01640
      DO 106 II=1,NVIB                                                  ANH01650
      IVIB=NFRQ(II)                                                     ANH01660
      IF(NDEG(IVIB).LE.1) GO TO 106                                     ANH01670
      FRQI=FREQ(IVIB)                                                   ANH01680
      VALUH=CLTYP(IVIB,1,1)                                             ANH01690
      VALUA=CLTYP(IVIB,2,1)                                             ANH01700
      VALUT=CLTYP(IVIB,3,1)                                             ANH01710
      VALUM=VALUT*CL*1.0D+04                                            ANH01720
      WRITE(6,6) II,IVIB,FRQI,VALUH,VALUA,VALUT,VALUM                   ANH01730
  106 CONTINUE                                                          ANH01740
C                                                                       ANH01750
C   HIGHER-ORDER L-TYPE DOUBLING CONSTANTS                              ANH01760
      VALU0=A00                                                         ANH01770
      DO 107 IVIB=1,NVIB                                                ANH01780
      IF(NDEG(IVIB).NE.1) GO TO 107                                     ANH01790
      FRQI=FREQ(IVIB)                                                   ANH01800
      FACT=CXY(IXX,IXX,IVIB)*CXY(IXX,IXX,IVIB)                          ANH01810
      VALU0=VALU0+FACT*FRQI                                             ANH01820
  107 CONTINUE                                                          ANH01830
      DE=VALU0*HALF                                                     ANH01840
      BE=ROTCM(IXX)                                                     ANH01850
      BE2=BE*BE                                                         ANH01860
      BE4=BE2*BE2                                                       ANH01870
      IF(IPRNT.GT.2) THEN                                               ANH01880
        WRITE(6,7) DE,BE                                                ANH01890
      END IF                                                            ANH01900
C                                                                       ANH01910
      DO 130 IVIB=1,NVIB                                                ANH01920
      IF(NDEG(IVIB).LE.1) GO TO 130                                     ANH01930
      FRQI=FREQ(IVIB)                                                   ANH01940
      FRQI2=FRQI*FRQI                                                   ANH01950
      FRQI3=FRQI2*FRQI                                                  ANH01960
      FRQI4=FRQI2*FRQI2                                                 ANH01970
C                                                                       ANH01980
C   CALCULATE QE(T) AND AE(T)                                           ANH01990
C   1ST , 2ND AND 3RD TERMS                                             ANH02000
      VAL1=A00                                                          ANH02010
      VAL2=A00                                                          ANH02020
      DO 111 JVIB=1,NVIB                                                ANH02030
      IF(NDEG(JVIB).NE.1) GO TO 111                                     ANH02040
      FRQJ=FREQ(JVIB)                                                   ANH02050
      FRQJ2=FRQJ*FRQJ                                                   ANH02060
      DNUM=(FRQI2*3.0D+00+FRQJ2)                                        ANH02070
      FACT=(ZETA(JVIB,IVIB,IXX)**2+ZETA(JVIB,IVIB,IYY)**2)              ANH02080
      DNOM=(FRQJ2-FRQI2)                                                ANH02090
      VAL1=VAL1+(FACT*DNUM)/DNOM                                        ANH02100
      VAL2=VAL2+CXY(IXX,IXX,JVIB)*F3Q(JVIB,IVIB,IVIB)                   ANH02110
  111 CONTINUE                                                          ANH02120
      QET=VAL1*BE2*TWO/FRQI                                             ANH02130
      AET=QET*HALF-VAL2*HALF                                            ANH02140
      IF(IPRNT.GT.2) THEN                                               ANH02150
        WRITE(6,8) QET,AET                                              ANH02160
      END IF                                                            ANH02170
C                                                                       ANH02180
      VALU1=-(4.0D+00*DE*QET)/BE                                        ANH02190
      VALU2=(TWO*BE4)/FRQI3                                             ANH02200
      VALU3=BE2*(AET*TWO-QET)/FRQI2                                     ANH02210
C                                                                       ANH02220
C   4TH , 5TH AND 6TH TERMS                                             ANH02230
      VAL4=A00                                                          ANH02240
      VAL5=A00                                                          ANH02250
      VAL6=A00                                                          ANH02260
      DO 112 JVIB=1,NVIB                                                ANH02270
      IF(NDEG(JVIB).NE.1) GO TO 112                                     ANH02280
      FRQJ=FREQ(JVIB)                                                   ANH02290
      FRQJ2=FRQJ*FRQJ                                                   ANH02300
      FRQJH=DSQRT(FRQJ)                                                 ANH02310
      FACT4=(ZETA(JVIB,IVIB,IXX)**2+ZETA(JVIB,IVIB,IYY)**2)             ANH02320
      DNUM4=FACT4*(FRQJ2+FRQI2)                                         ANH02330
      DNOM4=(FRQJ2-FRQI2)**2                                            ANH02340
      VAL4=VAL4+DNUM4/DNOM4                                             ANH02350
      DNUM5=CXY(IXX,IXX,JVIB)*(ZETA(JVIB,IVIB,IXX)+ZETA(JVIB,IVIB,IYY)) ANH02360
     1     *FRQJ**(1.50D+00)                                            ANH02370
      DNOM5=(FRQJ2-FRQI2)                                               ANH02380
      VAL5=VAL5+DNUM5/DNOM5                                             ANH02390
      VAL6=VAL6+CXY(IXX,IXX,JVIB)*(ZETA(JVIB,IVIB,IXX)                  ANH02400
     1                            +ZETA(JVIB,IVIB,IYY))/FRQJH           ANH02410
  112 CONTINUE                                                          ANH02420
      VALU4H=-(8.0D+00*BE2*QET)*VAL4                                    ANH02430
      VALU4A=(4.0D+00*BE2*VAL2)*VAL4                                    ANH02440
      VALU4=VALU4H+VALU4A                                               ANH02450
      VALU5=-(12.0D+00*BE*FRQI)*(VAL5*VAL5)                             ANH02460
      VALU6=-(8.0D+00*BE*FRQI)*(VAL5*VAL6)                              ANH02470
      IF(IPRNT.GT.2) THEN                                               ANH02480
        VALU4X=-(4.0D+00*BE2*(TWO*AET+QET))*VAL4                        ANH02490
        WRITE(6,9) VALU4H,VALU4A,VALU4,VALU4X                           ANH02500
      END IF                                                            ANH02510
C                                                                       ANH02520
C   7TH TERM                                                            ANH02530
      VAL7=A00                                                          ANH02540
      DO 115 JVIB=1,NVIB                                                ANH02550
      IF(NDEG(JVIB).NE.1) GO TO 115                                     ANH02560
      FRQJ=FREQ(JVIB)                                                   ANH02570
      FRQJ2=FRQJ*FRQJ                                                   ANH02580
      DO 114 KVIB=1,NVIB                                                ANH02590
      IF(NDEG(KVIB).NE.1) GO TO 114                                     ANH02600
      FRQK=FREQ(KVIB)                                                   ANH02610
      FRQK2=FRQK*FRQK                                                   ANH02620
      FRQKH=DSQRT(FRQK)                                                 ANH02630
      DO 113 LVIB=1,NVIB                                                ANH02640
      IF(NDEG(LVIB).NE.1) GO TO 113                                     ANH02650
      FRQL=FREQ(LVIB)                                                   ANH02660
      FRQL2=FRQL*FRQL                                                   ANH02670
      FRQLH=DSQRT(FRQL)                                                 ANH02680
      DNUM=CXY(IXX,IXX,JVIB)*F3Q(JVIB,KVIB,LVIB)                        ANH02690
     1    *(ZETA(KVIB,IVIB,IXX)+ZETA(KVIB,IVIB,IYY))                    ANH02700
     2    *(ZETA(LVIB,IVIB,IXX)+ZETA(LVIB,IVIB,IYY))* FRQKH*FRQLH       ANH02710
      DNOM=(FRQK2-FRQI2)*(FRQL2-FRQI2)                                  ANH02720
      VAL7=VAL7+DNUM/DNOM                                               ANH02730
  113 CONTINUE                                                          ANH02740
  114 CONTINUE                                                          ANH02750
  115 CONTINUE                                                          ANH02760
      VALU7=-(VAL7*BE2*FRQI)*8.0D+00                                    ANH02770
C                                                                       ANH02780
C   8TH TERM                                                            ANH02790
      VAL8T=A00                                                         ANH02800
      VAL8S=A00                                                         ANH02810
      DO 117 JVIB=1,NVIB                                                ANH02820
      IF(NDEG(JVIB).NE.2) GO TO 117                                     ANH02830
      FRQJ=FREQ(JVIB)                                                   ANH02840
      FRQJ2=FRQJ*FRQJ                                                   ANH02850
      IF(DABS(FRQI-FRQJ).LE.DLIMIT) GO TO 117                           ANH02860
      VAL8T=VAL8T+(FRQI2*3.0D+00+FRQJ2)/(FRQI2-FRQJ2)                   ANH02870
      DO 116 KVIB=1,NVIB                                                ANH02880
      IF(NDEG(KVIB).NE.1) GO TO 116                                     ANH02890
      FRQK=FREQ(KVIB)                                                   ANH02900
      FRQK2=FRQK*FRQK                                                   ANH02910
      DNOMS=(FRQK2-FRQI2)                                               ANH02920
      VAL8S=VAL8S+(ZETA(KVIB,IVIB,IXX)+ZETA(KVIB,IVIB,IYY))             ANH02930
     1           *(ZETA(KVIB,JVIB,IXX)+ZETA(KVIB,JVIB,IYY))/DNOMS       ANH02940
  116 CONTINUE                                                          ANH02950
  117 CONTINUE                                                          ANH02960
      VAL8T=VAL8T*HALF                                                  ANH02970
      VAL8S=VAL8S*HALF                                                  ANH02980
      VALU8=-(BE4*FRQI*VAL8T*VAL8S*VAL8S)*32.0D+00                      ANH02990
      IF(IPRNT.GT.2) THEN                                               ANH03000
        WRITE(6,10) VAL8T,VAL8S,VALU8                                   ANH03010
      END IF                                                            ANH03020
C                                                                       ANH03030
C   9TH TERM                                                            ANH03040
      VAL9T=A00                                                         ANH03050
      VAL9S=A00                                                         ANH03060
      VAL9SP=A00                                                        ANH03070
      DO 120 JVIB=1,NVIB                                                ANH03080
      IF(NDEG(JVIB).NE.2) GO TO 120                                     ANH03090
      FRQJ=FREQ(JVIB)                                                   ANH03100
      FRQJ2=FRQJ*FRQJ                                                   ANH03110
      IF(DABS(FRQI-FRQJ).LE.DLIMIT) GO TO 120                           ANH03120
      VAL9T=VAL9T+DSQRT(FRQI3*FRQJ)/(FRQI2-FRQJ2)                       ANH03130
      DO 118 KVIB=1,NVIB                                                ANH03140
      IF(NDEG(KVIB).NE.1) GO TO 118                                     ANH03150
      FRQK=FREQ(KVIB)                                                   ANH03160
      FRQK2=FRQK*FRQK                                                   ANH03170
      DNOMS=FRQK2-FRQI2                                                 ANH03180
      VAL9S=VAL9S+(ZETA(KVIB,IVIB,IXX)+ZETA(KVIB,IVIB,IYY))             ANH03190
     1           *(ZETA(KVIB,JVIB,IXX)+ZETA(KVIB,JVIB,IYY))/DNOMS       ANH03200
  118 CONTINUE                                                          ANH03210
      DO 119 LVIB=1,NVIB                                                ANH03220
      IF(NDEG(LVIB).NE.1) GO TO 119                                     ANH03230
      VAL9SP=VAL9SP+CXY(IXX,IXX,LVIB)*F3Q(LVIB,IVIB,JVIB)               ANH03240
  119 CONTINUE                                                          ANH03250
  120 CONTINUE                                                          ANH03260
      VAL9T=VAL9T*HALF                                                  ANH03270
      VAL9S=VAL9S*HALF                                                  ANH03280
      VAL9SP=VAL9SP*HALF                                                ANH03290
      VALU9=(BE2*VAL9T*VAL9S*VAL9SP)*16.0D+00                           ANH03300
      IF(IPRNT.GT.2) THEN                                               ANH03310
        WRITE(6,11) VAL9T,VAL9S,VAL9SP,VALU9                            ANH03320
      END IF                                                            ANH03330
C                                                                       ANH03340
      VAL10=A00                                                         ANH03350
      DO 121 JVIB=1,NVIB                                                ANH03360
      IF(NDEG(JVIB).NE.1) GO TO 121                                     ANH03370
      FRQJ=FREQ(JVIB)                                                   ANH03380
      FRQJ2=FRQJ*FRQJ                                                   ANH03390
      FRQJ4=FRQJ2*FRQJ2                                                 ANH03400
      DNOM10=FRQI*(FRQJ2-FRQI2)**3                                      ANH03410
      DNUM10=(ZETA(JVIB,IVIB,IXX)**2+ZETA(JVIB,IVIB,IYY)**2)            ANH03420
     1      *(FRQJ4+FRQJ2*FRQI2*10.0D+00+FRQI4*5.0D+00)                 ANH03430
      VAL10=VAL10+DNUM10/DNOM10                                         ANH03440
  121 CONTINUE                                                          ANH03450
      VALU10=-(VAL10*BE4)*8.0D+00                                       ANH03460
C                                                                       ANH03470
      VALUH=VALU1+VALU2+VALU4H+VALU5+VALU6+VALU8                        ANH03480
      VALUA=VALU3+VALU4A+VALU7+VALU9                                    ANH03490
      QJT=VALU1+VALU2+VALU3+VALU4+VALU5+VALU6+VALU7+VALU8+VALU9         ANH03500
CC    IF(IPRNT.GT.2) THEN                                               ANH03510
CC      VALU1=VALU1*1.0D+08                                             ANH03520
CC      VALU2=VALU2*1.0D+08                                             ANH03530
CC      VALU3=VALU3*1.0D+08                                             ANH03540
CC      VALU4=VALU4*1.0D+08                                             ANH03550
CC      VALU5=VALU5*1.0D+08                                             ANH03560
CC      VALU6=VALU6*1.0D+08                                             ANH03570
CC      VALU7=VALU7*1.0D+08                                             ANH03580
CC      VALU8=VALU8*1.0D+08                                             ANH03590
CC      VALU9=VALU9*1.0D+08                                             ANH03600
CC      VALU10=VALU10*1.0D+08                                           ANH03610
CC      VALUH=VALUH*1.0D+08                                             ANH03620
CC      VALUA=VALUA*1.0D+08                                             ANH03630
C                                                                       ANH03640
CC      VALUTX=VALUH+VALUA                                              ANH03650
CC      WRITE(6,12) IVIB,VALU1,VALU2,VALU3,VALU4,VALU5,VALU6,VALU7,     ANH03660
CC   1              VALU8,VALU9,VALU10,QJT,VALUH,VALUA,VALUTX           ANH03670
CC    END IF                                                            ANH03680
      CLTYP(IVIB,1,2)=VALUH                                             ANH03690
      CLTYP(IVIB,2,2)=VALUA                                             ANH03700
      CLTYP(IVIB,3,2)=QJT                                               ANH03710
      QKT=-QJT+VALU10                                                   ANH03720
      CLTYP(IVIB,1,3)=-VALUH+VALU10                                     ANH03730
      CLTYP(IVIB,2,3)=-VALUA                                            ANH03740
      CLTYP(IVIB,3,3)=QKT                                               ANH03750
C                                                                       ANH03760
  130 CONTINUE                                                          ANH03770
C                                                                       ANH03780
      WRITE(6,13)                                                       ANH03790
      WRITE(6,14)                                                       ANH03800
      DO 132 II=1,NVIB                                                  ANH03810
      IVIB=NFRQ(II)                                                     ANH03820
      IF(NDEG(IVIB).LE.1) GO TO 132                                     ANH03830
      FRQI=FREQ(IVIB)                                                   ANH03840
      VALUH=CLTYP(IVIB,1,2)*1.0D+07                                     ANH03850
      VALUA=CLTYP(IVIB,2,2)*1.0D+07                                     ANH03860
      VALUT=CLTYP(IVIB,3,2)*1.0D+07                                     ANH03870
      VALUM=VALUT*CL                                                    ANH03880
      WRITE(6,6) II,IVIB,FRQI,VALUH,VALUA,VALUT,VALUM                   ANH03890
  132 CONTINUE                                                          ANH03900
      WRITE(6,15)                                                       ANH03910
      WRITE(6,14)                                                       ANH03920
      DO 133 II=1,NVIB                                                  ANH03930
      IVIB=NFRQ(II)                                                     ANH03940
      IF(NDEG(IVIB).LE.1) GO TO 133                                     ANH03950
      FRQI=FREQ(IVIB)                                                   ANH03960
      VALUH=CLTYP(IVIB,1,3)*1.0D+07                                     ANH03970
      VALUA=CLTYP(IVIB,2,3)*1.0D+07                                     ANH03980
      VALUT=CLTYP(IVIB,3,3)*1.0D+07                                     ANH03990
      VALUM=VALUT*CL                                                    ANH04000
      WRITE(6,6) II,IVIB,FRQI,VALUH,VALUA,VALUT,VALUM                   ANH04010
  133 CONTINUE                                                          ANH04020
C                                                                       ANH04030
      RETURN                                                            ANH04040
      END                                                               ANH04050
      SUBROUTINE ROVLIR(ZETA,AXY,F3Q,EE,RVCNST)                         ANH04060
C   THE CALCULATION OF VIBRATION-ROTATION INTERACTION CONSTANTS         ANH04070
C   FOR A LINEAR POLYATOMIC MOLECULE                                    ANH04080
      IMPLICIT REAL*8 (A-H,O-Z)                                         ANH04090
      CHARACTER*1 RAXIS                                                 ANH04100
      CHARACTER*25 RTYPE                                                ANH04110
      DIMENSION F3Q(N3N,N3N,N3N),ZETA(N3N,N3N,3),AXY(3,3,N3N)           ANH04120
      DIMENSION EE(N3N,N3N),RVCNST(NVIB,3)                              ANH04130
      COMMON/VIB101/NATOM,N3N,NATRI,ILIN,NVIB                           ANH04140
      COMMON/VIB103/PARA,WAVE,CONST,CYCL,CONV                           ANH04150
      COMMON/VIB106/ROTAA(3),ROTGC(3),ROTCM(3),ROTMH(3)                 ANH04160
      COMMON/VIB108/IAXIS(3),NDEG(150),NDAB(150,5),IMAG(150)            ANH04170
      COMMON/VIB109/CLIMIT,FLIM1,FLIM2                                  ANH04180
      COMMON/VIB203/IOFF(150),IPRNT                                     ANH04190
      COMMON/VIB204/SQM(150),ROOT(150),FREQ(150)                        ANH04200
      COMMON/VIB205/IFREQ,NFRQ(150)                                     ANH04210
      COMMON/VIB206/RAXIS(3),RTYPE(6)                                   ANH04220
      COMMON/VIB209/SALP(3)                                             ANH04230
      DATA PH,AVN / 6.626176D+00 , 6.022045D+00 /                       ANH04240
      DATA PI,CL / 3.1415926536D+00 , 2.99792458D+00 /                  ANH04250
      DATA A00,HALF,ONE,TWO / 0.0D+00 , 0.5D+00 , 1.0D+00 , 2.0D+00 /   ANH04260
      DATA WLIMIT / 1.0D+00 /                                           ANH04270
      DATA ZLIMIT / 1.0D-05 /                                           ANH04280
    1 FORMAT(//,2X,' **************************************************'ANH04290
     1        /,2X,' ***THE VIBRATION-ROTATION INTERACTION CONSTANTS***'ANH04300
     2        /,2X,' **************************************************'ANH04310
     3        )                                                         ANH04320
    2 FORMAT(//,2X,' ZETA MATRIX, IABC = ',I5/)                         ANH04330
    3 FORMAT(//,2X,' F3(QI*QJ*QK) MATRIX (IN CM-1), IVIB = ',I5/)       ANH04340
    4 FORMAT(//,2X,' AXY MATRIX, IVIB = ',I5/)                          ANH04350
    5 FORMAT(//,2X,' QVIB',5X,' IVIB',10X,' A',21X,' B'/                ANH04360
     1          27X,' IN CM-1',15X,' IN CM-1'/)                         ANH04370
    6 FORMAT(2X,I5,5X,I5,3X,F20.10,3X,F20.10)                           ANH04380
    7 FORMAT(//,2X,' PRINCIPAL MOMENTS OF INERTIA WITH ZERO-POINT VIBRATANH04390
     1IONAL CORRECTION'/                                                ANH04400
     2 4X,5H AXIS,9X,11H IN AMU.A+2,9X,10H IN G.CM+2,10X,               ANH04410
     3 8H IN CM-1,15X,7H IN MHZ/38X,8H (*D+39)/)                        ANH04420
    8 FORMAT(8X,A1,3X,3F20.10,F20.5)                                    ANH04430
    9 FORMAT(//,2X,' ::::::::::::::::::::::::::::::::::::'/             ANH04440
     1          2X,' :::CORIOLIS RESONANCE IS OBSERVED:::'/             ANH04450
     2          2X,' ::::::::::::::::::::::::::::::::::::'/)            ANH04460
   10 FORMAT(2X,' RESONANCE LIMIT IS SET TO ',F10.3,' CM-1'/)           ANH04470
C                                                                       ANH04480
      WRITE(6,1)                                                        ANH04490
      RCONST=1.0D-01*PI*DSQRT(CL/(PH*AVN))                              ANH04500
      IAA=IAXIS(1)                                                      ANH04510
      IBB=IAXIS(2)                                                      ANH04520
      ICC=IAXIS(3)                                                      ANH04530
      NCORIO=0                                                          ANH04540
C                                                                       ANH04550
      IF(IPRNT.LE.2) GO TO 201                                          ANH04560
      DO 101 IABC=1,3                                                   ANH04570
      WRITE(6,2) IABC                                                   ANH04580
      CALL MATOUT(ZETA(1,1,IABC),N3N,N3N,N3N,N3N,6)                     ANH04590
  101 CONTINUE                                                          ANH04600
      DO 103 III=1,NVIB                                                 ANH04610
      IVIB=NFRQ(III)                                                    ANH04620
      WRITE(6,3) IVIB                                                   ANH04630
      DO 102 I=1,NVIB                                                   ANH04640
      II=NFRQ(I)                                                        ANH04650
      DO 102 J=1,NVIB                                                   ANH04660
      JJ=NFRQ(J)                                                        ANH04670
      EE(I,J)=F3Q(II,JJ,IVIB)                                           ANH04680
  102 CONTINUE                                                          ANH04690
      CALL MATOUT(EE,N3N,N3N,NVIB,NVIB,6)                               ANH04700
  103 CONTINUE                                                          ANH04710
      DO 104 IVIB=1,NVIB                                                ANH04720
      WRITE(6,4) IVIB                                                   ANH04730
      CALL MATOUT(AXY(1,1,IVIB),3,3,3,3,6)                              ANH04740
  104 CONTINUE                                                          ANH04750
C                                                                       ANH04760
C::::::::::::::::::::::::::::::                                         ANH04770
C:::FOR NON-DEGENERATE MODES:::                                         ANH04780
C::::::::::::::::::::::::::::::                                         ANH04790
  201 CONTINUE                                                          ANH04800
      DO 113 IVIB=1,NVIB                                                ANH04810
      IF(NDEG(IVIB).NE.1) GO TO 113                                     ANH04820
      FRQI=FREQ(IVIB)                                                   ANH04830
      FRQI2=FRQI*FRQI                                                   ANH04840
      VALSA1=(AXY(IAA,IAA,IVIB)**2/ROTAA(1))*0.75D+00                   ANH04850
      VALSB1=(AXY(IBB,IBB,IVIB)**2/ROTAA(2))*0.75D+00                   ANH04860
C   SUM OVER DEGENERATE MODES                                           ANH04870
      VALSA2=A00                                                        ANH04880
      VALSB2=A00                                                        ANH04890
      DO 111 JVIB=1,NVIB                                                ANH04900
      IF(NDEG(JVIB).EQ.1) GO TO 111                                     ANH04910
      FRQJ=FREQ(JVIB)                                                   ANH04920
      FRQJ2=FRQJ*FRQJ                                                   ANH04930
      IF(DABS(FRQI-FRQJ).LE.CLIMIT) GO TO 203                           ANH04940
      FACT2=0.5D+00*(FRQI2*3.0D+00+FRQJ2)/(FRQI2-FRQJ2)                 ANH04950
      GO TO 204                                                         ANH04960
C////////////////////////////////////                                   ANH04970
C///CORIOLIS RESONANCE CORRECTIONS///                                   ANH04980
C////////////////////////////////////                                   ANH04990
  203 CONTINUE                                                          ANH05000
      FACT2=-(0.25D+00*(FRQI-FRQJ)**2)/(FRQJ*(FRQI+FRQJ))               ANH05010
      IF(DABS(ZETA(IVIB,JVIB,IAA)).GT.ZLIMIT.OR.                        ANH05020
     *   DABS(ZETA(IVIB,JVIB,IBB)).GT.ZLIMIT) THEN                      ANH05030
        NCORIO=NCORIO+1                                                 ANH05040
      END IF                                                            ANH05050
C///////////////////////////////                                        ANH05060
C///END OF CORIOLIS RESONANCE///                                        ANH05070
C///////////////////////////////                                        ANH05080
  204 CONTINUE                                                          ANH05090
      ZFACT=ZETA(IVIB,JVIB,IAA)**2+ZETA(IVIB,JVIB,IBB)**2               ANH05100
      VALSA2=VALSA2+FACT2*ZFACT                                         ANH05110
      VALSB2=VALSB2+FACT2*ZFACT                                         ANH05120
  111 CONTINUE                                                          ANH05130
C   SUM OVER NON-DEGENERATE MODES                                       ANH05140
      VALSA3=A00                                                        ANH05150
      VALSB3=A00                                                        ANH05160
      DO 112 JVIB=1,NVIB                                                ANH05170
      IF(NDEG(JVIB).NE.1) GO TO 112                                     ANH05180
      FRQJ=FREQ(JVIB)                                                   ANH05190
      FACT3=(FRQI*F3Q(IVIB,IVIB,JVIB)*RCONST)/(FRQJ**(1.5D+00))         ANH05200
      VALSA3=VALSA3+AXY(IAA,IAA,JVIB)*FACT3                             ANH05210
      VALSB3=VALSB3+AXY(IBB,IBB,JVIB)*FACT3                             ANH05220
  112 CONTINUE                                                          ANH05230
      ROTSA=-(VALSA1+VALSA2+VALSA3)*(ROTCM(1)*ROTCM(1)*TWO)/FRQI        ANH05240
      ROTSB=-(VALSB1+VALSB2+VALSB3)*(ROTCM(2)*ROTCM(2)*TWO)/FRQI        ANH05250
      RVCNST(IVIB,1)=ROTSA                                              ANH05260
      RVCNST(IVIB,2)=ROTSB                                              ANH05270
  113 CONTINUE                                                          ANH05280
C                                                                       ANH05290
C::::::::::::::::::::::::::                                             ANH05300
C:::FOR DEGENERATE MODES:::                                             ANH05310
C::::::::::::::::::::::::::                                             ANH05320
      DO 115 IVIB=1,NVIB                                                ANH05330
      IF(NDEG(IVIB).EQ.1) GO TO 115                                     ANH05340
      FRQI=FREQ(IVIB)                                                   ANH05350
      FRQI2=FRQI*FRQI                                                   ANH05360
      VALTA1=A00                                                        ANH05370
      VALTB1=A00                                                        ANH05380
      VALTA2=A00                                                        ANH05390
      VALTB2=A00                                                        ANH05400
      VALTA3=A00                                                        ANH05410
      VALTB3=A00                                                        ANH05420
C   SUM OVER NON-DEGENERATE MODES                                       ANH05430
      DO 114 JVIB=1,NVIB                                                ANH05440
      IF(NDEG(JVIB).NE.1) GO TO 114                                     ANH05450
      FRQJ=FREQ(JVIB)                                                   ANH05460
      FRQJ2=FRQJ*FRQJ                                                   ANH05470
      IF(DABS(FRQI-FRQJ).LE.CLIMIT) GO TO 205                           ANH05480
      FACT2=0.5D+00*(FRQI2*3.0D+00+FRQJ2)/(FRQI2-FRQJ2)                 ANH05490
      GO TO 206                                                         ANH05500
C////////////////////////////////////                                   ANH05510
C///CORIOLIS RESONANCE CORRECTIONS///                                   ANH05520
C////////////////////////////////////                                   ANH05530
  205 CONTINUE                                                          ANH05540
      FACT2=-(0.25D+00*(FRQI-FRQJ)**2)/(FRQJ*(FRQI+FRQJ))               ANH05550
      IF(DABS(ZETA(IVIB,JVIB,IAA)).GT.ZLIMIT.OR.                        ANH05560
     *   DABS(ZETA(IVIB,JVIB,IBB)).GT.ZLIMIT) THEN                      ANH05570
        NCORIO=NCORIO+1                                                 ANH05580
      END IF                                                            ANH05590
C///////////////////////////////                                        ANH05600
C///END OF CORIOLIS RESONANCE///                                        ANH05610
C///////////////////////////////                                        ANH05620
  206 CONTINUE                                                          ANH05630
      ZFACT=ZETA(IVIB,JVIB,IAA)**2+ZETA(IVIB,JVIB,IBB)**2               ANH05640
      VALTA2=VALTA2+FACT2*ZFACT                                         ANH05650
      VALTB2=VALTB2+FACT2*ZFACT                                         ANH05660
      FACT3=(FRQI*F3Q(IVIB,IVIB,JVIB)*RCONST)/(FRQJ**(1.5D+00))         ANH05670
      VALTA3=VALTA3+AXY(IAA,IAA,JVIB)*FACT3                             ANH05680
      VALTB3=VALTB3+AXY(IBB,IBB,JVIB)*FACT3                             ANH05690
  114 CONTINUE                                                          ANH05700
      ROTTA=-(VALTA1+VALTA2+VALTA3)*(ROTCM(1)*ROTCM(1)*TWO)/FRQI        ANH05710
      ROTTB=-(VALTB1+VALTB2+VALTB3)*(ROTCM(2)*ROTCM(2)*TWO)/FRQI        ANH05720
      RVCNST(IVIB,1)=ROTTA                                              ANH05730
      RVCNST(IVIB,2)=ROTTB                                              ANH05740
  115 CONTINUE                                                          ANH05750
C                                                                       ANH05760
C   CALCULATE THE ROTATIONAL CONSTANTS                                  ANH05770
C   WITH ZERO-POINT VIBRATIONAL CORRECTION                              ANH05780
      CALL ZERO(SALP,3)                                                 ANH05790
      WRITE(6,5)                                                        ANH05800
      DO 120 III=1,NVIB                                                 ANH05810
      IVIB=NFRQ(III)                                                    ANH05820
      ROTA=RVCNST(IVIB,1)                                               ANH05830
      ROTB=RVCNST(IVIB,2)                                               ANH05840
      SALP(1)=SALP(1)+ROTA                                              ANH05850
      SALP(2)=SALP(2)+ROTB                                              ANH05860
      WRITE(6,6) III,IVIB,ROTA,ROTB                                     ANH05870
  120 CONTINUE                                                          ANH05880
      WRITE(6,7)                                                        ANH05890
      DO 122 I=1,2                                                      ANH05900
      PB=ROTCM(I)-SALP(I)*HALF                                          ANH05910
      PT=CONST/PB                                                       ANH05920
      PA=PT*PARA                                                        ANH05930
      PC=CYCL/PT                                                        ANH05940
      WRITE(6,8) RAXIS(I),PT,PA,PB,PC                                   ANH05950
  122 CONTINUE                                                          ANH05960
C                                                                       ANH05970
      IF(NCOR.NE.0) THEN                                                ANH05980
        WRITE(6,9)                                                      ANH05990
        WRITE(6,10) CLIMIT                                              ANH06000
      END IF                                                            ANH06010
C                                                                       ANH06020
      RETURN                                                            ANH06030
      END                                                               ANH06040
      SUBROUTINE CENTLI(BXY,CXY,F3Q,EE)                                 ANH06050
C   THE CALCULATION OF CENTRIFUGAL DISTORTION CONSTANTS                 ANH06060
C   FOR A LINEAR POLYATOMIC MOLECULE                                    ANH06070
      IMPLICIT REAL*8 (A-H,O-Z)                                         ANH06080
      DIMENSION BXY(3,3,N3N),CXY(3,3,N3N)                               ANH06090
      DIMENSION F3Q(N3N,N3N,N3N),EE(N3N,N3N)                            ANH06100
      COMMON/VIB101/NATOM,N3N,NATRI,ILIN,NVIB                           ANH06110
      COMMON/VIB106/ROTAA(3),ROTGC(3),ROTCM(3),ROTMH(3)                 ANH06120
      COMMON/VIB108/IAXIS(3),NDEG(150),NDAB(150,5),IMAG(150)            ANH06130
      COMMON/VIB203/IOFF(150),IPRNT                                     ANH06140
      COMMON/VIB204/SQM(150),ROOT(150),FREQ(150)                        ANH06150
      COMMON/VIB205/IFREQ,NFRQ(150)                                     ANH06160
      COMMON/VIB207/TABCD(3,3,3,3),PABC(3,3,3)                          ANH06170
      COMMON/VIB208/TAABB(3,3),TABAB(3,3)                               ANH06180
      DATA WLIMIT / 1.0D+00 /                                           ANH06190
      DATA CL / 2.99792458D+00 /                                        ANH06200
      DATA A00,HALF,ONE,TWO / 0.0D+00 , 0.5D+00 , 1.0D+00 , 2.0D+00 /   ANH06210
    1 FORMAT(//,2X,' ******************************************'/       ANH06220
     1          2X,' ***THE CENTRIFUGAL DISTORTION CONSTANTS***'/       ANH06230
     2          2X,' ******************************************')       ANH06240
    2 FORMAT(//,2X,' F3(QI*QJ*QK) MATRIX (IN CM-1), IVIB = ',I5/)       ANH06250
    3 FORMAT(//,2X,' TAABB MATRIX FROM BXY MATRICES (IN CM-1)'/)        ANH06260
    4 FORMAT(//,2X,' TABAB MATRIX FROM BXY MATRICES (IN CM-1)'/)        ANH06270
    5 FORMAT(//,2X,' FOR A LINEAR POLYATOMIC MOLECULE'/                 ANH06280
     1          2X,' DJ = ',F20.10,' D-10 CM-1',                        ANH06290
     2          8X,' DJ = ',F20.10,' MHZ'/)                             ANH06300
    6 FORMAT(2X,' HJ = ',F20.10,' D-10 CM-1',                           ANH06310
     1          8X,' HJ = ',F20.10,' HZ'/)                              ANH06320
C                                                                       ANH06330
      WRITE(6,1)                                                        ANH06340
C   ELEMENTS OF F3Q ARE IN CM-1                                         ANH06350
      IF(IPRNT.LE.2) GO TO 201                                          ANH06360
      DO 102 III=1,NVIB                                                 ANH06370
      IVIB=NFRQ(III)                                                    ANH06380
      WRITE(6,2) IVIB                                                   ANH06390
      DO 101 I=1,NVIB                                                   ANH06400
      II=NFRQ(I)                                                        ANH06410
      DO 101 J=1,NVIB                                                   ANH06420
      JJ=NFRQ(J)                                                        ANH06430
      EE(I,J)=F3Q(II,JJ,IVIB)                                           ANH06440
  101 CONTINUE                                                          ANH06450
      CALL MATOUT(EE,N3N,N3N,NVIB,NVIB,6)                               ANH06460
  102 CONTINUE                                                          ANH06470
C                                                                       ANH06480
  201 CONTINUE                                                          ANH06490
C   FOR A LINEAR POLYATOMIC MOLECULES                                   ANH06500
C   QUARTIC CENTRIFUGAL CONSTANTS                                       ANH06510
C                                                                       ANH06520
C***********************************************************************ANH06530
C***FOLLOWING CODES ARE FOR ALTERNATE EXPRESSIONS                      *ANH06540
C***THEY ARE TESTED TO GIVE THE SAME RESULTS                           *ANH06550
C***********************************************************************ANH06560
C*    DO 103 I=1,3                                                     *ANH06570
C*    DO 103 J=1,3                                                     *ANH06580
C*    TAABB(I,J)=TABCD(I,I,J,J)                                        *ANH06590
C*    TABAB(I,J)=TABCD(I,J,I,J)                                        *ANH06600
C*103 CONTINUE                                                         *ANH06610
C*    IF(IPRNT.LE.2) GO TO 202                                         *ANH06620
C*    WRITE(6,3)                                                       *ANH06630
C*    CALL MATOUT(TAABB,3,3,3,3,6)                                     *ANH06640
C*    WRITE(6,4)                                                       *ANH06650
C*    CALL MATOUT(TABAB,3,3,3,3,6)                                     *ANH06660
C                                                                      *ANH06670
C*202 CONTINUE                                                         *ANH06680
C*    IXX=1                                                            *ANH06690
C*    DJ=-TAABB(IXX,IXX)*0.25D+00                                      *ANH06700
C*    DJCM=DJ*1.0D+10                                                  *ANH06710
C*    DJMH=DJ*CL*1.0D+4                                                *ANH06720
C*    WRITE(6,5) DJCM,DJMH                                             *ANH06730
C***********************************************************************ANH06740
C*    DJ=A00                                                           *ANH06750
C*    DO 104 IVIB=1,NVIB                                               *ANH06760
C*    IF(NDEG(IVIB).NE.1) GO TO 104                                    *ANH06770
C*    FRQI=FREQ(IVIB)                                                  *ANH06780
C*    IF(DABS(FRQI).LE.WLIMIT) GO TO 104                               *ANH06790
C*    FACT=BXY(IXX,IXX,IVIB)*BXY(IXX,IXX,IVIB)/FRQI                    *ANH06800
C*    DJ=DJ+FACT/TWO                                                   *ANH06810
C*104 CONTINUE                                                         *ANH06820
C*    DJCM=DJ*1.0D+10                                                  *ANH06830
C*    DJMH=DJ*CL*1.0D+4                                                *ANH06840
C*    WRITE(6,5) DJCM,DJMH                                             *ANH06850
C***********************************************************************ANH06860
C                                                                       ANH06870
      IXX=1                                                             ANH06880
      DJ=A00                                                            ANH06890
      DO 105 IVIB=1,NVIB                                                ANH06900
      IF(NDEG(IVIB).NE.1) GO TO 105                                     ANH06910
      FRQI=FREQ(IVIB)                                                   ANH06920
      IF(DABS(FRQI).LE.WLIMIT) GO TO 105                                ANH06930
      FACT=CXY(IXX,IXX,IVIB)*CXY(IXX,IXX,IVIB)*FRQI                     ANH06940
      DJ=DJ+FACT/TWO                                                    ANH06950
  105 CONTINUE                                                          ANH06960
      DJCM=DJ*1.0D+10                                                   ANH06970
      DJMH=DJ*CL*1.0D+4                                                 ANH06980
      WRITE(6,5) DJCM,DJMH                                              ANH06990
C                                                                       ANH07000
C   SEXTIC CENTRIFUGAL DISTORTION CONSTANTS                             ANH07010
C***********************************************************************ANH07020
C***FOLLWOING CODES ARE FOR AN ALTERNATE EXPRESSION                    *ANH07030
C***THEY ARE TESTED TO GIVE THE SAME RESUTLS                           *ANH07040
C***********************************************************************ANH07050
C*    IXX=1                                                            *ANH07060
C*    VALU1=(DJ*DJ*4.0D+00)/ROTCM(IXX)                                 *ANH07070
C*    VALU2=A00                                                        *ANH07080
C*    VALU3=A00                                                        *ANH07090
C*    ROTSQ2=ROTCM(IXX)*ROTCM(IXX)*TWO                                 *ANH07100
C*    DO 110 IVIB=1,NVIB                                               *ANH07110
C*    IF(NDEG(IVIB).NE.1) GO TO 110                                    *ANH07120
C*    FRQI=FREQ(IVIB)                                                  *ANH07130
C*    FRQI3=FRQI**3                                                    *ANH07140
C*    BI=BXY(IXX,IXX,IVIB)                                             *ANH07150
C*    VALU2=VALU2-ROTSQ2*BXY(IXX,IXX,IVIB)*BXY(IXX,IXX,IVIB)/FRQI3     *ANH07160
C*    DO 109 JVIB=1,NVIB                                               *ANH07170
C*    IF(NDEG(JVIB).NE.1) GO TO 109                                    *ANH07180
C*    FRQJ=FREQ(JVIB)                                                  *ANH07190
C*    BJ=BXY(IXX,IXX,JVIB)                                             *ANH07200
C*    DO 108 KVIB=1,NVIB                                               *ANH07210
C*    IF(NDEG(KVIB).NE.1) GO TO 108                                    *ANH07220
C*    FRQK=FREQ(KVIB)                                                  *ANH07230
C*    BK=BXY(IXX,IXX,KVIB)                                             *ANH07240
C*    FACT=BI*BJ*BK                                                    *ANH07250
C*    VALU3=VALU3-FACT*F3Q(IVIB,JVIB,KVIB)/(FRQI*FRQJ*FRQK)            *ANH07260
C*108 CONTINUE                                                         *ANH07270
C*109 CONTINUE                                                         *ANH07280
C*110 CONTINUE                                                         *ANH07290
C*    HJ=VALU1+VALU2+VALU3/6.0D+00                                     *ANH07300
C*    HJCM=HJ*1.0D+10                                                  *ANH07310
C*    HJHZ=HJ*CL*1.0D+10                                               *ANH07320
C*    WRITE(6,6) HJCM,HJHZ                                             *ANH07330
C***********************************************************************ANH07340
C                                                                       ANH07350
      IXX=1                                                             ANH07360
      VALU1=(DJ*DJ*4.0D+00)/ROTCM(IXX)                                  ANH07370
      VALU2=A00                                                         ANH07380
      VALU3=A00                                                         ANH07390
      ROTSQ2=ROTCM(IXX)*ROTCM(IXX)*TWO                                  ANH07400
      DO 115 IVIB=1,NVIB                                                ANH07410
      IF(NDEG(IVIB).NE.1) GO TO 115                                     ANH07420
      FRQI=FREQ(IVIB)                                                   ANH07430
      CI=CXY(IXX,IXX,IVIB)                                              ANH07440
      VALU2=VALU2-ROTSQ2*CXY(IXX,IXX,IVIB)*CXY(IXX,IXX,IVIB)/FRQI       ANH07450
      DO 114 JVIB=1,NVIB                                                ANH07460
      IF(NDEG(JVIB).NE.1) GO TO 114                                     ANH07470
      CJ=CXY(IXX,IXX,JVIB)                                              ANH07480
      DO 113 KVIB=1,NVIB                                                ANH07490
      IF(NDEG(KVIB).NE.1) GO TO 113                                     ANH07500
      CK=CXY(IXX,IXX,KVIB)                                              ANH07510
      FACT=CI*CJ*CK                                                     ANH07520
      VALU3=VALU3+FACT*F3Q(IVIB,JVIB,KVIB)                              ANH07530
  113 CONTINUE                                                          ANH07540
  114 CONTINUE                                                          ANH07550
  115 CONTINUE                                                          ANH07560
      HJ=VALU1+VALU2+VALU3/6.0D+00                                      ANH07570
C                                                                       ANH07580
      HJCM=HJ*1.0D+10                                                   ANH07590
      HJHZ=HJ*CL*1.0D+10                                                ANH07600
      WRITE(6,6) HJCM,HJHZ                                              ANH07610
      RETURN                                                            ANH07620
      END                                                               ANH07630
      SUBROUTINE VLTYLI(F3Q,F4Q,CLTYP)                                  ANH07640
C   THE CALCULATION OF VIBRATIONAL L-TYPE DOUBLING CONSTANTS            ANH07650
C   FOR A LINEAR POLYATOMIC MOLECULE                                    ANH07660
      IMPLICIT REAL*8 (A-H,O-Z)                                         ANH07670
      DIMENSION CLTYP(N3N,N3N,3)                                        ANH07680
      DIMENSION F3Q(N3N,N3N,N3N),F4Q(N3N,N3N,N3N,N3N)                   ANH07690
      COMMON/VIB101/NATOM,N3N,NATRI,ILIN,NVIB                           ANH07700
      COMMON/VIB108/IAXIS(3),NDEG(150),NDAB(150,5),IMAG(150)            ANH07710
      COMMON/VIB109/CLIMIT,FLIM1,FLIM2                                  ANH07720
      COMMON/VIB203/IOFF(150),IPRNT                                     ANH07730
      COMMON/VIB204/SQM(150),ROOT(150),FREQ(150)                        ANH07740
      COMMON/VIB205/IFREQ,NFRQ(150)                                     ANH07750
      DATA A00,HALF,ONE,TWO / 0.0D+00 , 0.5D+00 , 1.0D+00 , 2.0D+00 /   ANH07760
      DATA WLIMIT / 1.0D+00 /                                           ANH07770
      DATA DLIMIT / 1.0D-01 /                                           ANH07780
      DATA PLIMIT / 1.0D-05 /                                           ANH07790
      DATA PI,CL / 3.1415926536D+00 , 2.99792458D+00 /                  ANH07800
    1 FORMAT(//,2X,' ***********************************************'/  ANH07810
     1          2X,' ***THE VIBRATIONAL L-TYPE DOUBLING CONSTANTS***'/  ANH07820
     2          2X,' ***********************************************')  ANH07830
    2 FORMAT(2X,3I5,3X,5F15.5)                                          ANH07840
    3 FORMAT(//,2X,' ::::::::::::::::::::::::::::::::::::::::::::::::'/ ANH07850
     1          2X,' :::ANHARMONIC RESONANCE (TYPE II) IS OBSERVED:::'/ ANH07860
     2          2X,' ::::::::::::::::::::::::::::::::::::::::::::::::'/)ANH07870
    4 FORMAT(2X,' RESONANCE LIMIT IS SET TO ',F10.3,' CM-1'/)           ANH07880
    5 FORMAT(//,2X,' L-TYPE DOUBLING CONSTANTS',62X,' TOTAL',9X,' TOTAL'ANH07890
     1        /,4X,' QI',5X,' QJ',8X,' FREQ(I)',5X,' FREQ(J)',          ANH07900
     2         10X,' CUBIC',12X,' QUARTIC',10X,' IN CM-1',7X,' IN MHZ'/)ANH07910
    6 FORMAT(2X,I5,3X,I5,3X,2(F13.3),3X,5(F15.5,3X))                    ANH07920
C                                                                       ANH07930
      WRITE(6,1)                                                        ANH07940
C                                                                       ANH07950
C   NOTE : THE MOLECULAR AXIS IS CHOSEN AS THE Z AXIS                   ANH07960
      IXX=1                                                             ANH07970
      IYY=2                                                             ANH07980
      IZZ=3                                                             ANH07990
C                                                                       ANH08000
      IANH2=0                                                           ANH08010
      IVLTYP=0                                                          ANH08020
      CALL ZERO(CLTYP,N3N*N3N*3)                                        ANH08030
C                                                                       ANH08040
C:::::::::::::::::::::::::::::::                                        ANH08050
C:::ONLY FOR DEGENERATE MODES:::                                        ANH08060
C:::::::::::::::::::::::::::::::                                        ANH08070
      DO 105 IVIB=1,NVIB                                                ANH08080
      IF(NDEG(IVIB).EQ.1) GO TO 105                                     ANH08090
      FRQI=FREQ(IVIB)                                                   ANH08100
      IF(DABS(FRQI).LE.WLIMIT) GO TO 105                                ANH08110
      FRQI2=FRQI*FRQI                                                   ANH08120
      IVA=NDAB(IVIB,1)                                                  ANH08130
      IVB=NDAB(IVIB,2)                                                  ANH08140
      DO 104 JVIB=1,NVIB                                                ANH08150
      IF(NDEG(JVIB).EQ.1) GO TO 104                                     ANH08160
      FRQJ=FREQ(JVIB)                                                   ANH08170
      IF(DABS(FRQJ).LE.WLIMIT) GO TO 104                                ANH08180
      FRQJ2=FRQJ*FRQJ                                                   ANH08190
      JVA=NDAB(JVIB,1)                                                  ANH08200
      JVB=NDAB(JVIB,2)                                                  ANH08210
      IF(IVIB.EQ.JVIB) GO TO 104                                        ANH08220
      IF(DABS(FRQI-FRQJ).LE.DLIMIT) GO TO 104                           ANH08230
      IVLTYP=IVLTYP+1                                                   ANH08240
C                                                                       ANH08250
      VALU1=(F4Q(IVA,IVA,JVA,JVA)-F4Q(IVA,IVA,JVB,JVB)                  ANH08260
     *      +F4Q(JVA,JVA,IVA,IVA)-F4Q(JVA,JVA,IVB,IVB))/8.0D+00         ANH08270
      VALU1=DABS(VALU1)                                                 ANH08280
C                                                                       ANH08290
C   SUM OVER NON-DEGENERATE MODES                                       ANH08300
      VALU2=A00                                                         ANH08310
      DO 103 KVIB=1,NVIB                                                ANH08320
      IF(NDEG(KVIB).NE.1) GO TO 103                                     ANH08330
      FRQK=FREQ(KVIB)                                                   ANH08340
      IF(DABS(FRQK).LE.WLIMIT) GO TO 103                                ANH08350
      FRQK2=FRQK*FRQK                                                   ANH08360
      DNUM=(FRQK2-FRQI2-FRQJ2)*FRQK                                     ANH08370
      FAC=(F3Q(IVA,JVA,KVIB)+F3Q(IVA,JVB,KVIB)                          ANH08380
     *    +F3Q(IVB,JVA,KVIB)+F3Q(IVB,JVB,KVIB))*0.5D+00                 ANH08390
      IF(DABS(FAC).LE.WLIMIT) GO TO 103                                 ANH08400
      FACT=FAC*FAC                                                      ANH08410
      WPPP=FRQI+FRQJ+FRQK                                               ANH08420
      WPPM=FRQI+FRQJ-FRQK                                               ANH08430
      WPMP=FRQI-FRQJ+FRQK                                               ANH08440
      WPMM=FRQI-FRQJ-FRQK                                               ANH08450
      DIJK=WPPP*WPPM*WPMP*WPMM                                          ANH08460
      IF(DABS(WPPM).LE.FLIM2) GO TO 201                                 ANH08470
      IF(DABS(WPMP).LE.FLIM2) GO TO 202                                 ANH08480
      IF(DABS(WPMM).LE.FLIM2) GO TO 203                                 ANH08490
      VALU2=VALU2+(FACT*DNUM)/(DIJK)                                    ANH08500
      GO TO 103                                                         ANH08510
C/////////////////////////////////                                      ANH08520
C///FERMI RESONANCE CORRECTIONS///                                      ANH08530
C/////////////////////////////////                                      ANH08540
C   WI+WJ=WK                                                            ANH08550
  201 CONTINUE                                                          ANH08560
      DNUM=ONE/WPPP+ONE/WPMP-ONE/WPMM                                   ANH08570
      IF(DABS(FAC).GT.PLIMIT) THEN                                      ANH08580
        IANH2=IANH2+1                                                   ANH08590
      END IF                                                            ANH08600
      VALU2=VALU2+(FACT*DNUM)/4.0D+00                                   ANH08610
      GO TO 103                                                         ANH08620
C   WI+WK=WJ                                                            ANH08630
  202 CONTINUE                                                          ANH08640
      DNUM=ONE/WPPP-ONE/WPPM-ONE/WPMM                                   ANH08650
      IF(DABS(FAC).GT.PLIMIT) THEN                                      ANH08660
        IANH2=IANH2+1                                                   ANH08670
      END IF                                                            ANH08680
      VALU2=VALU2+(FACT*DNUM)/4.0D+00                                   ANH08690
      GO TO 103                                                         ANH08700
C   WJ+WK=WI                                                            ANH08710
  203 CONTINUE                                                          ANH08720
      DNUM=ONE/WPPP-ONE/WPPM+ONE/WPMP                                   ANH08730
      IF(DABS(FAC).GT.PLIMIT) THEN                                      ANH08740
        IANH2=IANH2+1                                                   ANH08750
      END IF                                                            ANH08760
      VALU2=VALU2+(FACT*DNUM)/4.0D+00                                   ANH08770
C////////////////////////////                                           ANH08780
C///END OF FERMI RESONANCE///                                           ANH08790
C////////////////////////////                                           ANH08800
  103 CONTINUE                                                          ANH08810
      VALU2=-VALU2*0.5D+00                                              ANH08820
      VALUT=VALU1+VALU2                                                 ANH08830
      CLTYP(IVIB,JVIB,1)=VALU1                                          ANH08840
      CLTYP(IVIB,JVIB,2)=VALU2                                          ANH08850
      CLTYP(IVIB,JVIB,3)=VALUT                                          ANH08860
      IF(IPRNT.LE.2) GO TO 104                                          ANH08870
      WRITE(6,2) IVIB,JVIB,IVLTYP,FRQI,FRQJ,VALU1,VALU2,VALUT           ANH08880
  104 CONTINUE                                                          ANH08890
  105 CONTINUE                                                          ANH08900
C                                                                       ANH08910
      IF(IANH2.NE.0) THEN                                               ANH08920
        WRITE(6,3)                                                      ANH08930
        WRITE(6,4) FLIM2                                                ANH08940
      END IF                                                            ANH08950
C                                                                       ANH08960
      WRITE(6,5)                                                        ANH08970
      DO 107 II=1,NVIB                                                  ANH08980
      IVIB=NFRQ(II)                                                     ANH08990
      IF(NDEG(IVIB).EQ.1) GO TO 107                                     ANH09000
      FRQI=FREQ(IVIB)                                                   ANH09010
      DO 106 JJ=1,NVIB                                                  ANH09020
      JVIB=NFRQ(JJ)                                                     ANH09030
      IF(NDEG(JVIB).EQ.1) GO TO 106                                     ANH09040
      FRQJ=FREQ(JVIB)                                                   ANH09050
      IF(IVIB.EQ.JVIB) GO TO 106                                        ANH09060
      IF(DABS(FRQI-FRQJ).LE.DLIMIT) GO TO 106                           ANH09070
      VALU1=CLTYP(IVIB,JVIB,1)                                          ANH09080
      VALU2=CLTYP(IVIB,JVIB,2)                                          ANH09090
      VALUT=CLTYP(IVIB,JVIB,3)                                          ANH09100
      VALUM=VALUT*CL*1.0D+04                                            ANH09110
      WRITE(6,6) II,JJ,FRQI,FRQJ,VALU1,VALU2,VALUT,VALUM                ANH09120
  106 CONTINUE                                                          ANH09130
  107 CONTINUE                                                          ANH09140
C                                                                       ANH09150
      RETURN                                                            ANH09160
      END                                                               ANH09170
      SUBROUTINE LINEAR(ZETA,F3Q,F4Q,XIJ,GIJ)                           ANH09180
C   THE VIBRATIONAL ANHARMONICITY CONSTANTS                             ANH09190
C   FOR A LINEAR POLYATOMIC MOLECULE                                    ANH09200
      IMPLICIT REAL*8 (A-H,O-Z)                                         ANH09210
      DIMENSION F3Q(N3N,N3N,N3N),F4Q(N3N,N3N,N3N,N3N),ZETA(N3N,N3N,3)   ANH09220
      DIMENSION XIJ(N3N,N3N),GIJ(N3N,N3N)                               ANH09230
      COMMON/VIB101/NATOM,N3N,NATRI,ILIN,NVIB                           ANH09240
      COMMON/VIB106/ROTAA(3),ROTGC(3),ROTCM(3),ROTMH(3)                 ANH09250
      COMMON/VIB108/IAXIS(3),NDEG(150),NDAB(150,5),IMAG(150)            ANH09260
      COMMON/VIB109/CLIMIT,FLIM1,FLIM2                                  ANH09270
      COMMON/VIB111/ANHTOT,HZPVE,AZPVE
      COMMON/VIB203/IOFF(150),IPRNT                                     ANH09280
      COMMON/VIB204/SQM(150),ROOT(150),FREQ(150)                        ANH09290
      COMMON/VIB205/IFREQ,NFRQ(150)                                     ANH09300
      DATA A00,ONE / 0.0D+00 , 1.0D+00 /                                ANH09310
      DATA WLIMIT / 1.0D+00 /                                           ANH09320
      DATA DLIMIT / 1.0D-01 /                                           ANH09330
      DATA PLIMIT / 1.0D-05 /                                           ANH09340
    1 FORMAT(//,2X,' *********************************************'/    ANH09350
     1          2X,' ***THE VIBRATIONAL ANHARMONICITY CONSTANTS***'/    ANH09360
     2          2X,' *********************************************'/)   ANH09370
    2 FORMAT(//,2X,' ANHARMONIC CONSTANTS IN CM-1'/                     ANH09380
     1          2X,' IVIB',3X,' JVIB',8X,' FREQ(I)',5X,' FREQ(J)',      ANH09390
     2          7X,' VALU1',7X,' VALU2',7X,' VALU3',7X,' VALU4',        ANH09400
     3          7X,' XIJ',9X,' GIJ'/)                                   ANH09410
    3 FORMAT(2X,I5,3X,I5,3X,8(F13.3))                                   ANH09420
    4 FORMAT(//,2X,' ANHARMONIC CONSTANTS IN CM-1'/                     ANH09430
     1          4X,' QI',5X,' QJ',8X,' FREQ(I)',5X,' FREQ(J)',          ANH09440
     2          7X,' XIJ',9X,' GIJ'/)                                   ANH09450
    5 FORMAT(2X,I5,3X,I5,3X,4(F13.3))                                   ANH09460
    6 FORMAT(//,2X,' THEORETICAL FREQUENCIES IN CM-1'/                  ANH09470
     1          2X,' NO.',9X,' HARMONIC',14X,' DELTA',                  ANH09480
     2         13X,' FUNDAMENTAL'/)                                     ANH09490
    7 FORMAT(2X,I5,4X,F13.3,7X,F13.3,12X,F13.3)                         ANH09500
    8 FORMAT(//,2X,' :::::::::::::::::::::::::::::::::::::::::::::::'/  ANH09510
     1          2X,' :::ANHARMONIC RESONANCE (TYPE I) IS OBSERVED:::'/  ANH09520
     2          2X,' :::::::::::::::::::::::::::::::::::::::::::::::'/) ANH09530
    9 FORMAT(//,2X,' ::::::::::::::::::::::::::::::::::::::::::::::::'/ ANH09540
     1          2X,' :::ANHARMONIC RESONANCE (TYPE II) IS OBSERVED:::'/ ANH09550
     2          2X,' ::::::::::::::::::::::::::::::::::::::::::::::::'/)ANH09560
   10 FORMAT(2X,' RESONANCE LIMIT IS SET TO ',F10.3,' CM-1'/)           ANH09570
C                                                                       ANH09580
C:::::::::::::::::::::::::::::::::::::::::::::::::::                    ANH09590
C:::CALCULATE VIBRATIONAL ANHARMONICITY CONSTANTS:::                    ANH09600
C:::::::::::::::::::::::::::::::::::::::::::::::::::                    ANH09610
      IANH1=0                                                           ANH09620
      IANH2=0                                                           ANH09630
C                                                                       ANH09640
      WRITE(6,1)                                                        ANH09650
      IF(IPRNT.GT.2)                                                    ANH09660
     *WRITE(6,2)                                                        ANH09670
      IAA=IAXIS(1)                                                      ANH09680
      IBB=IAXIS(2)                                                      ANH09690
      ICC=IAXIS(3)                                                      ANH09700
C                                                                       ANH09710
      DO 120 IVIB=1,NVIB                                                ANH09720
      FRQI=FREQ(IVIB)                                                   ANH09730
      FRQI2=FRQI*FRQI                                                   ANH09740
      IVA=NDAB(IVIB,1)                                                  ANH09750
      IVB=NDAB(IVIB,2)                                                  ANH09760
      DO 120 JVIB=1,NVIB                                                ANH09770
      FRQJ=FREQ(JVIB)                                                   ANH09780
      FRQJ2=FRQJ*FRQJ                                                   ANH09790
      JVA=NDAB(JVIB,1)                                                  ANH09800
      JVB=NDAB(JVIB,2)                                                  ANH09810
C                                                                       ANH09820
      VALU1=A00                                                         ANH09830
      VALU2=A00                                                         ANH09840
      VALU3=A00                                                         ANH09850
      VALU4=A00                                                         ANH09860
      VALU5=A00                                                         ANH09870
      GVAL1=A00                                                         ANH09880
      GVAL2=A00                                                         ANH09890
      XRS=A00                                                           ANH09900
      GRS=A00                                                           ANH09910
      IF(DABS(FRQI).LE.WLIMIT) GO TO 230                                ANH09920
      IF(DABS(FRQJ).LE.WLIMIT) GO TO 230                                ANH09930
      IF(IVIB.NE.JVIB) GO TO 210                                        ANH09940
C                                                                       ANH09950
C+++++++++++++++++++++++++++++++++++++++++                              ANH09960
C+++DIAGONAL CONSTANTS (NON-DEGENERATE)+++                              ANH09970
C+++++++++++++++++++++++++++++++++++++++++                              ANH09980
      IF(NDEG(IVIB).NE.1) GO TO 205                                     ANH09990
      VALU1=F4Q(IVIB,IVIB,IVIB,IVIB)                                    ANH10000
      DO 101 KVIB=1,NVIB                                                ANH10010
C   SUM OVER ONLY NON-DEGENERATE MODES                                  ANH10020
      FRQK=FREQ(KVIB)                                                   ANH10030
      IF(DABS(FRQK).LE.WLIMIT) GO TO 101                                ANH10040
      IF(NDEG(KVIB).NE.1) GO TO 101                                     ANH10050
      FRQK2=FRQK*FRQK                                                   ANH10060
      FAC=F3Q(IVIB,IVIB,KVIB)                                           ANH10070
      FACT=FAC*FAC                                                      ANH10080
      IF(DABS(FRQI+FRQI-FRQK).LE.FLIM1) GO TO 201                       ANH10090
      DNUM=FRQI2*8.0D+00-FRQK2*3.0D+00                                  ANH10100
      DENM=FRQI2*4.0D+00-FRQK2                                          ANH10110
      VALU2=VALU2+(DNUM*FACT)/(FRQK*DENM)                               ANH10120
C***********************************************************************ANH10130
C***FOLLOWING CODES ARE FOR AN ALTERNATE EXPRESSION                    *ANH10140
C***THEY ARE TESTED TO GIVE THE SAME RESULTS                           *ANH10150
C***********************************************************************ANH10160
C*    DNUM=ONE/((FRQI+FRQI+FRQK)*2.0D+00)-ONE/((FRQI+FRQI-FRQK)*2.0D+00)ANH10170
C*   *    +2.0D+00/FRQK                                                *ANH10180
C     VALU2=VALU2+(DNUM*FACT)                                          *ANH10190
C***********************************************************************ANH10200
      GO TO 101                                                         ANH10210
C/////////////////////////////////                                      ANH10220
C///FERMI RESONANCE CORRECTIONS///                                      ANH10230
C/////////////////////////////////                                      ANH10240
  201 CONTINUE                                                          ANH10250
      DNUM=ONE/((FRQI+FRQI+FRQK)*2.0D+00)+2.0D+00/FRQK                  ANH10260
      IF(DABS(FAC).GT.PLIMIT) THEN                                      ANH10270
        IANH1=IANH1+1                                                   ANH10280
      END IF                                                            ANH10290
      VALU2=VALU2+FACT*DNUM                                             ANH10300
C////////////////////////////                                           ANH10310
C///END OF FERMI RESONANCE///                                           ANH10320
C////////////////////////////                                           ANH10330
  101 CONTINUE                                                          ANH10340
      XRS=(VALU1-VALU2)/16.0D+00                                        ANH10350
      GO TO 230                                                         ANH10360
C                                                                       ANH10370
C+++++++++++++++++++++++++++++++++++++                                  ANH10380
C+++DIAGONAL CONSTANTS (DEGENERATE)+++                                  ANH10390
C+++++++++++++++++++++++++++++++++++++                                  ANH10400
  205 CONTINUE                                                          ANH10410
      VALU1=F4Q(IVIB,IVIB,IVIB,IVIB)                                    ANH10420
      GVAL1=-VALU1                                                      ANH10430
      DO 102 KVIB=1,NVIB                                                ANH10440
C   SUM OVER NON-DEGENERATE MODES                                       ANH10450
      FRQK=FREQ(KVIB)                                                   ANH10460
      IF(DABS(FRQK).LE.WLIMIT) GO TO 102                                ANH10470
      IF(NDEG(KVIB).NE.1) GO TO 102                                     ANH10480
      FRQK2=FRQK*FRQK                                                   ANH10490
      FAC=F3Q(IVIB,IVIB,KVIB)                                           ANH10500
      FACT=FAC*FAC                                                      ANH10510
      IF(DABS(FRQI+FRQI-FRQK).LE.FLIM1) GO TO 206                       ANH10520
      DNUM=FRQI2*8.0D+00-FRQK2*3.0D+00                                  ANH10530
      DENM=FRQI2*4.0D+00-FRQK2                                          ANH10540
      VALU2=VALU2+(DNUM*FACT)/(FRQK*DENM)                               ANH10550
      GVAL2=GVAL2+(FACT*FRQK)/DENM                                      ANH10560
      GO TO 102                                                         ANH10570
C/////////////////////////////////                                      ANH10580
C///FERMI RESONANCE CORRECTIONS///                                      ANH10590
C/////////////////////////////////                                      ANH10600
  206 CONTINUE                                                          ANH10610
      DNUM=ONE/((FRQI+FRQI+FRQK)*2.0D+00)+2.0D+00/FRQK                  ANH10620
      IF(DABS(FAC).GT.PLIMIT) THEN                                      ANH10630
        IANH1=IANH1+1                                                   ANH10640
      END IF                                                            ANH10650
      VALU2=VALU2+FACT*DNUM                                             ANH10660
      GNUM=ONE/((FRQI+FRQI+FRQK)*2.0D+00)                               ANH10670
      GVAL2=GVAL2-FACT*GNUM                                             ANH10680
C////////////////////////////                                           ANH10690
C///END OF FERMI RESONANCE///                                           ANH10700
C////////////////////////////                                           ANH10710
  102 CONTINUE                                                          ANH10720
      XRS=(VALU1-VALU2)/16.0D+00                                        ANH10730
      GRS=(GVAL1/3.0D+00-GVAL2)/16.0D+00                                ANH10740
      GO TO 230                                                         ANH10750
C                                                                       ANH10760
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++         ANH10770
C+++OFF-DIAGONAL CONSTANTS (NON-DEGENERATE & NON-DEGENERATE)+++         ANH10780
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++         ANH10790
  210 CONTINUE                                                          ANH10800
      IF(NDEG(IVIB).NE.1) GO TO 215                                     ANH10810
      IF(NDEG(JVIB).NE.1) GO TO 215                                     ANH10820
      VALU1=F4Q(IVIB,IVIB,JVIB,JVIB)                                    ANH10830
      DO 103 KVIB=1,NVIB                                                ANH10840
C   SUM OVER ONLY NON-DEGENERATE MODES                                  ANH10850
      FRQK=FREQ(KVIB)                                                   ANH10860
      IF(DABS(FRQK).LE.WLIMIT) GO TO 103                                ANH10870
      IF(NDEG(KVIB).NE.1) GO TO 103                                     ANH10880
      FRQK2=FRQK*FRQK                                                   ANH10890
      WPPP=FRQI+FRQJ+FRQK                                               ANH10900
      WPPM=FRQI+FRQJ-FRQK                                               ANH10910
      WPMP=FRQI-FRQJ+FRQK                                               ANH10920
      WPMM=FRQI-FRQJ-FRQK                                               ANH10930
      DIJK=WPPP*WPPM*WPMP*WPMM                                          ANH10940
      DNUM2=F3Q(IVIB,IVIB,KVIB)*F3Q(KVIB,JVIB,JVIB)                     ANH10950
      VALU2=VALU2+DNUM2/FRQK                                            ANH10960
      FAC=F3Q(IVIB,JVIB,KVIB)                                           ANH10970
      FACT=FAC*FAC                                                      ANH10980
      IF(DABS(WPPM).LE.FLIM2) GO TO 211                                 ANH10990
      IF(DABS(WPMP).LE.FLIM2) GO TO 212                                 ANH11000
      IF(DABS(WPMM).LE.FLIM2) GO TO 213                                 ANH11010
      DNUM3=(FRQK2-FRQI2-FRQJ2)*FRQK                                    ANH11020
      VALU3=VALU3+(DNUM3*FACT)/DIJK                                     ANH11030
C***********************************************************************ANH11040
C***FOLLOWING CODES ARE FOR AN ALTERNATIVE EXPRESSION                  *ANH11050
C***THEY ARE TESTED TO GIVE THE SAME RESULTS                           *ANH11060
C***********************************************************************ANH11070
C*    DNUM3=ONE/WPPP-ONE/WPPM+ONE/WPMP-ONE/WPMM                        *ANH11080
C*    VALU3=VALU3+(DNUM3*FACT)/4.0D+00                                 *ANH11090
C***********************************************************************ANH11100
      GO TO 103                                                         ANH11110
C/////////////////////////////////                                      ANH11120
C///FERMI RESONANCE CORRECTIONS///                                      ANH11130
C/////////////////////////////////                                      ANH11140
C   WI+WJ=WK                                                            ANH11150
  211 CONTINUE                                                          ANH11160
      DNUM3=ONE/WPPP+ONE/WPMP-ONE/WPMM                                  ANH11170
      IF(DABS(FAC).GT.PLIMIT) THEN                                      ANH11180
        IANH2=IANH2+1                                                   ANH11190
      END IF                                                            ANH11200
      VALU3=VALU3+(DNUM3*FACT)/4.0D+00                                  ANH11210
      GO TO 103                                                         ANH11220
C   WI+WK=WJ                                                            ANH11230
  212 CONTINUE                                                          ANH11240
      DNUM3=ONE/WPPP-ONE/WPPM-ONE/WPMM                                  ANH11250
      IF(DABS(FAC).GT.PLIMIT) THEN                                      ANH11260
        IANH2=IANH2+1                                                   ANH11270
      END IF                                                            ANH11280
      VALU3=VALU3+(DNUM3*FACT)/4.0D+00                                  ANH11290
      GO TO 103                                                         ANH11300
C   WJ+WK=WI                                                            ANH11310
  213 CONTINUE                                                          ANH11320
      DNUM3=ONE/WPPP-ONE/WPPM+ONE/WPMP                                  ANH11330
      IF(DABS(FAC).GT.PLIMIT) THEN                                      ANH11340
        IANH2=IANH2+1                                                   ANH11350
      END IF                                                            ANH11360
      VALU3=VALU3+(DNUM3*FACT)/4.0D+00                                  ANH11370
C////////////////////////////                                           ANH11380
C///END OF FERMI RESONANCE///                                           ANH11390
C////////////////////////////                                           ANH11400
  103 CONTINUE                                                          ANH11410
      XRS=(VALU1-VALU2-VALU3*2.0D+00)/4.0D+00                           ANH11420
      GO TO 230                                                         ANH11430
C                                                                       ANH11440
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++             ANH11450
C+++OFF-DIAGONAL CONSTANTS (DEGENERATE & NON-DEGENERATE)+++             ANH11460
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++             ANH11470
C   NOTE : RELATIONSHIPS AMONG CUBIC FORCE CONSTANTS                    ANH11480
C     (S,TA,TA)=(S,TB,TB)                                               ANH11490
C     (S,TA,TB)=0                                                       ANH11500
C     (S,TA,TA')=(S,TB,TB')                                             ANH11510
C     (S,TA,TB')=0                                                      ANH11520
  215 CONTINUE                                                          ANH11530
      IF(NDEG(IVIB).NE.1.AND.NDEG(JVIB).NE.1) GO TO 225                 ANH11540
      VALU1=F4Q(IVIB,IVIB,JVIB,JVIB)                                    ANH11550
C   SUM OVER NON-DEGENERATE MODES                                       ANH11560
      DO 105 KVIB=1,NVIB                                                ANH11570
      FRQK=FREQ(KVIB)                                                   ANH11580
      IF(DABS(FRQK).LE.WLIMIT) GO TO 105                                ANH11590
      IF(NDEG(KVIB).NE.1) GO TO 105                                     ANH11600
      DNUM2=F3Q(IVIB,IVIB,KVIB)*F3Q(KVIB,JVIB,JVIB)                     ANH11610
      VALU2=VALU2+DNUM2/FRQK                                            ANH11620
  105 CONTINUE                                                          ANH11630
C   SUM OVER DEGENERATE MODES                                           ANH11640
      DO 106 KVIB=1,NVIB                                                ANH11650
      FRQK=FREQ(KVIB)                                                   ANH11660
      IF(DABS(FRQK).LE.WLIMIT) GO TO 106                                ANH11670
      FRQK2=FRQK*FRQK                                                   ANH11680
      IF(NDEG(KVIB).EQ.1) GO TO 106                                     ANH11690
      WPPP=FRQI+FRQJ+FRQK                                               ANH11700
      WPPM=FRQI+FRQJ-FRQK                                               ANH11710
      WPMP=FRQI-FRQJ+FRQK                                               ANH11720
      WPMM=FRQI-FRQJ-FRQK                                               ANH11730
      DIJK=WPPP*WPPM*WPMP*WPMM                                          ANH11740
      DNUM3=(FRQK2-FRQI2-FRQJ2)*FRQK                                    ANH11750
      FAC=F3Q(IVIB,JVIB,KVIB)                                           ANH11760
      FACT=FAC*FAC                                                      ANH11770
      IF(DABS(WPPM).LE.FLIM2) GO TO 217                                 ANH11780
      IF(DABS(WPMP).LE.FLIM2) GO TO 218                                 ANH11790
      IF(DABS(WPMM).LE.FLIM2) GO TO 219                                 ANH11800
      VALU3=VALU3+(DNUM3*FACT)/DIJK                                     ANH11810
      GO TO 106                                                         ANH11820
C/////////////////////////////////                                      ANH11830
C///FERMI RESONANCE CORRECTIONS///                                      ANH11840
C/////////////////////////////////                                      ANH11850
C   WI+WJ=WK                                                            ANH11860
  217 CONTINUE                                                          ANH11870
      DNUM3=ONE/WPPP+ONE/WPMP-ONE/WPMM                                  ANH11880
      IF(DABS(FAC).GT.PLIMIT) THEN                                      ANH11890
        IANH2=IANH2+1                                                   ANH11900
      END IF                                                            ANH11910
      VALU3=VALU3+(DNUM3*FACT)/4.0D+00                                  ANH11920
      GO TO 106                                                         ANH11930
C   WI+WK=WJ                                                            ANH11940
  218 CONTINUE                                                          ANH11950
      DNUM3=ONE/WPPP-ONE/WPPM-ONE/WPMM                                  ANH11960
      IF(DABS(FAC).GT.PLIMIT) THEN                                      ANH11970
        IANH2=IANH2+1                                                   ANH11980
      END IF                                                            ANH11990
      VALU3=VALU3+(DNUM3*FACT)/4.0D+00                                  ANH12000
      GO TO 106                                                         ANH12010
C   WJ+WK=WI                                                            ANH12020
  219 CONTINUE                                                          ANH12030
      DNUM3=ONE/WPPP-ONE/WPPM+ONE/WPMP                                  ANH12040
      IF(DABS(FAC).GT.PLIMIT) THEN                                      ANH12050
        IANH2=IANH2+1                                                   ANH12060
      END IF                                                            ANH12070
      VALU3=VALU3+(DNUM3*FACT)/4.0D+00                                  ANH12080
C////////////////////////////                                           ANH12090
C///END OF FERMI RESONANCE///                                           ANH12100
C////////////////////////////                                           ANH12110
  106 CONTINUE                                                          ANH12120
      ZFACT=(ZETA(IVIB,JVIB,IAA)**2+ZETA(IVIB,JVIB,IBB)**2)*ROTCM(IBB)  ANH12130
      VALU4=(FRQI/FRQJ+FRQJ/FRQI)*ZFACT                                 ANH12140
      XRS=(VALU1-VALU2-VALU3*2.0D+00)/4.0D+00+VALU4                     ANH12150
      GO TO 230                                                         ANH12160
C                                                                       ANH12170
C++++++++++++++++++++++++++++++++++++++++++++++++++++++                 ANH12180
C+++OFF-DIAGONAL CONSTANTS (DEGENERATE & DEGENERATE)+++                 ANH12190
C++++++++++++++++++++++++++++++++++++++++++++++++++++++                 ANH12200
  225 CONTINUE                                                          ANH12210
      IF(DABS(FRQI-FRQJ).LE.DLIMIT) GO TO 230                           ANH12220
      VALU1=(F4Q(IVA,IVA,JVA,JVA)+F4Q(IVA,IVA,JVB,JVB)                  ANH12230
     *      +F4Q(JVA,JVA,IVA,IVA)+F4Q(JVA,JVA,IVB,IVB))*0.25D+00        ANH12240
      DO 109 KVIB=1,NVIB                                                ANH12250
C   SUM OVER NON-DEGENERATE MODES                                       ANH12260
      FRQK=FREQ(KVIB)                                                   ANH12270
      IF(DABS(FRQK).LE.WLIMIT) GO TO 109                                ANH12280
      IF(NDEG(KVIB).NE.1) GO TO 109                                     ANH12290
      FRQK2=FRQK*FRQK                                                   ANH12300
      WPPP=FRQI+FRQJ+FRQK                                               ANH12310
      WPPM=FRQI+FRQJ-FRQK                                               ANH12320
      WPMP=FRQI-FRQJ+FRQK                                               ANH12330
      WPMM=FRQI-FRQJ-FRQK                                               ANH12340
      DIJK=WPPP*WPPM*WPMP*WPMM                                          ANH12350
      DNUM2=F3Q(IVIB,IVIB,KVIB)*F3Q(KVIB,JVIB,JVIB)                     ANH12360
      DNUM3=(FRQK2-FRQI2-FRQJ2)*FRQK                                    ANH12370
      FAC=(F3Q(IVA,JVA,KVIB)+F3Q(IVA,JVB,KVIB)                          ANH12380
     *    +F3Q(IVB,JVA,KVIB)+F3Q(IVB,JVB,KVIB))*0.5D+00                 ANH12390
      FACT=FAC*FAC                                                      ANH12400
      VALU2=VALU2+DNUM2/FRQK                                            ANH12410
      IF(DABS(WPPM).LE.FLIM2) GO TO 226                                 ANH12420
      IF(DABS(WPMP).LE.FLIM2) GO TO 227                                 ANH12430
      IF(DABS(WPMM).LE.FLIM2) GO TO 228                                 ANH12440
      VALU3=VALU3+(DNUM3*FACT)/DIJK                                     ANH12450
      GVAL1=GVAL1+(FRQI*FRQJ*FRQK*FACT)/DIJK                            ANH12460
      GO TO 109                                                         ANH12470
C/////////////////////////////////                                      ANH12480
C///FERMI RESONANCE CORRECTIONS///                                      ANH12490
C/////////////////////////////////                                      ANH12500
C   WI+WJ=WK                                                            ANH12510
  226 CONTINUE                                                          ANH12520
      DNUM3=ONE/WPPP+ONE/WPMP-ONE/WPMM                                  ANH12530
      IF(DABS(FAC).GT.PLIMIT) THEN                                      ANH12540
        IANH2=IANH2+1                                                   ANH12550
      END IF                                                            ANH12560
      VALU3=VALU3+(DNUM3*FACT)/4.0D+00                                  ANH12570
      GFACT=ONE/WPPP-ONE/WPMP+ONE/WPMM                                  ANH12580
      GVAL1=GVAL1+(GFACT*FACT)/8.0D+00                                  ANH12590
      GO TO 109                                                         ANH12600
C   WI+WK=WJ                                                            ANH12610
  227 CONTINUE                                                          ANH12620
      DNUM3=ONE/WPPP-ONE/WPPM-ONE/WPMM                                  ANH12630
      IF(DABS(FAC).GT.PLIMIT) THEN                                      ANH12640
        IANH2=IANH2+1                                                   ANH12650
      END IF                                                            ANH12660
      VALU3=VALU3+(DNUM3*FACT)/4.0D+00                                  ANH12670
      GFACT=ONE/WPPP-ONE/WPPM+ONE/WPMM                                  ANH12680
      GVAL1=GVAL1+(GFACT*FACT)/8.0D+00                                  ANH12690
      GO TO 109                                                         ANH12700
C   WJ+WK=WI                                                            ANH12710
  228 CONTINUE                                                          ANH12720
      DNUM3=ONE/WPPP-ONE/WPPM+ONE/WPMP                                  ANH12730
      IF(DABS(FAC).GT.PLIMIT) THEN                                      ANH12740
        IANH2=IANH2+1                                                   ANH12750
      END IF                                                            ANH12760
      VALU3=VALU3+(DNUM3*FACT)/4.0D+00                                  ANH12770
      GFACT=ONE/WPPP-ONE/WPPM-ONE/WPMP                                  ANH12780
      GVAL1=GVAL1+(GFACT*FACT)/8.0D+00                                  ANH12790
C////////////////////////////                                           ANH12800
C///END OF FERMI RESONANCE///                                           ANH12810
C////////////////////////////                                           ANH12820
  109 CONTINUE                                                          ANH12830
      XRS=(VALU1-VALU2-VALU3)/4.0D+00                                   ANH12840
      GRS=GVAL1*0.5D+00                                                 ANH12850
C                                                                       ANH12860
  230 CONTINUE                                                          ANH12870
      XIJ(IVIB,JVIB)=XRS                                                ANH12880
      GIJ(IVIB,JVIB)=GRS                                                ANH12890
      IF(IPRNT.LE.2) GO TO 120                                          ANH12900
      WRITE(6,3) IVIB,JVIB,FRQI,FRQJ,VALU1,VALU2,VALU3,VALU4,XRS,GRS    ANH12910
  120 CONTINUE                                                          ANH12920
C                                                                       ANH12930
      IF(IANH1.NE.0) THEN                                               ANH12940
        WRITE(6,8)                                                      ANH12950
        WRITE(6,10) FLIM1                                               ANH12960
      END IF                                                            ANH12970
      IF(IANH2.NE.0) THEN                                               ANH12980
        WRITE(6,9)                                                      ANH12990
        WRITE(6,10) FLIM2                                               ANH13000
      END IF                                                            ANH13010
C                                                                       ANH13020
      IF(IPRNT.GT.3) GO TO 231                                          ANH13030
      WRITE(6,4)                                                        ANH13040
      DO 122 II=1,NVIB                                                  ANH13050
      IVIB=NFRQ(II)                                                     ANH13060
      FRQI=FREQ(IVIB)                                                   ANH13070
C*    DO 122 JJ=1,II                                                    ANH13080
      DO 122 JJ=1,NVIB                                                  ANH13090
      JVIB=NFRQ(JJ)                                                     ANH13100
      FRQJ=FREQ(JVIB)                                                   ANH13110
      XRS=XIJ(IVIB,JVIB)                                                ANH13120
      GRS=GIJ(IVIB,JVIB)                                                ANH13130
      WRITE(6,5) II,JJ,FRQI,FRQJ,XRS,GRS                                ANH13140
  122 CONTINUE                                                          ANH13150
C                                                                       ANH13160
C:::::::::::::::::::::::::::::::::::::::                                ANH13170
C:::CALCULATE FUNDAMENTAL FREQUENCIES:::                                ANH13180
C:::::::::::::::::::::::::::::::::::::::                                ANH13190
  231 CONTINUE                                                          ANH13200
      WRITE(6,6)                                                        ANH13210
      DO 125 II=1,NVIB                                                  ANH13220
      IVIB=NFRQ(II)                                                     ANH13230
      DEGI=NDEG(IVIB)                                                   ANH13240
      VALU1=FREQ(IVIB)                                                  ANH13250
      VALU2=XIJ(IVIB,IVIB)*(ONE+DEGI)                                   ANH13260
      VALU3=A00                                                         ANH13270
      VALU4=GIJ(IVIB,IVIB)                                              ANH13280
      DO 124 JJ=1,NVIB                                                  ANH13290
      JVIB=NFRQ(JJ)                                                     ANH13300
      IF(IVIB.EQ.JVIB) GO TO 124                                        ANH13310
      VALU3=VALU3+XIJ(IVIB,JVIB)                                        ANH13320
  124 CONTINUE                                                          ANH13330
      DELTA=VALU2+VALU3*0.5D+00+VALU4                                   ANH13340
      VALUT=VALU1+DELTA                                                 ANH13350
      WRITE(6,7) II,VALU1,DELTA,VALUT                                   ANH13360
  125 CONTINUE                                                          ANH13370
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
      SUBROUTINE LINZPE(ZETA,F3Q,F4Q) 
C   THE ANHARMONIC ZPVE FOR A LINEAR MOLECULE
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
    2 FORMAT(//,2X,' FOR A LINEAR MOLECULE')           
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
      DO 107 IABC=1,2
      II=IAXIS(IABC)
      BE=ROTCM(IABC)
      TERM=A00
      DO 106 IVIB=1,NVIB
      FRQI=FREQ(IVIB)
      DO 105 JVIB=1,NVIB
      FRQJ=FREQ(JVIB)
      ZET=ZETA(IVIB,JVIB,II)
      ZETSQ=ZET*ZET
      TERM=TERM+ZETSQ*((FRQI/FRQJ)-ONE)
  105 CONTINUE
  106 CONTINUE
      VALU5=VALU5+TERM*BE
  107 CONTINUE
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
      SUBROUTINE LING00(ZETA,F3Q,F4Q)
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
    5 FORMAT(//,2X,' FOR A LINEAR MOLECULE')
    6 FORMAT(   2X,' VALU1      = ',F20.10,' CM-1',/
     1          2X,' VALU2      = ',F20.10,' CM-1',/
     2          2X,' VALU3      = ',F20.10,' CM-1',/
     3          2X,' VALU4      = ',F20.10,' CM-1',/
     4          2X,' VALU5      = ',F20.10,' CM-1',/ 
     5          2X,' VALU6      = ',F20.10,' CM-1',/ 
     6          2X,' VALU7      = ',F20.10,' CM-1',/ 
     7          2X,' VALU8      = ',F20.10,' CM-1',/ 
     8          2X,' VALU9      = ',F20.10,' CM-1',/)
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
C  Valu6 is always zero
C  (included for consistency with other loops)
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
      VALU1 =+VALU1/(64.0D+00)
      VALU2 =-VALU2/(36.0D+00)
      VALU3 =+(VALU3*3.0D+00)/(64.0D+00)
      VALU4 =-VALU4/(4.0D+00)
      VALU5 =+VALU5/(96.0D+00)
      VALU6 =-VALU6/(36.0D+00)
      VALU7 =+VALU7/(32.0D+00)
      VALU8 =-VALU8/(4.0D+00)
      BE=ROTCM(1)
      VALU9=-DFLOAT(NATOM-2)*BE
      GCUBIC= VALU2+VALU3+VALU4+VALU6+VALU7+VALU8
      GQUART= VALU1+VALU5
      GKINET= VALU9
      G00   = GCUBIC+GQUART+GKINET
      ANHG00=ANHTOT+G00
      WRITE(6,5)
      IF(IPRNT.GT.2)
     *WRITE(6,6) VALU1,VALU2,VALU3,VALU4,
     *           VALU5,VALU6,VALU7,VALU8,
     *           VALU9
      WRITE(6,7) GCUBIC,GQUART,GKINET,G00
      WRITE(6,8) HZPVE,AZPVE,ANHTOT,ANHG00
C
      RETURN
      END
