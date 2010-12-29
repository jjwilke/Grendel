c================================================================
      PROGRAM ANHARM                                   
c================================================================
c   This is a vibrational anharmonicity program
c   based on the second-order vibrational perturbation theory
c================================================================
c   Last updated on October 19, 2007 by Dr. Yukio Yamaguchi   ***
c   Modification for maxcor = 100,000 -> 1000,000             ***
c================================================================
c   Last updated on January 10, 2006 by Dr. Yukio Yamaguchi   ***
c   Modification for a stand-alone version                    ***
c================================================================
c   Last updated on April 20, 2004 by Dr. Yukio Yamaguchi     ***
c   Modification for anharmonic ZPVE and G_0 values           ***
c================================================================
c   Last updated on September 10, 1995 by Dr. Yukio Yamaguchi ***
c   Modification for Unix System                              ***
c================================================================
C***LAST UPDATED ON JANUARY 19, 1989 BY DR. YUKIO YAMAGUCHI   ***   
c================================================================
      IMPLICIT REAL*8 (A-H,O-Z)                                 
      parameter (maxcor=1000000)
cOLD     parameter (maxcor=100000)
      parameter (llength=80)
      dimension cc(maxcor),ia(maxcor*2)
      DIMENSION COORD(3,50)   
      LOGICAL ISOMOL                                                    ANH00090
      CHARACTER*1 RAXIS                                                 ANH00100
      CHARACTER*25 RTYPE                                                ANH00110
      DIMENSION ATM(36)                                                 ANH00120
      COMMON/VIB101/NATOM,N3N,NATRI,ILIN,NVIB                           ANH00130
      COMMON/VIB102/ITHREE,IFOUR,N3TOT,N4TOT                            ANH00140
      COMMON/VIB103/PARA,WAVE,CONST,CYCL,CONV                           ANH00150
      COMMON/VIB104/AIXX,AIYY,AIZZ,SUMW                                 ANH00160
      COMMON/VIB105/FACT3,FACT4                                         ANH00170
      COMMON/VIB106/ROTAA(3),ROTGC(3),ROTCM(3),ROTMH(3)                 ANH00180
      COMMON/VIB107/ITOP,IVSGL,IVDBL,IVTRP                              ANH00190
      COMMON/VIB108/IAXIS(3),NDEG(150),NDAB(150,5),IMAG(150)            ANH00200
      COMMON/VIB109/CLIMIT,FLIM1,FLIM2                                  ANH00210
      COMMON/VIB110/ISIGMA                                              ANH00220
      COMMON/VIB111/ANHTOT,HZPVE,AZPVE
      COMMON/VIB112/LINEWD
      COMMON/VIB201/CHG(50),X(50),Y(50),Z(50),W(50)                     ANH00230
      COMMON/VIB202/R(1275)                                             ANH00240
      COMMON/VIB203/IOFF(150),IPRNT                                     ANH00250
      COMMON/VIB204/SQM(150),ROOT(150),FREQ(150)                        ANH00260
      COMMON/VIB205/IFREQ,NFRQ(150)                                     ANH00270
      COMMON/VIB206/RAXIS(3),RTYPE(6)                                   ANH00280
      COMMON/VIB207/TABCD(3,3,3,3),PABC(3,3,3)                          ANH00290
      COMMON/VIB208/TAABB(3,3),TABAB(3,3)                               ANH00300
      COMMON/VIB209/SALP(3)                                             ANH00310
      DATA ATM / 1.007825D+00 ,  4.00260D+00 ,  7.01600D+00 ,
     *           9.01218D+00  , 11.00931D+00 , 12.00000D+00 ,           NOR00620
     *           14.00307D+00 , 15.99491D+00 , 18.99840D+00 ,           NOR00630
     *           19.99244D+00 , 22.98977D+00 , 23.98504D+00 ,           NOR00640
     *           26.98154D+00 , 27.97693D+00 , 30.97376D+00 ,           NOR00650
     *           31.97207D+00 , 34.96885D+00 , 39.96238D+00 ,           NOR00660
     *           38.96371D+00 , 39.96259D+00 , 44.95591D+00 ,           NOR00670
     *           47.94795D+00 , 50.94396D+00 , 51.94051D+00 ,           NOR00680
     *           54.93805D+00 , 55.93494D+00 , 58.93320D+00 ,           NOR00690
     *           57.93535D+00 , 62.92960D+00 , 63.92915D+00 ,           NOR00700
     *           68.92558D+00 , 73.92118D+00 , 74.92160D+00 ,           NOR00710
     *           79.91652D+00 , 78.91834D+00 , 83.91151D+00   /         NOR00720
      DATA PH,BK,AVN / 6.626176D+00 , 1.380662D+00 , 6.022045D+00 /     ANH00390
      DATA A0 / 0.529177249D+00 /
      DATA HE,EV / 4.359813653D+00 , 1.6021892D+00 /
      DATA CL,TEMP / 2.99792458D+00 , 300.0D+00 /                       ANH00410
      DATA PI,ALN / 3.1415926536D+00 , 2.3025850930D+00 /               ANH00420
    1 FORMAT(//,2X,' -----------------------------------------'/
     *          2X,'   THE VIBRATIONAL ANHARMONICITY PROGRAM  '/
     *          2X,'   Center for Computational Chemistry     '/
     *          2X,'          University of Georgia           '/
     *          2X,'         Athens, Georgia, 30602           '/
     *          2X,'                 U.S.A.                   '/
     *          2x,'        (Version of 10/19/07)             '/ 
     *          2x,'        (     15:20 PM      )             '/ 
     *          2X,' -----------------------------------------'/)
    2 FORMAT(10I5)                                                      ANH00440
    3 FORMAT(F20.10)                                                    ANH00450
    4 FORMAT(//,2X,' PARAMETERS FOR THE CALCULATION'//                  ANH00460
     * 2X,' ISOTOP = ',I5/                                              ANH00470
     * 2X,' IGEOMT = ',I5/                                              ANH00480
     * 2X,' ITHREE = ',I5/                                              ANH00490
     * 2X,' IFOUR  = ',I5/                                              ANH00500
     * 2X,' IFREQ  = ',I5/                                              ANH00510
     * 2X,' ICORIO = ',I5/                                              ANH00520
     * 2X,' IFERM1 = ',I5/                                              ANH00530
     * 2X,' IFERM2 = ',I5/                                              ANH00540
     * 2X,' ISIGMA = ',I5/                                              ANH00550
     * 2X,' IPRNT  = ',I5/                                              ANH00560
     * 2X,' NATOM  = ',I5/                                              ANH00570
     * 2X,' N3N    = ',I5/                                              ANH00580
     * 2X,' NATRI  = ',I5/                                              ANH00590
     * 2X,' N3TOT  = ',I5/                                              ANH00600
     * 2X,' N4TOT  = ',I5/)                                             ANH00610
    5 FORMAT(//,2X,' SCF ENERGY = ',F20.10/)                            ANH00620
    6 FORMAT(//,2X,' MOLECULAR GEOMETRY (IN BOHR)'/                     ANH00630
     1 5X,4H NO.,11X,2H X,18X,2H Y,18X,2H Z,18X,2H W/)                  ANH00640
    7 FORMAT(4F20.10)                                                   ANH00650
    8 FORMAT(2X,I7,5X,4F20.10)                                          ANH00660
    9 FORMAT(//,2X,' MOLECULAR GEOMETRY (IN A)'/                        ANH00670
     1 5X,4H NO.,11X,2H X,18X,2H Y,18X,2H Z,18X,2H W/)                  ANH00680
   10 FORMAT(F20.10)                                                    ANH00690
   11 FORMAT(//,2X,' ICMAX = ',I10,5X,' MAXCOR = ',I10/)                ANH00700
   12 FORMAT(//,2X,' -----------------------------------------'/
     *          2X,'   THE VIBRATIONAL ANHARMONICITY PROGRAM  '/
     *          2X,'   Center for Computational Chemistry     '/
     *          2X,'          University of Georgia           '/
     *          2X,' -----------------------------------------'/)
   13 FORMAT(/,2X,' -----------------------------------------'/
     *         2X,'      End of the anharm program           '/
     *         2X,' -----------------------------------------'/)
C                     
      input=5
      itape3=3
      itape6=6
      itap15=15
      itap20=20
      itap24=24
      linewd=llength
c
      call ffile(input,'input',0)
      call ffile(itape3,'tape3',0)
      call ffile(itape6,'output',0)
      call ffile(itap15,'file15',0)
      call ffile(itap20,'file20',0)
      call ffile(itap24,'file24',0)
C                                                                       ANH00750
C   PHYSICAL CONSTANTS                                                  ANH00760
      PARA=1.0D+00/AVN                                                  ANH00770
      WAVE=1.0D+04*DSQRT(AVN)/(2.0D+00*PI*CL)                           ANH00780
      CONST=1.0D+02*(PH*AVN)/(8.0D+00*PI*PI*CL)                         ANH00790
      CYCL=1.0D+06*(PH*AVN)/(8.0D+00*PI*PI)                             ANH00800
      CONV=1.0D-01*(PH*CL)/(2.0D+00*BK*TEMP)                            ANH00810
C   ALAM IS IN AMU-1.A-2                                                ANH00820
      ALAM=4.0D-02*(PI*PI*CL)/(PH*AVN)                                  ANH00830
C   SQLAM IS IN AMU-1/2.A-1                                             ANH00840
      SQLAM=DSQRT(ALAM)                                                 ANH00850
C   FACT3 IS IN CM-1                                                    ANH00860
      FACT3=1.0D+06/(SQLAM*SQLAM*SQLAM*PH*CL)                           ANH00870
C   FACT4 IS IN CM-1                                                    ANH00880
      FACT4=1.0D+06/(ALAM*ALAM*PH*CL)                                   ANH00890
C                                                                       ANH00900
      CALL LOCATE(INPUT,'# ANHARM #',IERR)                              ANH00990
C                                                                       ANH01000
      WRITE(6,1)                                                        ANH01010
      write(6,*) '  maxcor = ',maxcor
      write(6,*) '  what in the name of God is going on '
C                                                                       ANH01020
C***INPUT PARAMETERS FOR THE CALCULATION***                             ANH01030
C   ISOTOP IS A PARAMETER FOR ISOTOPIC MOLECULES                        ANH01040
C     ISOTOP =0 AND 1 :USE REGULAR ATOMIC MASS                          ANH01050
C     ISOTOP =2 OR MORE :NUMBER OF ISOTOPIC MOLECULES                   ANH01060
C   IGEOMT IS A PARAMETER FOR INPUT DATA                                ANH01070
C     IGEOMT =0 :READ IN GEOMETRY FROM TAPE 5                           ANH01080
C     IGEOMT =1 :READ IN GEOMETRY FROM TAPE 5                           ANH01080
C     IGEOMT =5 :READ IN GEOMETRY FROM TAPE 5                           ANH01090
C   ITHREE IS A PARAMETER FOR THIRD DERIVATIVES                         ANH01110
C     ITHREE =0 :DO NOT READ IN THIRD DERIVATIVES                       ANH01120
C     ITHREE =1 OR 20:READ IN THIRD DERIVATIVES FROM FILE20             ANH01130
C   IFOUR IS A PARAMETER FOR FOURTH DERIVATIVES                         ANH01150
C     IFOUR  =0 :DO NOT READ IN FOURTH DERIVATIVES                      ANH01160
C     IFOUR  =1 OR 24:READ IN FOURTH DERIVATIVES FROM FILE24            ANH01170
C   IFREQ IS A PARAMETER FOR ORDERING OF NORMAL MODES
C     IFREQ  =0 :CALCULATE 3N-6(OR 3N-5) SETS OF ANHARMONIC CONSTANTS   ANH01190
C     IFREQ  =N :REORDER IFREQ SETS OF ANHARMONIC CONSTANTS             ANH01200
C   ICORIO IS A PARAMETER FOR CORIOLIS RESONANCE
C     ICORIO =0 :USE A DEFAULT VALU FOR CORIOLIS RESONANCE              ANH01210
C     ICORIO =1 :READ IN A THRESHOLD VALU FOR CORIOLIS RESONANCE        ANH01220
C   IFERM1 IS A PARAMETER FOR TYPE 1 FERMI RESONANCE
C     IFERM1 =0 :USE A DEFAULT VALU FOR TYPE 1 FERMI RESONANCE          ANH01230
C                TYPE1----W(I)=W(J)                                     ANH01240
C     IFERM1 =1 :READ IN A THRESHOLD VALU FOR TYPE 1 FERMI RESONANCE    ANH01250
C   IFERM2 IS A PARAMETER FOR TYPE 1 FERMI RESONANCE
C     IFERM2 =0 :USE A DEFAULT VALU FOR TYPE 2 FERMI RESONANCE          ANH01260
C                TYPE2---W(I)+W(J)=W(K)                                 ANH01270
C     IFERM2 =1 :READ IN A THRESHOLD VALU FOR TYPE 2 FERMI RESONANCE    ANH01280
C   ISIGMA IS A PARAMETER FOR SIGMA VALUE
C     ISIGMA =0 :USE A(0') VALUE TO EVALUATE SIGMA VALUE THAT IS NEEDED ANH01290
C                FOR CENTRIFUGAL DISTORTION CONSTANTS                   ANH01300
C     ISIGMA =1 :USE A(E)  VALUE TO EVALUATE SIGMA VALUE THAT IS NEEDED ANH01310
C                FOR CENTRIFUGAL DISTORTION CONSTANTS                   ANH01320
C   IPRNT IS A PARAMETER FOR OUTPUT
C     IPRNT  =0 :MINIMUM OUTPUT
C     IPRNT  =N :MORE INTERMEDIATE INFORMATION
C***********************************************************************
      READ(5,2) ISOTOP,IGEOMT,ITHREE,IFOUR,IFREQ,ICORIO,IFERM1,IFERM2,  ANH01330
     1          ISIGMA,IPRNT                                            ANH01340

c     write(6,*)  ' After read(5,2)'
c     write(6,2) ISOTOP,IGEOMT,ITHREE,IFOUR,IFREQ,ICORIO,IFERM1,IFERM2,
c    1           ISIGMA,IPRNT                                          

      IF(IGEOMT.EQ.0) IGEOMT=5
      IF(IGEOMT.EQ.1) IGEOMT=5
      IF(ITHREE.EQ.1) ITHREE=20                                         ANH01360
      IF(IFOUR.EQ.1) IFOUR=24                                           ANH01370
C                                                                       ANH01380
      REWIND ITAP15                                                     ANH01390
      READ(ITAP15,2) NATOM,NSTORE                                       ANH01400

c     write(6,*) ' After read itap15 '
c     write(6,2) NATOM,NSTORE      
C                                                                       ANH01410
      IOFF(1)=0                                                         ANH01420
      DO 101 I=1,149                                                    ANH01430
      IOFF(I+1)=IOFF(I)+I                                               ANH01440
  101 CONTINUE                                                          ANH01450
C                                                                       ANH01460
      N3N=NATOM*3                                                       ANH01800
      NATRI=IOFF(N3N+1)                                                 ANH01810
      N3TOT=N3N*(N3N+1)*(N3N+2)/6                                       ANH01820
      N4TOT=N3N*(N3N+1)*(N3N+2)*(N3N+3)/24                              ANH01830
C                                                                       ANH01840
C   SET AN ARRAY FOR ANHARMONICITY CALCULATION                          ANH01850
      DO 102 I=1,N3N                                                    ANH01860
      NFRQ(I)=I                                                         ANH01870
  102 CONTINUE
      CALL LOCATE(INPUT,'# ANHARM #',IERR)                              ANH01880
      READ(5,2) ISOTOP,IGEOMT,ITHXXX,IFXXX,IFREQ,ICORIO,IFERM1,IFERM2,  ANH01890
     1          ISIGMA,IPRNT                                            ANH01900

c     write(6,*) ' After read(5,2)'
c     write(6,2) ISOTOP,IGEOMT,ITHXXX,IFXXX,IFREQ,ICORIO,IFERM1,IFERM2, 
c    1           ISIGMA,IPRNT                                          

      IF(ISOTOP.EQ.0) ISOTOP=1                                          ANH01910
      IF(IFREQ.NE.0) THEN                                               ANH01920
        NFREQ=IABS(IFREQ)                                               ANH01930
        READ(5,2) (NFRQ(I),I=1,NFREQ)                                   ANH01940
      END IF                                                            ANH01950
      IF(ICORIO.NE.0) THEN                                              ANH01960
        READ(5,3) CLIMIT                                                ANH01970
      ELSE                                                              ANH01980
        CLIMIT=100.0D+00                                                ANH01990
      END IF                                                            ANH02000
      IF(IFERM1.NE.0) THEN                                              ANH02010
        READ(5,3) FLIM1                                                 ANH02020
      ELSE                                                              ANH02030
        FLIM1=100.0D+00                                                 ANH02040
      END IF                                                            ANH02050
      IF(IFERM2.NE.0) THEN                                              ANH02060
        READ(5,3) FLIM2                                                 ANH02070
      ELSE                                                              ANH02080
        FLIM2=100.0D+00                                                 ANH02090
      END IF                                                            ANH02100
C                                                                       ANH02110
      WRITE(6,4) ISOTOP,IGEOMT,ITHREE,IFOUR,IFREQ,ICORIO,IFERM1,IFERM2, ANH02120
     1           ISIGMA,IPRNT,NATOM,N3N,NATRI,N3TOT,N4TOT               ANH02130
C                                                                       ANH02140
C  Commented out for a stand-alone version
c 12-19-05 yy  WRITE(6,5) ESCF   
      WRITE(6,6)                                                        ANH02160
      DO 103 I=1,NATOM                                                  ANH02170
      READ(5,7) COORD(1,I),COORD(2,I),COORD(3,I),W(I)                   ANH02270
      X(I)=COORD(1,I)                                                   ANH02280
      Y(I)=COORD(2,I)                                                   ANH02290
      Z(I)=COORD(3,I)                                                   ANH02300
      WRITE(6,8) I,X(I),Y(I),Z(I),W(I)                                  ANH02320
  103 CONTINUE                                                          ANH02330
      CALL DIST(1)                                                      ANH02340
      WRITE(6,9)                                                        ANH02350
      DO 104 I=1,NATOM                                                  ANH02360
      X(I)=X(I)*A0                                                      ANH02370
      Y(I)=Y(I)*A0                                                      ANH02380
      Z(I)=Z(I)*A0                                                      ANH02390
      WRITE(6,8) I,X(I),Y(I),Z(I),W(I)                                  ANH02400
  104 CONTINUE                                                          ANH02410
      CALL DIST(1)                                                      ANH02420
C                                                                       ANH02430
C   DYNAMIC ALLOCATION                                                  ANH02440
C   IC1  : FX   (N3N*N3N)                                               ANH02450
C   IC2  : ELX  (N3N*N3N)                                               ANH02460
C   IC3  : ELXM (N3N*N3N)                                               ANH02470
C   IC4  : UX   (N3N*N3N)                                               ANH02480
C   IC5  : EPX  (N3N*N3N)                                               ANH02490
C   IC6  : EPXM (N3N*N3N)                                               ANH02500
C   IC7  : ZETA (N3N*N3N*3)                                             ANH02510
C   IC8  : AXY  (3*3*N3N)                                               ANH02520
C   IC9  : BXY  (3*3*N3N)                                               ANH02530
C   IC10 : CXY  (3*3*N3N)                                               ANH02540
C   IC11 : BXYZ (3*3*3*N3N)                                             ANH02550
C   IC12 : CXYZ (3*3*3*N3N)                                             ANH02560
C   IC13 : F33  (N3TOT)                                                 ANH02570
C   IC14 : F3X  (N3N*N3N*N3N)                                           ANH02580
C   IC15 : F3Q  (N3N*N3N*N3N)                                           ANH02590
C   IC16 : F44  (N4TOT)                                                 ANH02600
C   IC17 : F4X  (N3N*N3N*N3N*N3N)                                       ANH02610
C   IC18 : F4Q  (N3N*N3N*N3N*N3N)                                       ANH02620
C                                                                       ANH02630
C   DYNAMIC ALLOCATION OF CORE MEMORY AND CLEARANCE OF IT               ANH02640
      IC1=1                                                             ANH02650
      IC2=IC1+N3N*N3N                                                   ANH02660
      IC3=IC2+N3N*N3N                                                   ANH02670
      IC4=IC3+N3N*N3N                                                   ANH02680
      IC5=IC4+N3N*N3N                                                   ANH02690
      IC6=IC5+N3N*N3N                                                   ANH02700
      IC7=IC6+N3N*N3N                                                   ANH02710
      IC8=IC7+N3N*N3N*3                                                 ANH02720
      IC9=IC8+3*3*N3N                                                   ANH02730
      IC10=IC9+3*3*N3N                                                  ANH02740
      IC11=IC10+3*3*N3N                                                 ANH02750
      IC12=IC11+3*3*3*N3N                                               ANH02760
      IC13=IC12+3*3*3*N3N                                               ANH02770
      IC14=IC13+N3TOT                                                   ANH02780
      IC15=IC14+N3N*N3N*N3N                                             ANH02790
      IC16=IC15+N3N*N3N*N3N                                             ANH02800
      IC17=IC16+N4TOT                                                   ANH02810
      IC18=IC17+N3N*N3N*N3N*N3N                                         ANH02820
      ICKEEP=IC18+N3N*N3N*N3N*N3N                                       ANH02830
      if(ickeep.gt.maxcor) then
        write(6,*) '  '
        write(6,*) ' ickeep is larger than maxcor ---> STOP!'
        write(6,*) ' increase maxcor'
        write(6,*) ' ickeep  = ',ickeep
        write(6,*) ' maxcor  = ',maxcor
        write(6,*) '  '
        stop
      end if
C
      if(iprnt.gt.5) then
        write(6,*) '  ic1    = ',ic1
        write(6,*) '  ic2    = ',ic2
        write(6,*) '  ic3    = ',ic3
        write(6,*) '  ic4    = ',ic4
        write(6,*) '  ic5    = ',ic5
        write(6,*) '  ic6    = ',ic6
        write(6,*) '  ic7    = ',ic7
        write(6,*) '  ic8    = ',ic8
        write(6,*) '  ic9    = ',ic9
        write(6,*) '  ic10   = ',ic10
        write(6,*) '  ic11   = ',ic11
        write(6,*) '  ic12   = ',ic12
        write(6,*) '  ic13   = ',ic13
        write(6,*) '  ic14   = ',ic14
        write(6,*) '  ic15   = ',ic15
        write(6,*) '  ic16   = ',ic16
        write(6,*) '  ic17   = ',ic17
        write(6,*) '  ic18   = ',ic18
        write(6,*) '  ickeep = ',ickeep
      end if
C
      CALL ZERO(CC(1),ICKEEP)                                           ANH02840
C                                                                       ANH02850
C   FORM THE FORCE CONSTANT MATRIX                                      ANH02860
      WRITE(3,21)                                                       ANH02870
   21 FORMAT(//,2X,' NOW YOU ARE IN FCONST'/)                           ANH02880
      ICMAX=ICKEEP                                                      ANH02890
      WRITE(3,11) ICMAX,MAXCOR                                          ANH02900
C.................FX......                                              ANH02910
      CALL FCONST(CC(IC1))                                              ANH02920
C                                                                       ANH02930
C    CALCULATE THE CENTER OF MASS AND MOMENTS OF INERTIA                ANH02940
C    FOR THE PARENT MOLECULE                                            ANH02950
      WRITE(3,22)                                                       ANH02960
   22 FORMAT(//,2X,' NOW YOU ARE IN MOMENT'/)                           ANH02970
      IC19=ICKEEP                                                       ANH02980
      ICMAX=IC19+N3N*N3N                                                ANH02990
      WRITE(3,11) ICMAX,MAXCOR                                          ANH03000
C.................UX      EE.......                                     ANH03010
      CALL MOMENT(CC(IC4),CC(IC19))                                     ANH03020
C                                                                       ANH03030
      ISOMOL=.FALSE.                                                    ANH03040
      NISO=0                                                            ANH03050
  400 CONTINUE
      NISO=NISO+1                                                       ANH03060
      IF(NISO.EQ.1) GO TO 401                                           ANH03070
      IF(NISO.GT.ISOTOP) GO TO 220                                      ANH03080
      WRITE(6,12)                                                       ANH03090
      ISOMOL=.TRUE.                                                     ANH03100
C                                                                       ANH03110
C   READ IN ISOTOPIC ATOMIC MASS                                        ANH03120
      WRITE(6,6)                                                        ANH03130
      DO 105 I=1,NATOM                                                  ANH03140
      READ(5,10) W(I)                                                   ANH03150
      X(I)=COORD(1,I)                                                   ANH03160
      Y(I)=COORD(2,I)                                                   ANH03170
      Z(I)=COORD(3,I)                                                   ANH03180
      WRITE(6,8) I,X(I),Y(I),Z(I),W(I)                                  ANH03190
  105 CONTINUE                                                          ANH03200
      CALL DIST(1)                                                      ANH03210
      WRITE(6,9)                                                        ANH03220
      DO 106 I=1,NATOM                                                  ANH03230
      X(I)=X(I)*A0                                                      ANH03240
      Y(I)=Y(I)*A0                                                      ANH03250
      Z(I)=Z(I)*A0                                                      ANH03260
      WRITE(6,8) I,X(I),Y(I),Z(I),W(I)                                  ANH03270
  106 CONTINUE                                                          ANH03280
      CALL DIST(1)                                                      ANH03290
C                                                                       ANH03300
C    CALCULATE THE CENTER OF MASS AND MOMENTS OF INERTIA                ANH03310
C    FOR AN ISOTOPIC MOLECULE                                           ANH03320
      CALL ZERO(CC(IC2),ICKEEP)                                         ANH03330
      WRITE(3,22)                                                       ANH03340
      IC19=ICKEEP                                                       ANH03350
      ICMAX=IC19+N3N*N3N                                                ANH03360
      WRITE(3,11) ICMAX,MAXCOR                                          ANH03370
C.................UX      EE.......                                     ANH03380
      CALL MOMENT(CC(IC4),CC(IC19))                                     ANH03390
C                                                                       ANH03400
C   THE NORMAL COORDINATE ANALYSIS                                      ANH03410
C                                                                       ANH03420
C   THE FX MATRIX METHOD                                                ANH03430
  401 CONTINUE                                                          ANH03440
      WRITE(3,23)                                                       ANH03450
   23 FORMAT(//,2X,' NOW YOU ARE IN VIBFX'/)                            ANH03460
      IC19=ICKEEP                                                       ANH03470
      IC20=IC19+N3N*N3N                                                 ANH03480
      IC21=IC20+N3N                                                     ANH03490
      IC22=IC21+N3N                                                     ANH03500
      IC23=IC22+N3N                                                     ANH03510
      IC24=IC23+NATRI                                                   ANH03520
      ICMAX=IC24+N3N*N3N                                                ANH03530
      WRITE(3,11) ICMAX,MAXCOR                                          ANH03540
C................FX      ELX     ELXM    UX      EPX     EPXM....       ANH03550
      CALL VIBFX(CC(IC1),CC(IC2),CC(IC3),CC(IC4),CC(IC5),CC(IC6),       ANH03560
C................EE       VALU     FV1      FV2      FXM      DD....... ANH03570
     1           CC(IC19),CC(IC20),CC(IC21),CC(IC22),CC(IC23),CC(IC24), ANH03580
C................CC.......                                              ANH03590
     2           CC(IC24))                                              ANH03600
C                                                                       ANH03610
C   CORIOLIS COUPLING CONSTANTS                                         ANH03620
      WRITE(3,24)                                                       ANH03630
   24 FORMAT(//,2X,' NOW YOU ARE IN CORIOL'/)                           ANH03640
      IC19=ICKEEP                                                       ANH03650
      IC20=IC19+N3N*N3N                                                 ANH03660
      IC21=IC20+N3N*N3N                                                 ANH03670
      ICMAX=IC21+N3N*N3N*3                                              ANH03680
      WRITE(3,11) ICMAX,MAXCOR                                          ANH03690
C.................EPXM    ZETA    EE       FF       AM.......           ANH03700
      CALL CORIOL(CC(IC6),CC(IC7),CC(IC19),CC(IC20),CC(IC21))           ANH03710
C                                                                       ANH03720
C   CALCULATE AXY MATRICES                                              ANH03730
      WRITE(3,25)                                                       ANH03740
   25 FORMAT(//,2X,' NOW YOU ARE IN ABCMAT'/)                           ANH03750
      IC19=ICKEEP                                                       ANH03760
      ICMAX=IC19+N3N*N3N                                                ANH03770
      WRITE(3,11) ICMAX,MAXCOR                                          ANH03780
C.................EPXM    AXY     BXY     CXY      EE.......            ANH03790
      CALL ABCMAT(CC(IC6),CC(IC8),CC(IC9),CC(IC10),CC(IC19))            ANH03800
C                                                                       ANH03810
C   CALCULATE TAU MATRICES                                              ANH03820
      WRITE(3,27)                                                       ANH03830
   27 FORMAT(//,2X,' NOW YOU ARE IN TAUMAT'/)                           ANH03840
      IC19=ICKEEP                                                       ANH03850
      ICMAX=IC19+N3N*N3N                                                ANH03860
      WRITE(3,11) ICMAX,MAXCOR                                          ANH03870
C.................AXY     BXY.....                                      ANH03880
      CALL TAUMAT(CC(IC8),CC(IC9))                                      ANH03890
C                                                                       ANH03900
C   CALCULATE BXYZ MATRICES                                             ANH03910
      WRITE(3,28)                                                       ANH03920
   28 FORMAT(//,2X,' NOW YOU ARE IN BCXYZ'/)                            ANH03930
      ICMAX=ICKEEP                                                      ANH03940
      WRITE(3,11) ICMAX,MAXCOR                                          ANH03950
C................ZETA    BXY     BXYZ     CXYZ.....                     ANH03960
      CALL BCXYZ(CC(IC7),CC(IC9),CC(IC11),CC(IC12))                     ANH03970
C                                                                       ANH03980
C   CALCULATE CUBIC FORCE CONSTANTS                                     ANH04090
      IF(ITHREE.EQ.0) GO TO 303                                         ANH04110
      WRITE(3,30)                                                       ANH04130
   30 FORMAT(//,2X,' NOW YOU ARE IN FORCE3'/)                           ANH04140
      IC19=ICKEEP                                                       ANH04150
      IC20=IC19+N3N*N3N                                                 ANH04160
      IC21=IC20+N3N*N3N                                                 ANH04170
      ICMAX=IC21+N3N*N3N                                                ANH04180
      WRITE(3,11) ICMAX,MAXCOR                                          ANH04190
C.................ELX     F33      F3X      F3Q      EE       DD....... ANH04200
      CALL FORCE3(CC(IC2),CC(IC13),CC(IC14),CC(IC15),CC(IC19),CC(IC20), ANH04210
C.................GG.......                                             ANH04220
     1            CC(IC21))                                             ANH04230
C                                                                       ANH04240
C   CALCULATE PHI MATRICES                                              ANH04250
  303 CONTINUE                                                          ANH04260
      WRITE(3,31)                                                       ANH04270
   31 FORMAT(//,2X,' NOW YOU ARE IN PHIMAT'/)                           ANH04280
      IC19=ICKEEP                                                       ANH04290
      ICMAX=IC19+N3N*N3N                                                ANH04300
      WRITE(3,11) ICMAX,MAXCOR                                          ANH04310
C.................BXY     CXY      BXYZ     CXYZ     F3Q      EE....... ANH04320
      CALL PHIMAT(CC(IC9),CC(IC10),CC(IC11),CC(IC12),CC(IC15),CC(IC19)) ANH04330
C                                                                       ANH04340
C   CALCULATE QUARTIC FORCE CONSTANTS                                   ANH04350
      IF(IFOUR.EQ.0) GO TO 304                                          ANH04360
      WRITE(3,32)                                                       ANH04380
   32 FORMAT(//,2X,' NOW YOU ARE IN FORCE4'/)                           ANH04390
      IC19=ICKEEP                                                       ANH04400
      IC20=IC19+N3N*N3N                                                 ANH04410
      IC21=IC20+N3N*N3N                                                 ANH04420
      ICMAX=IC21+N3N*N3N                                                ANH04430
      WRITE(3,11) ICMAX,MAXCOR                                          ANH04440
C.................ELX     F44      F4X      F4Q      EE       DD....... ANH04450
      CALL FORCE4(CC(IC2),CC(IC16),CC(IC17),CC(IC18),CC(IC19),CC(IC20), ANH04460
C.................GG.......                                             ANH04470
     1            CC(IC21))                                             ANH04480
C                                                                       ANH04490
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                          ANH04500
C%%%ANHARMONICITY ANALYSIS %%%                                          ANH04510
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                          ANH04520
C                                                                       ANH04530
C   SELECT THE ROTATIONAL TYPE OF THE MOLECULE                          ANH04540
  304 CONTINUE                                                          ANH04550
      GO TO (501,502,503,504,504,505),ITOP                              ANH04560
C                                                                       ANH04570
C:::::::::::::::::::::::                                                ANH04580
C:::DIATOMIC MOLECULE:::                                                ANH04590
C:::::::::::::::::::::::                                                ANH04600
  501 CONTINUE                                                          ANH04610
      WRITE(3,33)                                                       ANH04620
   33 FORMAT(//,2X,' NOW YOU ARE IN DIATOM'/)                           ANH04630
      IC19=ICKEEP                                                       ANH04640
      ICMAX=IC19+N3N*N3N                                                ANH04650
      WRITE(3,11) ICMAX,MAXCOR                                          ANH04660
C.................F3Q      F4Q      EE.......                           ANH04670
      CALL DIATOM(CC(IC15),CC(IC18),CC(IC19))                           ANH04680
      GO TO 215                                                         ANH04690
C                                                                       ANH04700
C::::::::::::::::::::::::::::::::                                       ANH04710
C:::LINEAR POLYATOMIC MOLECULE:::                                       ANH04720
C::::::::::::::::::::::::::::::::                                       ANH04730
  502 CONTINUE                                                          ANH04740
      WRITE(3,34)                                                       ANH04750
   34 FORMAT(//,2X,' NOW YOU ARE IN LINEAR'/)                           ANH04760
      IC19=ICKEEP                                                       ANH04770
      IC20=IC19+N3N*N3N                                                 ANH04780
      ICMAX=IC20+N3N*N3N                                                ANH04790
      WRITE(3,11) ICMAX,MAXCOR                                          ANH04800
C.................ZETA    CXY      F3Q      CLTYP....                   ANH04810
      CALL LTYPLI(CC(IC7),CC(IC10),CC(IC15),CC(IC19))                   ANH04820
C.................ZETA    AXY     F3Q      EE       RVCNST...           ANH04830
      CALL ROVLIR(CC(IC7),CC(IC8),CC(IC15),CC(IC19),CC(IC20))           ANH04840
C.................BXY     CXY      F3Q      EE.......                   ANH04850
      CALL CENTLI(CC(IC9),CC(IC10),CC(IC15),CC(IC19))                   ANH04860
      IF(IVDBL.LE.1) GO TO 205                                          ANH04870
C.................F3Q      F4Q      CLTYP....                           ANH04880
      CALL VLTYLI(CC(IC15),CC(IC18),CC(IC19))                           ANH04890
  205 CONTINUE                                                          ANH04900
C.................ZETA    F3Q      F4Q      XIJ      GIJ......          ANH04910
      CALL LINEAR(CC(IC7),CC(IC15),CC(IC18),CC(IC19),CC(IC20))          ANH04920
C.................ZETA    F3Q      F4Q......
      CALL LINZPE(CC(IC7),CC(IC15),CC(IC18))
C.................ZETA    F3Q      F4Q......
      CALL LING00(CC(IC7),CC(IC15),CC(IC18))
      GO TO 215                                                         ANH04930
C                                                                       ANH04940
C:::::::::::::::::::::::::::::                                          ANH04950
C:::ASYMMETRIC TOP MOLECULE:::                                          ANH04960
C:::::::::::::::::::::::::::::                                          ANH04970
  503 CONTINUE                                                          ANH04980
      WRITE(3,35)                                                       ANH04990
   35 FORMAT(//,2X,' NOW YOU ARE IN ASYTOP'/)                           ANH05000
      IC19=ICKEEP                                                       ANH05010
      IC20=IC19+N3N*N3N                                                 ANH05020
      ICMAX=IC20+N3N*N3N                                                ANH05030
      WRITE(3,11) ICMAX,MAXCOR                                          ANH05040
C.................ZETA    AXY     F3Q      EE       RVCNST...           ANH05050
      CALL ROVASY(CC(IC7),CC(IC8),CC(IC15),CC(IC19),CC(IC20))           ANH05060
C.................F3Q      EE.......                                    ANH05070
      CALL CENTAS(CC(IC15),CC(IC19))                                    ANH05080
C.................ZETA    F3Q      F4Q      XIJ......                   ANH05090
      CALL ASYTOP(CC(IC7),CC(IC15),CC(IC18),CC(IC19))                   ANH05100
C.................ZETA    F3Q      F4Q...... 
      CALL ASYZPE(CC(IC7),CC(IC15),CC(IC18))            
C.................ZETA    F3Q      F4Q...... 
      CALL ASYG00(CC(IC7),CC(IC15),CC(IC18))            
      GO TO 215                                                         ANH05110
C                                                                       ANH05120
C::::::::::::::::::::::::::::                                           ANH05130
C:::SYMMETRIC TOP MOLECULE:::                                           ANH05140
C::::::::::::::::::::::::::::                                           ANH05150
  504 CONTINUE                                                          ANH05160
      WRITE(3,36)                                                       ANH05170
   36 FORMAT(//,2X,' NOW YOU ARE IN SYMTOP'/)                           ANH05180
      IC19=ICKEEP                                                       ANH05190
      IC20=IC19+N3N*N3N                                                 ANH05200
      ICMAX=IC20+N3N*N3N                                                ANH05210
      WRITE(3,11) ICMAX,MAXCOR                                          ANH05220
C.................ZETA    AXY     F3Q      EE       CLTYP....           ANH05230
      CALL LTYPSY(CC(IC7),CC(IC8),CC(IC15),CC(IC19),CC(IC20))           ANH05240
C.................ZETA    AXY     F3Q      EE       RVCNST...           ANH05250
      CALL ROVSYM(CC(IC7),CC(IC8),CC(IC15),CC(IC19),CC(IC20))           ANH05260
C.................F3Q      EE.......                                    ANH05270
      CALL CENTSY(CC(IC15),CC(IC19))                                    ANH05280
      IF(IVDBL.LE.1) GO TO 206                                          ANH05290
C.................ZETA    AXY     F3Q      F4Q      CLTYP....           ANH05300
      CALL VLTYSY(CC(IC7),CC(IC8),CC(IC15),CC(IC18),CC(IC19))           ANH05310
  206 CONTINUE                                                          ANH05320
C.................ZETA    F3Q      F4Q      XIJ      GIJ......          ANH05330
      CALL SYMTOP(CC(IC7),CC(IC15),CC(IC18),CC(IC19),CC(IC20))          ANH05340
C.................ZETA    F3Q      F4Q......
      CALL SYMZPE(CC(IC7),CC(IC15),CC(IC18))
C.................ZETA    F3Q      F4Q......
      CALL SYMG00(CC(IC7),CC(IC15),CC(IC18))
      GO TO 215                                                         ANH05350
C                                                                       ANH05360
C::::::::::::::::::::::::::::                                           ANH05370
C:::SPHERICAL TOP MOLECULE:::                                           ANH05380
C::::::::::::::::::::::::::::                                           ANH05390
  505 CONTINUE                                                          ANH05400
      WRITE(3,37)                                                       ANH05410
   37 FORMAT(//,2X,' NOW YOU ARE IN SPHTOP'/)                           ANH05420
      IC19=ICKEEP                                                       ANH05430
      IC20=IC19+N3N*N3N                                                 ANH05440
      ICMAX=IC20+N3N*N3N                                                ANH05450
      WRITE(3,11) ICMAX,MAXCOR                                          ANH05460
C.................ZETA    AXY     F3Q      EE       RVCNST...           ANH05470
      CALL ROVSPH(CC(IC7),CC(IC8),CC(IC15),CC(IC19),CC(IC20))           ANH05480
C.................F3Q      EE.......                                    ANH05490
      CALL CENTSP(CC(IC15),CC(IC19))                                    ANH05500
C.................ZETA    F3Q      F4Q      XIJ      GIJ......          ANH05510
      CALL SPHTOP(CC(IC7),CC(IC15),CC(IC18),CC(IC19),CC(IC20))          ANH05520
C                                                                       ANH05530
C##########################################                             ANH05540
C   RETURN FOR ANOTHER ISOTOPIC MOLECULE                                ANH05550
C##########################################                             ANH05560
  215 CONTINUE                                                          ANH05570
      GO TO 400                                                         ANH05580
C                                                                       ANH05590
  220 CONTINUE                                                          ANH05600
      WRITE(6,13)
C
      STOP                                                              ANH05620
      END                                                               ANH05630
      SUBROUTINE DIST(IPRNT)                                            ANH05640
      IMPLICIT REAL*8 (A-H,O-Z)                                         ANH05650
      COMMON/VIB101/NATOM,N3N,NATRI,ILIN,NVIB                           ANH05660
      COMMON/VIB201/CHG(50),X(50),Y(50),Z(50),W(50)                     ANH05670
      COMMON/VIB202/R(1275)                                             ANH05680
      DATA A00 / 0.0D+00 /                                              ANH05690
    1 FORMAT(//,2X,' INTERATOMIC DISTANCE MATRIX'/)                     ANH05700
C                                                                       ANH05710
      IJ=0                                                              ANH05720
      DO 101 I=1,NATOM                                                  ANH05730
      DO 101 J=1,I                                                      ANH05740
      IJ=IJ+1                                                           ANH05750
      R(IJ)=A00                                                         ANH05760
      IF(I.EQ.J) GO TO 101                                              ANH05770
      XD=X(I)-X(J)                                                      ANH05780
      YD=Y(I)-Y(J)                                                      ANH05790
      ZD=Z(I)-Z(J)                                                      ANH05800
      R(IJ)=DSQRT(XD*XD+YD*YD+ZD*ZD)                                    ANH05810
  101 CONTINUE                                                          ANH05820
      IF(IPRNT.EQ.0) GO TO 201                                          ANH05830
      WRITE(6,1)                                                        ANH05840
      CALL PRINT(R,1275,NATOM,6)                                        ANH05850
  201 CONTINUE
C
      RETURN                                                            ANH05860
      END                                                               ANH05870
      SUBROUTINE FCONST(FX)                                             ANH05880
      IMPLICIT REAL*8 (A-H,O-Z)                                         ANH05890
      DIMENSION FX(N3N,N3N)                                             ANH05900
      COMMON/VIB101/NATOM,N3N,NATRI,ILIN,NVIB                           ANH05910
      DATA HE,A0 / 4.359813653D+00 , 0.529177249D+00 /  
      DATA ITAP15 / 15 /                                                ANH05930
    1 FORMAT(2I5)                                                       ANH05940
    2 FORMAT(3F20.10)                                                   ANH05950
    3 FORMAT(//,2X,' FX MATRIX (IN HARTREE/BOHR+2)'/)                   ANH05960
    4 FORMAT(//,2X,' FX MATRIX (IN MDYN/A)'/)                           ANH05970
C                                                                       ANH05980
C   FUNIT IS IN MDYN/A                                                  ANH05990
      FUNIT=HE/(A0*A0)                                                  ANH06000
C                                                                       ANH06010
C   READ IN ANALYTICAL SECOND DERIVATIVES                               ANH06020
C   WITH RESPECT TO NUCLEAR COORDINATE                                  ANH06030

c     write(6,*) ' Now you are in FCONST'
c     write(6,*) ' NATOM in fconst = ',natom
c     write(6,*) ' N3N in fconst   = ',n3n

      REWIND ITAP15                                                     ANH06040
      READ(ITAP15,1) NNATOM,NNSTOR                                      ANH06050

c     write(6,1) NNATOM,NNSTOR          

      READ(ITAP15,2) ((FX(I,J),J=1,N3N),I=1,N3N)                        ANH06060

c     write(6,2) ((FX(I,J),J=1,N3N),I=1,N3N)      

      REWIND ITAP15                                                     ANH06070
      WRITE(6,3)                                                        ANH06080
      CALL MATOUT(FX,N3N,N3N,N3N,N3N,6)                                 ANH06090
      DO 101 I=1,N3N                                                    ANH06100
      DO 101 J=1,N3N                                                    ANH06110
      FX(I,J)=FX(I,J)*FUNIT                                             ANH06120
  101 CONTINUE                                                          ANH06130
      WRITE(6,4)                                                        ANH06140
      CALL MATOUT(FX,N3N,N3N,N3N,N3N,6)                                 ANH06150
C                                                                       ANH06160
      RETURN                                                            ANH06170
      END                                                               ANH06180
      SUBROUTINE MOMENT(UX,EE)                                          ANH06190
      IMPLICIT REAL*8 (A-H,O-Z)                                         ANH06200
      DIMENSION UX(N3N,N3N),EE(N3N,N3N)                                 ANH06210
      DIMENSION T(6),PT(3),EG(3,3),FV1(3),FV2(3)                        ANH06220
      CHARACTER*1 RAXIS                                                 ANH06230
      CHARACTER*25 RTYPE                                                ANH06240
      COMMON/VIB101/NATOM,N3N,NATRI,ILIN,NVIB                           ANH06250
      COMMON/VIB103/PARA,WAVE,CONST,CYCL,CONV                           ANH06260
      COMMON/VIB104/AIXX,AIYY,AIZZ,SUMW                                 ANH06270
      COMMON/VIB106/ROTAA(3),ROTGC(3),ROTCM(3),ROTMH(3)                 ANH06280
      COMMON/VIB107/ITOP,IVSGL,IVDBL,IVTRP                              ANH06290
      COMMON/VIB108/IAXIS(3),NDEG(150),NDAB(150,5),IMAG(150)            ANH06300
      COMMON/VIB201/CHG(50),X(50),Y(50),Z(50),W(50)                     ANH06310
      COMMON/VIB203/IOFF(150),IPRNT
      COMMON/VIB206/RAXIS(3),RTYPE(6)                                   ANH06320
      DATA A00,ONE,TWO / 0.0D+00 , 1.0D+00 , 2.0D+00 /                  ANH06330
      DATA A0 / 0.529177249D+00 /
      DATA DLIMIT / 1.0D-02 /                                           ANH06350
      DATA RLIMIT / 1.0D-05 /                                           ANH06360
    1 FORMAT(//,2X,' CENTER OF MASS'/2X,' CMX',11X,' CMY',11X,' CMZ'//  ANH06370
     1 2X,3(5X,F10.7)/)                                                 ANH06380
    2 FORMAT(//,2X,' CARTESIAN COORDINATES W.R.T. CENTER OF MASS (IN A)'ANH06390
     1/5X,4H NO.,11X,2H X,18X,2H Y,18X,2H Z,18X,2H W/)                  ANH06400
    3 FORMAT(2X,I7,5X,4F20.10)                                          ANH06410
    4 FORMAT(//,2X,' INERTIA TENSOR'/)                                  ANH06420
    5 FORMAT(//,2X,' PRINCIPAL MOMENTS OF INERTIA AT EQUILIBRIUM GEOMETRANH06430
     1Y'/                                                               ANH06440
     2 4X,5H AXIS,9X,11H IN AMU.A+2,9X,10H IN G.CM+2,10X,               ANH06450
     3 8H IN CM-1,15X,7H IN MHZ/38X,8H (*D+39)/)                        ANH06460
    6 FORMAT(8X,A1,3X,3F20.10,F20.5)                                    ANH06470
    7 FORMAT(//,2X,' EIGENVECTOR OF INERTIA TENSOR'/)                   ANH06480
    8 FORMAT(//,2X,' THE ROTATIONAL TYPE OF THIS ISOTOPIC MOLECULE  = ',ANH06490
     1 2X,A25/)                                                         ANH06500
    9 FORMAT(/,2X,' ASYMMETRY PARAMETERS'/                              ANH06510
     1 2X,' B(P)    = ',F10.5/                                          ANH06520
     2 2X,' B(O)    = ',F10.5/                                          ANH06530
     3 2X,' KAPPA   = ',F10.5/)                                         ANH06540
   10 FORMAT(//,2X,' UX MATRIX'/)                                       ANH06550
   11 FORMAT(//,2X,' UX * UX(T) = E '/)                                 ANH06560
   12 FORMAT(//,2X,' CARTESIAN COORDINATES W.R.T. PRINCIPAL AXIS (IN BOHANH06570
     1R)'/5X,4H NO.,11X,2H X,18X,2H Y,18X,2H Z,18X,2H W/)               ANH06580
   13 FORMAT(//,2X,' CARTESIAN COORDINATES W.R.T. PRINCIPAL AXIS (IN A)'ANH06590
     1/5X,4H NO.,11X,2H X,18X,2H Y,18X,2H Z,18X,2H W/)                  ANH06600
C                                                                       ANH06610
C   THE CALCULATION OF CENTER OF MASS                                   ANH06620
      SUMW=A00                                                          ANH06630
      SUMWX=A00                                                         ANH06640
      SUMWY=A00                                                         ANH06650
      SUMWZ=A00                                                         ANH06660
      DO 101 I=1,NATOM                                                  ANH06670
      SUMW=SUMW+W(I)                                                    ANH06680
      SUMWX=SUMWX+W(I)*X(I)                                             ANH06690
      SUMWY=SUMWY+W(I)*Y(I)                                             ANH06700
      SUMWZ=SUMWZ+W(I)*Z(I)                                             ANH06710
  101 CONTINUE                                                          ANH06720
      CMX=SUMWX/SUMW                                                    ANH06730
      CMY=SUMWY/SUMW                                                    ANH06740
      CMZ=SUMWZ/SUMW                                                    ANH06750
      WRITE(6,1) CMX,CMY,CMZ                                            ANH06760
C   THE CARTESIAN COORDINATES W.R.T. CENTER OF MASS                     ANH06770
      WRITE(6,2)                                                        ANH06780
      DO 102 I=1,NATOM                                                  ANH06790
      X(I)=X(I)-CMX                                                     ANH06800
      Y(I)=Y(I)-CMY                                                     ANH06810
      Z(I)=Z(I)-CMZ                                                     ANH06820
      WRITE(6,3) I,X(I),Y(I),Z(I),W(I)                                  ANH06830
  102 CONTINUE                                                          ANH06840
C                                                                       ANH06850
C   THE CALCULATION OF INERTIA TENSOR                                   ANH06860
      CALL ZERO(T,6)                                                    ANH06870
      DO 104 I=1,NATOM                                                  ANH06880
      T(1)=T(1)+W(I)*(Y(I)**2+Z(I)**2)                                  ANH06890
      T(3)=T(3)+W(I)*(Z(I)**2+X(I)**2)                                  ANH06900
      T(6)=T(6)+W(I)*(X(I)**2+Y(I)**2)                                  ANH06910
      T(2)=T(2)-W(I)*X(I)*Y(I)                                          ANH06920
      T(4)=T(4)-W(I)*Z(I)*X(I)                                          ANH06930
      T(5)=T(5)-W(I)*Y(I)*Z(I)                                          ANH06940
  104 CONTINUE                                                          ANH06950
      AIXX=T(1)                                                         ANH06960
      AIYY=T(3)                                                         ANH06970
      AIZZ=T(6)                                                         ANH06980
      ILIN=0                                                            ANH06990
      IF(AIXX.LE.DLIMIT) ILIN=1                                         ANH07000
      IF(AIYY.LE.DLIMIT) ILIN=1                                         ANH07010
      IF(AIZZ.LE.DLIMIT) ILIN=1                                         ANH07020
      WRITE(6,4)                                                        ANH07030
      CALL PRINT(T,6,3,6)                                               ANH07040
      NVIB=N3N-6+ILIN                                                   ANH07050
C                                                                       ANH07060
C   THE CALCULATION OF PRINCIPAL MOMENTS OF INERTIA                     ANH07070
      CALL RSP(3,3,6,T,PT,1,EG,FV1,FV2)                                 ANH07080
      CALL ZERO(ROTAA,3)                                                ANH07090
      CALL ZERO(ROTGC,3)                                                ANH07100
      CALL ZERO(ROTCM,3)                                                ANH07110
      CALL ZERO(ROTMH,3)                                                ANH07120
      WRITE(6,5)                                                        ANH07130
      II=0                                                              ANH07140
      DO 106 I=1,3                                                      ANH07150
      PA=A00                                                            ANH07160
      PB=A00                                                            ANH07170
      PC=A00                                                            ANH07180
      IF(DABS(PT(I)).LE.DLIMIT) GO TO 201                               ANH07190
      II=II+1                                                           ANH07200
      PA=PT(I)*PARA                                                     ANH07210
      PB=CONST/PT(I)                                                    ANH07220
      PC=CYCL/PT(I)                                                     ANH07230
      ROTAA(II)=PT(I)                                                   ANH07240
      ROTGC(II)=PA                                                      ANH07250
      ROTCM(II)=PB                                                      ANH07260
      ROTMH(II)=PC                                                      ANH07270
  201 CONTINUE                                                          ANH07280
      WRITE(6,6) RAXIS(I),PT(I),PA,PB,PC                                ANH07290
  106 CONTINUE                                                          ANH07300
      WRITE(6,7)                                                        ANH07310
      CALL MATOUT(EG,3,3,3,3,6)                                         ANH07320
C                                                                       ANH07330
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%              ANH07340
C%  ASSIGN THE ROTATIONAL TYPE OF THE ISOTOPIC MOLECULE  %              ANH07350
C%  ITOP IS A PARAMETER FOR ROTATIONAL TOP               %              ANH07360
C%  ITOP  = 1  : DIATOMIC MOLECULE                       %              ANH07370
C%  ITOP  = 2  : LINEAR POLYATOMIC MOLECULE              %              ANH07380
C%  ITOP  = 3  : ASYMMETRIC TOP                          %              ANH07390
C%  ITOP  = 4  : SYMMETRIC TOP (PROLATE)                 %              ANH07400
C%  ITOP  = 5  : SYMMETRIC TOP (OBLATE)                  %              ANH07410
C%  ITOP  = 6  : SPHERICAL TOP                           %              ANH07420
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%              ANH07430
      ITOP=3                                                            ANH07440
      IF(DABS(ROTAA(2)-ROTAA(3)).LE.RLIMIT) ITOP=4                      ANH07450
      IF(DABS(ROTAA(1)-ROTAA(2)).LE.RLIMIT) ITOP=5                      ANH07460
      IF(DABS(ROTAA(1)-ROTAA(2)).LE.RLIMIT.AND.                         ANH07470
     *   DABS(ROTAA(2)-ROTAA(3)).LE.RLIMIT) ITOP=6                      ANH07480
      IF(ILIN.NE.0) ITOP=2                                              ANH07490
      IF(NATOM.EQ.2) ITOP=1                                             ANH07500
      WRITE(6,8) RTYPE(ITOP)                                            ANH07510
C                                                                       ANH07520
C   ASSIGN PRINCIPAL AXIS                                               ANH07530
      IAXIS(1)=1                                                        ANH07540
      IAXIS(2)=2                                                        ANH07550
      IAXIS(3)=3                                                        ANH07560
C                                                                       ANH07570
C   FORM UX MATRIX                                                      ANH07580
      CALL ZERO(UX,N3N*N3N)                                             ANH07590
C                                                                       ANH07600
C%%%%%%%%%%%%%%%%%%                                                     ANH07610
C%  FOR DIATOMIC  %                                                     ANH07620
C%%%%%%%%%%%%%%%%%%                                                     ANH07630
      IF(ITOP.EQ.1) THEN                                                ANH07640
        CALL ZERO(EG,9)                                                 ANH07650
        EG(1,1)=ONE                                                     ANH07660
        EG(2,2)=ONE                                                     ANH07670
        EG(3,3)=ONE                                                     ANH07680
      END IF                                                            ANH07690
C                                                                       ANH07700
C   FOR LINEAR POLYATOMIC                                               ANH07710
      IF(ITOP.EQ.2) THEN                                                ANH07720
        DO 108 I=1,3                                                    ANH07730
        TG=EG(I,1)                                                      ANH07740
        EG(I,1)=EG(I,3)                                                 ANH07750
        EG(I,3)=TG                                                      ANH07760
  108 CONTINUE                                                          ANH07770
      END IF                                                            ANH07780
C                                                                       ANH07790
C%%%%%%%%%%%%%%%%%%%%%%%%                                               ANH07800
C%  FOR ASYMMETRIC TOP  %                                               ANH07810
C%%%%%%%%%%%%%%%%%%%%%%%%                                               ANH07820
C     FOR PROLATE TOP RAY'S KAPPA = -1                                  ANH07830
C     FOR OBLATE  TOP RAY'S KAPPA = +1                                  ANH07840
      IF(ITOP.EQ.3) THEN                                                ANH07850
        BP=(ROTCM(3)-ROTCM(2))/(ROTCM(1)*TWO-ROTCM(2)-ROTCM(3))         ANH07860
        BO=(ROTCM(1)-ROTCM(2))/(ROTCM(3)*TWO-ROTCM(2)-ROTCM(1))         ANH07870
        RAYK=(ROTCM(2)*TWO-ROTCM(1)-ROTCM(3))/(ROTCM(1)-ROTCM(3))       ANH07880
        WRITE(6,9) BP,BO,RAYK                                           ANH07890
C     FOR AN ASYMMETRIC TOP THE COORDINATE SYSTEM WILL BE ROTATED       ANH07900
C     TO SATISFY THE FOLLOWING RELATIONSHIPS                            ANH07910
C       I(X)<I(Y)<I(Z)                                                  ANH07920
C       I(A)<I(B)<I(C)                                                  ANH07930
C         A >  B >  C                                                   ANH07940
C                                                                       ANH07950
        IAXIS(1)=1                                                      ANH07960
        IAXIS(2)=2                                                      ANH07970
        IAXIS(3)=3                                                      ANH07980
      END IF                                                            ANH07990
C                                                                       ANH08000
C%%%%%%%%%%%%%%%%%%%%%%%                                                ANH08010
C%  FOR SYMMETRIC TOP  %                                                ANH08020
C%%%%%%%%%%%%%%%%%%%%%%%                                                ANH08030
C     FOR A SYMMETRIC TOP THE MOLECULAR AXIS IS ALONG THE Z AXIS        ANH08040
      IF(ITOP.EQ.4.OR.ITOP.EQ.5) THEN                                   ANH08050
        CALL ZERO(EG,9)                                                 ANH08060
        EG(1,1)=ONE                                                     ANH08070
        EG(2,2)=ONE                                                     ANH08080
        EG(3,3)=ONE                                                     ANH08090
      END IF                                                            ANH08100
C     FOR A PROLATE TOP                                                 ANH08110
C        I(X)<I(Y)=I(Z)                                                 ANH08120
C        I(A)<I(B)=I(C)                                                 ANH08130
C          A >  B =  C                                                  ANH08140
      IF(ITOP.EQ.4) THEN                                                ANH08150
        IAXIS(1)=3                                                      ANH08160
        IAXIS(2)=2                                                      ANH08170
        IAXIS(3)=1                                                      ANH08180
      END IF                                                            ANH08190
C                                                                       ANH08200
C     FOR AN OBLATE TOP                                                 ANH08210
C        I(X)=I(Y)<I(Z)                                                 ANH08220
C        I(A)=I(B)<I(C)                                                 ANH08230
C          A =  B >  C                                                  ANH08240
      IF(ITOP.EQ.5) THEN                                                ANH08250
        IAXIS(1)=1                                                      ANH08260
        IAXIS(2)=2                                                      ANH08270
        IAXIS(3)=3                                                      ANH08280
      END IF                                                            ANH08290
C                                                                       ANH08300
C%%%%%%%%%%%%%%%%%%%%%%%                                                ANH08310
C%  FOR SPHERICAL TOP  %                                                ANH08320
C%%%%%%%%%%%%%%%%%%%%%%%                                                ANH08330
      IF(ITOP.EQ.6) THEN                                                ANH08340
        CALL ZERO(EG,9)                                                 ANH08350
        EG(1,1)=ONE                                                     ANH08360
        EG(2,2)=ONE                                                     ANH08370
        EG(3,3)=ONE                                                     ANH08380
        IAXIS(1)=1                                                      ANH08390
        IAXIS(2)=2                                                      ANH08400
        IAXIS(3)=3                                                      ANH08410
      END IF                                                            ANH08420
C                                                                       ANH08430
C   FORM THE UX MATRIX                                                  ANH08440
      DO 111 IABC=1,NATOM                                               ANH08450
      IPOS=(IABC-1)*3                                                   ANH08460
      DO 110 I=1,3                                                      ANH08470
      DO 110 J=1,3                                                      ANH08480
      II=IPOS+I                                                         ANH08490
      JJ=IPOS+J                                                         ANH08500
      UX(II,JJ)=EG(I,J)                                                 ANH08510
  110 CONTINUE                                                          ANH08520
  111 CONTINUE                                                          ANH08530
      WRITE(6,10)                                                       ANH08540
      CALL MATOUT(UX,N3N,N3N,N3N,N3N,6)                                 ANH08550
C                                                                       ANH08560
C   CHECK UNITARITY OF UX MATRIX                                        ANH08570
      IF(IPRNT.GE.2) THEN
        CALL MTXMPY(UX,N3N,UX,N3N,EE,N3N,EE,N3N,N3N,3)    
        WRITE(6,11)                                      
        CALL MATOUT(EE,N3N,N3N,N3N,N3N,6)               
      END IF
C                                                                       ANH08610
C   THE CARTESIAN COORDINATES W.R.T. PRINCIPAL AXIS                     ANH08620
      WRITE(6,12)                                                       ANH08630
      DO 112 I=1,NATOM                                                  ANH08640
      XT=EG(1,1)*X(I)+EG(2,1)*Y(I)+EG(3,1)*Z(I)                         ANH08650
      YT=EG(1,2)*X(I)+EG(2,2)*Y(I)+EG(3,2)*Z(I)                         ANH08660
      ZT=EG(1,3)*X(I)+EG(2,3)*Y(I)+EG(3,3)*Z(I)                         ANH08670
      X(I)=XT                                                           ANH08680
      Y(I)=YT                                                           ANH08690
      Z(I)=ZT                                                           ANH08700
      XT=XT/A0                                                          ANH08710
      YT=YT/A0                                                          ANH08720
      ZT=ZT/A0                                                          ANH08730
      WRITE(6,3) I,XT,YT,ZT,W(I)                                        ANH08740
  112 CONTINUE                                                          ANH08750
      WRITE(6,13)                                                       ANH08760
      DO 113 I=1,NATOM                                                  ANH08770
      WRITE(6,3) I,X(I),Y(I),Z(I),W(I)                                  ANH08780
  113 CONTINUE
C                                                                       ANH08790
      RETURN                                                            ANH08800
      END                                                               ANH08810
      SUBROUTINE VIBFX(FX,ELX,ELXM,UX,EPX,EPXM,EE,VALU,FV1,FV2,FXM,DD,  ANH08820
     1                 CC)                                              ANH08830
C   THE NORMAL COORDINATE ANALYSIS FOR CARTESIAN COORDINATE SYSTEM      ANH08840
      IMPLICIT REAL*8 (A-H,O-Z)                                         ANH08850
      DIMENSION FX(N3N,N3N),ELX(N3N,N3N),ELXM(N3N,N3N)                  ANH08860
      DIMENSION UX(N3N,N3N),EPX(N3N,N3N),EPXM(N3N,N3N)                  ANH08870
      DIMENSION EE(N3N,N3N),DD(N3N,N3N),CC(1)                           ANH08880
      DIMENSION VALU(N3N),FV1(N3N),FV2(N3N),FXM(NATRI)                  ANH08890
      DIMENSION IATOMX(10)                                              ANH08900
      COMMON/VIB101/NATOM,N3N,NATRI,ILIN,NVIB                           ANH08910
      COMMON/VIB107/ITOP,IVSGL,IVDBL,IVTRP                              ANH08920
      COMMON/VIB108/IAXIS(3),NDEG(150),NDAB(150,5),IMAG(150)            ANH08930
      COMMON/VIB201/CHG(50),X(50),Y(50),Z(50),W(50)                     ANH08940
      COMMON/VIB203/IOFF(150),IPRNT                                     ANH08950
      COMMON/VIB204/SQM(150),ROOT(150),FREQ(150)                        ANH08960
      DATA A00,ONE / 0.0D+00 , 1.0D+00 /                                ANH08970
      DATA DLIMIT / 1.0D-01 /                                           ANH08980
      DATA ELIMIT / 1.0D-04 /                                           ANH08990
      DATA WLIMIT / 3.0D+00 /                                           ANH09000
      DATA XLIMIT / 1.0D-05 /                                           ANH09010
    1 FORMAT(//,2X,' THE NORMAL COORDINATE ANALYSIS BY FX MATRIX METHOD'ANH09020
     1/)                                                                ANH09030
    2 FORMAT(//,2X,' SECULAR EQUATIONS'/2X,' FXM MATRIX'/)              ANH09040
    3 FORMAT(//,2X,' *****************************'/                    ANH09050
     1          2X,' ***VIBRATIONAL FREQUENCIES***'/                    ANH09060
     2          2X,' *****************************')                    ANH09070
    4 FORMAT(//,2X,' LX MATRIX'/)                                       ANH09080
    5 FORMAT(//,2X,' LX * LX(T) = M(-1)'/)                              ANH09090
    6 FORMAT(//,2X,' LAMDA = LX(T) * FX * LX'/)                         ANH09100
    7 FORMAT(//,2X,' LX'' MATRIX'/)                                     ANH09110
    8 FORMAT(//,2X,' LX'' * LX''(T) = M(-1)'/)                          ANH09120
    9 FORMAT(//,2X,' LXM'' MATRIX'/)                                    ANH09130
   10 FORMAT(//,2X,' LXM'' * LXM''(T) = E'/)                            ANH09140
   11 FORMAT(//,2X,' ANALYSIS OF VIBRATIONAL MODES'//                   ANH09150
     1          2X,' NO.OF NON-DEGENERATE MODES    = ',I5/              ANH09160
     1          2X,' NO.OF DOUBLY-DEGENERATE MODES = ',I5/              ANH09170
     2          2X,' NO.OF TRIPLY-DEGENERATE MODES = ',I5/              ANH09180
     3          2X,' NO.OF TOTAL VIBRATIONAL MODES = ',I5/)             ANH09190
   12 FORMAT(//,2X,' ROTATED LXM'' MATRIX'/)                            ANH09200
C                                                                       ANH09210
      IF(IPRNT.NE.0)                                                    ANH09220
     *WRITE(6,1)                                                        ANH09230
C                                                                       ANH09240
C   THE FORMATION OF THE SECULAR EQUATION                               ANH09250
      DO 101 I=1,NATOM                                                  ANH09260
      SQM(I)=DSQRT(W(I))                                                ANH09270
  101 CONTINUE
      DO 102 I=1,N3N                                                    ANH09280
      II=I/3                                                            ANH09290
      IF(I.GT.3*II) II=II+1                                             ANH09300
      DO 102 J=I,N3N                                                    ANH09310
      JJ=J/3                                                            ANH09320
      IF(J.GT.3*JJ) JJ=JJ+1                                             ANH09330
      IJ=I+IOFF(J)                                                      ANH09340
      FXM(IJ)=(ONE/SQM(II))*FX(I,J)*(ONE/SQM(JJ))                       ANH09350
  102 CONTINUE                                                          ANH09360
      WRITE(6,2)                                                        ANH09370
      CALL PRINT(FXM,NATRI,N3N,6)                                       ANH09380
C                                                                       ANH09390
C   THE NORMAL COORDINATE ANALYSIS                                      ANH09400
      WRITE(6,3)                                                        ANH09410
      IF(ITOP.EQ.2) THEN                                                ANH09420
        NAA=IOFF(NATOM+1)                                               ANH09430
        IC1=1                                                           ANH09440
        IC2=IC1+NAA                                                     ANH09450
        IC3=IC2+NATOM                                                   ANH09460
        IC4=IC3+NATOM*NATOM                                             ANH09470
        ICMAX=IC4+N3N                                                   ANH09480
C                                                                       ANH09490
C...................                         FXMX    VALUX   ELXMX...   ANH09500
        CALL NORMLI(ELXM,EE,VALU,FV1,FV2,FXM,CC(IC1),CC(IC2),CC(IC3),   ANH09510
C...................NORD........                                        ANH09520
     1              CC(IC4),NAA)                                        ANH09530
      ELSE                                                              ANH09540
        CALL NORMFX(ELXM,EE,VALU,FV1,FV2,FXM)                           ANH09550
      END IF                                                            ANH09560
C                                                                       ANH09570
C   THE CALCULATION OF LX MATRIX                                        ANH09580
      DO 103 I=1,N3N                                                    ANH09590
      II=I/3                                                            ANH09600
      IF(I.GT.3*II) II=II+1                                             ANH09610
      DO 103 J=1,N3N                                                    ANH09620
      ELX(I,J)=(ONE/SQM(II))*ELXM(I,J)                                  ANH09630
  103 CONTINUE                                                          ANH09640
      WRITE(6,4)                                                        ANH09650
      CALL FRQOUT(ELX,FREQ,N3N,N3N,N3N,N3N,6)                           ANH09660
      IF(IPRNT.LE.2) GO TO 201                                          ANH09670
      CALL MTXMPY(ELX,N3N,ELX,N3N,EE,N3N,EE,N3N,N3N,3)                  ANH09680
      WRITE(6,5)                                                        ANH09690
      CALL MATOUT(EE,N3N,N3N,N3N,N3N,6)                                 ANH09700
C                                                                       ANH09710
C   TRANSFORM FX MATRIX FOR A TEST                                      ANH09720
      CALL MTXMPY(ELX,N3N,FX,N3N,EE,N3N,DD,N3N,N3N,5)                   ANH09730
      WRITE(6,6)                                                        ANH09740
      CALL MATOUT(EE,N3N,N3N,N3N,N3N,6)                                 ANH09750
C                                                                       ANH09760
C   THE CALCULATION OF LX' MATRIX                                       ANH09770
C==============================================                         ANH09780
C===FOR ASYMMETRIC TOP AND  LINEAR MOLECULES===                         ANH09790
C==============================================                         ANH09800
  201 CONTINUE                                                          ANH09810
      CALL MTXMPY(UX,N3N,ELX,N3N,EPX,N3N,33,N3N,N3N,2)                  ANH09820
      WRITE(6,7)                                                        ANH09830
      CALL FRQOUT(EPX,FREQ,N3N,N3N,N3N,N3N,6)                           ANH09840
      IF(IPRNT.GE.2) THEN
        CALL MTXMPY(EPX,N3N,EPX,N3N,EE,N3N,EE,N3N,N3N,3)    
        WRITE(6,8)                                         
        CALL MATOUT(EE,N3N,N3N,N3N,N3N,6)                 
      END IF
C                                                                       ANH09880
C   THE CALCULATION OF LXM' MATRIX                                      ANH09890
      DO 104 I=1,N3N                                                    ANH09900
      II=I/3                                                            ANH09910
      IF(I.GT.3*II) II=II+1                                             ANH09920
      DO 104 J=1,N3N                                                    ANH09930
      EPXM(I,J)=EPX(I,J)*SQM(II)                                        ANH09940
  104 CONTINUE                                                          ANH09950
      WRITE(6,9)                                                        ANH09960
      CALL FRQOUT(EPXM,FREQ,N3N,N3N,N3N,N3N,6)                          ANH09970
      IF(IPRNT.GE.2) THEN
        CALL MTXMPY(EPXM,N3N,EPXM,N3N,EE,N3N,EE,N3N,N3N,3)    
        WRITE(6,10)                                          
        CALL MATOUT(EE,N3N,N3N,N3N,N3N,6)                     
      END IF
C                                                                       ANH10010
C   FIND DEGENERACY OF NORMAL MODES                                     ANH10020
  202 CONTINUE                                                          ANH10030
      DO 105 I=1,N3N                                                    ANH10040
      NDEG(I)=0                                                         ANH10050
      DO 105 J=1,5                                                      ANH10060
      NDAB(I,J)=0                                                       ANH10070
  105 CONTINUE                                                          ANH10080
      DO 107 I=1,NVIB                                                   ANH10090
      FRQI=FREQ(I)                                                      ANH10100
      IF(DABS(FRQI).LE.WLIMIT) GO TO 107                                ANH10110
      II=0                                                              ANH10120
      DO 106 J=1,NVIB                                                   ANH10130
      FRQJ=FREQ(J)                                                      ANH10140
      IF(DABS(FRQJ).LE.WLIMIT) GO TO 106                                ANH10150
      IF(DABS(FRQI-FRQJ).GT.DLIMIT) GO TO 106                           ANH10160
      II=II+1                                                           ANH10170
      NDAB(I,II)=J                                                      ANH10180
  106 CONTINUE                                                          ANH10190
      NDEG(I)=II                                                        ANH10200
  107 CONTINUE                                                          ANH10210
C                                                                       ANH10220
C   FIND NUMBERS OF DEGENERATE PAIRS                                    ANH10230
      IVDBL=0                                                           ANH10240
      IVTRP=0                                                           ANH10250
      DO 108 I=1,N3N                                                    ANH10260
      IF(NDEG(I).EQ.2) IVDBL=IVDBL+1                                    ANH10270
      IF(NDEG(I).EQ.3) IVTRP=IVTRP+1                                    ANH10280
  108 CONTINUE                                                          ANH10290
      IVDBL=IVDBL/2                                                     ANH10300
      IVTRP=IVTRP/3                                                     ANH10310
      IVSGL=NVIB-IVDBL*2-IVTRP*3                                        ANH10320
      WRITE(6,11) IVSGL,IVDBL,IVTRP,NVIB                                ANH10330
C                                                                       ANH10340
      IF(ITOP.EQ.2) GO TO 205                                           ANH10350
      IF(IVDBL.EQ.0) GO TO 205                                          ANH10360
C                                                                       ANH10370
C   ROTATE EIGENVECTORS FOR DEGENERATE MODES                            ANH10380
C==================================                                     ANH10390
C===FOR A SYMMETRIC TOP MOLECULE===                                     ANH10400
C==================================                                     ANH10410
C   FIND ALL ATOMS ON THE C3V AXIS                                      ANH10420
      II=0                                                              ANH10430
      DO 110 I=1,NATOM                                                  ANH10440
      XYSUM=X(I)*X(I)+Y(I)*Y(I)                                         ANH10450
      IF(XYSUM.GT.XLIMIT) GO TO 110                                     ANH10460
      II=II+1                                                           ANH10470
      IATOMX(II)=I                                                      ANH10480
  110 CONTINUE                                                          ANH10500
      IAXTOT=II                                                         ANH10510
C                                                                       ANH10530
      IASAV=0                                                           ANH10540
      DO 120 IVIB=1,NVIB                                                ANH10550
      NN=NDEG(IVIB)                                                     ANH10560
      IF(NN.NE.2) GO TO 120                                             ANH10570
      IA=NDAB(IVIB,1)                                                   ANH10580
      IB=NDAB(IVIB,2)                                                   ANH10590
      IF(IASAV.EQ.IA) GO TO 120                                         ANH10600
      IASAV=IA                                                          ANH10610
      DO 111 I=1,IAXTOT                                                 ANH10620
      IAX=(IATOMX(I)-1)*3+1                                             ANH10630
      EP1=EPXM(IAX,IA)                                                  ANH10640
      EP2=EPXM(IAX,IB)                                                  ANH10650
      DNORM=DSQRT(EP1*EP1+EP2*EP2)                                      ANH10660
      IF(DNORM.LT.ELIMIT) GO TO 111                                     ANH10670
      XNORM=ONE/DNORM                                                   ANH10680
      EP1=-EP1*XNORM                                                    ANH10690
      EP2= EP2*XNORM                                                    ANH10700
      GO TO 203                                                         ANH10710
  111 CONTINUE                                                          ANH10720
  203 CONTINUE                                                          ANH10730
      DO 113 I=1,N3N                                                    ANH10740
      FV1(I)=EPXM(I,IA)                                                 ANH10750
      FV2(I)=EPXM(I,IB)                                                 ANH10760
  113 CONTINUE                                                          ANH10770
      DO 114 I=1,N3N                                                    ANH10780
      EPXM(I,IA)=EP1*FV1(I)-EP2*FV2(I)                                  ANH10790
      EPXM(I,IB)=EP2*FV1(I)+EP1*FV2(I)                                  ANH10800
  114 CONTINUE                                                          ANH10810
      TESTX=EPXM(IAX,IA)+EPXM(IAX+1,IA)                                 ANH10820
      IF(TESTX.LT.A00) THEN                                             ANH10830
        DO 115 I=1,N3N                                                  ANH10840
        EPXM(I,IA)=-EPXM(I,IA)                                          ANH10850
  115   CONTINUE                                                        ANH10860
      END IF                                                            ANH10870
      TESTX=EPXM(IAX,IB)+EPXM(IAX+1,IB)                                 ANH10880
      IF(TESTX.LT.A00) THEN                                             ANH10890
        DO 116 I=1,N3N                                                  ANH10900
        EPXM(I,IB)=-EPXM(I,IB)                                          ANH10910
  116   CONTINUE                                                        ANH10920
      END IF                                                            ANH10930
  120 CONTINUE                                                          ANH10940
      WRITE(6,12)                                                       ANH10950
      CALL FRQOUT(EPXM,FREQ,N3N,N3N,N3N,N3N,6)                          ANH10960
      IF(IPRNT.GE.2) THEN
        CALL MTXMPY(EPXM,N3N,EPXM,N3N,EE,N3N,EE,N3N,N3N,3)      
        WRITE(6,10)                                            
        CALL MATOUT(EE,N3N,N3N,N3N,N3N,6)                     
      END IF
C                                                                       ANH11000
C   THE CALCULATION OF LX' MATRIX                                       ANH11010
      DO 121 I=1,N3N                                                    ANH11020
      II=I/3                                                            ANH11030
      IF(I.GT.3*II) II=II+1                                             ANH11040
      DO 121 J=1,N3N                                                    ANH11050
      EPX(I,J)=(ONE/SQM(II))*EPXM(I,J)                                  ANH11060
  121 CONTINUE                                                          ANH11070
      WRITE(6,7)                                                        ANH11080
      CALL FRQOUT(EPX,FREQ,N3N,N3N,N3N,N3N,6)                           ANH11090
      IF(IPRNT.GE.2) THEN
        CALL MTXMPY(EPX,N3N,EPX,N3N,EE,N3N,EE,N3N,N3N,3)    
        WRITE(6,8)                                         
        CALL MATOUT(EE,N3N,N3N,N3N,N3N,6)                 
      END IF
      DO 122 I=1,N3N                                                    ANH11130
      DO 122 J=1,N3N                                                    ANH11140
      ELX(I,J)=EPX(I,J)                                                 ANH11150
  122 CONTINUE                                                          ANH11160
C                                                                       ANH11170
  205 CONTINUE                                                          ANH11180
      RETURN                                                            ANH11190
      END                                                               ANH11200
      SUBROUTINE NORMLI(ELXM,EE,VALU,FV1,FV2,FXM,FXMX,VALUX,ELXMX,NORD, ANH11210
     1                  NAA)                                            ANH11220
      IMPLICIT REAL*8 (A-H,O-Z)                                         ANH11230
      DIMENSION ELXM(N3N,N3N),EE(N3N,N3N)                               ANH11240
      DIMENSION VALU(N3N),FV1(N3N),FV2(N3N),FXM(NATRI)                  ANH11250
      DIMENSION FXMX(NAA),VALUX(NATOM),ELXMX(NATOM,NATOM)               ANH11260
      DIMENSION NORD(1)                                                 ANH11270
      COMMON/VIB101/NATOM,N3N,NATRI,ILIN,NVIB                           ANH11280
      COMMON/VIB103/PARA,WAVE,CONST,CYCL,CONV                           ANH11290
      COMMON/VIB108/IAXIS(3),NDEG(150),NDAB(150,5),IMAG(150)            ANH11300
      COMMON/VIB203/IOFF(150),IPRNT                                     ANH11310
      COMMON/VIB204/SQM(150),ROOT(150),FREQ(150)                        ANH11320
      DATA A00 / 0.0D+00 /                                              ANH11330
    1 FORMAT(//,2X,' FXMX MATRIX, IXYZ = ',I5/)                         ANH11340
    2 FORMAT(//,2X,' ELXM MATRIX BEFORE SORTING'/)                      ANH11350
    3 FORMAT(/,2X,' EIGENVALUES AND FREQUENCIES'/                       ANH11360
     1 13X,8H IN A.U.,9X,8H IN CM-1/)                                   ANH11370
    4 FORMAT(2X,I5,F15.6,2X,F15.6)                                      ANH11380
    5 FORMAT(//,2X,' LXM MATRIX'/)                                      ANH11390
    6 FORMAT(//,2X,' LXM * LXM(T) = E '/)                               ANH11400
C                                                                       ANH11410
C   THE CALCULATION OF EIGENVALUE AND EIGENVECTOR                       ANH11420
      DO 105 IXYZ=1,3                                                   ANH11430
C                                                                       ANH11440
      IIJJ=0                                                            ANH11450
      DO 101 IABC=1,NATOM                                               ANH11460
      IADD=(IABC-1)*3                                                   ANH11470
      IX=IADD+IXYZ                                                      ANH11480
      DO 101 JABC=1,IABC                                                ANH11490
      JADD=(JABC-1)*3                                                   ANH11500
      JX=JADD+IXYZ                                                      ANH11510
      IJX=IOFF(IX)+JX                                                   ANH11520
      IIJJ=IIJJ+1                                                       ANH11530
      FXMX(IIJJ)=FXM(IJX)                                               ANH11540
  101 CONTINUE                                                          ANH11550
      IF(IPRNT.LE.2) GO TO 201                                          ANH11560
      WRITE(6,1) IXYZ                                                   ANH11570
      CALL PRINT(FXMX,NAA,NATOM,6)                                      ANH11580
  201 CONTINUE                                                          ANH11590
      CALL RSP(NATOM,NATOM,NAA,FXMX,VALUX,1,ELXMX,FV1,FV2)              ANH11600
C                                                                       ANH11610
      DO 103 IABC=1,NATOM                                               ANH11620
      IADD=(IABC-1)*3                                                   ANH11630
      IX=IADD+IXYZ                                                      ANH11640
      VALU(IX)=VALUX(IABC)                                              ANH11650
      DO 103 JABC=1,NATOM                                               ANH11660
      JADD=(JABC-1)*3                                                   ANH11670
      JX=JADD+IXYZ                                                      ANH11680
      ELXM(IX,JX)=ELXMX(IABC,JABC)                                      ANH11690
  103 CONTINUE                                                          ANH11700
  105 CONTINUE                                                          ANH11710
      IF(IPRNT.LE.2) GO TO 202                                          ANH11720
      WRITE(6,2)                                                        ANH11730
      CALL MATOUT(ELXM,N3N,N3N,N3N,N3N,6)                               ANH11740
C                                                                       ANH11750
C   SORT EIGENVALUES AND EIGENVECTORS                                   ANH11760
  202 CONTINUE                                                          ANH11770
      CALL REORDR(VALU,NORD,N3N)                                        ANH11780
      WRITE(6,3)                                                        ANH11790
      DO 106 I=1,N3N                                                    ANH11800
      ROOT(I)=VALU(I)                                                   ANH11810
      IMAG(I)=0                                                         ANH11820
      IF(VALU(I).LT.A00) IMAG(I)=1                                      ANH11830
      IF(VALU(I).LT.A00) VALU(I)=-VALU(I)                               ANH11840
      AVAL=VALU(I)                                                      ANH11850
      FREQ(I)=WAVE*DSQRT(AVAL)                                          ANH11860
      WRITE(6,4) I,ROOT(I),FREQ(I)                                      ANH11870
  106 CONTINUE
C                                                                       ANH11880
C   EIGENVECTOR CORRESPONDS TO LXM VECTOR                               ANH11890
      DO 108 II=1,N3N                                                   ANH11900
      NN=NORD(II)                                                       ANH11910
      DO 107 I=1,N3N                                                    ANH11920
      EE(I,II)=ELXM(I,NN)                                               ANH11930
  107 CONTINUE
  108 CONTINUE                                                          ANH11940
      DO 109 I=1,N3N                                                    ANH11950
      DO 109 J=1,N3N                                                    ANH11960
      ELXM(I,J)=EE(I,J)                                                 ANH11970
  109 CONTINUE                                                          ANH11980
      WRITE(6,5)                                                        ANH11990
      CALL FRQOUT(ELXM,FREQ,N3N,N3N,N3N,N3N,6)                          ANH12000
C                                                                       ANH12010
C   THE EXAMINATION OF ABOVE CALCULATION                                ANH12020
      IF(IPRNT.LE.2) GO TO 203                                          ANH12030
      CALL MTXMPY(ELXM,N3N,ELXM,N3N,EE,N3N,EE,N3N,N3N,3)                ANH12040
      WRITE(6,6)                                                        ANH12050
      CALL MATOUT(EE,N3N,N3N,N3N,N3N,6)                                 ANH12060
C                                                                       ANH12070
  203 CONTINUE                                                          ANH12080
      RETURN                                                            ANH12090
      END                                                               ANH12100
      SUBROUTINE ORDER(A,NORD,N)                                        ANH12110
      IMPLICIT REAL*8 (A-H,O-Z)                                         ANH12120
      DIMENSION A(N),NORD(N)                                            ANH12130
C                                                                       ANH12140
      DO 101 I=1,N                                                      ANH12150
      NORD(I)=I                                                         ANH12160
  101 CONTINUE                                                          ANH12170
      IF(N.LE.1) GO TO 201                                              ANH12180
      N1=N-1                                                            ANH12190
      DO 103 I=1,N1                                                     ANH12200
      IR=I+1                                                            ANH12210
      DO 102 J=IR,N                                                     ANH12220
      IF(A(I).LT.A(J)) GO TO 102                                        ANH12230
      AA=A(I)                                                           ANH12240
      A(I)=A(J)                                                         ANH12250
      A(J)=AA                                                           ANH12260
      NN=NORD(I)                                                        ANH12270
      NORD(I)=NORD(J)                                                   ANH12280
      NORD(J)=NN                                                        ANH12290
  102 CONTINUE                                                          ANH12300
  103 CONTINUE                                                          ANH12310
  201 CONTINUE
      RETURN                                                            ANH12320
C                                                                       ANH12330
      ENTRY REORDR(A,NORD,N)                                            ANH12340
      DO 105 I=1,N                                                      ANH12350
      NORD(I)=I                                                         ANH12360
  105 CONTINUE
      IF(N.LE.1) GO TO 205                                              ANH12370
      N1=N-1                                                            ANH12380
      DO 107 I=1,N1                                                     ANH12390
      IR=I+1                                                            ANH12400
      DO 106 J=IR,N                                                     ANH12410
      IF(A(I).GE.A(J)) GO TO 106                                        ANH12420
      AA=A(I)                                                           ANH12430
      A(I)=A(J)                                                         ANH12440
      A(J)=AA                                                           ANH12450
      NN=NORD(I)                                                        ANH12460
      NORD(I)=NORD(J)                                                   ANH12470
      NORD(J)=NN                                                        ANH12480
  106 CONTINUE                                                          ANH12490
  107 CONTINUE                                                          ANH12500
  205 CONTINUE                                                          ANH12510
C                                                                       ANH12520
      RETURN                                                            ANH12530
      END                                                               ANH12540
      SUBROUTINE NORMFX(ELXM,EE,VALU,FV1,FV2,FXM)                       ANH12550
      IMPLICIT REAL*8 (A-H,O-Z)                                         ANH12560
      DIMENSION ELXM(N3N,N3N),EE(N3N,N3N)                               ANH12570
      DIMENSION VALU(N3N),FV1(N3N),FV2(N3N),FXM(NATRI)                  ANH12580
      COMMON/VIB101/NATOM,N3N,NATRI,ILIN,NVIB                           ANH12590
      COMMON/VIB103/PARA,WAVE,CONST,CYCL,CONV                           ANH12600
      COMMON/VIB108/IAXIS(3),NDEG(150),NDAB(150,5),IMAG(150)            ANH12610
      COMMON/VIB203/IOFF(150),IPRNT
      COMMON/VIB204/SQM(150),ROOT(150),FREQ(150)                        ANH12620
      DATA A00 / 0.0D+00 /                                              ANH12630
    1 FORMAT(//,2X,28H EIGENVALUES AND EIGENVECTOR)                     ANH12640
    2 FORMAT(/,2X,28H EIGENVALUES AND FREQUENCIES/                      ANH12650
     1 13X,8H IN A.U.,9X,8H IN CM-1/)                                   ANH12660
    3 FORMAT(2X,I5,F15.6,2X,F15.6)                                      ANH12670
    4 FORMAT(//,2X,' LXM MATRIX'/)                                      ANH12680
    5 FORMAT(//,2X,' LXM * LXM(T) = E '/)                               ANH12690
C                                                                       ANH12700
C   THE CALCULATION OF EIGENVALUE AND EIGENVECTOR                       ANH12710
      CALL RSP(N3N,N3N,NATRI,FXM,VALU,1,ELXM,FV1,FV2)                   ANH12720
C     WRITE(6,1)                                                        ANH12730
      WRITE(6,2)                                                        ANH12740
      DO 101 I=1,N3N                                                    ANH12750
      II=N3N-I+1                                                        ANH12760
      ROOT(I)=VALU(II)                                                  ANH12770
      IMAG(I)=0                                                         ANH12780
      IF(VALU(II).LT.A00) IMAG(I)=1                                     ANH12790
      IF(VALU(II).LT.A00) VALU(II)=-VALU(II)                            ANH12800
      AVAL=VALU(II)                                                     ANH12810
      FREQ(I)=WAVE*DSQRT(AVAL)                                          ANH12820
      WRITE(6,3) I,ROOT(I),FREQ(I)                                      ANH12830
  101 CONTINUE
C                                                                       ANH12840
C   EIGENVECTOR CORRESPONDS TO LXM VECTOR                               ANH12850
      NN=N3N/2                                                          ANH12860
      DO 102 I=1,N3N                                                    ANH12870
      DO 102 J=1,NN                                                     ANH12880
      JJ=N3N-J+1                                                        ANH12890
      EJ=ELXM(I,J)                                                      ANH12900
      EJJ=ELXM(I,JJ)                                                    ANH12910
      ELXM(I,J)=EJJ                                                     ANH12920
      ELXM(I,JJ)=EJ                                                     ANH12930
  102 CONTINUE                                                          ANH12940
      WRITE(6,4)                                                        ANH12950
      CALL FRQOUT(ELXM,FREQ,N3N,N3N,N3N,N3N,6)                          ANH12960
C                                                                       ANH12970
C   THE EXAMINATION OF ABOVE CALCULATION                                ANH12980
      IF(IPRNT.LE.2) GO TO 201                                          ANH12990
      CALL MTXMPY(ELXM,N3N,ELXM,N3N,EE,N3N,EE,N3N,N3N,3)                ANH13000
      WRITE(6,5)                                                        ANH13010
      CALL MATOUT(EE,N3N,N3N,N3N,N3N,6)                                 ANH13020
C                                                                       ANH13030
  201 CONTINUE                                                          ANH13040
      RETURN                                                            ANH13050
      END                                                               ANH13060
      SUBROUTINE CORIOL(EPXM,ZETA,EE,FF,AM)                             ANH13070
C   THE CORIOLIS COUPLING CONSTANTS                                     ANH13080
      IMPLICIT REAL*8 (A-H,O-Z)                                         ANH13090
      CHARACTER*3 AXIS(3)                                               ANH13100
      DIMENSION EPXM(N3N,N3N),ZETA(N3N,N3N,3)                           ANH13110
      DIMENSION EE(N3N,N3N),FF(N3N,N3N),AM(N3N,N3N,3)                   ANH13120
      COMMON/VIB101/NATOM,N3N,NATRI,ILIN,NVIB                           ANH13130
      COMMON/VIB203/IOFF(150),IPRNT                                     ANH13140
      DATA A00,ONE / 0.0D+00 , 1.0D+00 /                                ANH13150
      DATA AXIS / '(X)' , '(Y)' , '(Z)' /                               ANH13160
    1 FORMAT(//,2X,' *************************************'/            ANH13170
     1          2X,' ***THE CORIOLIS COUPLING CONSTANTS***'/            ANH13180
     2          2X,' *************************************')            ANH13190
    2 FORMAT(//,2X,' LXM'' MATRIX'/)                                    ANH13200
    3 FORMAT(//,2X,' LXM'' * LXM''(T) = E'/)                            ANH13210
    4 FORMAT(//,2X,' MALPH MATRIX',A3/)                                 ANH13220
    5 FORMAT(//,2X,' ZETA MATRIX ',A3/)                                 ANH13230
C                                                                       ANH13240
      WRITE(6,1)                                                        ANH13250
      IF(IPRNT.LE.2) GO TO 201                                          ANH13260
      WRITE(6,2)                                                        ANH13270
      CALL MATOUT(EPXM,N3N,N3N,N3N,N3N,6)                               ANH13280
      CALL MTXMPY(EPXM,N3N,EPXM,N3N,EE,N3N,EE,N3N,N3N,3)                ANH13290
      WRITE(6,3)                                                        ANH13300
      CALL MATOUT(EE,N3N,N3N,N3N,N3N,6)                                 ANH13310
C                                                                       ANH13320
  201 CONTINUE                                                          ANH13330
C****************************************************                   ANH13340
C***FORMATION OF ZETA MATRICES                      *                   ANH13350
C***AN ALTERNATE EXPRESSION (IN A SIMPLER MANNER)   *                   ANH13360
C***GIVES THE SAME RESULTS                          *                   ANH13370
C****************************************************                   ANH13380
C*    CALL ZERO(AM,N3N*N3N*3)                                           ANH13390
C*    DO 102 I=1,NATOM                                                  ANH13400
C*    IC=(I-1)*3                                                        ANH13410
C*    IROW=IC+2                                                         ANH13420
C*    ICOL=IC+3                                                         ANH13430
C*    AM(IROW,ICOL,1)=ONE                                               ANH13440
C*    AM(ICOL,IROW,1)=-ONE                                              ANH13450
C*    IROW=IC+3                                                         ANH13460
C*    ICOL=IC+1                                                         ANH13470
C*    AM(IROW,ICOL,2)=ONE                                               ANH13480
C*    AM(ICOL,IROW,2)=-ONE                                              ANH13490
C*    IROW=IC+1                                                         ANH13500
C*    ICOL=IC+2                                                         ANH13510
C*    AM(IROW,ICOL,3)=ONE                                               ANH13520
C*    AM(ICOL,IROW,3)=-ONE                                              ANH13530
C*102 CONTINUE                                                          ANH13540
C*    DO 105 KABC=1,3                                                   ANH13550
C*    IF(IPRNT.LE.2) GO TO 202                                          ANH13560
C*    WRITE(6,4) AXIS(KABC)                                             ANH13570
C*    CALL MATOUT(AM(1,1,KABC),N3N,N3N,N3N,N3N,6)                       ANH13580
C*202 CONTINUE                                                          ANH13590
C*    CALL MTXMPY(EPXM,N3N,AM(1,1,KABC),N3N,FF,N3N,EE,N3N,N3N,5)        ANH13600
C*    WRITE(6,5) AXIS(KABC)                                             ANH13610
C*    CALL MATOUT(FF,N3N,N3N,N3N,N3N,6)                                 ANH13620
C*    DO 103 I=1,N3N                                                    ANH13630
C*    DO 103 J=1,N3N                                                    ANH13640
C*    ZETA(I,J,KABC)=FF(I,J)                                            ANH13650
C*103 CONTINUE
C*105 CONTINUE                                                          ANH13660
C                                                                       ANH13670
C   FORM ZETA MATRICES USING AN ALTERNATE EXPRESSION                    ANH13680
      CALL ZERO(ZETA,N3N*N3N*3)                                         ANH13690
      DO 107 IVIB=1,NVIB                                                ANH13700
      DO 107 JVIB=1,NVIB                                                ANH13710
      VALUX=A00                                                         ANH13720
      VALUY=A00                                                         ANH13730
      VALUZ=A00                                                         ANH13740
      DO 106 K=1,NATOM                                                  ANH13750
      IX=(K-1)*3+1                                                      ANH13760
      IY=IX+1                                                           ANH13770
      IZ=IX+2                                                           ANH13780
      VALUX=VALUX+EPXM(IY,IVIB)*EPXM(IZ,JVIB)                           ANH13790
     1           -EPXM(IZ,IVIB)*EPXM(IY,JVIB)                           ANH13800
      VALUY=VALUY+EPXM(IZ,IVIB)*EPXM(IX,JVIB)                           ANH13810
     1           -EPXM(IX,IVIB)*EPXM(IZ,JVIB)                           ANH13820
      VALUZ=VALUZ+EPXM(IX,IVIB)*EPXM(IY,JVIB)                           ANH13830
     1           -EPXM(IY,IVIB)*EPXM(IX,JVIB)                           ANH13840
  106 CONTINUE                                                          ANH13850
      ZETA(IVIB,JVIB,1)=VALUX                                           ANH13860
      ZETA(IVIB,JVIB,2)=VALUY                                           ANH13870
      ZETA(IVIB,JVIB,3)=VALUZ                                           ANH13880
  107 CONTINUE                                                          ANH13890
      DO 108 IXYZ=1,3                                                   ANH13900
      WRITE(6,5) AXIS(IXYZ)                                             ANH13910
      CALL MATOUT(ZETA(1,1,IXYZ),N3N,N3N,NVIB,NVIB,6)                   ANH13920
  108 CONTINUE                                                          ANH13930
C                                                                       ANH13940
      RETURN                                                            ANH13950
      END                                                               ANH13960
      SUBROUTINE ABCMAT(EPXM,AXY,BXY,CXY,EE)                            ANH13970
C   FORM A (XY) MATRICES                                                ANH13980
      IMPLICIT REAL*8 (A-H,O-Z)                                         ANH13990
      DIMENSION EPXM(N3N,N3N),AXY(3,3,N3N),BXY(3,3,N3N),CXY(3,3,N3N)    ANH14000
      DIMENSION EE(N3N,N3N)                                             ANH14010
      COMMON/VIB101/NATOM,N3N,NATRI,ILIN,NVIB                           ANH14020
      COMMON/VIB106/ROTAA(3),ROTGC(3),ROTCM(3),ROTMH(3)                 ANH14030
      COMMON/VIB108/IAXIS(3),NDEG(150),NDAB(150,5),IMAG(150)            ANH14040
      COMMON/VIB201/CHG(50),X(50),Y(50),Z(50),W(50)                     ANH14050
      COMMON/VIB203/IOFF(150),IPRNT                                     ANH14060
      COMMON/VIB204/SQM(150),ROOT(150),FREQ(150)                        ANH14070
      DATA PH,AVN / 6.626176D+00 , 6.022045D+00 /                       ANH14080
      DATA PI,CL / 3.1415926536D+00 , 2.99792458D+00 /                  ANH14090
      DATA A00,HALF,ONE,TWO / 0.0D+00 , 0.5D+00 , 1.0D+00 , 2.0D+00 /   ANH14100
      DATA WLIMIT / 1.0D-2 /                                            ANH14110
    1 FORMAT(//,2X,' ********************************************'/     ANH14120
     1          2X,' ***THE A (XY) , B (XY) AND C(XY) MATRICES***'/     ANH14130
     2          2X,' ********************************************')     ANH14140
    2 FORMAT(//,2X,' LXM'' MATRIX'/)                                    ANH14150
    3 FORMAT(//,2X,' LXM'' * LXM''(T) = E'/)                            ANH14160
    4 FORMAT(//,2X,' CARTESIAN COORDINATES W.R.T. CENTER OF MASS (IN A)'ANH14170
     1/5X,4H NO.,11X,2H X,18X,2H Y,18X,2H Z,18X,2H W/)                  ANH14180
    5 FORMAT(2X,I7,5X,4F15.7)                                           ANH14190
    6 FORMAT(//,2X,' AXY MATRIX, IVIB = ',I5/)                          ANH14200
    7 FORMAT(//,2X,' BXY MATRIX, IVIB = ',I5/)                          ANH14210
    8 FORMAT(//,2X,' CXY MATRIX, IVIB = ',I5/)                          ANH14220
    9 FORMAT(//,2X,' AXY MATRIX, IVIB = ',I5,18X,' BXY MATRIX, IVIB = ',ANH14230
     1       I5,18X,' CXY MATRIX, IVIB = ',I5/)                         ANH14240
C                                                                       ANH14250
      IF(IPRNT.GE.1)                                                    ANH14260
     *WRITE(6,1)                                                        ANH14270
      IF(IPRNT.LE.2) GO TO 201                                          ANH14280
      WRITE(6,2)                                                        ANH14290
      CALL MATOUT(EPXM,N3N,N3N,N3N,N3N,6)                               ANH14300
      CALL MTXMPY(EPXM,N3N,EPXM,N3N,EE,N3N,EE,N3N,N3N,3)                ANH14310
      WRITE(6,3)                                                        ANH14320
      CALL MATOUT(EE,N3N,N3N,N3N,N3N,6)                                 ANH14330
      WRITE(6,4)                                                        ANH14340
      DO 101 I=1,NATOM                                                  ANH14350
      WRITE(6,5) I,X(I),Y(I),Z(I),W(I)                                  ANH14360
  101 CONTINUE                                                          ANH14370
C                                                                       ANH14380
C   FORM AXY MATRICES                                                   ANH14390
C   ELEMENTS OF AXY MATRIX ARE IN AMU(1/2).A                            ANH14400
  201 CONTINUE                                                          ANH14410
      CALL ZERO(AXY,3*3*N3N)                                            ANH14420
      DO 105 IVIB=1,NVIB                                                ANH14430
      CALL ZERO(EE,N3N*N3N)                                             ANH14440
      DO 103 I=1,NATOM                                                  ANH14450
      IX=(I-1)*3+1                                                      ANH14460
      IY=IX+1                                                           ANH14470
      IZ=IX+2                                                           ANH14480
      EE(1,1)=EE(1,1)+SQM(I)*(Y(I)*EPXM(IY,IVIB)+Z(I)*EPXM(IZ,IVIB))*TWOANH14490
      EE(2,2)=EE(2,2)+SQM(I)*(Z(I)*EPXM(IZ,IVIB)+X(I)*EPXM(IX,IVIB))*TWOANH14500
      EE(3,3)=EE(3,3)+SQM(I)*(X(I)*EPXM(IX,IVIB)+Y(I)*EPXM(IY,IVIB))*TWOANH14510
      EE(1,2)=EE(1,2)-SQM(I)*(X(I)*EPXM(IY,IVIB)+Y(I)*EPXM(IX,IVIB))    ANH14520
      EE(1,3)=EE(1,3)-SQM(I)*(Z(I)*EPXM(IX,IVIB)+X(I)*EPXM(IZ,IVIB))    ANH14530
      EE(2,3)=EE(2,3)-SQM(I)*(Y(I)*EPXM(IZ,IVIB)+Z(I)*EPXM(IY,IVIB))    ANH14540
  103 CONTINUE                                                          ANH14550
      EE(2,1)=EE(1,2)                                                   ANH14560
      EE(3,1)=EE(1,3)                                                   ANH14570
      EE(3,2)=EE(2,3)                                                   ANH14580
      DO 104 I=1,3                                                      ANH14590
      DO 104 J=1,3                                                      ANH14600
      AXY(I,J,IVIB)=EE(I,J)                                             ANH14610
  104 CONTINUE                                                          ANH14620
      IF(IPRNT.LE.2) GO TO 105                                          ANH14630
      IF(DABS(FREQ(IVIB)).LE.WLIMIT) GO TO 105                          ANH14640
      WRITE(6,6) IVIB                                                   ANH14650
      CALL MATOUT(AXY(1,1,IVIB),3,3,3,3,6)                              ANH14660
  105 CONTINUE                                                          ANH14670
C                                                                       ANH14680
C   BCONST IS IN CM-1                                                   ANH14690
      BCONST=1.0D+03*((PH*AVN/CL)**(1.5D+00))/(16.0D+00*(PI**3))        ANH14700
C                                                                       ANH14710
C   FORM BXY AND CXY MATRICES                                           ANH14720
C   ELEMENTS OF BXY MATRIX ARE IN CM-1                                  ANH14730
C   ELEMENTS OF CXY MATRIX DO NOT HAVE UNITS                            ANH14740
      CALL ZERO(BXY,3*3*N3N)                                            ANH14750
      CALL ZERO(CXY,3*3*N3N)                                            ANH14760
      DO 108 IVIB=1,NVIB                                                ANH14770
      FRQI=FREQ(IVIB)                                                   ANH14780
      IF(DABS(FRQI).LE.WLIMIT) GO TO 108                                ANH14790
      FRQIH=DSQRT(FRQI)                                                 ANH14800
      DO 107 IA=1,3                                                     ANH14810
      DO 106 IB=1,3                                                     ANH14820
      IAA=IAXIS(IA)                                                     ANH14830
      IBB=IAXIS(IB)                                                     ANH14840
      DENAB=ROTAA(IAA)*ROTAA(IBB)                                       ANH14850
      IF(DABS(DENAB).LE.WLIMIT) GO TO 106                               ANH14860
      DENAB=-DENAB*FRQIH                                                ANH14870
      BXY(IA,IB,IVIB)=BCONST*AXY(IA,IB,IVIB)/DENAB                      ANH14880
      CXY(IA,IB,IVIB)=-BXY(IA,IB,IVIB)/FRQI                             ANH14890
  106 CONTINUE                                                          ANH14900
  107 CONTINUE                                                          ANH14910
      IF(IPRNT.LE.2) GO TO 108                                          ANH14920
      WRITE(6,7) IVIB                                                   ANH14930
      CALL MATOUT(BXY(1,1,IVIB),3,3,3,3,6)                              ANH14940
      WRITE(6,8) IVIB                                                   ANH14950
      CALL MATOUT(CXY(1,1,IVIB),3,3,3,3,6)                              ANH14960
  108 CONTINUE                                                          ANH14970
C                                                                       ANH14980
      IF(IPRNT.LE.1) GO TO 202                                          ANH14990
      DO 110 IVIB=1,NVIB                                                ANH15000
      WRITE(6,9) IVIB,IVIB,IVIB                                         ANH15010
      CALL MX3OUT(AXY(1,1,IVIB),3,BXY(1,1,IVIB),3,CXY(1,1,IVIB),3,3,6)  ANH15020
  110 CONTINUE                                                          ANH15030
C                                                                       ANH15040
  202 CONTINUE                                                          ANH15050
      RETURN                                                            ANH15060
      END                                                               ANH15070
      SUBROUTINE TAUMAT(AXY,BXY)                                        ANH15080
C   THE CALCULATION OF CENTRIFUGAL DISTORTION CONSTANTS                 ANH15090
      IMPLICIT REAL*8 (A-H,O-Z)                                         ANH15100
      DIMENSION AXY(3,3,N3N),BXY(3,3,N3N)                               ANH15110
      DIMENSION TAABB(3,3),TABAB(3,3)                                   ANH15120
      COMMON/VIB101/NATOM,N3N,NATRI,ILIN,NVIB                           ANH15130
      COMMON/VIB102/ITHREE,IFOUR,N3TOT,N4TOT                            ANH15140
      COMMON/VIB106/ROTAA(3),ROTGC(3),ROTCM(3),ROTMH(3)                 ANH15150
      COMMON/VIB108/IAXIS(3),NDEG(150),NDAB(150,5),IMAG(150)            ANH15160
      COMMON/VIB203/IOFF(150),IPRNT                                     ANH15170
      COMMON/VIB204/SQM(150),ROOT(150),FREQ(150)                        ANH15180
      COMMON/VIB205/IFREQ,NFRQ(150)                                     ANH15190
      COMMON/VIB207/TABCD(3,3,3,3),PABC(3,3,3)                          ANH15200
      DATA WLIMIT / 1.0D-2 /                                            ANH15210
      DATA PH,AVN / 6.626176D+00 , 6.022045D+00 /                       ANH15220
      DATA PI,CL / 3.1415926536D+00 , 2.99792458D+00 /                  ANH15230
      DATA A00,HALF,ONE,TWO / 0.0D+00 , 0.5D+00 , 1.0D+00 , 2.0D+00 /   ANH15240
    1 FORMAT(//,2X,' *****************************'/                    ANH15250
     1          2X,' ***THE TAU (ABCD) MATRICES***'/                    ANH15260
     2          2X,' *****************************')                    ANH15270
    2 FORMAT(//,2X,' AXY MATRIX, IVIB = ',I5/)                          ANH15280
    3 FORMAT(//,2X,' BXY MATRIX, IVIB = ',I5/)                          ANH15290
    4 FORMAT(//,2X,' TAABB MATRIX FROM AXY MATRICES (IN CM-1)',         ANH15300
     1          2X,' TABAB MATRIX FROM AXY MATRICES (IN CM-1)'/)        ANH15310
    5 FORMAT(//,2X,' TAABB MATRIX FROM BXY MATRICES (IN CM-1)',         ANH15320
     1          2X,' TABAB MATRIX FROM BXY MATRICES (IN CM-1)'/)        ANH15330
C                                                                       ANH15340
      WRITE(6,1)                                                        ANH15350
C   TCONST IS IN CM-1                                                   ANH15360
      TCONST=1.0D+06*((PH*AVN/CL)**3)/(128.0D+00*(PI**6))               ANH15370
C                                                                       ANH15380
C   ELEMENTS OF AXY MATRIX ARE IN AMU(1/2).A                            ANH15390
      IF(IPRNT.LE.2) GO TO 201                                          ANH15400
      DO 101 IVIB=1,NVIB                                                ANH15410
      WRITE(6,2) IVIB                                                   ANH15420
      CALL MATOUT(AXY(1,1,IVIB),3,3,3,3,6)                              ANH15430
  101 CONTINUE                                                          ANH15440
C                                                                       ANH15450
C   ELEMENTS OF BXY MATRIX ARE IN CM-1                                  ANH15460
      DO 102 IVIB=1,NVIB                                                ANH15470
      WRITE(6,3) IVIB                                                   ANH15480
      CALL MATOUT(BXY(1,1,IVIB),3,3,3,3,6)                              ANH15490
  102 CONTINUE                                                          ANH15500
C                                                                       ANH15510
C   FORM TABCD MATRIX USING AXY MATRICES                                ANH15520
C   ELEMENTS OF TABCD MATRIX ARE IN CM-1                                ANH15530
  201 CONTINUE                                                          ANH15540
      CALL ZERO(TABCD,81)                                               ANH15550
      DO 105 IVIB=1,NVIB                                                ANH15560
      FRQI=FREQ(IVIB)                                                   ANH15570
      IF(DABS(FRQI).LE.WLIMIT) GO TO 105                                ANH15580
      FRQI2=FRQI*FRQI                                                   ANH15590
      DO 103 IA=1,3                                                     ANH15600
      DO 103 IB=1,3                                                     ANH15610
      DO 103 IC=1,3                                                     ANH15620
      DO 103 ID=1,3                                                     ANH15630
      IAA=IAXIS(IA)                                                     ANH15640
      IBB=IAXIS(IB)                                                     ANH15650
      ICC=IAXIS(IC)                                                     ANH15660
      IDD=IAXIS(ID)                                                     ANH15670
      DABCD=ROTAA(IAA)*ROTAA(IBB)*ROTAA(ICC)*ROTAA(IDD)                 ANH15680
      IF(DABS(DABCD).LE.WLIMIT) GO TO 103                               ANH15690
      FABCD=TCONST*(AXY(IA,IB,IVIB)*AXY(IC,ID,IVIB))/DABCD              ANH15700
      TABCD(IA,IB,IC,ID)=TABCD(IA,IB,IC,ID)-FABCD/FRQI2                 ANH15710
  103 CONTINUE                                                          ANH15720
  105 CONTINUE                                                          ANH15730
      DO 106 IA=1,3                                                     ANH15740
      DO 106 IB=1,3                                                     ANH15750
      TAABB(IA,IB)=TABCD(IA,IA,IB,IB)                                   ANH15760
      TABAB(IA,IB)=TABCD(IA,IB,IA,IB)                                   ANH15770
  106 CONTINUE                                                          ANH15780
      WRITE(6,4)                                                        ANH15790
      CALL MX2OUT(TAABB,3,TABAB,3,3,6)                                  ANH15800
C                                                                       ANH15810
C   FORM TABCD USING BXY MATRICES                                       ANH15820
C   ELEMENTS OF TABCD MATRIX ARE IN CM-1                                ANH15830
      CALL ZERO(TABCD,81)                                               ANH15840
      DO 108 IVIB=1,NVIB                                                ANH15850
      FRQI=FREQ(IVIB)                                                   ANH15860
      IF(DABS(FRQI).LE.WLIMIT) GO TO 108                                ANH15870
      DO 107 IA=1,3                                                     ANH15880
      DO 107 IB=1,3                                                     ANH15890
      DO 107 IC=1,3                                                     ANH15900
      DO 107 ID=1,3                                                     ANH15910
      FABCD=-TWO*BXY(IA,IB,IVIB)*BXY(IC,ID,IVIB)/FRQI                   ANH15920
      TABCD(IA,IB,IC,ID)=TABCD(IA,IB,IC,ID)+FABCD                       ANH15930
  107 CONTINUE                                                          ANH15940
  108 CONTINUE                                                          ANH15950
      DO 109 IA=1,3                                                     ANH15960
      DO 109 IB=1,3                                                     ANH15970
      TAABB(IA,IB)=TABCD(IA,IA,IB,IB)                                   ANH15980
      TABAB(IA,IB)=TABCD(IA,IB,IA,IB)                                   ANH15990
  109 CONTINUE                                                          ANH16000
      WRITE(6,5)                                                        ANH16010
      CALL MX2OUT(TAABB,3,TABAB,3,3,6)                                  ANH16020
C                                                                       ANH16030
      RETURN                                                            ANH16040
      END                                                               ANH16050
      SUBROUTINE BCXYZ(ZETA,BXY,BXYZ,CXYZ)                              ANH16060
C   THE CALCULATION OF CENTRIFUGAL DISTORTION CONSTANTS                 ANH16070
      IMPLICIT REAL*8 (A-H,O-Z)                                         ANH16080
      DIMENSION ZETA(N3N,N3N,3),BXY(3,3,N3N)                            ANH16090
      DIMENSION BXYZ(3,3,3,N3N),CXYZ(3,3,3,N3N)                         ANH16100
      COMMON/VIB101/NATOM,N3N,NATRI,ILIN,NVIB                           ANH16110
      COMMON/VIB102/ITHREE,IFOUR,N3TOT,N4TOT                            ANH16120
      COMMON/VIB106/ROTAA(3),ROTGC(3),ROTCM(3),ROTMH(3)                 ANH16130
      COMMON/VIB108/IAXIS(3),NDEG(150),NDAB(150,5),IMAG(150)            ANH16140
      COMMON/VIB203/IOFF(150),IPRNT                                     ANH16150
      COMMON/VIB204/SQM(150),ROOT(150),FREQ(150)                        ANH16160
      COMMON/VIB205/IFREQ,NFRQ(150)                                     ANH16170
      COMMON/VIB207/TABCD(3,3,3,3),PABC(3,3,3)                          ANH16180
      DATA WLIMIT / 1.0D-2 /                                            ANH16190
      DATA PH,AVN / 6.626176D+00 , 6.022045D+00 /                       ANH16200
      DATA PI,CL / 3.1415926536D+00 , 2.99792458D+00 /                  ANH16210
      DATA A00,HALF,ONE,TWO / 0.0D+00 , 0.5D+00 , 1.0D+00 , 2.0D+00 /   ANH16220
    1 FORMAT(//,2X,' **************************************'/           ANH16230
     1          2X,' ***THE B (XYZ) AND C (XYZ) MATRICES***'/           ANH16240
     2          2X,' **************************************')           ANH16250
    2 FORMAT(//,2X,' BXY MATRIX, IVIB = ',I5/)                          ANH16260
    3 FORMAT(//,2X,' BXYZ (X) MATRIX (IN D-06 CM-1), IVIB = ',I2,       ANH16270
     1          1X,' BXYZ (Y) MATRIX (IN D-06 CM-1), IVIB = ',I2,       ANH16280
     2          1X,' BXYZ (Z) MATRIX (IN D-06 CM-1), IVIB = ',I2/)      ANH16290
    4 FORMAT(//,2X,' CXYZ (X) MATRIX (IN D-08 CM-1), IVIB = ',I2,       ANH16300
     1          1X,' CXYZ (Y) MATRIX (IN D-08 CM-1), IVIB = ',I2,       ANH16310
     2          1X,' CXYZ (Z) MATRIX (IN D-08 CM-1), IVIB = ',I2/)      ANH16320
C                                                                       ANH16330
      IF(IPRNT.GE.1)                                                    ANH16340
     *WRITE(6,1)                                                        ANH16350
      IXX=1                                                             ANH16360
      IYY=2                                                             ANH16370
      IZZ=3                                                             ANH16380
C                                                                       ANH16390
C   ELEMENTS OF BXY MATRIX ARE IN CM-1                                  ANH16400
C   ELEMENTS OF CXY MATRIX ARE IN NO UNITS                              ANH16410
      IF(IPRNT.LE.2) GO TO 201                                          ANH16420
      DO 101 IVIB=1,NVIB                                                ANH16430
      WRITE(6,2) IVIB                                                   ANH16440
      CALL MATOUT(BXY(1,1,IVIB),3,3,3,3,6)                              ANH16450
  101 CONTINUE                                                          ANH16460
C                                                                       ANH16470
C   FORM BXYZ MATRICES                                                  ANH16480
C   ELEMENTS OF BXYZ MATRIX ARE IN CM(-1/2)                             ANH16490
  201 CONTINUE                                                          ANH16500
      CALL ZERO(BXYZ,3*3*3*N3N)                                         ANH16510
      CALL ZERO(CXYZ,3*3*3*N3N)                                         ANH16520
      DO 110 IVIB=1,NVIB                                                ANH16530
      FRQI=FREQ(IVIB)                                                   ANH16540
      FRQI2=FRQI*FRQI                                                   ANH16550
      FRQIH=DSQRT(FRQI)                                                 ANH16560
C                                                                       ANH16570
C   B(IA,IA,IA,IVIB) AND CXYZ(IA,IA,IA,IVIB) TERMS                      ANH16580
      DO 103 IA=1,3                                                     ANH16590
      VALU=A00                                                          ANH16600
      DO 102 JVIB=1,NVIB                                                ANH16610
      FRQJ=FREQ(JVIB)                                                   ANH16620
      FRQJ2=FRQJ*FRQJ                                                   ANH16630
      FRQJH=DSQRT(FRQJ)                                                 ANH16640
      IAA=IAXIS(IA)                                                     ANH16650
      FACT=BXY(IA,IA,JVIB)*ROTCM(IAA)*ZETA(JVIB,IVIB,IA)                ANH16660
      VALU=VALU+(FACT*FRQJH)*(ONE/FRQJ2+TWO/FRQI2)                      ANH16670
  102 CONTINUE                                                          ANH16680
      BXYZ(IA,IA,IA,IVIB)=VALU                                          ANH16690
      CXYZ(IA,IA,IA,IVIB)=-VALU/FRQIH                                   ANH16700
  103 CONTINUE                                                          ANH16710
C                                                                       ANH16720
C   BXYZ(IA,IA,IB,IVIB) AND CXYZ(IA,IA,IB,IVIB) TERMS                   ANH16730
      DO 105 IA=1,3                                                     ANH16740
      DO 105 IB=1,3                                                     ANH16750
      IF(IA.EQ.IB) GO TO 105                                            ANH16760
      VALU=A00                                                          ANH16770
      DO 104 JVIB=1,NVIB                                                ANH16780
      FRQJ=FREQ(JVIB)                                                   ANH16790
      FRQJ2=FRQJ*FRQJ                                                   ANH16800
      FRQJH=DSQRT(FRQJ)                                                 ANH16810
      IAA=IAXIS(IA)                                                     ANH16820
      IBB=IAXIS(IB)                                                     ANH16830
      FACT=BXY(IA,IA,JVIB)*ROTCM(IBB)*ZETA(JVIB,IVIB,IB)                ANH16840
     1    +BXY(IA,IB,JVIB)*ROTCM(IAA)*ZETA(JVIB,IVIB,IA)*TWO            ANH16850
      VALU=VALU+(FACT*FRQJH)*(ONE/FRQJ2+TWO/FRQI2)                      ANH16860
  104 CONTINUE                                                          ANH16870
      BXYZ(IA,IA,IB,IVIB)=VALU                                          ANH16880
      BXYZ(IA,IB,IA,IVIB)=VALU                                          ANH16890
      BXYZ(IB,IA,IA,IVIB)=VALU                                          ANH16900
      CXYZ(IA,IA,IB,IVIB)=-VALU/FRQIH                                   ANH16910
      CXYZ(IA,IB,IA,IVIB)=-VALU/FRQIH                                   ANH16920
      CXYZ(IB,IA,IA,IVIB)=-VALU/FRQIH                                   ANH16930
  105 CONTINUE                                                          ANH16940
C                                                                       ANH16950
C   BXYZ(IA,IB,IC,IVIB) AND CXYZ(IA,IB,IC,IVIB) TERMS                   ANH16960
      VALU=A00                                                          ANH16970
      IAA=IAXIS(1)                                                      ANH16980
      IBB=IAXIS(2)                                                      ANH16990
      ICC=IAXIS(3)                                                      ANH17000
      DO 106 JVIB=1,NVIB                                                ANH17010
      FRQJ=FREQ(JVIB)                                                   ANH17020
      FRQJ2=FRQJ*FRQJ                                                   ANH17030
      FRQJH=DSQRT(FRQJ)                                                 ANH17040
      FACT=BXY(IYY,IZZ,JVIB)*ROTCM(IAA)*ZETA(JVIB,IVIB,IXX)             ANH17050
     1    +BXY(IZZ,IXX,JVIB)*ROTCM(IBB)*ZETA(JVIB,IVIB,IYY)             ANH17060
     2    +BXY(IXX,IYY,JVIB)*ROTCM(ICC)*ZETA(JVIB,IVIB,IZZ)             ANH17070
      VALU=VALU+(FACT*FRQJH)*(ONE/FRQJ2+TWO/FRQI2)                      ANH17080
  106 CONTINUE                                                          ANH17090
      BXYZ(IXX,IYY,IZZ,IVIB)=VALU                                       ANH17100
      BXYZ(IXX,IZZ,IYY,IVIB)=VALU                                       ANH17110
      BXYZ(IYY,IXX,IZZ,IVIB)=VALU                                       ANH17120
      BXYZ(IYY,IZZ,IXX,IVIB)=VALU                                       ANH17130
      BXYZ(IZZ,IXX,IYY,IVIB)=VALU                                       ANH17140
      BXYZ(IZZ,IYY,IXX,IVIB)=VALU                                       ANH17150
      CXYZ(IXX,IYY,IZZ,IVIB)=-VALU/FRQIH                                ANH17160
      CXYZ(IXX,IZZ,IYY,IVIB)=-VALU/FRQIH                                ANH17170
      CXYZ(IYY,IXX,IZZ,IVIB)=-VALU/FRQIH                                ANH17180
      CXYZ(IYY,IZZ,IXX,IVIB)=-VALU/FRQIH                                ANH17190
      CXYZ(IZZ,IXX,IYY,IVIB)=-VALU/FRQIH                                ANH17200
      CXYZ(IZZ,IYY,IXX,IVIB)=-VALU/FRQIH                                ANH17210
C                                                                       ANH17220
  110 CONTINUE                                                          ANH17230
      IF(IPRNT.LE.1) GO TO 202                                          ANH17240
      DO 115 IVIB=1,NVIB                                                ANH17250
      DO 113 IA=1,3                                                     ANH17260
      DO 113 IB=1,3                                                     ANH17270
      DO 113 IC=1,3                                                     ANH17280
      BXYZ(IA,IB,IC,IVIB)=BXYZ(IA,IB,IC,IVIB)*1.0D+06                   ANH17290
      CXYZ(IA,IB,IC,IVIB)=CXYZ(IA,IB,IC,IVIB)*1.0D+08                   ANH17300
  113 CONTINUE                                                          ANH17310
      WRITE(6,3) IVIB,IVIB,IVIB                                         ANH17320
      CALL MX3OUT(BXYZ(1,1,1,IVIB),3,BXYZ(1,1,2,IVIB),3,                ANH17330
     1            BXYZ(1,1,3,IVIB),3,3,6)                               ANH17340
      WRITE(6,4) IVIB,IVIB,IVIB                                         ANH17350
      CALL MX3OUT(CXYZ(1,1,1,IVIB),3,CXYZ(1,1,2,IVIB),3,                ANH17360
     1            CXYZ(1,1,3,IVIB),3,3,6)                               ANH17370
      DO 114 IA=1,3                                                     ANH17380
      DO 114 IB=1,3                                                     ANH17390
      DO 114 IC=1,3                                                     ANH17400
      BXYZ(IA,IB,IC,IVIB)=BXYZ(IA,IB,IC,IVIB)/1.0D+06                   ANH17410
      CXYZ(IA,IB,IC,IVIB)=CXYZ(IA,IB,IC,IVIB)/1.0D+08                   ANH17420
  114 CONTINUE                                                          ANH17430
  115 CONTINUE                                                          ANH17440
C                                                                       ANH17450
  202 CONTINUE                                                          ANH17460
      RETURN                                                            ANH17470
      END                                                               ANH17480
      SUBROUTINE FORCE3(ELX,F33,F3X,F3Q,EE,DD,GG)                       ANH18400
C   THE THIRD ORDER ANHARMONICITY ANALYSIS                              ANH18410
      IMPLICIT REAL*8 (A-H,O-Z)                                         ANH18420
      DIMENSION ELX(N3N,N3N)                                            ANH18430
      DIMENSION F33(N3TOT),F3X(N3N,N3N,N3N),F3Q(N3N,N3N,N3N)            ANH18440
      DIMENSION EE(N3N,N3N),DD(N3N,N3N),GG(N3N,N3N)                     ANH18450
      COMMON/VIB101/NATOM,N3N,NATRI,ILIN,NVIB                           ANH18460
      COMMON/VIB102/ITHREE,IFOUR,N3TOT,N4TOT                            ANH18470
      COMMON/VIB105/FACT3,FACT4                                         ANH18480
      COMMON/VIB203/IOFF(150),IPRNT                                     ANH18490
      COMMON/VIB204/SQM(150),ROOT(150),FREQ(150)                        ANH18500
      COMMON/VIB205/IFREQ,NFRQ(150)                                     ANH18510
      DATA A00 / 0.0D+00 /                                              ANH18520
      DATA ITAP20 / 20 /                                                ANH18530
      DATA A0,HE / 0.529177249D+00 , 4.359813653D+00 /    
      DATA WLIMIT / 3.0D+00 /                                           ANH18550
    1 FORMAT(//,2X,' THE THIRD ORDER ANHARMONICITY ANALYSIS'/)          ANH18560
    2 FORMAT(2I5)                                                       ANH18570
    3 FORMAT(3F20.10)                                                   ANH18580
    4 FORMAT(/,2X,' *************************************************'/ ANH18590
     1         2X,' ***THIRD DERIVATIVES W.R.T. NORMAL COORDINATES***'/ ANH18600
     2         2X,' *************************************************') ANH18610
    5 FORMAT(//,2X,' F3X MATRIX (IN HARTREE/BOHR+3), KABC = ',I5/)      ANH18620
    6 FORMAT(//,2X,' F3X MATRIX (IN MDYN/A+2), KABC = ',I5/)            ANH18630
    7 FORMAT(//,2X,' LX(T) * F3X * LX MATRIX, KABC = ',I5/)             ANH18640
    8 FORMAT(//,2X,' F3(QI*QJ*QK) MATRIX (IN MDYN/(A+2.AMU-3/2)),  QVIB ANH18650
     1=',I3,2X,', IVIB = ',I3/)                                         ANH18660
    9 FORMAT(//,2X,' F3(QI*QJ*QK) MATRIX (IN CM-1),  QVIB = ',I3,2X,    ANH18670
     1      ', IVIB = ',I3/)                                            ANH18680
C                                                                       ANH18690
      IF(IPRNT.NE.0)                                                    ANH18700
     *WRITE(6,1)                                                        ANH18710
C                                                                       ANH18720
C   FIND CONVERSION FACTORS                                             ANH18730
C   FUNIT3 IS IN MDYN/A+2                                               ANH18740
      FUNIT3=HE/(A0*A0*A0)                                              ANH18750
C                                                                       ANH18760
C   READ IN THIRD DERIVATIVES                                           ANH18770
      REWIND ITAP20                                                     ANH18780
      READ(ITAP20,2) NNATOM,NF3                                         ANH18790
      READ(ITAP20,3) (F33(I),I=1,N3TOT)                                 ANH18800
      REWIND ITAP20                                                     ANH18810
C                                                                       ANH18820
C   TRANSFORM THIRD DERIVATIVES                                         ANH18830
C   FROM CARTESIAN TO NORMAL COORDINATE SYSTEM                          ANH18840
      WRITE(6,4)                                                        ANH18850
      LABC=0                                                            ANH18860
      DO 101 IABC=1,N3N                                                 ANH18870
      DO 101 JABC=1,IABC                                                ANH18880
      DO 101 KABC=1,JABC                                                ANH18890
      LABC=LABC+1                                                       ANH18900
      VAL=F33(LABC)                                                     ANH18910
      F3X(IABC,JABC,KABC)=VAL                                           ANH18920
      F3X(IABC,KABC,JABC)=VAL                                           ANH18930
      F3X(JABC,IABC,KABC)=VAL                                           ANH18940
      F3X(JABC,KABC,IABC)=VAL                                           ANH18950
      F3X(KABC,IABC,JABC)=VAL                                           ANH18960
      F3X(KABC,JABC,IABC)=VAL                                           ANH18970
  101 CONTINUE                                                          ANH18980
C                                                                       ANH18990
      DO 105 KABC=1,N3N                                                 ANH19000
      DO 102 I=1,N3N                                                    ANH19010
      DO 102 J=1,N3N                                                    ANH19020
      DD(I,J)=F3X(I,J,KABC)                                             ANH19030
  102 CONTINUE
      IF(IPRNT.LE.2) GO TO 201                                          ANH19040
      WRITE(6,5) KABC                                                   ANH19050
      CALL MATOUT(DD,N3N,N3N,N3N,N3N,6)                                 ANH19060
  201 CONTINUE                                                          ANH19070
      DO 103 I=1,N3N                                                    ANH19080
      DO 103 J=1,N3N                                                    ANH19090
      DD(I,J)=DD(I,J)*FUNIT3                                            ANH19100
  103 CONTINUE
      IF(IPRNT.LE.2) GO TO 202                                          ANH19110
      WRITE(6,6) KABC                                                   ANH19120
      CALL MATOUT(DD,N3N,N3N,N3N,N3N,6)                                 ANH19130
C                                                                       ANH19140
  202 CONTINUE                                                          ANH19150
      CALL MTXMPY(ELX,N3N,DD,N3N,EE,N3N,GG,N3N,N3N,5)                   ANH19160
      IF(IPRNT.LE.2) GO TO 203                                          ANH19170
      WRITE(6,7) KABC                                                   ANH19180
      CALL MATOUT(EE,N3N,N3N,N3N,N3N,6)                                 ANH19190
  203 CONTINUE                                                          ANH19200
      DO 104 I=1,N3N                                                    ANH19210
      DO 104 J=1,N3N                                                    ANH19220
      F3X(I,J,KABC)=EE(I,J)                                             ANH19230
  104 CONTINUE
  105 CONTINUE                                                          ANH19240
      DO 107 I=1,N3N                                                    ANH19250
      DO 107 J=1,N3N                                                    ANH19260
      DO 107 K=1,N3N                                                    ANH19270
      VAL=A00                                                           ANH19280
      DO 106 L=1,N3N                                                    ANH19290
      VAL=VAL+F3X(I,J,L)*ELX(L,K)                                       ANH19300
  106 CONTINUE
      F3Q(I,J,K)=VAL                                                    ANH19310
  107 CONTINUE                                                          ANH19320
C                                                                       ANH19330
C   SET NUMBER OF ANHARMONICITY SETS                                    ANH19340
C   REGULARLY NVIB=3N-6(OR 3N-5)                                        ANH19350
C   IF IFREQ IS NOT EQUAL TO 0 THEN NVIB=IFREQ                          ANH19360
      IF(IFREQ.NE.0) NVIB=IABS(IFREQ)                                   ANH19370
C                                                                       ANH19380
      IF(IPRNT.LE.1) GO TO 204                                          ANH19390
      DO 109 II=1,NVIB                                                  ANH19400
      IVIB=NFRQ(II)                                                     ANH19410
      WRITE(6,8) II,IVIB                                                ANH19420
      DO 108 I=1,N3N                                                    ANH19430
      DO 108 J=1,N3N                                                    ANH19440
      DD(I,J)=F3Q(I,J,IVIB)                                             ANH19450
  108 CONTINUE
      CALL MATOUT(DD,N3N,N3N,N3N,N3N,6)                                 ANH19460
  109 CONTINUE                                                          ANH19470
C                                                                       ANH19480
  204 CONTINUE                                                          ANH19490
      DO 112 III=1,N3N                                                  ANH19500
      WK=FREQ(III)                                                      ANH19510
      DO 111 I=1,N3N                                                    ANH19520
      WI=FREQ(I)                                                        ANH19530
      DO 111 J=1,N3N                                                    ANH19540
      WJ=FREQ(J)                                                        ANH19550
      WIJK=WI*WJ*WK                                                     ANH19560
      IF(DABS(WK).LE.WLIMIT) GO TO 205                                  ANH19570
      IF(DABS(WI).LE.WLIMIT) GO TO 205                                  ANH19580
      IF(DABS(WJ).LE.WLIMIT) GO TO 205                                  ANH19590
      IF(DABS(WIJK).LE.WLIMIT) GO TO 205                                ANH19600
      SQWS=DSQRT(WIJK)                                                  ANH19610
      FACT=FACT3/SQWS                                                   ANH19620
      GO TO 206                                                         ANH19630
  205 CONTINUE                                                          ANH19640
      FACT=A00                                                          ANH19650
  206 CONTINUE                                                          ANH19660
      F3Q(I,J,III)=F3Q(I,J,III)*FACT                                    ANH19670
  111 CONTINUE                                                          ANH19680
  112 CONTINUE                                                          ANH19690
      DO 114 III=1,NVIB                                                 ANH19700
      IVIB=NFRQ(III)                                                    ANH19710
      WRITE(6,9) III,IVIB                                               ANH19720
      DO 113 I=1,NVIB                                                   ANH19730
      II=NFRQ(I)                                                        ANH19740
      DO 113 J=1,NVIB                                                   ANH19750
      JJ=NFRQ(J)                                                        ANH19760
      DD(I,J)=F3Q(II,JJ,IVIB)                                           ANH19770
  113 CONTINUE                                                          ANH19780
      CALL ANHOUT(DD,N3N,N3N,NVIB,NVIB,6)                               ANH19790
  114 CONTINUE                                                          ANH19800
C                                                                       ANH19810
      RETURN                                                            ANH19820
      END                                                               ANH19830
      SUBROUTINE PHIMAT(BXY,CXY,BXYZ,CXYZ,F3Q,EE)                       ANH19840
C   THE CALCULATION OF PHI MATRICES                                     ANH19850
      IMPLICIT REAL*8 (A-H,O-Z)                                         ANH19860
      DIMENSION BXY(3,3,N3N),CXY(3,3,N3N)                               ANH19870
      DIMENSION BXYZ(3,3,3,N3N),CXYZ(3,3,3,N3N)                         ANH19880
      DIMENSION F3Q(N3N,N3N,N3N),EE(N3N,N3N)                            ANH19890
      COMMON/VIB101/NATOM,N3N,NATRI,ILIN,NVIB                           ANH19900
      COMMON/VIB102/ITHREE,IFOUR,N3TOT,N4TOT                            ANH19910
      COMMON/VIB106/ROTAA(3),ROTGC(3),ROTCM(3),ROTMH(3)                 ANH19920
      COMMON/VIB108/IAXIS(3),NDEG(150),NDAB(150,5),IMAG(150)            ANH19930
      COMMON/VIB203/IOFF(150),IPRNT                                     ANH19940
      COMMON/VIB204/SQM(150),ROOT(150),FREQ(150)                        ANH19950
      COMMON/VIB205/IFREQ,NFRQ(150)                                     ANH19960
      COMMON/VIB207/TABCD(3,3,3,3),PABC(3,3,3)                          ANH19970
      DATA WLIMIT / 1.0D-02 /                                           ANH19980
      DATA RLIMIT / 1.0D-05 /                                           ANH19990
      DATA PH,AVN / 6.626176D+00 , 6.022045D+00 /                       ANH20000
      DATA PI,CL / 3.1415926536D+00 , 2.99792458D+00 /                  ANH20010
      DATA A00,HALF,ONE,TWO / 0.0D+00 , 0.5D+00 , 1.0D+00 , 2.0D+00 /   ANH20020
    1 FORMAT(//,2X,' ****************************'/                     ANH20030
     1          2X,' ***THE PHI (ABC) MATRICES***'/                     ANH20040
     2          2X,' ****************************')                     ANH20050
    2 FORMAT(//,2X,' PHI (A) MATIRX (IN D-08 CM-1)',                    ANH20060
     1         13X,' PHI (B) MATRIX (IN D-08 CM-1)',                    ANH20070
     2         13X,' PHI (C) MATRIX (IN D-08 CM-1)'/)                   ANH20080
    3 FORMAT(//,2X,' PHI (A) MATIRX (IN D-08 CM-1)'/)                   ANH20060
    4 FORMAT(//,2X,' PHI (B) MATIRX (IN D-08 CM-1)'/)                   ANH20060
    5 FORMAT(//,2X,' PHI (C) MATIRX (IN D-08 CM-1)'/)                   ANH20060
C                                                                       ANH20090
      WRITE(6,1)                                                        ANH20100
      IXX=1                                                             ANH20110
      IYY=2                                                             ANH20120
      IZZ=3                                                             ANH20130
C                                                                       ANH20140
C   FORM PHI MATRIX                                                     ANH20150
      CALL ZERO(PABC,3*3*3)                                             ANH20160
C                                                                       ANH20170
C   PHI(IA,IA,IA) ELEMENTS                                              ANH20180
      DO 105 IA=1,3                                                     ANH20190
      IAA=IAXIS(IA)                                                     ANH20200
      VALU1=A00                                                         ANH20210
      VALU2=A00                                                         ANH20220
      VALU3=A00                                                         ANH20230
      VALU4=A00                                                         ANH20240
      DO 101 IB=1,3                                                     ANH20250
      IBB=IAXIS(IB)                                                     ANH20260
      IF(ROTCM(IBB).LE.RLIMIT) GO TO 101                                ANH20270
      FACT1=TABCD(IA,IA,IA,IB)**2                                       ANH20280
      VALU1=VALU1+FACT1/ROTCM(IBB)                                      ANH20290
  101 CONTINUE                                                          ANH20300
      DO 102 IVIB=1,NVIB                                                ANH20310
      FRQI=FREQ(IVIB)                                                   ANH20320
      FACT2=CXYZ(IA,IA,IA,IVIB)**2                                      ANH20330
      VALU2=VALU2-FRQI*FACT2                                            ANH20340
  102 CONTINUE                                                          ANH20350
      DO 103 IVIB=1,NVIB                                                ANH20360
      DO 103 JVIB=1,NVIB                                                ANH20370
      DO 103 KVIB=1,NVIB                                                ANH20380
      FACT3=CXY(IA,IA,IVIB)*CXY(IA,IA,JVIB)*CXY(IA,IA,KVIB)             ANH20390
      VALU3=VALU3+FACT3*F3Q(IVIB,JVIB,KVIB)                             ANH20400
  103 CONTINUE                                                          ANH20410
      DO 104 IB=1,3                                                     ANH20420
      IBB=IAXIS(IB)                                                     ANH20430
      IF(IA.EQ.IB) GO TO 104                                            ANH20440
      FACT4=TABCD(IA,IA,IA,IB)**2                                       ANH20450
      ROTAB=ROTCM(IAA)-ROTCM(IBB)                                       ANH20460
      IF(DABS(ROTAB).LE.RLIMIT) GO TO 104                               ANH20470
      VALU4=VALU4+FACT4/ROTAB                                           ANH20480
  104 CONTINUE                                                          ANH20490
      VALUT=(VALU1*3.0D+00/16.0D+00)+VALU2*TWO                          ANH20500
     1      +VALU3/6.0D+00+VALU4/4.0D+00                                ANH20510
      PABC(IA,IA,IA)=VALUT                                              ANH20520
  105 CONTINUE                                                          ANH20530
C                                                                       ANH20540
C   PHI(IA,IA,IB) ELEMENTS                                              ANH20550
      DO 110 IA=1,3                                                     ANH20560
      DO 110 IB=1,3                                                     ANH20570
      IAA=IAXIS(IA)                                                     ANH20580
      IBB=IAXIS(IB)                                                     ANH20590
      IF(IA.EQ.IB) GO TO 110                                            ANH20600
      VALU1=A00                                                         ANH20610
      VALU2=A00                                                         ANH20620
      VALU3=A00                                                         ANH20630
      VALU4=A00                                                         ANH20640
      VALU5=A00                                                         ANH20650
      VALU6=A00                                                         ANH20660
      DO 106 IC=1,3                                                     ANH20670
      ICC=IAXIS(IC)                                                     ANH20680
      IF(ROTCM(ICC).LE.RLIMIT) GO TO 106                                ANH20690
      FAC1=TABCD(IA,IA,IB,IC)+TABCD(IA,IB,IA,IC)*TWO                    ANH20700
      FAC2=TABCD(IB,IB,IA,IC)+TABCD(IA,IB,IB,IC)*TWO                    ANH20710
      FACT1=FAC1*FAC1+TABCD(IA,IA,IA,IC)*FAC2*TWO                       ANH20720
      VALU1=VALU1+FACT1/ROTCM(ICC)                                      ANH20730
  106 CONTINUE                                                          ANH20740
      DO 107 IVIB=1,NVIB                                                ANH20750
      FRQI=FREQ(IVIB)                                                   ANH20760
      FACT2=CXYZ(IA,IA,IB,IVIB)**2                                      ANH20770
     1     +CXYZ(IA,IA,IA,IVIB)*CXYZ(IA,IB,IB,IVIB)*TWO                 ANH20780
      VALU2=VALU2-FACT2*FRQI                                            ANH20790
  107 CONTINUE                                                          ANH20800
      DO 108 IVIB=1,NVIB                                                ANH20810
      DO 108 JVIB=1,NVIB                                                ANH20820
      DO 108 KVIB=1,NVIB                                                ANH20830
      FACT3=CXY(IA,IA,IVIB)*(CXY(IA,IA,JVIB)*CXY(IB,IB,KVIB)            ANH20840
     1                      +CXY(IA,IB,JVIB)*CXY(IA,IB,KVIB)*4.0D+00)   ANH20850
      VALU3=VALU3+FACT3*F3Q(IVIB,JVIB,KVIB)                             ANH20860
  108 CONTINUE                                                          ANH20870
      ROTAB=ROTCM(IAA)-ROTCM(IBB)                                       ANH20880
      IF(DABS(ROTAB).LE.RLIMIT) GO TO 201                               ANH20890
      FACT4=TABCD(IA,IA,IA,IB)*(TABCD(IB,IB,IB,IA)*4.0D+00              ANH20900
     1                         -TABCD(IA,IA,IA,IB)*3.0D+00)             ANH20910
      VALU4=FACT4/ROTAB                                                 ANH20920
  201 CONTINUE                                                          ANH20930
      DO 109 IC=1,3                                                     ANH20940
      ICC=IAXIS(IC)                                                     ANH20950
      IF(IC.EQ.IA) GO TO 109                                            ANH20960
      IF(IC.EQ.IB) GO TO 109                                            ANH20970
      ROTAC=ROTCM(IAA)-ROTCM(ICC)                                       ANH20980
      ROTAC2=ROTAC*ROTAC                                                ANH20990
      ROTAB=ROTCM(IAA)-ROTCM(IBB)                                       ANH21000
      IF(DABS(ROTAC).LE.RLIMIT) GO TO 202                               ANH21010
      FACT5=ROTAC*(TABCD(IB,IB,IA,IC)+TABCD(IB,IA,IB,IC)*TWO)           ANH21020
     1     +ROTAB*(TABCD(IC,IC,IC,IA)-TABCD(IA,IA,IA,IC))               ANH21030
      VALU5=VALU5+(FACT5*TABCD(IA,IA,IA,IC))/ROTAC2                     ANH21040
  202 CONTINUE                                                          ANH21050
      ROTBC=ROTCM(IBB)-ROTCM(ICC)                                       ANH21060
      ROTBC2=ROTBC*ROTBC                                                ANH21070
      IF(DABS(ROTBC).LE.RLIMIT) GO TO 109                               ANH21080
      FACT6=ROTBC*(TABCD(IA,IA,IB,IC)+TABCD(IA,IB,IA,IC)*TWO)           ANH21090
     1     +ROTAB*(TABCD(IB,IB,IB,IC)-TABCD(IC,IC,IC,IB))*TWO           ANH21100
      VALU6=VALU6+(FACT6*(TABCD(IA,IA,IB,IC)                            ANH21110
     1                   +TABCD(IA,IB,IA,IC)*TWO))/ROTBC2               ANH21120
  109 CONTINUE                                                          ANH21130
      VALUT=(VALU1*3.0D+00/32.0D+00)+VALU2+VALU3/4.0D+00                ANH21140
     1      +VALU4/8.0D+00+VALU5/4.0D+00+VALU6/8.0D+00                  ANH21150
      PABC(IA,IA,IB)=VALUT                                              ANH21160
      PABC(IA,IB,IA)=VALUT                                              ANH21170
      PABC(IB,IA,IA)=VALUT                                              ANH21180
  110 CONTINUE                                                          ANH21190
C                                                                       ANH21200
C   PHI(IA,IB,IC) ELEMENTS                                              ANH21210
      VALU1=A00                                                         ANH21220
      VALU2=A00                                                         ANH21230
      VALU3=A00                                                         ANH21240
      VALU4=A00                                                         ANH21250
      VALU5=A00                                                         ANH21260
      VALU6=A00                                                         ANH21270
      DO 111 IA=1,3                                                     ANH21280
      IAA=IAXIS(IA)                                                     ANH21290
      IF(ROTCM(IAA).LE.RLIMIT) GO TO 111                                ANH21300
      FACT=TABCD(IA,IXX,IYY,IZZ)+TABCD(IA,IYY,IZZ,IXX)                  ANH21310
     1    +TABCD(IA,IZZ,IXX,IYY)                                        ANH21320
      FAC1=FACT*FACT*TWO                                                ANH21330
      FAC2=(TABCD(IA,IXX,IYY,IYY)+TABCD(IA,IYY,IXX,IYY)*TWO)            ANH21340
     1    *(TABCD(IA,IXX,IZZ,IZZ)+TABCD(IA,IZZ,IXX,IZZ)*TWO)            ANH21350
      FAC3=(TABCD(IA,IYY,IZZ,IZZ)+TABCD(IA,IZZ,IYY,IZZ)*TWO)            ANH21360
     1    *(TABCD(IA,IYY,IXX,IXX)+TABCD(IA,IXX,IYY,IXX)*TWO)            ANH21370
      FAC4=(TABCD(IA,IZZ,IXX,IXX)+TABCD(IA,IXX,IZZ,IXX)*TWO)            ANH21380
     1    *(TABCD(IA,IZZ,IYY,IYY)+TABCD(IA,IYY,IZZ,IYY)*TWO)            ANH21390
      VALU1=VALU1+(FAC1+FAC2+FAC3+FAC4)/ROTCM(IAA)                      ANH21400
  111 CONTINUE                                                          ANH21410
      DO 112 IVIB=1,NVIB                                                ANH21420
      FRQI=FREQ(IVIB)                                                   ANH21430
      FACT2=CXYZ(IXX,IYY,IZZ,IVIB)**2*TWO                               ANH21440
     1     +CXYZ(IYY,IYY,IXX,IVIB)*CXYZ(IZZ,IZZ,IXX,IVIB)               ANH21450
     2     +CXYZ(IZZ,IZZ,IYY,IVIB)*CXYZ(IXX,IXX,IYY,IVIB)               ANH21460
     3     +CXYZ(IXX,IXX,IZZ,IVIB)*CXYZ(IYY,IYY,IZZ,IVIB)               ANH21470
      VALU2=VALU2-FRQI*FACT2                                            ANH21480
  112 CONTINUE                                                          ANH21490
      DO 113 IVIB=1,NVIB                                                ANH21500
      DO 113 JVIB=1,NVIB                                                ANH21510
      DO 113 KVIB=1,NVIB                                                ANH21520
      FACT3=CXY(IXX,IXX,IVIB)*CXY(IYY,IYY,JVIB)*CXY(IZZ,IZZ,KVIB)       ANH21530
     1     +CXY(IXX,IXX,IVIB)*CXY(IYY,IZZ,JVIB)*CXY(IYY,IZZ,KVIB)*TWO   ANH21540
     2     +CXY(IYY,IYY,IVIB)*CXY(IZZ,IXX,JVIB)*CXY(IZZ,IXX,KVIB)*TWO   ANH21550
     3     +CXY(IZZ,IZZ,IVIB)*CXY(IXX,IYY,JVIB)*CXY(IXX,IYY,KVIB)*TWO   ANH21560
     4     +CXY(IYY,IZZ,IVIB)*CXY(IZZ,IXX,JVIB)*CXY(IXX,IYY,KVIB)       ANH21570
     5     *8.0D+00                                                     ANH21580
      VALU3=VALU3+FACT3*F3Q(IVIB,JVIB,KVIB)                             ANH21590
  113 CONTINUE                                                          ANH21600
      IAA=IAXIS(1)                                                      ANH21610
      IBB=IAXIS(2)                                                      ANH21620
      ICC=IAXIS(3)                                                      ANH21630
      ROTYZ=ROTCM(IBB)-ROTCM(ICC)                                       ANH21640
      ROTYZ2=ROTYZ*ROTYZ                                                ANH21650
      ROTZX=ROTCM(ICC)-ROTCM(IAA)                                       ANH21660
      ROTZX2=ROTZX*ROTZX                                                ANH21670
      ROTXY=ROTCM(IAA)-ROTCM(IBB)                                       ANH21680
      ROTXY2=ROTXY*ROTXY                                                ANH21690
      IF(DABS(ROTYZ).LE.RLIMIT) GO TO 203                               ANH21700
      FACT4=ROTZX*TABCD(IYY,IYY,IYY,IZZ)                                ANH21710
     1     +ROTXY*TABCD(IZZ,IZZ,IZZ,IYY)                                ANH21720
      VALU4=(TABCD(IYY,IYY,IYY,IZZ)-TABCD(IZZ,IZZ,IZZ,IYY))*FACT4       ANH21730
     1      /ROTYZ2                                                     ANH21740
  203 CONTINUE                                                          ANH21750
      IF(DABS(ROTZX).LE.RLIMIT) GO TO 204                               ANH21760
      FACT5=ROTXY*TABCD(IZZ,IZZ,IZZ,IXX)                                ANH21770
     1     +ROTYZ*TABCD(IXX,IXX,IXX,IZZ)                                ANH21780
      VALU5=(TABCD(IZZ,IZZ,IZZ,IXX)-TABCD(IXX,IXX,IXX,IZZ))*FACT5       ANH21790
     1      /ROTZX2                                                     ANH21800
  204 CONTINUE                                                          ANH21810
      IF(DABS(ROTXY).LE.RLIMIT) GO TO 205                               ANH21820
      FACT6=ROTYZ*TABCD(IXX,IXX,IXX,IYY)                                ANH21830
     1     +ROTZX*TABCD(IYY,IYY,IYY,IXX)                                ANH21840
      VALU6=(TABCD(IXX,IXX,IXX,IYY)-TABCD(IYY,IYY,IYY,IXX))*FACT6       ANH21850
     1      /ROTXY2                                                     ANH21860
  205 CONTINUE                                                          ANH21870
      VALUT=(VALU1*3.0D+00/16.0D+00)+VALU2*TWO+VALU3*HALF               ANH21880
     1     +(VALU4*3.0D+00)/4.0D+00+(VALU5*3.0D+00)/4.0D+00             ANH21890
     2     +(VALU6*3.0D+00)/4.0D+00                                     ANH21900
      PABC(IXX,IYY,IZZ)=VALUT                                           ANH21910
      PABC(IXX,IZZ,IYY)=VALUT                                           ANH21920
      PABC(IYY,IXX,IZZ)=VALUT                                           ANH21930
      PABC(IYY,IZZ,IXX)=VALUT                                           ANH21940
      PABC(IZZ,IXX,IYY)=VALUT                                           ANH21950
      PABC(IZZ,IYY,IXX)=VALUT                                           ANH21960
C                                                                       ANH21970
      DO 115 IA=1,3                                                     ANH21980
      DO 115 IB=1,3                                                     ANH21990
      DO 115 IC=1,3                                                     ANH22000
      PABC(IA,IB,IC)=PABC(IA,IB,IC)*1.0D+08                             ANH22010
  115 CONTINUE                                                          ANH22020
c  Commented out for a stand-alone version
c     WRITE(6,2)                                                        ANH22030
c     CALL MX3OUT(PABC(1,1,1),3,PABC(1,1,2),3,PABC(1,1,3),3,3,6)        ANH22040
      WRITE(6,3)                                                        ANH22030
      CALL MX1OUT(PABC(1,1,1),3,3,6) 
      WRITE(6,4)                                                        ANH22030
      CALL MX1OUT(PABC(1,1,2),3,3,6) 
      WRITE(6,5)                                                        ANH22030
      CALL MX1OUT(PABC(1,1,3),3,3,6) 
      DO 116 IA=1,3                                                     ANH22050
      DO 116 IB=1,3                                                     ANH22060
      DO 116 IC=1,3                                                     ANH22070
      PABC(IA,IB,IC)=PABC(IA,IB,IC)/1.0D+08                             ANH22080
  116 CONTINUE                                                          ANH22090
C                                                                       ANH22100
      RETURN                                                            ANH22110
      END                                                               ANH22120
      SUBROUTINE FORCE4(ELX,F44,F4X,F4Q,EE,DD,GG)                       ANH22130
C   THE FOURTH ORDER ANHARMONICITY ANALYSIS                             ANH22140
      IMPLICIT REAL*8 (A-H,O-Z)                                         ANH22150
      DIMENSION ELX(N3N,N3N)                                            ANH22160
      DIMENSION F44(N4TOT),F4X(N3N,N3N,N3N,N3N),F4Q(N3N,N3N,N3N,N3N)    ANH22170
      DIMENSION EE(N3N,N3N),DD(N3N,N3N),GG(N3N,N3N)                     ANH22180
      COMMON/VIB101/NATOM,N3N,NATRI,ILIN,NVIB                           ANH22190
      COMMON/VIB102/ITHREE,IFOUR,N3TOT,N4TOT                            ANH22200
      COMMON/VIB105/FACT3,FACT4                                         ANH22210
      COMMON/VIB203/IOFF(150),IPRNT                                     ANH22220
      COMMON/VIB204/SQM(150),ROOT(150),FREQ(150)                        ANH22230
      COMMON/VIB205/IFREQ,NFRQ(150)                                     ANH22240
      DATA A00,ONE / 0.0D+00 , 1.0D+00 /                                ANH22250
      DATA WLIMIT / 3.0D+00 /                                           ANH22260
      DATA ITAP24 / 24 /                                                ANH22270
      DATA A0,HE / 0.529177249D+00 , 4.359813653D+00 /   
    1 FORMAT(//,2X,' THE FOURTH ORDER ANHARMONICITY ANALYSIS'/)         ANH22290
    2 FORMAT(2I5)                                                       ANH22300
    3 FORMAT(3F20.10)                                                   ANH22310
    4 FORMAT(//,2X,' **************************************************'ANH22320
     1        /,2X,' ***FOURTH DERIVATIVES W.R.T. NORMAL COORDINATES***'ANH22330
     2        /,2X,' **************************************************'ANH22340
     3        /)                                                        ANH22350
    5 FORMAT(//,2X,' F4X MATRIX (IN HARTREE/BOHR+4), KABC = ',I5/)      ANH22360
    6 FORMAT(//,2X,' F4X MATRIX (IN MDYN/A+3), KABC = ',I5,             ANH22370
     1          5X,' LABC = ',I5/)                                      ANH22380
    7 FORMAT(//,2X,' LX(T) * F4X * LX MATRIX, KABC = ',I5,              ANH22390
     1          5X,' LABC = ',I5/)                                      ANH22400
    8 FORMAT(//,2X,' LX(T)*LX(T) * F4X *LX*LX MATRIX, KABC = ',I5,      ANH22410
     1          5X,' LABC = ',I5/)                                      ANH22420
    9 FORMAT(//,2X,' F4(QI*QJ*QK*QL) MATRIX (IN MDYN/(A+2.AMU-3/2))',   ANH22430
     1 1X,' QVIB (',I3,',',I3,')',2X,',',' IVIB (',I3,',',I3,')'/)      ANH22440
   10 FORMAT(//,2X,' F4(QI*QJ*QK*QL) MATRIX (IN CM-1),',                ANH22450
     1 1X,' QVIB (',I3,',',I3,')',2X,',',' IVIB (',I3,',',I3,')'/)      ANH22460
C                                                                       ANH22470
      IF(IPRNT.NE.0)                                                    ANH22480
     *WRITE(6,1)                                                        ANH22490
C                                                                       ANH22500
C   FIND CONVERSION FACTORS                                             ANH22510
C   FUNIT4 IS IN MDYN/A+3                                               ANH22520
      FUNIT4=HE/(A0*A0*A0*A0)                                           ANH22530
C                                                                       ANH22540
C   READ IN FOURTH DERIVATIVES                                          ANH22550
      REWIND ITAP24                                                     ANH22560
      READ(ITAP24,2) NNATOM,NF3                                         ANH22570
      READ(ITAP24,3) (F44(I),I=1,N4TOT)                                 ANH22580
      REWIND ITAP24                                                     ANH22590
C                                                                       ANH22600
C   TRANSFORM FOURTH DERIVATIVES                                        ANH22610
C   FROM CARTESIAN TO NORMAL COORDINATE SYSTEM                          ANH22620
      WRITE(6,4)                                                        ANH22630
      MABC=0                                                            ANH22640
      DO 101 IABC=1,N3N                                                 ANH22650
      DO 101 JABC=1,IABC                                                ANH22660
      DO 101 KABC=1,JABC                                                ANH22670
      DO 101 LABC=1,KABC                                                ANH22680
      MABC=MABC+1                                                       ANH22690
      VAL=F44(MABC)                                                     ANH22700
C                                                                       ANH22710
      F4X(IABC,JABC,KABC,LABC)=VAL                                      ANH22720
      F4X(IABC,JABC,LABC,KABC)=VAL                                      ANH22730
      F4X(IABC,KABC,JABC,LABC)=VAL                                      ANH22740
      F4X(IABC,KABC,LABC,JABC)=VAL                                      ANH22750
      F4X(IABC,LABC,JABC,KABC)=VAL                                      ANH22760
      F4X(IABC,LABC,KABC,JABC)=VAL                                      ANH22770
C                                                                       ANH22780
      F4X(JABC,IABC,KABC,LABC)=VAL                                      ANH22790
      F4X(JABC,IABC,LABC,KABC)=VAL                                      ANH22800
      F4X(JABC,KABC,IABC,LABC)=VAL                                      ANH22810
      F4X(JABC,KABC,LABC,IABC)=VAL                                      ANH22820
      F4X(JABC,LABC,IABC,KABC)=VAL                                      ANH22830
      F4X(JABC,LABC,KABC,IABC)=VAL                                      ANH22840
C                                                                       ANH22850
      F4X(KABC,IABC,JABC,LABC)=VAL                                      ANH22860
      F4X(KABC,IABC,LABC,JABC)=VAL                                      ANH22870
      F4X(KABC,JABC,IABC,LABC)=VAL                                      ANH22880
      F4X(KABC,JABC,LABC,IABC)=VAL                                      ANH22890
      F4X(KABC,LABC,IABC,JABC)=VAL                                      ANH22900
      F4X(KABC,LABC,JABC,IABC)=VAL                                      ANH22910
C                                                                       ANH22920
      F4X(LABC,IABC,JABC,KABC)=VAL                                      ANH22930
      F4X(LABC,IABC,KABC,JABC)=VAL                                      ANH22940
      F4X(LABC,JABC,IABC,KABC)=VAL                                      ANH22950
      F4X(LABC,JABC,KABC,IABC)=VAL                                      ANH22960
      F4X(LABC,KABC,IABC,JABC)=VAL                                      ANH22970
      F4X(LABC,KABC,JABC,IABC)=VAL                                      ANH22980
  101 CONTINUE                                                          ANH22990
C                                                                       ANH23000
      DO 105 KABC=1,N3N                                                 ANH23010
      DO 105 LABC=1,N3N                                                 ANH23020
      DO 102 I=1,N3N                                                    ANH23030
      DO 102 J=1,N3N                                                    ANH23040
      DD(I,J)=F4X(I,J,KABC,LABC)                                        ANH23050
  102 CONTINUE
      IF(IPRNT.LE.2) GO TO 201                                          ANH23060
      WRITE(6,5) KABC,LABC                                              ANH23070
      CALL MATOUT(DD,N3N,N3N,N3N,N3N,6)                                 ANH23080
  201 CONTINUE                                                          ANH23090
      DO 103 I=1,N3N                                                    ANH23100
      DO 103 J=1,N3N                                                    ANH23110
      DD(I,J)=DD(I,J)*FUNIT4                                            ANH23120
  103 CONTINUE
      IF(IPRNT.LE.2) GO TO 202                                          ANH23130
      WRITE(6,6) KABC,LABC                                              ANH23140
      CALL MATOUT(DD,N3N,N3N,N3N,N3N,6)                                 ANH23150
C                                                                       ANH23160
  202 CONTINUE                                                          ANH23170
      CALL MTXMPY(ELX,N3N,DD,N3N,EE,N3N,GG,N3N,N3N,5)                   ANH23180
      IF(IPRNT.LE.2) GO TO 203                                          ANH23190
      WRITE(6,7) KABC,LABC                                              ANH23200
      CALL MATOUT(EE,N3N,N3N,N3N,N3N,6)                                 ANH23210
  203 CONTINUE                                                          ANH23220
      DO 104 I=1,N3N                                                    ANH23230
      DO 104 J=1,N3N                                                    ANH23240
      F4X(I,J,KABC,LABC)=EE(I,J)                                        ANH23250
  104 CONTINUE
  105 CONTINUE                                                          ANH23260
C                                                                       ANH23270
      DO 110 KABC=1,N3N                                                 ANH23280
      DO 110 LABC=1,N3N                                                 ANH23290
      DO 106 I=1,N3N                                                    ANH23300
      DO 106 J=1,N3N                                                    ANH23310
      DD(I,J)=F4X(KABC,LABC,I,J)                                        ANH23320
  106 CONTINUE
      IF(IPRNT.LE.2) GO TO 204                                          ANH23330
      WRITE(6,7) KABC,LABC                                              ANH23340
      CALL MATOUT(DD,N3N,N3N,N3N,N3N,6)                                 ANH23350
  204 CONTINUE                                                          ANH23360
C                                                                       ANH23370
      CALL MTXMPY(ELX,N3N,DD,N3N,EE,N3N,GG,N3N,N3N,5)                   ANH23380
      IF(IPRNT.LE.2) GO TO 206                                          ANH23390
      WRITE(6,8) KABC,LABC                                              ANH23400
      CALL MATOUT(EE,N3N,N3N,N3N,N3N,6)                                 ANH23410
  206 CONTINUE                                                          ANH23420
      DO 108 I=1,N3N                                                    ANH23430
      DO 108 J=1,N3N                                                    ANH23440
      F4Q(KABC,LABC,I,J)=EE(I,J)                                        ANH23450
  108 CONTINUE
  110 CONTINUE                                                          ANH23460
C                                                                       ANH23470
C   SET NUMBER OF ANHARMONICITY SETS                                    ANH23480
C   REGULARLY NVIB=3N-6(OR 3N-5)                                        ANH23490
C   IF IFREQ IS NOT EQUAL TO 0 THEN NVIB=IFREQ                          ANH23500
      IF(IFREQ.NE.0) NVIB=IABS(IFREQ)                                   ANH23510
C                                                                       ANH23520
      IF(IPRNT.LE.2) GO TO 207                                          ANH23530
      DO 112 II=1,NVIB                                                  ANH23540
      IVIB=NFRQ(II)                                                     ANH23550
      DO 112 JJ=1,NVIB                                                  ANH23560
      JVIB=NFRQ(JJ)                                                     ANH23570
      WRITE(6,9) II,JJ,IVIB,JVIB                                        ANH23580
      DO 111 I=1,N3N                                                    ANH23590
      DO 111 J=1,N3N                                                    ANH23600
      DD(I,J)=F4Q(I,J,IVIB,JVIB)                                        ANH23610
  111 CONTINUE
      CALL MATOUT(DD,N3N,N3N,N3N,N3N,6)                                 ANH23620
  112 CONTINUE                                                          ANH23630
C                                                                       ANH23640
  207 CONTINUE                                                          ANH23650
      DO 115 III=1,N3N                                                  ANH23660
      WK=FREQ(III)                                                      ANH23670
      DO 115 JJJ=1,N3N                                                  ANH23680
      WL=FREQ(JJJ)                                                      ANH23690
      DO 114 I=1,N3N                                                    ANH23700
      WI=FREQ(I)                                                        ANH23710
      DO 114 J=1,N3N                                                    ANH23720
      WJ=FREQ(J)                                                        ANH23730
      WIJKL=WI*WJ*WK*WL                                                 ANH23740
      IF(DABS(WK).LE.WLIMIT) GO TO 208                                  ANH23750
      IF(DABS(WL).LE.WLIMIT) GO TO 208                                  ANH23760
      IF(DABS(WI).LE.WLIMIT) GO TO 208                                  ANH23770
      IF(DABS(WJ).LE.WLIMIT) GO TO 208                                  ANH23780
      IF(DABS(WIJKL).LE.WLIMIT) GO TO 208                               ANH23790
      SQWS=DSQRT(WIJKL)                                                 ANH23800
      FACT=FACT4/SQWS                                                   ANH23810
      GO TO 209                                                         ANH23820
  208 CONTINUE                                                          ANH23830
      FACT=A00                                                          ANH23840
  209 CONTINUE                                                          ANH23850
      F4Q(I,J,III,JJJ)=F4Q(I,J,III,JJJ)*FACT                            ANH23860
  114 CONTINUE                                                          ANH23870
  115 CONTINUE                                                          ANH23880
      DO 120 III=1,NVIB                                                 ANH23890
      IVIB=NFRQ(III)                                                    ANH23900
      DO 120 JJJ=1,NVIB                                                 ANH23910
      JVIB=NFRQ(JJJ)                                                    ANH23920
      WRITE(6,10) III,JJJ,IVIB,JVIB                                     ANH23930
      DO 116 I=1,NVIB                                                   ANH23940
      II=NFRQ(I)                                                        ANH23950
      DO 116 J=1,NVIB                                                   ANH23960
      JJ=NFRQ(J)                                                        ANH23970
      DD(I,J)=F4Q(II,JJ,IVIB,JVIB)                                      ANH23980
  116 CONTINUE                                                          ANH23990
      CALL ANHOUT(DD,N3N,N3N,NVIB,NVIB,6)                               ANH24000
  120 CONTINUE                                                          ANH24010
C                                                                       ANH24020
      RETURN                                                            ANH24030
      END                                                               ANH24040
      BLOCK DATA                                                        ANH25240
      CHARACTER*1  RAXIS                                                ANH25250
      CHARACTER*25 RTYPE                                                ANH25260
      COMMON/VIB206/RAXIS(3),RTYPE(6)                                   ANH25270
C                                                                       ANH25280
      DATA RAXIS(1) / 'A' /                                             ANH25290
      DATA RAXIS(2) / 'B' /                                             ANH25300
      DATA RAXIS(3) / 'C' /                                             ANH25310
      DATA RTYPE(1) / ' DIATOMIC                ' /                     ANH25320
      DATA RTYPE(2) / ' LINEAR POLYATOMIC       ' /                     ANH25330
      DATA RTYPE(3) / ' ASYMMETRIC TOP          ' /                     ANH25340
      DATA RTYPE(4) / ' SYMMETRIC TOP (PROLATE) ' /                     ANH25350
      DATA RTYPE(5) / ' SYMMETRIC TOP (OBLATE)  ' /                     ANH25360
      DATA RTYPE(6) / ' SPHERICAL TOP           ' /                     ANH25370
      END                                                               ANH25380
