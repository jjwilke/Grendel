C     PROGRAM ANHDIA                                                    ANH00010
C******************************************************                 ANH00020
C***LAST UPDATED ON JUNE 04, 1986 BY YUKIO YAMAGUCHI***                 ANH00030
C******************************************************                 ANH00040
      SUBROUTINE DIATOM(F3Q,F4Q,EE)                                     ANH00050
C   THE ANHARMONICITY CONTSTANTS FOR A DIATOMIC MOLECULE                ANH00060
      IMPLICIT REAL*8 (A-H,O-Z)                                         ANH00070
      DIMENSION F3Q(N3N,N3N,N3N),F4Q(N3N,N3N,N3N,N3N),EE(N3N,N3N)       ANH00080
      COMMON/VIB101/NATOM,N3N,NATRI,ILIN,NVIB                           ANH00090
      COMMON/VIB106/ROTAA(3),ROTGC(3),ROTCM(3),ROTMH(3)                 ANH00100
      COMMON/VIB203/IOFF(150),IPRNT                                     ANH00110
      COMMON/VIB204/SQM(150),ROOT(150),FREQ(150)                        ANH00120
      DATA CL / 2.99792458D+00 /                                        ANH00130
    1 FORMAT(//,2X,' *************************************'/            ANH00140
     1          2X,' ***CENTRIFUGAL DISTORTION CONSTANT***'/            ANH00150
     2          2X,' *************************************'//           ANH00160
     3          2X,' DJ = ',F20.10,' D-06 CM-1',                        ANH00170
     4          5X,' DJ = ',F20.10,' MHZ'/)                             ANH00180
    2 FORMAT(//,2X,' *********************************************'/    ANH00190
     1          2X,' ***VIBRATION-ROTATION INTERACTION CONSTANT***'/    ANH00200
     2          2X,' *********************************************'//   ANH00210
     3          2X,' ALPHA = ',F20.10,' CM-1',                          ANH00220
     4          7X,' ALPHA = ',F20.10,' MHZ'/)                          ANH00230
    3 FORMAT(//,2X,' *************************************'/            ANH00240
     1          2X,' ***THEORETICAL FREQUENCIES IN CM-1***'/            ANH00250
     2          2X,' *************************************'//           ANH00260
     3          2X,' HARMONIC',14X,' DELTA',13X,' FUNDAMENTAL'//        ANH00270
     4          1X,F13.3,8X,F13.3,8X,F13.3/)                            ANH00280
C                                                                       ANH00290
      BE=ROTCM(1)                                                       ANH00300
      BE2=BE*BE                                                         ANH00310
      BE3=BE2*BE                                                        ANH00320
      WE=FREQ(1)                                                        ANH00330
      WE2=WE*WE                                                         ANH00340
      P3=F3Q(1,1,1)                                                     ANH00350
      P4=F4Q(1,1,1,1)                                                   ANH00360
C                                                                       ANH00370
      DE=(BE3/WE2)*4.0D+00                                              ANH00380
      DECM=DE*1.0D+06                                                   ANH00390
      DEMH=DE*CL*1.0D+04                                                ANH00400
      WRITE(6,1) DECM,DEMH                                              ANH00410
C                                                                       ANH00420
      ALPHA=-DSQRT(BE3/WE)*(DSQRT(BE/WE)+(DSQRT(2.0D+00)*P3)/WE)        ANH00430
      ALPMH=ALPHA*CL*1.0D+04                                            ANH00440
      WRITE(6,2) ALPHA,ALPMH                                            ANH00450
C                                                                       ANH00460
      X11=P4/16.0D+00-(P3*P3*5.0D+00/WE)/48.0D+00                       ANH00470
      FW=WE+X11                                                         ANH00480
      WRITE(6,3) WE,X11,FW                                              ANH00490
C                                                                       ANH00500
      RETURN                                                            ANH00510
      END                                                               ANH00520
