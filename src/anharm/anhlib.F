      SUBROUTINE FFILE(ITAPE,name,code)
C
C USE THIS TO OPEN ALL FORMATTED SEQUENTAL ACCESS FILES
C
      IMPLICIT INTEGER (A-Z)
C
      integer code
      character*1 digit
      character*12 file
      character*(*) name
      character*1 junk

      if (name .eq. ' ') then
        file='file'//digit(itape/10)//digit(itape-10*(itape/10))
     &       //'.dat'
      else
        file=name//'.dat'
        endif

      open (unit=itape,access='SEQUENTIAL',
     &      form='FORMATTED',
     &      iostat=ierr,status='UNKNOWN',file=file)

      IF (IERR .ne. 0) then
        WRITE(*,*) ' ERROR ENCOUNTERED OPENING FILE ',ITAPE,' IERR =',
     &            ierr
        stop
        endif

      if (iand(code,1) .eq. 1) then
c       Go to the end of the file
 800    continue
        read(itape,'(A1)',end=900) junk
        goto 800
 900    continue
        endif
C
      RETURN
      END
      character*1 function digit(i)
      character*1 digits(0:9)
      data digits/'0','1','2','3','4','5','6','7','8','9'/
      digit = digits(i)
      return
      end
      SUBROUTINE ZERO(A,N)
C
      REAL*8 A(N)
C
      DO 101 I=1,N
           A(I)=0.0D+00
  101 CONTINUE
C
      RETURN
      END
      SUBROUTINE LOCATE(INPUT,TOKEN,IERROR)
C
C     ----- SEARCH THROUGH INPUT FILE FOR TOKEN BEGINNING
C           WITH # TO LOCATE INPUT FOR PROGRAM.  IERROR IS
C           SET TO 0 IF NO ERRORS, 1 IF ANY ERROR OCCURS.
C
C
      CHARACTER*10 TOKEN,LINE
C
      REWIND (UNIT=INPUT,ERR=99)
C
    1 CONTINUE
      READ (UNIT=INPUT,FMT='(A10)',END=99,ERR=99) LINE
      IF (LINE.NE.TOKEN) GO TO 1
C
      IERROR=0
      RETURN
C
C
   99 CONTINUE
      IERROR=1
      RETURN
      END
      SUBROUTINE PRINT(A,NAD,M,IOUT)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*6 LINE
      DIMENSION A(NAD)
      COMMON/VIB112/LINEWD
c     The number of data per line
      integer ndpl
    1 FORMAT(1X,10(7X,I5))
    2 FORMAT(1X,21A6)
    3 FORMAT(1X,I4,2X,10F12.7)
    4 FORMAT(/)
      DATA LINE / '------' /
C
c The width of the output page.
      ndpl = (linewd-7)/12
      if ((ndpl.le.0).or.(ndpl.gt.200)) then
        write(6,*) 'ERROR: bogus ndpl in print: check ''LINEWIDTH'''
        stop
      endif
C
      II=0
      JJ=0
  200 CONTINUE
      II=II+1
      JJ=JJ+1
      KK=ndpl*JJ
      NN=KK+KK*(KK-1)/2
      MM=M
      IF(M.GT.KK) MM=KK
      LL=2*(MM-II+1)+1
      WRITE (IOUT,1) (I,I=II,MM)
      WRITE (IOUT,2) (LINE,I=1,LL)
      DO 101 I=II,M
      I1=I*(I-1)/2+II
      I2=I+I*(I-1)/2
      IF(I2.GT.NN) I2=I1+ndpl-1
      WRITE (IOUT,3) I,(A(J),J=I1,I2)
  101 CONTINUE
      IF(M.LE.KK) GO TO 201
      WRITE (IOUT,4)
      II=KK
      GO TO 200
C
  201 CONTINUE
      RETURN
      END
      SUBROUTINE MATOUT(A,NAD,NBD,M,N,IOUT)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*6 LINE
      DIMENSION A(NAD,NBD)
      COMMON/VIB112/LINEWD
c     The number of data per line
      integer ndpl
    1 FORMAT(3X,10(7X,I5))
    2 FORMAT(1X,21A6)
    3 FORMAT(2X,I3,2X,10F12.7)
    4 FORMAT(/)
      DATA LINE / '------' /
C
      ndpl = (linewd-7)/12
      if ((ndpl.le.0).or.(ndpl.gt.200)) then
        write(6,*) 'ERROR: bogus ndpl in matout: check ''LINEWIDTH'''
        stop
      endif
C
      II=0
      JJ=0
  200 CONTINUE
      II=II+1
      JJ=JJ+1
      KK=ndpl*JJ
      NN=N
      IF(N.GT.KK) NN=KK
      LL=2*(NN-II+1)+1
      WRITE (IOUT,1) (I,I=II,NN)
      WRITE (IOUT,2) (LINE,I=1,LL)
      DO 101 I=1,M
      WRITE (IOUT,3) I,(A(I,J),J=II,NN)
  101 CONTINUE
      IF(N.LE.KK) GO TO 201
      WRITE (IOUT,4)
      II=KK
      GO TO 200
C
  201 CONTINUE
      RETURN
      END
      SUBROUTINE MTXMPY(A,NAD,B,NBD,C,NCD,D,NDD,N,IMPY)                 ANH24050
      IMPLICIT REAL*8 (A-H,O-Z)                                         ANH24060
      DIMENSION A(NAD,NAD),B(NBD,NBD),C(NCD,NCD),D(NDD,NDD)             ANH24070
      DATA A00 / 0.0D+00 /                                              ANH24080
C                                                                       ANH24090
      CALL ZERO(C,NCD*NCD)                                              ANH24100
      GO TO (201,202,203,204,205,206),IMPY                              ANH24110
C
C     C=A*B                                                             ANH24120
  201 CONTINUE
      DO 102 I=1,N                                                      ANH24130
      DO 102 J=1,N                                                      ANH24140
      DO 102 K=1,N                                                      ANH24150
      C(I,J)=C(I,J)+A(I,K)*B(K,J)                                       ANH24160
  102 CONTINUE
      RETURN                                                            ANH24170
C
C     C=TA*B                                                            ANH24180
  202 CONTINUE
      DO 103 I=1,N                                                      ANH24190
      DO 103 J=1,N                                                      ANH24200
      DO 103 K=1,N                                                      ANH24210
      C(I,J)=C(I,J)+A(K,I)*B(K,J)                                       ANH24220
  103 CONTINUE
      RETURN                                                            ANH24170
C
C     C=A*TB                                                            ANH24240
  203 CONTINUE
      DO 104 I=1,N                                                      ANH24250
      DO 104 J=1,N                                                      ANH24260
      DO 104 K=1,N                                                      ANH24270
      C(I,J)=C(I,J)+A(I,K)*B(J,K)                                       ANH24280
  104 CONTINUE
      RETURN                                                            ANH24290
C
C     C=TA*TB                                                           ANH24300
  204 CONTINUE
      DO 105 I=1,N                                                      ANH24310
      DO 105 J=1,N                                                      ANH24320
      DO 105 K=1,N                                                      ANH24330
      C(I,J)=C(I,J)+A(K,I)*B(J,K)                                       ANH24340
  105 CONTINUE
      RETURN                                                            ANH24350
C
C     C=TA*B*A                                                          ANH24360
  205 CONTINUE
      DO 106 I=1,N                                                      ANH24370
      DO 106 J=1,N                                                      ANH24380
      D(I,J)=A00                                                        ANH24390
      DO 106 K=1,N                                                      ANH24400
      D(I,J)=D(I,J)+A(K,I)*B(K,J)                                       ANH24410
  106 CONTINUE
      DO 107 I=1,N                                                      ANH24420
      DO 107 J=1,N                                                      ANH24430
      DO 107 K=1,N                                                      ANH24440
      C(I,J)=C(I,J)+D(I,K)*A(K,J)                                       ANH24450
  107 CONTINUE
      RETURN                                                            ANH24460
C
C     C=A*B*TA                                                          ANH24470
  206 CONTINUE
      DO 108 I=1,N                                                      ANH24480
      DO 108 J=1,N                                                      ANH24490
      D(I,J)=A00                                                        ANH24500
      DO 108 K=1,N                                                      ANH24510
      D(I,J)=D(I,J)+A(I,K)*B(K,J)                                       ANH24520
  108 CONTINUE
      DO 109 I=1,N                                                      ANH24530
      DO 109 J=1,N                                                      ANH24540
      DO 109 K=1,N                                                      ANH24550
      C(I,J)=C(I,J)+D(I,K)*A(J,K)                                       ANH24560
  109 CONTINUE
      RETURN                                                            ANH24570
      END                                                               ANH24580
      SUBROUTINE ANHOUT(A,NAD,NBD,M,N,IOUT)                             ANH24590
      IMPLICIT REAL*8 (A-H,O-Z)                                         ANH24600
      CHARACTER*6 LINE(21)                                              ANH24610
      DIMENSION A(NAD,NBD)                                              ANH24620
      DATA LINE / 21*'------' /                                         ANH24630
    1 FORMAT(2X,10(7X,I5))                                              ANH24640
    2 FORMAT(2X,21A6)                                                   ANH24650
    3 FORMAT(2X,I2,2X,10F12.3)                                          ANH24660
    4 FORMAT(/)                                                         ANH24670
C                                                                       ANH24680
      II=0                                                              ANH24690
      JJ=0                                                              ANH24700
  200 CONTINUE
      II=II+1                                                           ANH24710
      JJ=JJ+1                                                           ANH24720
      KK=10*JJ                                                          ANH24730
      NN=N                                                              ANH24740
      IF(N.GT.KK) NN=KK                                                 ANH24750
      LL=2*(NN-II+1)+1                                                  ANH24760
      WRITE (IOUT,1) (I,I=II,NN)                                        ANH24770
      WRITE (IOUT,2) (LINE(I),I=1,LL)                                   ANH24780
      DO 101 I=1,M                                                      ANH24790
      WRITE (IOUT,3) I,(A(I,J),J=II,NN)                                 ANH24800
  101 CONTINUE                                                          ANH24810
      IF(N.LE.KK) GO TO 201                                             ANH24820
      WRITE (IOUT,4)                                                    ANH24830
      II=KK                                                             ANH24840
      GO TO 200                                                         ANH24850
C
  201 CONTINUE
      RETURN                                                            ANH24860
      END                                                               ANH24870
      SUBROUTINE MX1OUT(A,NAD,M,IOUT)  
      IMPLICIT REAL*8 (A-H,O-Z)                                         ANH24890
      CHARACTER*6 LINE(7)                                               ANH24900
      DIMENSION A(NAD,NAD)
      DATA LINE / 7*'------' /                                          ANH24920
    1 FORMAT(3X,3(7X,I5))
    2 FORMAT(3X,7A6) 
    3 FORMAT(3X,I2,2X,3F12.7) 
    4 FORMAT(/)                                                         ANH24960
C                                                                       ANH24970
      WRITE (IOUT,1) (I,I=1,3)
      WRITE (IOUT,2) (LINE(I),I=1,7)
      DO 101 I=1,M                                                      ANH25000
      WRITE (IOUT,3) I,(A(I,J),J=1,3)
  101 CONTINUE                                                          ANH25020
c     WRITE (IOUT,4)                                                    ANH25030
C
      RETURN                                                            ANH25040
      END                                                               ANH25050
      SUBROUTINE MX2OUT(A,NAD,B,NBD,M,IOUT)                             ANH24880
      IMPLICIT REAL*8 (A-H,O-Z)                                         ANH24890
      CHARACTER*6 LINE(7)                                               ANH24900
      DIMENSION A(NAD,NAD),B(NBD,NBD)                                   ANH24910
      DATA LINE / 7*'------' /                                          ANH24920
    1 FORMAT(3X,3(7X,I5),7X,3(7X,I5))                                   ANH24930
    2 FORMAT(3X,7A6,1X,7A6,1X,7A6)                                      ANH24940
    3 FORMAT(3X,I2,2X,3F12.7,3X,I2,2X,3F12.7)                           ANH24950
    4 FORMAT(/)                                                         ANH24960
C                                                                       ANH24970
      WRITE (IOUT,1) (I,I=1,3),(I,I=1,3)                                ANH24980
      WRITE (IOUT,2) (LINE(I),I=1,7),(LINE(I),I=1,7)                    ANH24990
      DO 101 I=1,M                                                      ANH25000
      WRITE (IOUT,3) I,(A(I,J),J=1,3),I,(B(I,J),J=1,3)                  ANH25010
  101 CONTINUE                                                          ANH25020
      WRITE (IOUT,4)                                                    ANH25030
C
      RETURN                                                            ANH25040
      END                                                               ANH25050
      SUBROUTINE MX3OUT(A,NAD,B,NBD,C,NCD,M,IOUT)                       ANH25060
      IMPLICIT REAL*8 (A-H,O-Z)                                         ANH25070
      CHARACTER*6 LINE(7)                                               ANH25080
      DIMENSION A(NAD,NAD),B(NBD,NBD),C(NCD,NCD)                        ANH25090
      DATA LINE / 7*'------' /                                          ANH25100
    1 FORMAT(3X,3(7X,I5),7X,3(7X,I5),7X,3(7X,I5))                       ANH25110
    2 FORMAT(3X,7A6,1X,7A6,1X,7A6)                                      ANH25120
    3 FORMAT(3X,I2,2X,3F12.7,3X,I2,2X,3F12.7,3X,I2,2X,3F12.7)           ANH25130
    4 FORMAT(/)                                                         ANH25140
C                                                                       ANH25150
      WRITE (IOUT,1) (I,I=1,3),(I,I=1,3),(I,I=1,3)                      ANH25160
      WRITE (IOUT,2) (LINE(I),I=1,7),(LINE(I),I=1,7),(LINE(I),I=1,7)    ANH25170
      DO 101 I=1,M                                                      ANH25180
      WRITE (IOUT,3) I,(A(I,J),J=1,3),I,(B(I,J),J=1,3),I,(C(I,J),J=1,3) ANH25190
  101 CONTINUE                                                          ANH25200
      WRITE (IOUT,4)                                                    ANH25210
C
      RETURN                                                            ANH25220
      END                                                               ANH25230
      SUBROUTINE FRQOUT(A,B,NAD,NBD,M,N,IOUT)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*6 LINE
      DIMENSION A(NAD,NBD),B(NBD)
      COMMON/VIB112/LINEWD
c     The number of data per line
      integer ndpl
    1 FORMAT(2X,10(7X,I5))
    2 FORMAT(6X,10(4X,F8.2))
    3 FORMAT(2X,21A6)
    4 FORMAT(2X,I2,2X,10F12.7)
    5 FORMAT(/)
C
      DATA LINE / '------' /
C
      ndpl = (linewd-7)/12
      if ((ndpl.le.0).or.(ndpl.gt.200)) then
        write(6,*) 'ERROR: bogus ndpl in matout: check ''LINEWIDTH'''
        stop
      endif
C
      II=0
      JJ=0
  200 CONTINUE
      II=II+1
      JJ=JJ+1
      KK=ndpl*JJ
      NN=N
      IF(N.GT.KK) NN=KK
      LL=2*(NN-II+1)+1
      WRITE(IOUT,1) (I,I=II,NN)
      WRITE(IOUT,2) (B(I),I=II,NN)
      WRITE(IOUT,3) (LINE,I=1,LL)
      DO 101 I=1,M
      WRITE(IOUT,4) I,(A(I,J),J=II,NN)
  101 CONTINUE
      IF(N.LE.KK) GO TO 201
      WRITE(IOUT,5)
      II=KK
      GO TO 200
C
  201 CONTINUE
      RETURN
      END
      SUBROUTINE RSP(NM,N,NV,A,W,MATZ,Z,FV1,FV2)                        
C                                                                       
      IMPLICIT REAL*8 (A-H,O-Z)                                         
C                                                                       
C                                                                       
      INTEGER I,J,N,NM,NV,IERR,MATZ                                     
      REAL*8 A(NV),W(N),Z(NM,N),FV1(N),FV2(N)                           
      DATA A00,ONE / 0.0D+00,1.0D+00 /                                 
    1 FORMAT(//,2X,' IERR = ',I5//)                                     
C                                                                       
C     THIS SUBROUTINE CALLS THE RECOMMENDED SEQUENCE OF                 
C     SUBROUTINES FROM THE EIGENSYSTEM SUBROUTINE PACKAGE (EISPACK)     
C     TO FIND THE EIGENVALUES AND EIGENVECTORS (IF DESIRED)             
C     OF A REAL SYMMETRIC PACKED MATRIX.                                
C                                                                       
C     ON INPUT-                                                         
C                                                                       
C        NM  MUST BE SET TO THE ROW DIMENSION OF THE TWO-DIMENSIONAL    
C        ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM            
C        DIMENSION STATEMENT,                                           
C                                                                       
C        N  IS THE ORDER OF THE MATRIX  A,                              
C                                                                       
C        NV  IS AN INTEGER VARIABLE SET EQUAL TO THE                    
C        DIMENSION OF THE ARRAY  A  AS SPECIFIED FOR                    
C        A  IN THE CALLING PROGRAM.  NV  MUST NOT BE                    
C        LESS THAN  N*(N+1)/2,                                          
C                                                                       
C        A  CONTAINS THE LOWER TRIANGLE OF THE REAL SYMMETRIC           
C        PACKED MATRIX STORED ROW-WISE,                                 
C                                                                       
C        MATZ  IS AN INTEGER VARIABLE SET EQUAL TO ZERO IF              
C        ONLY EIGENVALUES ARE DESIRED,  OTHERWISE IT IS SET TO          
C        ANY NON-ZERO INTEGER FOR BOTH EIGENVALUES AND EIGENVECTORS.    
C                                                                       
C     ON OUTPUT-                                                        
C                                                                       
C        W  CONTAINS THE EIGENVALUES IN ASCENDING ORDER,                
C                                                                       
C        Z  CONTAINS THE EIGENVECTORS IF MATZ IS NOT ZERO,              
C                                                                       
C        IERR  IS AN INTEGER OUTPUT VARIABLE SET EQUAL TO AN            
C        ERROR COMPLETION CODE DESCRIBED IN SECTION 2B OF THE           
C        DOCUMENTATION.  THE NORMAL COMPLETION CODE IS ZERO,            
C                                                                       
C        FV1  AND  FV2  ARE TEMPORARY STORAGE ARRAYS.                   
C                                                                       
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,        
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY         
C                                                                       
C     ------------------------------------------------------------------
C                                                                       
      IF (N .LE. NM) GO TO 5                                            
      IERR = 10 * N                                                     
      GO TO 50                                                          
    5 IF (NV .GE. (N * (N + 1)) / 2) GO TO 10                           
      IERR = 20 * N                                                     
      GO TO 50                                                          
C                                                                       
   10 CALL  TRED3(N,NV,A,W,FV1,FV2)                                     
      IF (MATZ .NE. 0) GO TO 20                                         
C     ********** FIND EIGENVALUES ONLY **********                       
      CALL  TQLRAT(N,W,FV2,IERR)                                        
      GO TO 50                                                          
C     ********** FIND BOTH EIGENVALUES AND EIGENVECTORS **********      
   20 DO 40 I = 1, N                                                    
C                                                                       
           DO 30 J = 1, N                                               
                Z(J,I) =A00
   30      CONTINUE                                                     
C                                                                       
           Z(I,I) =ONE                                                  
   40 CONTINUE                                                          
C                                                                       
      CALL  TQL2(NM,N,W,FV1,Z,IERR)                                     
      IF (IERR .NE. 0) GO TO 50                                         
      CALL  TRBAK3(NM,N,NV,A,N,Z)                                       
   50 IF(IERR.NE.0) GO TO 60                                            
      RETURN                                                            
   60 WRITE(6,1) IERR                                                   
      RETURN                                                            
C     ********** LAST CARD OF RSP **********                            
      END                                                               
C                                                                       
C     ------------------------------------------------------------------
C                                                                       
      SUBROUTINE TRED3(N,NV,A,D,E,E2)                                   
C                                                                       
      IMPLICIT REAL*8 (A-H,O-Z)                                         
C                                                                       
C                                                                       
      INTEGER I,J,K,L,N,II,IZ,JK,NV                                     
      REAL*8 A(NV),D(N),E(N),E2(N)                                      
      REAL*8 F,G,H,HH,SCALE                                             
      DATA A00 /0.0D+00/      
C     REAL*8 SQRT,ABS,SIGN                                              
C                                                                       
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TRED3,    
C     NUM. MATH. 11, 181-195(1968) BY MARTIN, REINSCH, AND WILKINSON.   
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).   
C                                                                       
C     THIS SUBROUTINE REDUCES A REAL SYMMETRIC MATRIX, STORED AS        
C     A ONE-DIMENSIONAL ARRAY, TO A SYMMETRIC TRIDIAGONAL MATRIX        
C     USING ORTHOGONAL SIMILARITY TRANSFORMATIONS.                      
C                                                                       
C     ON INPUT-                                                         
C                                                                       
C        N IS THE ORDER OF THE MATRIX,                                  
C                                                                       
C        NV MUST BE SET TO THE DIMENSION OF THE ARRAY PARAMETER A       
C          AS DECLARED IN THE CALLING PROGRAM DIMENSION STATEMENT,      
C                                                                       
C        A CONTAINS THE LOWER TRIANGLE OF THE REAL SYMMETRIC            
C          INPUT MATRIX, STORED ROW-WISE AS A ONE-DIMENSIONAL           
C          ARRAY, IN ITS FIRST N*(N+1)/2 POSITIONS.                     
C                                                                       
C     ON OUTPUT-                                                        
C                                                                       
C        A CONTAINS INFORMATION ABOUT THE ORTHOGONAL                    
C          TRANSFORMATIONS USED IN THE REDUCTION,                       
C                                                                       
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE TRIDIAGONAL MATRIX,    
C                                                                       
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE TRIDIAGONAL         
C          MATRIX IN ITS LAST N-1 POSITIONS.  E(1) IS SET TO ZERO,      
C                                                                       
C        E2 CONTAINS THE SQUARES OF THE CORRESPONDING ELEMENTS OF E.    
C          E2 MAY COINCIDE WITH E IF THE SQUARES ARE NOT NEEDED.        
C                                                                       
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,        
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY         
C                                                                       
C     ------------------------------------------------------------------
C                                                                       
C     ********** FOR I=N STEP -1 UNTIL 1 DO -- **********               
      DO  300 II = 1, N                                                 
           I = N + 1 - II                                               
           L = I - 1                                                    
           IZ = (I * L) / 2                                             
           H =A00
           SCALE = A00
           IF (L .LT. 1) GO TO 130                                      
C     ********** SCALE ROW (ALGOL TOL THEN NOT NEEDED) **********       
           DO 120 K = 1, L                                              
                IZ = IZ + 1                                             
                D(K) = A(IZ)                                            
                SCALE = SCALE + DABS(D(K))                              
  120      CONTINUE                                                     
C                                                                       
           IF (SCALE .NE.A00) GO TO 140         
  130      E(I) =A00
           E2(I) = A00
           GO TO 290                                                    
C                                                                       
  140      DO 150 K = 1, L                                              
                D(K) = D(K) / SCALE                                     
                H = H + D(K) * D(K)                                     
  150      CONTINUE                                                     
C                                                                       
           E2(I) = SCALE * SCALE * H                                    
           F = D(L)                                                     
           G = -DSIGN(DSQRT(H),F)                                       
           E(I) = SCALE * G                                             
           H = H - F * G                                                
           D(L) = F - G                                                 
           A(IZ) = SCALE * D(L)                                         
           IF (L .EQ. 1) GO TO 290                                      
           F =A00
C                                                                       
           DO 240 J = 1, L                                              
                G =A00
                JK = (J * (J-1)) / 2                                    
C     ********** FORM ELEMENT OF A*U **********                         
                DO 180 K = 1, L                                         
                     JK = JK + 1                                        
                     IF (K .GT. J) JK = JK + K - 2                      
                     G = G + A(JK) * D(K)                               
  180           CONTINUE                                                
C     ********** FORM ELEMENT OF P **********                           
                E(J) = G / H                                            
                F = F + E(J) * D(J)                                     
  240      CONTINUE                                                     
C                                                                       
           HH = F / (H + H)                                             
           JK = 0                                                       
C     ********** FORM REDUCED A **********                              
           DO 260 J = 1, L                                              
                F = D(J)                                                
                G = E(J) - HH * F                                       
                E(J) = G                                                
C                                                                       
                DO 260 K = 1, J                                         
                     JK = JK + 1                                        
                     A(JK) = A(JK) - F * E(K) - G * D(K)                
  260      CONTINUE                                                     
C                                                                       
  290      D(I) = A(IZ+1)                                               
           A(IZ+1) = SCALE * DSQRT(H)                                   
  300 CONTINUE                                                          
C                                                                       
      RETURN                                                            
C     ********** LAST CARD OF TRED3 **********                          
      END                                                               
C                                                                       
C     ------------------------------------------------------------------
C                                                                       
      SUBROUTINE TQLRAT(N,D,E2,IERR)                                    
C                                                                       
      IMPLICIT REAL*8 (A-H,O-Z)                                         
C                                                                       
C                                                                       
      INTEGER I,J,L,M,N,II,L1,MML,IERR                                  
      REAL*8 D(N),E2(N)                                                 
      REAL*8 B,C,F,G,H,P,R,S,MACHEP                                     
      DATA A00,ONE,TWO /0.0D+00,1.0D+00,2.0D+00/      
C     REAL*8 SQRT,ABS,SIGN                                              
C                                                                       
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TQLRAT,   
C     ALGORITHM 464, COMM. ACM 16, 689(1973) BY REINSCH.                
C                                                                       
C     THIS SUBROUTINE FINDS THE EIGENVALUES OF A SYMMETRIC              
C     TRIDIAGONAL MATRIX BY THE RATIONAL QL METHOD.                     
C                                                                       
C     ON INPUT-                                                         
C                                                                       
C        N IS THE ORDER OF THE MATRIX,                                  
C                                                                       
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX,          
C                                                                       
C        E2 CONTAINS THE SQUARES OF THE SUBDIAGONAL ELEMENTS OF THE     
C          INPUT MATRIX IN ITS LAST N-1 POSITIONS.  E2(1) IS ARBITRARY. 
C                                                                       
C      ON OUTPUT-                                                       
C                                                                       
C        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN          
C          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT AND          
C          ORDERED FOR INDICES 1,2,...IERR-1, BUT MAY NOT BE            
C          THE SMALLEST EIGENVALUES,                                    
C                                                                       
C        E2 HAS BEEN DESTROYED,                                         
C                                                                       
C        IERR IS SET TO                                                 
C          ZERO       FOR NORMAL RETURN,                                
C          J          IF THE J-TH EIGENVALUE HAS NOT BEEN               
C                     DETERMINED AFTER 30 ITERATIONS.                   
C                                                                       
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,        
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY         
C                                                                       
C     ------------------------------------------------------------------
C                                                                       
C     ********** MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING     
C                THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.   
C                                                                       
C                **********                                             
      MACHEP = TWO**(-47)                                               
C                                                                       
      IERR = 0                                                          
      IF (N .EQ. 1) GO TO 1001                                          
C                                                                       
      DO 100 I = 2, N                                                   
  100 E2(I-1) = E2(I)                                                   
C                                                                       
      F =A00
      B =A00
      E2(N) =A00
C                                                                       
      DO 290 L = 1, N                                                   
           J = 0                                                        
           H = MACHEP * (DABS(D(L)) + DSQRT(E2(L)))                     
           IF (B .GT. H) GO TO 105                                      
           B = H                                                        
           C = B * B                                                    
C     ********** LOOK FOR SMALL SQUARED SUB-DIAGONAL ELEMENT ********** 
  105      DO 110 M = L, N                                              
                IF (E2(M) .LE. C) GO TO 120                             
C     ********** E2(N) IS ALWAYS ZERO, SO THERE IS NO EXIT              
C                THROUGH THE BOTTOM OF THE LOOP **********              
  110      CONTINUE                                                     
C                                                                       
  120      IF (M .EQ. L) GO TO 210                                      
  130      IF (J .EQ. 30) GO TO 1000                                    
           J = J + 1                                                    
C     ********** FORM SHIFT **********                                  
           L1 = L + 1                                                   
           S = DSQRT(E2(L))                                             
           G = D(L)                                                     
           P = (D(L1) - G) / (TWO * S)                                  
           R = DSQRT(P*P+ONE)                                           
           D(L) = S / (P + DSIGN(R,P))                                  
           H = G - D(L)                                                 
C                                                                       
           DO 140 I = L1, N                                             
  140      D(I) = D(I) - H                                              
C                                                                       
           F = F + H                                                    
C     ********** RATIONAL QL TRANSFORMATION **********                  
           G = D(M)                                                     
           IF (G .EQ.A00) G = B                                        
           H = G                                                        
           S =A00
           MML = M - L                                                  
C     ********** FOR I=M-1 STEP -1 UNTIL L DO -- **********             
           DO 200 II = 1, MML                                           
                I = M - II                                              
                P = G * H                                               
                R = P + E2(I)                                           
                E2(I+1) = S * R                                         
                S = E2(I) / R                                           
                D(I+1) = H + S * (H + D(I))                             
                G = D(I) - E2(I) / G                                    
                IF (G .EQ.A00) G = B                                   
                H = G * P / R                                           
  200      CONTINUE                                                     
C                                                                       
           E2(L) = S * G                                                
           D(L) = H                                                     
C     ********** GUARD AGAINST UNDERFLOW IN CONVERGENCE TEST ********** 
           IF (H .EQ.A00) GO TO 210                                    
           IF (DABS(E2(L)) .LE. DABS(C/H)) GO TO 210                    
           E2(L) = H * E2(L)                                            
           IF (E2(L) .NE.A00) GO TO 130                                
  210      P = D(L) + F                                                 
C     ********** ORDER EIGENVALUES **********                           
           IF (L .EQ. 1) GO TO 250                                      
C     ********** FOR I=L STEP -1 UNTIL 2 DO -- **********               
           DO 230 II = 2, L                                             
                I = L + 2 - II                                          
                IF (P .GE. D(I-1)) GO TO 270                            
                D(I) = D(I-1)                                           
  230      CONTINUE                                                     
C                                                                       
  250      I = 1                                                        
  270      D(I) = P                                                     
  290 CONTINUE                                                          
C                                                                       
      GO TO 1001                                                        
C     ********** SET ERROR -- NO CONVERGENCE TO AN                      
C                EIGENVALUE AFTER 30 ITERATIONS **********              
 1000 IERR = L                                                          
 1001 RETURN                                                            
C     ********** LAST CARD OF TQLRAT **********                         
      END                                                               
C                                                                       
C     ------------------------------------------------------------------
C                                                                       
      SUBROUTINE TQL2(NM,N,D,E,Z,IERR)                                  
C                                                                       
      IMPLICIT REAL*8 (A-H,O-Z)                                         
C                                                                       
C                                                                       
      INTEGER I,J,K,L,M,N,II,L1,NM,MML,IERR                             
      REAL*8 D(N),E(N),Z(NM,N)                                          
      REAL*8 B,C,F,G,H,P,R,S,MACHEP                                     
      DATA A00,ONE,TWO /0.0D+00,1.0D+00,2.0D+00/                       
C     REAL*8 SQRT,ABS,SIGN                                              
C                                                                       
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TQL2,     
C     NUM. MATH. 11, 293-306(1968) BY BOWDLER, MARTIN, REINSCH, AND     
C     WILKINSON.                                                        
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 227-240(1971).   
C                                                                       
C     THIS SUBROUTINE FINDS THE EIGENVALUES AND EIGENVECTORS            
C     OF A SYMMETRIC TRIDIAGONAL MATRIX BY THE QL METHOD.               
C     THE EIGENVECTORS OF A FULL SYMMETRIC MATRIX CAN ALSO              
C     BE FOUND IF  TRED2  HAS BEEN USED TO REDUCE THIS                  
C     FULL MATRIX TO TRIDIAGONAL FORM.                                  
C                                                                       
C     ON INPUT-                                                         
C                                                                       
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL         
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM          
C          DIMENSION STATEMENT,                                         
C                                                                       
C        N IS THE ORDER OF THE MATRIX,                                  
C                                                                       
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX,          
C                                                                       
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX        
C          IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY,               
C                                                                       
C        Z CONTAINS THE TRANSFORMATION MATRIX PRODUCED IN THE           
C          REDUCTION BY  TRED2, IF PERFORMED.  IF THE EIGENVECTORS      
C          OF THE TRIDIAGONAL MATRIX ARE DESIRED, Z MUST CONTAIN        
C          THE IDENTITY MATRIX.                                         
C                                                                       
C      ON OUTPUT-                                                       
C                                                                       
C        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN          
C          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT BUT          
C          UNORDERED FOR INDICES 1,2,...,IERR-1,                        
C                                                                       
C        E HAS BEEN DESTROYED,                                          
C                                                                       
C        Z CONTAINS ORTHONORMAL EIGENVECTORS OF THE SYMMETRIC           
C          TRIDIAGONAL (OR FULL) MATRIX.  IF AN ERROR EXIT IS MADE,     
C          Z CONTAINS THE EIGENVECTORS ASSOCIATED WITH THE STORED       
C          EIGENVALUES,                                                 
C                                                                       
C        IERR IS SET TO                                                 
C          ZERO       FOR NORMAL RETURN,                                
C          J          IF THE J-TH EIGENVALUE HAS NOT BEEN               
C                     DETERMINED AFTER 30 ITERATIONS.                   
C                                                                       
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,        
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY         
C                                                                       
C     ------------------------------------------------------------------
C                                                                       
C     ********** MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING     
C                THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.   
C                                                                       
C                **********                                             
      MACHEP =TWO**(-47)                                                
C                                                                       
      IERR = 0                                                          
      IF (N .EQ. 1) GO TO 1001                                          
C                                                                       
      DO 100 I = 2, N                                                   
  100 E(I-1) = E(I)                                                     
C                                                                       
      F =A00                                                           
      B =A00
      E(N) = A00
C                                                                       
      DO 240 L = 1, N                                                   
           J = 0                                                        
           H = MACHEP * (DABS(D(L)) + DABS(E(L)))                       
           IF (B .LT. H) B = H                                          
C     ********** LOOK FOR SMALL SUB-DIAGONAL ELEMENT **********         
           DO 110 M = L, N                                              
                IF (DABS(E(M)) .LE. B) GO TO 120                        
C     ********** E(N) IS ALWAYS ZERO, SO THERE IS NO EXIT               
C                THROUGH THE BOTTOM OF THE LOOP **********              
  110      CONTINUE                                                     
C                                                                       
  120      IF (M .EQ. L) GO TO 220                                      
  130      IF (J .EQ. 30) GO TO 1000                                    
           J = J + 1                                                    
C     ********** FORM SHIFT **********                                  
           L1 = L + 1                                                   
           G = D(L)                                                     
           P = (D(L1) - G) / (TWO * E(L))                               
           R = DSQRT(P*P+ONE)                                           
           D(L) = E(L) / (P + DSIGN(R,P))                               
           H = G - D(L)                                                 
C                                                                       
           DO 140 I = L1, N                                             
  140      D(I) = D(I) - H                                              
C                                                                       
           F = F + H                                                    
C     ********** QL TRANSFORMATION **********                           
           P = D(M)                                                     
           C = ONE                                                      
           S = A00
           MML = M - L                                                  
C     ********** FOR I=M-1 STEP -1 UNTIL L DO -- **********             
           DO 200 II = 1, MML                                           
                I = M - II                                              
                G = C * E(I)                                            
                H = C * P                                               
                IF (DABS(P) .LT. DABS(E(I))) GO TO 150                  
                C = E(I) / P                                            
                R = DSQRT(C*C+ONE)                                      
                E(I+1) = S * P * R                                      
                S = C / R                                               
                C = ONE / R                                             
                GO TO 160                                               
  150           C = P / E(I)                                            
                R = DSQRT(C*C+ONE)                                      
                E(I+1) = S * E(I) * R                                   
                S = ONE / R                                             
                C = C * S                                               
  160           P = C * D(I) - S * G                                    
                D(I+1) = H + S * (C * G + S * D(I))                     
C     ********** FORM VECTOR **********                                 
                DO 180 K = 1, N                                         
                     H = Z(K,I+1)                                       
                     Z(K,I+1) = S * Z(K,I) + C * H                      
                     Z(K,I) = C * Z(K,I) - S * H                        
  180           CONTINUE                                                
C                                                                       
  200      CONTINUE                                                     
C                                                                       
           E(L) = S * P                                                 
           D(L) = C * P                                                 
           IF (DABS(E(L)) .GT. B) GO TO 130                             
  220      D(L) = D(L) + F                                              
  240 CONTINUE                                                          
C     ********** ORDER EIGENVALUES AND EIGENVECTORS **********          
      DO 300 II = 2, N                                                  
           I = II - 1                                                   
           K = I                                                        
           P = D(I)                                                     
C                                                                       
           DO 260 J = II, N                                             
                IF (D(J) .GE. P) GO TO 260                              
                K = J                                                   
                P = D(J)                                                
  260      CONTINUE                                                     
C                                                                       
           IF (K .EQ. I) GO TO 300                                      
           D(K) = D(I)                                                  
           D(I) = P                                                     
C                                                                       
           DO 280 J = 1, N                                              
                P = Z(J,I)                                              
                Z(J,I) = Z(J,K)                                         
                Z(J,K) = P                                              
  280      CONTINUE                                                     
C                                                                       
  300 CONTINUE                                                          
C                                                                       
      GO TO 1001                                                        
C     ********** SET ERROR -- NO CONVERGENCE TO AN                      
C                EIGENVALUE AFTER 30 ITERATIONS **********              
 1000 IERR = L                                                          
 1001 RETURN                                                            
C     ********** LAST CARD OF TQL2 **********                           
      END                                                               
C                                                                       
C     ------------------------------------------------------------------
C                                                                       
      SUBROUTINE TRBAK3(NM,N,NV,A,M,Z)                                  
C                                                                       
      IMPLICIT REAL*8 (A-H,O-Z)                                         
C                                                                       
C                                                                       
      INTEGER I,J,K,L,M,N,IK,IZ,NM,NV                                   
      REAL*8 A(NV),Z(NM,M)                                              
      REAL*8 H,S                                                        
      DATA A00 /0.0D+00/                                               
C                                                                       
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TRBAK3,   
C     NUM. MATH. 11, 181-195(1968) BY MARTIN, REINSCH, AND WILKINSON.   
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).   
C                                                                       
C     THIS SUBROUTINE FORMS THE EIGENVECTORS OF A REAL SYMMETRIC        
C     MATRIX BY BACK TRANSFORMING THOSE OF THE CORRESPONDING            
C     SYMMETRIC TRIDIAGONAL MATRIX DETERMINED BY  TRED3.                
C                                                                       
C     ON INPUT-                                                         
C                                                                       
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL         
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM          
C          DIMENSION STATEMENT,                                         
C                                                                       
C        N IS THE ORDER OF THE MATRIX,                                  
C                                                                       
C        NV MUST BE SET TO THE DIMENSION OF THE ARRAY PARAMETER A       
C          AS DECLARED IN THE CALLING PROGRAM DIMENSION STATEMENT,      
C                                                                       
C        A CONTAINS INFORMATION ABOUT THE ORTHOGONAL TRANSFORMATIONS    
C          USED IN THE REDUCTION BY  TRED3  IN ITS FIRST                
C          N*(N+1)/2 POSITIONS,                                         
C                                                                       
C        M IS THE NUMBER OF EIGENVECTORS TO BE BACK TRANSFORMED,        
C                                                                       
C        Z CONTAINS THE EIGENVECTORS TO BE BACK TRANSFORMED             
C          IN ITS FIRST M COLUMNS.                                      
C                                                                       
C     ON OUTPUT-                                                        
C                                                                       
C        Z CONTAINS THE TRANSFORMED EIGENVECTORS                        
C          IN ITS FIRST M COLUMNS.                                      
C                                                                       
C     NOTE THAT TRBAK3 PRESERVES VECTOR EUCLIDEAN NORMS.                
C                                                                       
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,        
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY         
C                                                                       
C     ------------------------------------------------------------------
C                                                                       
      IF (M .EQ. 0) GO TO 200                                           
      IF (N .EQ. 1) GO TO 200                                           
C                                                                       
      DO 140 I = 2, N                                                   
           L = I - 1                                                    
           IZ = (I * L) / 2                                             
           IK = IZ + I                                                  
           H = A(IK)                                                    
           IF (H .EQ.A00) GO TO 140                                    
C                                                                       
           DO 130 J = 1, M                                              
                S =A00
                IK = IZ                                                 
C                                                                       
                DO 110 K = 1, L                                         
                     IK = IK + 1                                        
                     S = S + A(IK) * Z(K,J)                             
  110           CONTINUE                                                
C     ********** DOUBLE DIVISION AVOIDS POSSIBLE UNDERFLOW **********   
                S = (S / H) / H                                         
                IK = IZ                                                 
C                                                                       
                DO 120 K = 1, L                                         
                     IK = IK + 1                                        
                     Z(K,J) = Z(K,J) - S * A(IK)                        
  120           CONTINUE                                                
C                                                                       
  130      CONTINUE                                                     
C                                                                       
  140 CONTINUE                                                          
C                                                                       
  200 RETURN                                                            
C     ********** LAST CARD OF TRBAK3 **********                         
      END                                                               
