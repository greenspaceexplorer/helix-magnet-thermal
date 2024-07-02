C        OPTIONS/CHECK=ALL/EXTEND_SOURCE
        SUBROUTINE MOMENTM    
C   
c.        Momentum fit using the quintic spline model.  
c.            documentation = nucl. instr. meth.115,    
c.            431 (1974)  (h.wind), report on omega 
c.            experience dd 75.2 (d.townsend, j.wilson) 
c.            origin =  cern library x510,  
c.            amended d.townsend(1974)  
c.            status = cleaned-up from field proven code    
C   
C         INPUT PARAMETERS  
C             NXM             NUMBER OF MEASURED POINTS 
C             XM,YM,ZM        MEASURED POINT COORDINATES IN CM
C             WEIGHT          TO BE USED IN FIT ( 1/SIG^2) constructed from
C                              DXM, DYM, DZM sigmas
C             INSERT          NUMBER OF MAG.FIELD POINTS TO BE INSERTED 
C                             (LAST ELEMENT NOT USED, BUT SET = 0)  
C
        IMPLICIT NONE
C Include machine limits
        INCLUDE 'machine.inc'

        REAL*8 DET
        REAL*8 T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T31,T13,T33
        REAL*8 SX,SY,SC,SZ,SXX,SXY,SXC,SXZ,SYC,SCC,SCZ
        REAL*8 SZCZ,SXCZ,SCCZ,SCCSCZ,SYCZCZ

        INTEGER NXM,INSERT
        REAL*8 XM,YM,ZM,DXM,DYM,DZM,WEIGHT(99)
        COMMON /MOM_INPUT/XM(99),YM(99),ZM(99),DXM(99),
     +                    DYM(99),DZM(99),NXM,INSERT(99)
   
C         OUPUT PARAMETERS  
C             PMAG            FITTED MOMENTUM (APART FROM CLIGHT FACTOR) IN MV
C             PVEC(3)         COMPONENTS OF PMAG    
C             A1,A2,D,B1,B2   THE FITTED PARAMETERS 
C             ERR(5,5)        ERROR MATRIX FOR THESE    
C             IERR            =0 IF O.K., ELSE SETUP ERROR  

        INTEGER NXOUT,IERR
        REAL*8 PMAG,PVECIN,PVECOUT,A1,A2,D,B1,B2,ERR,XOUT,YOUT,ZOUT
        REAL*8 DEVX,DEVY,DEVZ,ERRX,ERRY,ERRZ
        COMMON /MOM_OUTPUT/PMAG,PVECIN(3),PVECOUT(3),
     +                     A1,A2,D,B1,B2,ERR(5,5),
     +                     XOUT(100),YOUT(100),ZOUT(100),
     +           DEVX,DEVY,DEVZ,ERRX(100),ERRY(100),ERRZ(100),
     +           NXOUT,IERR
   
C         LOCAL COMMUNICATION   
C             NTOT            TOTAL NUMBER OF POINTS WITH MAG. FIELD    
C                             ( = NXM + SUM(INSERT(I))   )  
C             X,Y,Z           COORDINATES OF ALL POINTS IN FIT SYSTEM in CM
C             YI,ZI           DY/DX AND DZ/DX ALL POINTS    
C             BX,BY,BZ        FIELD COMPONENTS OF ALL POINTS,FIT SYSTEM IN GAUSS
C             BXM,BYM,BZM     FIELD COMPONENTS,MEASURED POINTS, 
C                             MEASURED SYSTEM  in GAUSS
C                             (MAY BE MADE AN INPUT PARAMETER)  
C             W11,W12...W33   ORTHONORMAL ROTATION MATRIX BETWEEN   
C                             MEASURED AND FIT COORD. SYSTEMS   
C             PYII,PZII       SECOND DERIVATIVES ALL POINTS 
C             YII,ZII         X AND Y INTEGRATED, FROM DTSPL    
C
        REAL*8 B(3)

        INTEGER NTOT
        REAL*8 X,Y,Z,YI,ZI
        COMMON /LOCAL1/ X(100),Y(100),Z(100),YI(100),ZI(100), NTOT 
        REAL*8 BXM,BYM,BZM
        COMMON /LOCAL2/BXM(100),BYM(100),BZM(100) 
        REAL*8 W11,W12,W13,W21,W22,W23,W31,W32,W33    
        COMMON /LOCAL3/W11,W12,W13,W21,W22,W23,W31,W32,W33    
        REAL*8 BX,BY,BZ
        COMMON /LOCAL4/BX(100),BY(100),BZ(100)    
        REAL*8 PYII,PZII,YII,ZII
        COMMON /LOCAL5/PYII(100),PZII(100),YII(100),ZII(100)  
C Keep track of measured point locations
        INTEGER IMEAS(100)
   
c Subroutines called    

C  CALL GETB(X,Y,Z,BX,BY,BZ,IERROR);USER SUPPLIED FOR FIELD 

C  CALL SPLFIT     INTERNAL FOR INTERPOLATION    
C  CALL DTSPL      INTERNAL FOR INTEGRATION  
   
c Local variable definitions
        INTEGER NX,ITER,J,KL,KK,IERROR,INCR,I,NTO,NTOTX,ITEM,K
        INTEGER NITERX,IJK
        REAL*8 SQR,DX,DY,DZ,A,D1,CLIGHT,XX,YY,ZZ,PV1,PV2,PV3,XERR,WX,WY,WZ
        REAL*8 SQX,SQY,SQZ,YERR,ZERR,XTMP,YTMP,ZTMP

C Local parameters for iterating and dimensioning   
        DATA NITERX/3/,NTOTX/100/,CLIGHT/0.299E-3/  !CLIGHT is for cm, gauss and MV        
        
C We define an inproduct.   
        REAL*8 AIN,XP,YP,ZP,AP,BP,CP
        AIN(XP,YP,ZP,AP,BP,CP) = XP*AP+YP*BP+ZP*CP    
   
        IERR = 0  

C NX is the number of measured points...
        NX = NXM  

C We need at least three points 
        IF (NX.LT.3)    GO TO 140 
C and no more than ntotx    
        IF (NX.GT.NTOTX)    GO TO 140 
C Also check ntot   
        NTO = NX  

        DO I = 2,NX    
         NTO = NTO+INSERT(I-1) 
        ENDDO

        NTOT = NTO    

cc        print *,' MOMENTM: NX:',NX,' NTOT:',NTO

        IF (NTO.GT.NTOTX)    GO TO 140    

C Get the weights to use...~1/sigma**2
        DO I = 1, NX
           WEIGHT(I) = 1.0/(DXM(I)**2+DYM(I)**2+DZM(I)**2)
        ENDDO
   
C We now change the coordinate system. in the new system, we may    
c expect x to be increasing along the track, y' to be small,    
c and z' to be very small.  
c The new coordinate system will have the x-axis parallel to the    
c vector joining the first and the last point on the track. 
c The y-axis will have the direction of the vector joining the  
c line segment between the first and the last point on the track    
c with a point far from this segment. The z-axis will be    
c orthogonal to the x- and y-axes, and the system will be right-    
c handed.   
   
C Find the new x-axis.  
        W11 = XM(NX)-XM(1)    
        W21 = YM(NX)-YM(1)    
        W31 = ZM(NX)-ZM(1)    
        D = 1.0/SQRT(W11**2+W21**2+W31**2)    
        W11 = W11*D   
        W21 = W21*D   
        W31 = W31*D   

C Find the new y-axis   
        INCR = NX/7+1 
        D = 0.0   
        DO 20 I = 2,NX,INCR   
              DX = XM(I)-XM(1)  
              DY = YM(I)-YM(1)  
              DZ = ZM(I)-ZM(1)  
              A = AIN(W11,W21,W31,DX,DY,DZ) 
C Here we use x,y,z as temporary storage space. 
              X(I) = DX-A*W11   
              Y(I) = DY-A*W21   
              Z(I) = DZ-A*W31   
              D1 = X(I)**2+Y(I)**2+Z(I)**2  
              IF (D1.LT.D)    GO TO 20  
              ITEM = I  
              D = D1    
20      CONTINUE  
        D = 1.0/SQRT(D)   
        W12 = X(ITEM)*D   
        W22 = Y(ITEM)*D   
        W32 = Z(ITEM)*D   
   
C Find the new z-axis   

        W13 = W21*W32-W22*W31 
        W23 = W31*W12-W11*W32 
        W33 = W11*W22-W12*W21 

C We now have a matrix to convert from the measured coordinate system
C to the working system:
C
C            | W11 W21 W31 |
C    X' =    | W12 W23 W32 | * X
C            | W13 W23 W33 |
C
C  OR  to convert back to the measured system:
C
C            | W11 W12 W13 |
C    X =     | W21 W23 W23 | * X'
C            | W31 W32 W33 |

cc        print *, 'Momentm: Coord axis conversion matrix elements: ', W13, W23, W33, W12, W22, W32, W11, W21, W31

C Compute new coordinates x,y,z, new field-components bx,by,bz  
   
        INSERT(NX) = 0    
        K = 1 
        DO 30 I = 1,NX    
          IMEAS(I) = K
cc          print *,'Old coords: ', K, xm(i), ym(i), zm(i)
          X(K) = AIN(W11,W21,W31,XM(I),YM(I),ZM(I)) 
          Y(K) = AIN(W12,W22,W32,XM(I),YM(I),ZM(I)) 
          Z(K) = AIN(W13,W23,W33,XM(I),YM(I),ZM(I)) 
cc          print *,'New coords: ', K, x(k), y(k), z(k)
          CALL GETB(XM(I),YM(I),ZM(I),B(1),B(2),B(3),IERROR)
c            print *,'MOMENTM: getb returns',B(1),B(2),B(3)
          IF (IERROR.NE.0) THEN
            PRINT *,' Error in MOMENTM: GETB IERROR = ',IERROR
          ENDIF
c          print *, 'MOMENTM: calling AIN for BX(K)'
          BX(K) = AIN(W11,W21,W31,B(1),B(2),B(3)) 
          BY(K) = AIN(W12,W22,W32,B(1),B(2),B(3)) 
          BZ(K) = AIN(W13,W23,W33,B(1),B(2),B(3)) 
cc          print *,'New bfield: ', K, bx(k), by(k), bz(k)
C Increasing order?
          IF (I.NE.1.AND.X(K).LT.X(KK))    GO TO 150    
          KK = K    
          K = K+INSERT(I)+1 
cc          print *, 'K, KK, I, INSERT(I):', K, KK, I, INSERT(I)
   30 CONTINUE  
   
c Fill in the x-values of inserted points   
c the corresponding y and z values and the three components of  
c magnetic field will be set by splfit  
   
      D = 0 
      K = 1 
      DO 50 I = 2,NX    
          KK = K+INSERT(I-1)    
          IF (INSERT(I-1).EQ.0)    GO TO 50 
C Distribute points uniformly at distance D along X
          D = (X(KK+1)-X(K))/FLOAT(INSERT(I-1)+1)   
          KL = K+1  
          DO J = KL,KK   
            X(J) = X(J-1)+D   
          ENDDO
50    K = KK+1  
   
c A call to splfit follows. 
c Splfit must interpolate to find values for y and z at the 
c intermediate points. It must set the derivatives yi and zi at 
c all points, bx,by and bz at the interpolated points.  
c a different interpolation routine may be substituted. 
   
cc      print *, 'MOMENTM: calling splfit'
      CALL SPLFIT   
  
      ITER = 0  
C Beginning of iterations   
   60 ITER = ITER+1 
   
C Find p.d2y/dx2 and p.d2z/dx2  (pyii and pzii) at all points.  
C from equations of motion  
C [Eq (1) and (2) in section 5 of Winds paper...]
   
      DO 70 I = 1,NTO   
          SQR = SQRT(1.0+YI(I)**2+ZI(I)**2) 
          PYII(I) = 
     +     SQR*(BX(I)*ZI(I)+BY(I)*YI(I)*ZI(I)-BZ(I)*(1.0+YI(I)**2))    
          PZII(I) = 
     +     SQR*(BY(I)*(1.0+ZI(I)**2)-BX(I)*YI(I)-BZ(I)*YI(I)*ZI(I))    
70    CONTINUE  
   
c We make a fit through pyii and pzii, to estimate yii and zii. 
c To do this we call dtspl. 
c DTSPL must provide single and double integrals of p.d2y/dx and    
c p.d2z/dx2 (choosing the integration constants so as to make   
c the integrals zero at the first point. This choice is 
c arbitrary)    
   
cc      print *, 'MOMENTM: calling DTSPL'
      CALL DTSPL    
   
c We now compare our integrals (yii = p*(y-a2*x-a1) and zii =   
c p*(z-b2*x-b1)) with known values of y and z. This gives pmag  
c (the quotient of momentum and charge) and the integration 
c constants.    
   
C We make a least squares fit to find p, a1, a2, b1 and b2. 
      SCC = 0.0 
      SYC = 0.0 
      SXC = 0.0 
      SXY = 0.0 
      SXX = 0.0 
      SC = 0.0  
      SY = 0.0  
      SX = 0.0  
      A = 0.0   
      SXZ = 0.  
      SZ = 0.   
      SZCZ = 0. 
      SXCZ = 0. 
      SCZ = 0.  
      SCCZ = 0. 

      DO 80 I = 1,NX    
          K = IMEAS(I)
          A = A+WEIGHT(I)   
          SX = SX+X(K)*WEIGHT(I)    
          SY = SY+Y(K)*WEIGHT(I)    
          SC = SC+YII(K)*WEIGHT(I)  
          SXX = SXX+X(K)*X(K)*WEIGHT(I) 
          SXY = SXY+X(K)*Y(K)*WEIGHT(I) 
          SXC = SXC+X(K)*YII(K)*WEIGHT(I)   
          SYC = SYC+Y(K)*YII(K)*WEIGHT(I)   
          SCC = SCC+YII(K)*YII(K)*WEIGHT(I) 
          SCCZ = SCCZ+ZII(K)*ZII(K)*WEIGHT(I)   
          SCZ = SCZ+ZII(K)*WEIGHT(I)    
          SXCZ = SXCZ+X(K)*ZII(K)*WEIGHT(I) 
          SZCZ = SZCZ+Z(K)*ZII(K)*WEIGHT(I) 
          SZ = SZ+Z(K)*WEIGHT(I)    
          SXZ = SXZ+X(K)*Z(K)*WEIGHT(I) 
c          print *, 'MOMENTM: weight: ', weight(i)
   80 CONTINUE
      SCCSCZ = SCC+SCCZ 
      SYCZCZ = SYC+SZCZ 
   
C  The matrix to be inverted to make the least squares fit is    
   
C             A   SX   SC    0    0 
C            SX  SXX  SXC    0    0 
C            SC  SXC SCCSCCZ SCZ SXCZ   
C             0    0  SCZ    A   SX 
C             0    0 SXCZ   SX  SXX 
   
   
c  The unknowns are a1,a2,d,b1,b2   
c  The right hand sides are sy,sxy,syczcz,sz,sxz    
   
c  The meaning of t1 to t10 is clear from the expressions a1=...    
      T1 = SXX*SY-SX*SXY    
      T2 = SC*SXX-SXC*SX    
      T3 = A*SXX-SX*SX  
      T4 = -SX*SY+A*SXY 
      T5 = -SC*SX+A*SXC 
      T7 = SXX*SZ-SX*SXZ    
      T8 = SXX*SCZ-SX*SXCZ  
      T10 = -SX*SZ+A*SXZ    
      T11 = -SCZ*SX+A*SXCZ  
      T31 = 1./T3   
      T33 = T3*T3   
      pmag=.00101
      D = (SYCZCZ*T33-(SC*T1+SXC*T4)*T3-(SCZ*T7+SXCZ*T10)*T3)   
cc      print *, 'Momentm: LSQ fit matrix elements and more', SYCZCZ,T33,SC,T1,SXC,T4,SCZ,T7,SXCZ,T10,T3
c      print *,' Momemtm; D=',D
      DET = (SCCSCZ*T33-(SC*T2+SXC*T5)*T3-(SCZ*T8+SXCZ*T11)*T3) 
      IF (DET.NE.0.0) THEN
        DET = 1./DET  
      ELSE
        IERR = 4
        DET = XMAX
      ENDIF
cc      print *,' Momemtm: DET=',DET
      D = D*DET 
c      print *,' Momemtm; D=',D
      A1 = (T1-D*T2)*T31    
      A2 = (T4-D*T5)*T31    
      B1 = (T7-D*T8)*T31    
      B2 = (T10-D*T11)*T31
c      print *,' Momemtm; CLIGHT=',CLIGHT
      IF (D.NE.0.0) THEN
        PMAG = CLIGHT/D
      ELSE
        PMAG = XMAX
        IERR = 5
      ENDIF    
c      print *,' Momemtm; PMAG=',PMAG
c The parameters a1, a2, b1, b2, d, define a track that does not    
c neccesarily go exactly through all of the measured points.    

      IF (ITER.GE.NITERX)    GO TO 120  

C Prepare next iteration    
      K = 1 
      DO 100 I = 2,NX   
          IF (INSERT(I-1).EQ.0)    GO TO 100    
          IJK = INSERT(I-1) 
          DO 90 J = 1,IJK   
              K = K+1   
              Y(K) = D*YII(K)+A2*X(K)+A1    
              Z(K) = D*ZII(K)+B2*X(K)+B1    
              XX = AIN(W11,W12,W13,X(K),Y(K),Z(K))  
              YY = AIN(W21,W22,W23,X(K),Y(K),Z(K))  
              ZZ = AIN(W31,W32,W33,X(K),Y(K),Z(K))  
              CALL GETB(XX,YY,ZZ,B(1),B(2),B(3),IERROR)
              BX(K) = AIN(W11,W21,W31,B(1),B(2),B(3))  
              BY(K) = AIN(W12,W22,W32,B(1),B(2),B(3))  
              BZ(K) = AIN(W13,W23,W33,B(1),B(2),B(3))  
90        CONTINUE  
100   K = K+1   
      DO K = 1,NTO  
          YI(K) = D*YI(K)+A2    
          ZI(K) = D*ZI(K)+B2    
      ENDDO
      GO TO 60 

C After last iteration continue here    
  120 CONTINUE  
      T1 = D*YI(1)+A2   
      T2 = D*ZI(1)+B2   
      PV1 = ABS(PMAG)/SQRT(1.0+T1**2+T2**2) 
      PV2 = PV1*T1  
      PV3 = PV1*T2  
      PVECIN(1) = AIN(W11,W12,W13,PV1,PV2,PV3)    
      PVECIN(2) = AIN(W21,W22,W23,PV1,PV2,PV3)    
      PVECIN(3) = AIN(W31,W32,W33,PV1,PV2,PV3)    

      T1 = D*YI(IMEAS(NXM))+A2   
      T2 = D*ZI(IMEAS(NXM))+B2   
      PV1 = ABS(PMAG)/SQRT(1.0+T1**2+T2**2) 
      PV2 = PV1*T1  
      PV3 = PV1*T2  
      PVECOUT(1) = AIN(W11,W12,W13,PV1,PV2,PV3)    
      PVECOUT(2) = AIN(W21,W22,W23,PV1,PV2,PV3)    
      PVECOUT(3) = AIN(W31,W32,W33,PV1,PV2,PV3)    
   
C Set up error matrix of parameters a1,a2,d,b1 and b2   
      ERR(1,1) = SXX*T31+T2*T2*DET  
      ERR(1,2) = -SX*T31+T2*T5*DET  
      ERR(2,1) = -SX*T31+T2*T5*DET  
      ERR(1,3) = -T2*T3*DET 
      ERR(3,1) = -T2*T3*DET 
      ERR(1,4) = +T2*T8*DET 
      ERR(4,1) = +T2*T8*DET 
      ERR(1,5) = +T2*T11*DET    
      ERR(5,1) = +T2*T11*DET    
      ERR(2,2) = A*T31+T5*T5*DET    
      ERR(2,3) = -T5*T3*DET 
      ERR(3,2) = -T5*T3*DET 
      ERR(2,4) = +T5*T8*DET 
      ERR(4,2) = +T5*T8*DET 
      ERR(2,5) = +T5*T11*DET    
      ERR(5,2) = +T5*T11*DET    
      ERR(3,3) = T33*DET    
      ERR(3,4) = -T3*T8*DET 
      ERR(4,3) = -T3*T8*DET 
      ERR(3,5) = -T3*T11*DET    
      ERR(5,3) = -T3*T11*DET    
      ERR(4,4) = SXX*T31+T8*T8*DET  
      ERR(4,5) = -SX*T31+T8*T11*DET 
      ERR(5,4) = -SX*T31+T8*T11*DET 
      ERR(5,5) = A*T31+T11*T11*DET  
   
  130 continue

cc      print *,'Momentm: passing label 130'
      
C Get errors and something like a normalized CHI-SQUARE in HEAT coordinates
        WX = 0.0
        WY = 0.0
        WZ = 0.0
        SQX = 0.0 
        SQY = 0.0 
        SQZ = 0.0
        DO I = 1,NX   ! Loop over measured points only
            K = IMEAS(I)

C Get trajectory points based upon hit points (Same X in ROTATED system)
            XTMP = X(K)
            YTMP = A1+A2*X(K)+D*YII(K)
            ZTMP = B1+B2*X(K)+D*ZII(K)

C Get trajectory points in HEAT system
            XOUT(I) = AIN(W11,W12,W13,XTMP,YTMP,ZTMP)    
            YOUT(I) = AIN(W21,W22,W23,XTMP,YTMP,ZTMP)    
            ZOUT(I) = AIN(W31,W32,W33,XTMP,YTMP,ZTMP)    

C Get error in MOMENTM system
            XERR = AIN(W11,W21,W31,DXM(1),DYM(1),DZM(1)) 
            YERR = Y(K) - A1 - A2*X(K) - D*YII(K)
            ZERR = Z(K) - B1 - B2*X(K) - D*ZII(K)

C Convert errors to HEAT system
            ERRX(I) = AIN(W11,W12,W13,XERR,YERR,ZERR)
            ERRY(I) = AIN(W21,W22,W23,XERR,YERR,ZERR)
            ERRZ(I) = AIN(W31,W32,W33,XERR,YERR,ZERR)

C Find CHI-2 in HEAT system...
            WX = WX + 1.0/DXM(I)**2
            WY = WY + 1.0/DYM(I)**2
            WZ = WZ + 1.0/DZM(I)**2
            SQX = SQX + (ERRX(I)/DXM(I))**2    
            SQY = SQY + (ERRY(I)/DYM(I))**2    
            SQZ = SQZ + (ERRZ(I)/DZM(I))**2    
        ENDDO

        DEVX = SQRT(SQX/WX)   
        DEVY = SQRT(SQY/WY)   
        DEVZ = SQRT(SQZ/WZ)   
C Insure NXOUT is set correctly...
        NXOUT = NX

        RETURN    

C   Wrong setup conditions lead into this error section   
C     NTOT not in range 3....NTOTX  
  140 IERR = 1  
cc      print *,'Momentm: IERR=1'
      GO TO  130    
C   X-coordinates not in ascending order after rotation   
  150 IERR = 3  
cc      print *,'Momentm: IERR=3'
      GO TO  130    
      END
