C      OPTIONS/CHECK=ALL/EXTEND_SOURCE
      SUBROUTINE SPLFIT 
C   
C Performs a third-order spline fit through the measured    
C points to interpolate y,z,yi (=dy/dx) and zi (=dz/dx) 
C for the inserted points.  
C   
C         INPUT PARAMETERS  
C             NXM,NTOT        AS DEFINED IN MOMENTM 
C             X               FOR ALL POINTS    
C             Y,Z             FOR THE MEASURED POINTS   
C             INSERT          AS DEFINED IN MOMENTM 
C             EXM,EYM,EZM     Estimate of errors on measured points
C   
C         OUTPUT PARAMETERS 
C             Y,Z             FOR THE INSERTED POINTS   
C             YI,ZI           FOR ALL POINTS    
C             BX,BY,BZ        FOR THE INSERTED POINTS   
      IMPLICIT NONE

      INTEGER NXM,INSERT
      REAL*8 XM,YM,ZM,DXM,DYM,DZM
      COMMON /MOM_INPUT/XM(99),YM(99),ZM(99),
     +              DXM(99),DYM(99),DZM(99),NXM,INSERT(99)
      INTEGER NTOT
      REAL*8 X,Y,Z,YI,ZI
      COMMON /LOCAL1/ X(100),Y(100),Z(100),YI(100),ZI(100),NTOT
      REAL*8 W11,W12,W13,W21,W22,W23,W31,W32,W33    
      COMMON /LOCAL3/W11,W12,W13,W21,W22,W23,W31,W32,W33    
      REAL*8 BX,BY,BZ
      COMMON /LOCAL4/BX(100),BY(100),BZ(100)    
   
C  Locally used second derivatives   
      REAL*8 PYII,PZII,YII,ZII
      COMMON /LOCAL5/PYII(100),PZII(100),YII(100),ZII(100)  
   
      REAL*8 B(3)

      REAL*8 THIRD,SIXTH
      DATA THIRD,SIXTH/.333333333333333,.166666666666666/   

C Local variables:
      INTEGER N1,K2,N2,I,I1,K,I2,NTOT1,J,J1,K3,II,JJ,H,IERROR
      REAL*8 D2,DZ2,DY2,DD,DY1,DZ1,RY,RZ,T1,D1,T2,DX,DY,DZ,D1Y,D1Z
      REAL*8 DX1,DX2,D2Y,D2Z,YPPP,ZPPP,XX,YY,ZZ
      REAL*8 AIN,XP,YP,ZP,AP,BP,CP
      AIN(XP,YP,ZP,AP,BP,CP) = XP*AP+YP*BP+ZP*CP    
   
C Force the natural spline condition that second derivatives
C ZII (=D2Z/DX2), YII (=D2Y/DX2)
C be zero at the end points...
      YII(NTOT) = 0.0   
      YII(1) = 0.0  
      ZII(NTOT) = 0.0   
      ZII(1) = 0.0  

C Set first derivatives to 1 for ZI (=DZ/DX) and YI (=DY/DZ)
      ZI(1) = 1.0   
      YI(1) = 1.0   

C Get information for measured points (YI, ZI, YII, ZII):
      N1 = NXM-1    
      IF (NXM.EQ.3)    GO TO 40 
C Index of first inserted point?
      K2 = 2+INSERT(1)  
      D2 = X(K2)-X(1)   
      DY2 = Y(K2)-Y(1)  
      DZ2 = Z(K2)-Z(1)  
      N2 = NXM-2    
      I = 1 
      I1 = 2+INSERT(1)  
C Loop around measured points to get yi,zi,yii,zii  
      DO 10 K = 3,NXM   
        I2 = I1+INSERT(K-1)+1 
        D1 = D2   
        D2 = X(I2)-X(I1)  
        DD = D1+D2    
        DY1 = DY2 
        DZ1 = DZ2 
        DY2 = Y(I2)-Y(I1) 
        DZ2 = Z(I2)-Z(I1) 
        RY = 3.0*(DY2/D2-DY1/D1)/DD   
        RZ = 3.0*(DZ2/D2-DZ1/D1)/DD   
        T1 = 0.5*D1/DD    
        T2 = 0.5*D2/DD    
        YI(I1) = -T2/(1.0+T1*YI(I))   
        ZI(I1) = -T2/(1.0+T1*ZI(I))   
        YII(I1) = (RY-T1*YII(I))/(1.0+T1*YI(I))   
        ZII(I1) = (RZ-T1*ZII(I))/(1.0+T1*ZI(I))   
        I = I1    
        I1 = I2   
   10 CONTINUE
   
c      print *,' In SPLFIT...NXM = ',NXM
C Find index of second to last measured point...
      NTOT1 = NTOT-1-INSERT(NXM-1)  
C Extrapolate for YII, ZII based upon last measured point...
      YII(NTOT) = YII(NTOT1)/(1.0-YI(NTOT1))    
      ZII(NTOT) = ZII(NTOT1)/(1.0-ZI(NTOT1))    
C Set 2nd derivatives equal for last 2 points...
      YII(NTOT1) = YII(NTOT)    
      ZII(NTOT1) = ZII(NTOT)    

C Start with second to last measured point...
      J = NTOT1 
C Loop over all points but the last 2 (set above)
      DO 20 I = 2,N2    
        J1 = J    
C Get index for prior measured point...
        J = J-1-INSERT(NXM-I) 
C Find 2nd derivatives
        YII(J) = YI(J)*YII(J1)+YII(J) 
   20   ZII(J) = ZI(J)*ZII(J1)+ZII(J) 
C Set first measured point to 2nd derivates of first inserted point?
      YII(1) = YII(K2)  
      ZII(1) = ZII(K2)  

C Build first derivatives based upon second derivates 
C (Eq. (7) in Appendix in Winds paper (NIM 115 (1974) pg 431-434))
C Measured points first
      I = 1 
      I1 = INSERT(1)+2  
      DO 30 K = 2,NXM   
        DX = X(I1)-X(I)   
        DY = Y(I1)-Y(I)   
        DZ = Z(I1)-Z(I)   
        D1Y = DY/DX   
        D1Z = DZ/DX   
        YI(I) = D1Y-DX*(THIRD*YII(I)+SIXTH*YII(I1))   
        ZI(I) = D1Z-DX*(THIRD*ZII(I)+SIXTH*ZII(I1))   
        I = I1    
        I1 = I1+INSERT(K)+1   
   30 CONTINUE  
C Get values for ZI, YI at end...
      DX = X(NTOT)-X(NTOT1) 
      DY = Y(NTOT)-Y(NTOT1) 
      DZ = Z(NTOT)-Z(NTOT1) 
      D1Y = DY/DX   
      D1Z = DZ/DX   
      YI(NTOT) = D1Y+DX*0.5*YII(NTOT)   
      ZI(NTOT) = D1Z+DX*0.5*ZII(NTOT)   
      I = 1 
      I1 = 2+INSERT(1)  

c      print *,' In SPLFIT...before 50 line'

      GO TO  50 

C  For only three measured points short section  
   40 CONTINUE  
      K2 = 2+INSERT(1)  
      K3 = K2+1+INSERT(2)   
      DX1 = X(K2)-X(1)  
      DX2 = X(K3)-X(K2) 
      D1Y = (Y(K2)-Y(1))/DX1    
      D1Z = (Z(K2)-Z(1))/DX1    
      D2Y = (Y(K3)-Y(K2))/DX2   
      D2Z = (Z(K3)-Z(K2))/DX2   
      YII(K2) = 3.0*(D2Y-D1Y)/(DX1+DX2) 
      ZII(K2) = 3.0*(D2Z-D1Z)/(DX1+DX2) 
      YI(1) = D1Y-SIXTH*DX1*YII(K2) 
      ZI(1) = D1Z-SIXTH*DX1*ZII(K2) 
      YI(K3) = D2Y+SIXTH*DX2*YII(K2)    
      ZI(K3) = D2Z+SIXTH*DX2*ZII(K2)    
      YI(K2) = YI(1)+.5*DX1*YII(K2) 
      ZI(K2) = ZI(1)+.5*DX1*ZII(K2) 

c      print *,' In SPLFIT...Starting taylor series'

C Expand in taylor series for interpolated points   
   50 CONTINUE  
      II = 1    
      DO 70 K = 1,N1    
        I = II    
        II = II+INSERT(K)+1   
        IF (INSERT(K).EQ.0)    GO TO 70   
        YPPP = (YII(II)-YII(I))/(X(II)-X(I))  
        ZPPP = (ZII(II)-ZII(I))/(X(II)-X(I))  
        JJ = INSERT(K)    
        DO 60 J = 1,JJ    
            H = X(I+J)-X(I)   
            YI(I+J) = YI(I)+H*YII(I)+0.5*H**2*YPPP    
            ZI(I+J) = ZI(I)+H*ZII(I)+0.5*H**2*ZPPP    
            Y(I+J) = Y(I)+H*YI(I)+0.5*H**2*YII(I)+SIXTH*H**3*YPPP 
            Z(I+J) = Z(I)+H*ZI(I)+0.5*H**2*ZII(I)+SIXTH*H**3*ZPPP 
C   For mag. field have to go through measurement coord.system    
            XX = AIN(W11,W12,W13,X(I+J),Y(I+J),Z(I+J))    
            YY = AIN(W21,W22,W23,X(I+J),Y(I+J),Z(I+J))    
            ZZ = AIN(W31,W32,W33,X(I+J),Y(I+J),Z(I+J))    
            CALL GETB(XX,YY,ZZ,B(1),B(2),B(3),IERROR)
c            print *,'SPLFIT: getb returns',B(1),B(2),B(3)
            IF (IERROR.NE.0) THEN
              PRINT *,' Error in SPLFIT: GETB IERROR = ',IERROR
            ENDIF
            BX(I+J) = AIN(W11,W21,W31,B(1),B(2),B(3))    
            BY(I+J) = AIN(W12,W22,W32,B(1),B(2),B(3))    
            BZ(I+J) = AIN(W13,W23,W33,B(1),B(2),B(3))    
   60   CONTINUE  
   70 CONTINUE  
c      print *,' Finished SPLFIT'

      RETURN    
      END   
