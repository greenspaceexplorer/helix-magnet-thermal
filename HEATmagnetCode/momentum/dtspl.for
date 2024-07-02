C      OPTIONS/CHECK=ALL/EXTEND_SOURCE
      SUBROUTINE DTSPL  
   
      implicit none
      
c A subroutine to simultaneously fit two third-order splines,   
c and evaluate their single and double integrals.   
C   
C METHOD  SEE APPENDIX TO CERN-NP-DHG-73/5 (WIND)   
C   
C         INPUT PARAMETERS  
C             NTOT            AS DEFINED IN MOMENTM 
C             X               FOR ALL NTOT POINTS   
C             U,V             PYII, PZII USED AS FUNCTION OF X  
C   
C         OUTPUT PARAMETERS 
C             YPI,ZPI         INTEGRATED DY/DX AND DZ/DX, ALL POINTS    
C             YP,ZP           INTEGRATED Y AND Z, ALL POINTS    
C             UI,VI           THIRD DERIVATIVES 
C             UII,VII         FOURTH DERIVATIVES    
   
      INTEGER NTOT
      REAL*8 X,Y,Z,YPI,ZPI
      COMMON /LOCAL1/X(100),Y(100),Z(100),YPI(100),ZPI(100), NTOT    
      REAL*8 U,V,YP,ZP
      COMMON /LOCAL5/U(100),V(100),YP(100),ZP(100)  
      REAL*8 UI,VI,UII,VII
      COMMON /LOCAL6/UI(100),VI(100),UII(100),VII(100)  
C  Local data    
      REAL*8 THIRD,SIXTH,P7,P8,T24TH
      DATA THIRD,SIXTH,P7,P8/.333333333333333,.166666666666666,.7,.8/   
      DATA T24TH/.041666666666666/  
   
C Local variables
      INTEGER N,N1,J,I,N2
      REAL*8 DX,DY,T1,T2,RY,RZ,DZ1,DZ2,DY1,DY2,DD,D1,D2,D1Y,D1Z
      REAL*8 DXH,DXSQ,DX2,D2Y,DX1,D2X,DX22,DX12,D2Z,DZ

      N = NTOT  
      YP(1) = 0.0   
      YPI(1) = 0.0  
      UII(N) = 0.0  
      UII(1) = 0.0  
      ZP(1) = 0.0   
      ZPI(1) = 0.0  
      VII(N) = 0.0  
      VII(1) = 0.0  
      VI(1) = 1.0   
      UI(1) = 1.0   

      IF (N.EQ.3)    GO TO 50   
C For more than three points    
      D2 = X(2)-X(1)    
      DY2 = U(2)-U(1)   
      DZ2 = V(2)-V(1)   
      N2 = N-2  
C Recursion for p and q 
      DO 10 I = 1,N2    
        D1 = D2   
        D2 = X(I+2)-X(I+1)    
        DD = D1+D2    
        DY1 = DY2 
        DZ1 = DZ2 
        DY2 = U(I+2)-U(I+1)   
        DZ2 = V(I+2)-V(I+1)   
        RY = 3.0*(DY2/D2-DY1/D1)/DD   
        RZ = 3.0*(DZ2/D2-DZ1/D1)/DD   
        T2 = 0.5*D2/DD    
        T1 = 0.5*D1/DD    
        UI(I+1) = -T2/(1.0+T1*UI(I))  
        VI(I+1) = -T2/(1.0+T1*VI(I))  
        UII(I+1) = (RY-T1*UII(I))/(1.0+T1*UI(I))  
        VII(I+1) = (RZ-T1*VII(I))/(1.0+T1*VI(I))  
   10 CONTINUE  
   
      UII(N) = UII(N-1)/(1.0-UI(N-1))   
      VII(N) = VII(N-1)/(1.0-VI(N-1))   
      UII(N-1) = UII(N) 
      VII(N-1) = VII(N) 
C Recursion for uii,vii = y'''',z''''   
      DO 20 I = 2,N2    
        J = N-I   
        UII(J) = UI(J)*UII(J+1)+UII(J)    
   20   VII(J) = VI(J)*VII(J+1)+VII(J)    
      UII(1) = UII(2)   
      VII(1) = VII(2)   
      N1 = N-1  
C Recursion for ui,vi = y''',z'''   
      DO 30 I = 1,N1    
        DX = X(I+1)-X(I)  
        DY = U(I+1)-U(I)  
        DZ = V(I+1)-V(I)  
        D1Y = DY/DX   
        D1Z = DZ/DX   
        UI(I) = D1Y-DX*(THIRD*UII(I)+SIXTH*UII(I+1))  
        VI(I) = D1Z-DX*(THIRD*VII(I)+SIXTH*VII(I+1))  
        DXH = DX*.5   
        YPI(I+1) = DXH*(U(I+1)+U(I)-THIRD*(UII(I+1)+UII(I))*DXH**2)+YPI(I)    
   30   ZPI(I+1) = DXH*(V(I+1)+V(I)-THIRD*(VII(I+1)+VII(I))*DXH**2)+ZPI(I)    
      DX = X(N)-X(N-1)  
      DY = U(N)-U(N-1)  
      DZ = V(N)-V(N-1)  
      D1Y = DY/DX   
      D1Z = DZ/DX   
      UI(N) = D1Y+DX*0.5*UII(N) 
      VI(N) = D1Z+DX*0.5*VII(N) 

C Recursion for y'',z'' 
      DO 40 I = 1,N1    
        DX = X(I+1)-X(I)  
        DXSQ = DX*DX*SIXTH    
        DXH = .5*DX   
        YP(I+1) = DXSQ*((U(I)+U(I)+U(I+1))-DXSQ*(P8*UII(I)+P7*UII(I+  
     +  1)))+YP(I)+YPI(I)*DX  
        ZP(I+1) = DXSQ*((V(I)+V(I)+V(I+1))-DXSQ*(P8*VII(I)+P7*VII(I+  
     +  1)))+ZP(I)+ZPI(I)*DX  
   40 CONTINUE  

      GO TO  60 

C Three-point procedure 
   50 CONTINUE  
      DX1 = X(2)-X(1)   
      DX2 = X(3)-X(2)   
      D1Y = (U(2)-U(1))/DX1 
      D1Z = (V(2)-V(1))/DX1 
      D2Y = (U(3)-U(2))/DX2 
      D2Z = (V(3)-V(2))/DX2 
      UII(2) = 3.0*(D2Y-D1Y)/(DX1+DX2)  
      VII(2) = 3.0*(D2Z-D1Z)/(DX1+DX2)  
      UI(1) = D1Y-SIXTH*DX1*UII(2)  
      VI(1) = D1Z-SIXTH*DX1*VII(2)  
      UI(3) = D2Y+SIXTH*DX2*UII(2)  
      VI(3) = D2Z+SIXTH*DX2*VII(2)  
      UI(2) = UI(1)+.5*DX1*UII(2)   
      VI(2) = VI(1)+.5*DX1*VII(2)   
      DX12 = DX1**2 
      YPI(2) = DX1*(.5*(U(2)+U(1))-T24TH*DX12*UII(2))   
      ZPI(2) = DX1*(.5*(V(2)+V(1))-T24TH*DX12*VII(2))   
      DX22 = DX2**2 
      YPI(3) = DX2*(.5*(U(3)+U(2))-T24TH*DX22*UII(2))+YPI(2)    
      ZPI(3) = DX2*(.5*(V(3)+V(2))-T24TH*DX22*VII(2))+ZPI(2)    
      DXSQ = SIXTH*DX1**2   
      YP(2) = DXSQ*(U(1)+U(1)+U(2)-DXSQ*(P7*UII(2)))    
      ZP(2) = DXSQ*(V(1)+V(1)+V(2)-DXSQ*(P7*VII(2)))    
      DXSQ = SIXTH*DX2**2   
      YP(3) = DXSQ*(U(2)+U(2)+U(3)-DXSQ*(P8*UII(2)))+YP(2)+YPI(2)*DX2   
      ZP(3) = DXSQ*(V(2)+V(2)+V(3)-DXSQ*(P8*VII(2)))+ZP(2)+ZPI(2)*DX2   
   60 RETURN    

      END   
