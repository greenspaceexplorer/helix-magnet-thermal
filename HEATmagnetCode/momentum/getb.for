        subroutine getb(xx,yy,zz,bxx,byy,bzz,ierror)
        implicit none
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Gets HELIX magnetic field value at any xx,yy,zz 
c  where coordinate origin is hodoscope center, zz is up, and xx is 
c  towards stack (+B).
c
c  ierror=  0  Normal return
c           1  xx out of range
c           2  rho=sqrt(zz*zz+yy*yy) out of range
c
c  12/16 S. L. Nutter -NKU
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        real*8 xx,yy,zz,bxx,byy,bzz
        integer ierror


c define map size and spatial limits:
c------------------------------------
c compute bfield
c------------------------------------

C Uniform field along x axis of 1.0 T -- normal helix coords
c      bxx = 10000.00
c      byy = 0.0
c      bzz = 0.0
      
c rotated helix coords with B along z axis
      bzz = 10000.00
      byy = 0.0
      bxx = 0.0




      ierror = 0
      
        return
        end
c end of GETB
