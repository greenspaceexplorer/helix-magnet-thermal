      program heat3d
c Link with HEAT_MAG
      implicit none

      real r, x, y, z, Bx, By, Bz, Br, B 

      print 13
 13   format(' Welcome to the HEAT magnet 3D Bfield finder.')
      print 14
 14   format('z is up, x is towards stack (North pole)') 
      print 16
 16   format('See heat_mag.f code for detailed description of coordinate system and methodology')
 11   print 12
 12   format(' ')
 15   write(*, 20)
 20   format(' Input point coordinates x, y, z (in centimeters): ', $)
      read(*, *) x, y, z
      Br=0.0
      Bx=0.0
      B=0.0
      r=sqrt(y*y+z*z)
c      call heat_mag (r, 0, x, Br, Bx, B)
c Try this. SLN 10/22/15
      call heat_mag (x, y, z, Br, Bx, B)
      if(r.ne.0.0) then
        Bz=z/r*Br
        By=y/r*Br
      else
        Bz=Br
        By=Br
      endif
c      write(*, 60)
c 60   format(7x,'x(cm)',4x,'y(cm)',4x,'z(cm)',5x,'Br',10x,'Bx',7x,'Btot(gauss)')
c      write(*, 70) x, y, z, Br, Bz, B
c 70   format(3f10.1,4(f12.5))
      write(*, 40)
 40   format(7x,'x(cm)',4x,'y(cm)',4x,'z(cm)',5x,'Bx',10x,'By',10x,'Bz',7x,'Btot(gauss)')
      write(*, 50) x, y, z, Bx, By, Bz, B
 50   format(3f10.1,4(f12.5))
      goto 11
      stop
      end
