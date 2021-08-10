c**************************************************************************
c*                                                                        *
c*                        program tube generator                          *
c* (based on the scheme of PRB vol. 47. p5485 1993 C.T. White,...         *
c*                     written by Kim Kwang won                           *
c*                                                                        *
c**************************************************************************
c
c
      INTEGER n1,n2,N
      INTEGER p1,p2
      real angle,theta,phi
      REAL ro(3,3) 
      DIMENSION r1(3), r2(3), r(3),x(100,1000,3),tx(100,1000,3) 
c
      phi = 4.0*atan(1.0)
c     
c     open files
c
      open(unit=1,file='tube.dat',form='formatted',status='old')
      open(unit=2,file='tube.out',form='formatted')
c
c=========================================================================
c
c     part of reading input parameter
c     You must input heliticity vector (with respect to the standard primi
c     tive basis of graphite) and tube length
c
c=========================================================================
c
      thr = sqrt(3.0)
      two = 2.0
      unit = 0.0
      ang = 0.0
c
c     To obtain largest common dividor======>N
c
      read (1,*)
      read (1,*) n1,n2
      read (1,*) istep
      read (1,*) adis
      write(*,*)'C-C distance input=',adis
      if ( n1 .lt. n2 ) goto 190
      if ( n1 .eq. n2 ) then
         N = n1
         goto 15
        else if ( n1 .lt. n2) then
           mo = n1
           ja = n2
           else 
              mo = n2
              ja = n1
             endif
       if (n2 .eq. 0) then
           N = n1
           else
 10    mok = ja/mo
       nam = ja - mok*mo
       if (nam .eq. 0)then
          N = mo
          goto 15
       endif
       ja = mo
       mo = nam
      goto 10
      endif
 15    continue
c
c     ending of part
c============================================================
c
c
 20   format(3x,I2,4x,I2)
 22   format(6x,I2)
 24   format(7x,f4.2)
      close(unit=1)
c      
c
c     First Motif ( first and second atom )
c
      r1(1) = thr*adis
      r1(2) = 0.0
      r1(3) = 0.0
      r2(1) = thr*0.5*adis
      r2(2) = 1.5*adis
      r2(3) = 0.0
      r(1) = thr*adis*(n1+0.5*n2)
      r(2) = 1.5*n2*adis
      r(3) = 0.0
      tubr = (3.0*adis**2)*((n1+n2)**2-n1*n2)
      r_t = (sqrt(tubr))/(two*phi)
      theta = phi*(n1+n2)/((n1+n2)**2-n1*n2)
      i = 1
      j = 1
      call translation(n1,n2,i,j,unit)
      x(1,1,1) = r_t
      x(1,1,2) = 0.0
      x(1,1,3) = 0.0
      x(1,2,1) = r_t*cos(theta)
      x(1,2,2) = r_t*sin(theta)
      x(1,2,3) = adis*unit
c
c     ending of first motif and begining of ist floor
c
      ang = two*phi/float(N)
      call rot(ang,ro)
      do 30 i = 1,N-1
       k = 2*i - 1
       bx = x(1,k,1)
       by = x(1,k,2)
       bz = x(1,k,3)
       fx = x(1,k+1,1)
       fy = x(1,k+1,2)
       fz = x(1,k+1,3)
       ang=ang*i
       x(1,k+2,1) = ro(1,1)*bx + ro(1,2)*by + ro(1,3)*bz
       x(1,k+3,1) = ro(1,1)*fx + ro(1,2)*fy + ro(1,3)*fz
       x(1,k+2,2) = ro(2,1)*bx + ro(2,2)*by + ro(2,3)*bz
       x(1,k+3,2) = ro(2,1)*fx + ro(2,2)*fy + ro(2,3)*fz
       x(1,k+2,3) = ro(3,1)*bx + ro(3,2)*by + ro(3,3)*bz
       x(1,k+3,3) = ro(3,1)*fx + ro(3,2)*fy + ro(3,3)*fz
 30   continue
 40   format(3f10.5)
c 40   format(f10.5)
c
c     ending of 1st floor and beginning of next floor
c     using operation S(h,a)
c
      call screw(n1,n2,N,p1,p2)
      call translation(n1,n2,p1,p2,unit_p)
      unit = 3.0*unit_p*adis
      u_theta = (p1*n1 + p2*n2 + 0.5*(p1*n2 + n1*p2))
      d_theta = (n1+n2)**2 - n1*n2
      theta = 2.0*phi*(u_theta/d_theta)
      call rot (theta,ro)
c
c     ifloor indicates its locating block
c
      do ifloor = 2,istep
        m = 2*N
        do i = 1,m
           rx = x(ifloor-1,i,1)
           ry = x(ifloor-1,i,2)
           rz = x(ifloor-1,i,3)
           x(ifloor,i,1) = ro(1,1)*rx + ro(1,2)*ry + ro(1,3)*rz
           x(ifloor,i,2) = ro(2,1)*rx + ro(2,2)*ry + ro(2,3)*rz
           x(ifloor,i,3) = x(ifloor-1,i,3)+unit
        enddo
      enddo  
 190  continue
      inumb = 2*N*istep 
c
c     print-out
c
      write(2,*)'Results of tube-generator'
      write(2,200) n1,n2
      write(2,210) N
      write(2,220) adis
      write(2,230) p1,p2
      write(2,240) unit
      write(2,250) theta
      write(2,260) inumb
      write(2,*) 'Atom position'
      do ifloor = 1,istep
         m = 2*N
         do i = 1,m
            write(2,40)(x(ifloor,i,j),j=1,3)
            write(9,*)(x(ifloor,i,j),j=1,3)
         enddo
      enddo   
 200  format(2x,'Heliticity vector = >',I2,'R1 +',I2,'R2')
 210  format(2x,'Rotational symmetry =>  C',i2,'group')
 220  format(2x,'Size of n.n distance = ',f10.5)
 230  format(2x,'H vector =>',i2,'R1 +',i2,'R2')
 240  format(2x,'h-translation =',f8.5)
 250  format(2x,'h-rotation =',f8.5)
 260  format(2x,'number of atoms =',i3)
c
      close (unit=2)
      if (n1 .lt. n2) write (*,*)'n1 must be lager than n2!!'
      stop
      end
c
      subroutine rot(theta,R)
c        
      REAL R(3,3)
      REAL theta
c
         R(1,1) = cos(theta)
         R(1,2) = -sin(theta)
         R(1,3) = 0.0
         R(2,1) = -R(1,2)
         R(2,2) = R(1,1)
         R(2,3) = 0.0
         R(3,1) = 0.0
         R(3,2) = 0.0
         R(3,3) = 1.0
         return 
         end
c
      subroutine translation(i,j,k,m,unit)
	u_arg = abs( k*j - m*i )
        d_arg = sqrt(float((i+j)**2-i*j))
        unit = abs(0.5*(u_arg/d_arg))
        return
        end
c
      subroutine screw(n1,n2,N,ip1,ip2)
      if (( n1 .eq. 0) .or. (n2 .eq.0)) then
        if ( n1 .eq. 0 ) then
           write(*,*) 'You must enter data corectly.'
           else
              if ( n2 .eq. 0) then
                 ip1 = 0
                 ip2 = 1
                 goto 20
                 endif
        endif         
        else 
        do i = 1,n1
           do j = 1,n2
              ieq = (j*n1-i*n2) - N
              if ( ieq .eq. 0 ) then
                 ip1 = i
                 ip2 = j
                 goto 20
              endif   
           enddo
        enddo
 20     continue
        return
        endif
        end

































