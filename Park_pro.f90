!123456     
      real, allocatable :: eigup(:,:),eigdn(:,:),feigup(:,:)
      real, allocatable :: feigdn(:,:),xrk(:,:), xk(:,:), kval(:)
      integer :: nband, nks, nline, istart, nrem, check, line, next_line
      integer :: nion
      real :: k, kpt, dis_k
      real :: fermi,spin
      real :: b11, b12, b13, b21, b22, b23, b31, b32, b33
      real :: c11, c12, c13, c21, c22, c23, c31, c32, c33
      real :: a11, a12, a13, a21, a22, a23, a31, a32, a33
      real, allocatable :: ts(:,:,:), tp(:,:,:), td(:,:,:), tt(:,:,:)
      real, allocatable :: xs(:,:,:), xp(:,:,:), xd(:,:,:), xt(:,:,:)
      real, allocatable :: ys(:,:,:), yp(:,:,:), yd(:,:,:), yt(:,:,:)
      real, allocatable :: zs(:,:,:), zp(:,:,:), zd(:,:,:), zt(:,:,:)
      real, allocatable :: tts(:,:), ttp(:,:), ttd(:,:), ttt(:,:)
      real, allocatable :: txs(:,:), txp(:,:), txd(:,:), txt(:,:)
      real, allocatable :: tys(:,:), typ(:,:), tyd(:,:), tyt(:,:)
      real, allocatable :: tzs(:,:), tzp(:,:), tzd(:,:), tzt(:,:)
      character (len=20) :: tot
      real, allocatable :: xtot(:,:), ytot(:,:), ztot(:,:)
      real, allocatable :: fspinup(:,:),fspindown(:,:)
      real :: Pi
!
!   please remember the number of writting field in the format should be larger than the number of bands
!123456
      OPEN(1,file='band.dat')
      read(1,*) fermi
      write(6,*) "spin unpolar(type '0') or spin polar(type '1')"
      write(6,*) "fermi level"
      write(6,*) "b vector"

      open(17,file='CONTCAR')
      read(17,*)
      read(17,*)
      read(17,*) a11, a12, a13
      read(17,*) a21, a22, a23
      read(17,*) a31, a32, a33
      
      c11 = a22*a33 - a23*a32
      c12 = a23*a31 - a21*a33
      c13 = a21*a32 - a22*a31
      c21 = a32*a13 - a33*a12
      c22 = a33*a11 - a31*a13
      c23 = a31*a12 - a32*a11
      c31 = a12*a23 - a13*a22
      c32 = a13*a21 - a11*a23
      c33 = a11*a22 - a12*a21

      Pi = Atan(1.0)*4

      b11 = 2*Pi*c11/(a11*c11 + a12*c12 + a13*c13)
      b12 = 2*Pi*c12/(a11*c11 + a12*c12 + a13*c13)
      b13 = 2*Pi*c13/(a11*c11 + a12*c12 + a13*c13)
      b21 = 2*Pi*c21/(a21*c21 + a22*c22 + a23*c23)
      b22 = 2*Pi*c22/(a21*c21 + a22*c22 + a23*c23)
      b23 = 2*Pi*c23/(a21*c21 + a22*c22 + a23*c23)
      b31 = 2*Pi*c31/(a31*c31 + a32*c32 + a33*c33)
      b32 = 2*Pi*c32/(a31*c31 + a32*c32 + a33*c33)
      b33 = 2*Pi*c33/(a31*c31 + a32*c32 + a33*c33)

      open(15,file='EIGENVAL')
      read(15,*) nion, c2, c3, spin

       If(spin.eq.1) then
      do line=1,4
        read(15,*)
        next_line=next_line+line
      enddo
      read(15,*)Nelectron,nk,nband
      allocate (feigup(nk,nband),eigup(nk,nband))
      allocate (xrk(nk,3),xk(nk,3),kval(nk))

      do i =1,nk
         read(15,*)
         read(15,*)xrk(i,1),xrk(i,2),xrk(i,3),tmp123
         do j = 1,nband 
            read(15,*)itmp,eigup(i,j) 
         enddo
      enddo

      xk(:,1) = b11*xrk(:,1) + b21*xrk(:,2) + b31*xrk(:,3)
      xk(:,2) = b12*xrk(:,1) + b22*xrk(:,2) + b32*xrk(:,3)
      xk(:,3) = b13*xrk(:,1) + b23*xrk(:,2) + b33*xrk(:,3)
      
      k = 0.0
      kpt = 0.0
      kval(1) = 0.0
      do i=2,nk
         k=sqrt((xk(i,1)-xk(i-1,1))**2 + (xk(i,2)-xk(i-1,2))**2 + &
                 (xk(i,3)-xk(i-1,3))**2)
         kpt = kpt + k
         kval(i) = kpt
      enddo

      feigup(:,:) = eigup(:,:) - fermi

      open(30,file='PROCAR')
      allocate (ts(nk,nband,nion), tp(nk,nband,nion), td(nk,nband,nion))
      allocate (tt(nk,nband,nion), xs(nk,nband,nion), xp(nk,nband,nion))
      allocate (xd(nk,nband,nion), xt(nk,nband,nion), ys(nk,nband,nion))
      allocate (yp(nk,nband,nion), yd(nk,nband,nion), yt(nk,nband,nion))
      allocate (zs(nk,nband,nion), zp(nk,nband,nion), zd(nk,nband,nion))
      allocate (zt(nk,nband,nion))
      allocate (tts(nk,nband), ttp(nk,nband), ttd(nk,nband))
      allocate (ttt(nk,nband), txs(nk,nband), txp(nk,nband))
      allocate (txd(nk,nband), txt(nk,nband), tys(nk,nband))
      allocate (typ(nk,nband), tyd(nk,nband), tyt(nk,nband))
      allocate (tzs(nk,nband), tzp(nk,nband), tzd(nk,nband))
      allocate (tzt(nk,nband))

      read(30,*)
      read(30,*)
      read(30,*)
      read(30,*)
      read(30,*)
      read(30,*)
      read(30,*)
      read(30,*)
      
      do i = 1,nk-1
        do j = 1,nband
          do k = 1,nion
            read(30,*)n,ts(i,j,k),tp(i,j,k),td(i,j,k),tt(i,j,k)
          enddo

          read(30,*)tot,tts(i,j),ttp(i,j),ttd(i,j),ttt(i,j)
          
          do k = 1,nion
            read(30,*)n,xs(i,j,k),xp(i,j,k),xd(i,j,k),xt(i,j,k)
          enddo

          read(30,*)tot,txs(i,j),txp(i,j),txd(i,j),txt(i,j)

          do k = 1,nion
            read(30,*)n,ys(i,j,k),yp(i,j,k),yd(i,j,k),yt(i,j,k)
          enddo

          read(30,*)tot,tys(i,j),typ(i,j),tyd(i,j),tyt(i,j)

          do k = 1,nion
            read(30,*)n,zs(i,j,k),zp(i,j,k),zd(i,j,k),zt(i,j,k)

          enddo

          read(30,*)tot,tzs(i,j),tzp(i,j),tzd(i,j),tzt(i,j)

        
          read(30,*)
          read(30,*)
          read(30,*)
          read(30,*)
      enddo

          read(30,*)
          read(30,*)
          read(30,*)
      enddo

       i = nk
        do j = 1,nband-1
          do k = 1,nion
            read(30,*)n,ts(i,j,k),tp(i,j,k),td(i,j,k),tt(i,j,k)
          enddo

          read(30,*)tot,tts(i,j),ttp(i,j),ttd(i,j),ttt(i,j)

          do k = 1,nion
            read(30,*)n,xs(i,j,k),xp(i,j,k),xd(i,j,k),xt(i,j,k)
          enddo

          read(30,*)tot,txs(i,j),txp(i,j),txd(i,j),txt(i,j)

          do k = 1,nion
            read(30,*)n,ys(i,j,k),yp(i,j,k),yd(i,j,k),yt(i,j,k)
          enddo
          
          read(30,*)tot,tys(i,j),typ(i,j),tyd(i,j),tyt(i,j)
          
          do k = 1,nion
            read(30,*)n,zs(i,j,k),zp(i,j,k),zd(i,j,k),zt(i,j,k)
          
          enddo
          
          read(30,*)tot,tzs(i,j),tzp(i,j),tzd(i,j),tzt(i,j)

        
          read(30,*)
          read(30,*)
          read(30,*)
          read(30,*)
        enddo

      i = nk
      j = nband
          do k = 1,nion
            read(30,*)n,ts(i,j,k),tp(i,j,k),td(i,j,k),tt(i,j,k)
          enddo

          read(30,*)tot,tts(i,j),ttp(i,j),ttd(i,j),ttt(i,j)

          do k = 1,nion
            read(30,*)n,xs(i,j,k),xp(i,j,k),xd(i,j,k),xt(i,j,k)
          enddo

          read(30,*)tot,txs(i,j),txp(i,j),txd(i,j),txt(i,j)

          do k = 1,nion
            read(30,*)n,ys(i,j,k),yp(i,j,k),yd(i,j,k),yt(i,j,k)
          enddo

          read(30,*)tot,tys(i,j),typ(i,j),tyd(i,j),tyt(i,j)

          do k = 1,nion
            read(30,*)n,zs(i,j,k),zp(i,j,k),zd(i,j,k),zt(i,j,k)

          enddo

          read(30,*)tot,tzs(i,j),tzp(i,j),tzd(i,j),tzt(i,j)


      open(24,file='band.txt')
 118  format(f14.7,3x,f14.7)
        write(24,*) '#kval, total'
      do j = 1, nband
      do i = 1, nk
          write(24,118)kval(i), feigup(i,j)
      enddo
          write(24,*)
      enddo
        
        allocate(xtot(nk,nband),ytot(nk,nband),ztot(nk,nband))

        do i = 1, nk
          do j = 1, nband
        if(txt(i,j)**2.GT.tyt(i,j)**2.and.txt(i,j)**2.GT.tzt(i,j)**2)  then
              xtot(i,j) = txt(i,j)/ttt(i,j)
            endif
            
         if(tyt(i,j)**2.GT.tzt(i,j)**2.and.tyt(i,j)**2.GT.txt(i,j)**2) then
              ytot(i,j) = tyt(i,j)/ttt(i,j)
            endif
            
         if(tzt(i,j)**2.GT.txt(i,j)**2.and.tzt(i,j)**2.GT.tyt(i,j)**2) then
              ztot(i,j) = tzt(i,j)/ttt(i,j)
            endif
          enddo
        enddo

        allocate(fspinup(nk,nband), fspindown(nk,nband))

        do i = 1, nk
          do j = 1, nband
            if(xtot(i,j) .GT. 0) then
              fspinup(i,j) = feigup(i,j)
             endif
         
            if(xtot(i,j) .LT. 0 ) then  
              fspindown(i,j) = feigup(i,j)
            endif
          
            if(ytot(i,j) .GT. 0) then
              fspinup(i,j) = feigup(i,j)
           endif 
               
            if(ytot(i,j) .LT. 0) then
              fspindown(i,j) = feigup(i,j)
            endif

            if(ztot(i,j) .GT. 0) then
              fspinup(i,j) = feigup(i,j)
            endif 
               
            if(ztot(i,j) .LT. 0) then
              fspindown(i,j) = feigup(i,j)
            endif
           
          enddo
        enddo

        do i = 1,nk
        do j = 1, nband
          if(fspinup(i,j).eq.0) then
        fspinup(i,j) = 1000000
          endif
          if(fspindown(i,j).eq.0) then
        fspindown(i,j) = 1000000
         endif
        enddo   
        enddo

     open(40,file='band2.txt')
 140  format(f14.7,3x,f14.7,3x,f14.7,3x,f14.7)
        write(40,*) '#kval, total'
      do j = 1, nband
      do i = 1, nk
        write(40,140)kval(i), feigup(i,j), fspinup(i,j), fspindown(i,j)
      enddo
          write(40,*)
      enddo


       else

      read(15,*)
      read(15,*)
      read(15,*)
      read(15,*)
      read(15,*)Nelectron,nk,nband
      allocate (feigup(nk,nband),feigdn(nk,nband),eigup(nk,nband))
      allocate (eigdn(nk,nband),xk(nk,3),xrk(nk,3))

      do i =1,nk
         read(15,*)
         read(15,*)xrk(i,1),xrk(i,2),xrk(i,3),tmp123
         do j = 1,nband 
            read(15,*)itmp,eigup(i,j),eigdn(i,j)
         enddo
      enddo

      xk(:,1) = b11*xrk(:,1) + b21*xrk(:,2) + b31*xrk(:,3)
      xk(:,2) = b12*xrk(:,1) + b22*xrk(:,2) + b32*xrk(:,3)
      xk(:,3) = b13*xrk(:,1) + b23*xrk(:,2) + b33*xrk(:,3)
      
      k = 0.0
      kpt = 0.0
      kval(1) = 0.0
      do i=2,nks
         k=sqrt((xk(i,1)-xk(i-1,1))**2 + (xk(i,2)-xk(i-1,2))**2 + &
                (xk(i,3)-xk(i-1,3))**2)
         kpt = kpt + k
         kval(i) = kpt
      enddo

      feigup(:,:) = eigup(:,:) - fermi
      feigdn(:,:) = eigdn(:,:) - fermi

      open(30,file='PROCAR')
      allocate (ts(nk,nband,nion), tp(nk,nband,nion), td(nk,nband,nion))
      allocate (tt(nk,nband,nion), xs(nk,nband,nion), xp(nk,nband,nion))
      allocate (xd(nk,nband,nion), xt(nk,nband,nion), ys(nk,nband,nion))
      allocate (yp(nk,nband,nion), yd(nk,nband,nion), yt(nk,nband,nion))
      allocate (zs(nk,nband,nion), zp(nk,nband,nion), zd(nk,nband,nion))
      allocate (zt(nk,nband,nion))
      allocate (tts(nk,nband), ttp(nk,nband), ttd(nk,nband))
      allocate (ttt(nk,nband), txs(nk,nband), txp(nk,nband))
      allocate (txd(nk,nband), txt(nk,nband), tys(nk,nband))
      allocate (typ(nk,nband), tyd(nk,nband), tyt(nk,nband))
      allocate (tzs(nk,nband), tzp(nk,nband), tzd(nk,nband))
      allocate (tzt(nk,nband))

      read(30,*)
      read(30,*)
      read(30,*)
      read(30,*)
      read(30,*)
      read(30,*)
      read(30,*)
      read(30,*)

      do i = 1,nk-1
        do j = 1,nband
          do k = 1,nion
            read(30,*)n,ts(i,j,k),tp(i,j,k),td(i,j,k),tt(i,j,k)
          enddo

          read(30,*)tot,tts(i,j),ttp(i,j),ttd(i,j),ttt(i,j)

          do k = 1,nion
            read(30,*)n,xs(i,j,k),xp(i,j,k),xd(i,j,k),xt(i,j,k)
          enddo

          read(30,*)tot,txs(i,j),txp(i,j),txd(i,j),txt(i,j)

          do k = 1,nion
            read(30,*)n,ys(i,j,k),yp(i,j,k),yd(i,j,k),yt(i,j,k)
          enddo

          read(30,*)tot,tys(i,j),typ(i,j),tyd(i,j),tyt(i,j)

          do k = 1,nion
            read(30,*)n,zs(i,j,k),zp(i,j,k),zd(i,j,k),zt(i,j,k)

          enddo

          read(30,*)tot,tzs(i,j),tzp(i,j),tzd(i,j),tzt(i,j)


          read(30,*)
          read(30,*)
          read(30,*)
          read(30,*)
      enddo

          read(30,*)
          read(30,*)
          read(30,*)
      enddo
      
       i = nk
        do j = 1,nband-1
          do k = 1,nion
            read(30,*)n,ts(i,j,k),tp(i,j,k),td(i,j,k),tt(i,j,k)
          enddo

          read(30,*)tot,tts(i,j),ttp(i,j),ttd(i,j),ttt(i,j)

          do k = 1,nion
            read(30,*)n,xs(i,j,k),xp(i,j,k),xd(i,j,k),xt(i,j,k)
          enddo

          read(30,*)tot,txs(i,j),txp(i,j),txd(i,j),txt(i,j)

          do k = 1,nion
            read(30,*)n,ys(i,j,k),yp(i,j,k),yd(i,j,k),yt(i,j,k)
          enddo

          read(30,*)tot,tys(i,j),typ(i,j),tyd(i,j),tyt(i,j)

          do k = 1,nion
            read(30,*)n,zs(i,j,k),zp(i,j,k),zd(i,j,k),zt(i,j,k)

          enddo

          read(30,*)tot,tzs(i,j),tzp(i,j),tzd(i,j),tzt(i,j)


          read(30,*)
          read(30,*)
          read(30,*)
          read(30,*)
        enddo

      i = nk
      j = nband
          do k = 1,nion
            read(30,*)n,ts(i,j,k),tp(i,j,k),td(i,j,k),tt(i,j,k)
          enddo

          read(30,*)tot,tts(i,j),ttp(i,j),ttd(i,j),ttt(i,j)

          do k = 1,nion
            read(30,*)n,xs(i,j,k),xp(i,j,k),xd(i,j,k),xt(i,j,k)
          enddo

          read(30,*)tot,txs(i,j),txp(i,j),txd(i,j),txt(i,j)

          do k = 1,nion
            read(30,*)n,ys(i,j,k),yp(i,j,k),yd(i,j,k),yt(i,j,k)
          enddo

          read(30,*)tot,tys(i,j),typ(i,j),tyd(i,j),tyt(i,j)

          do k = 1,nion
            read(30,*)n,zs(i,j,k),zp(i,j,k),zd(i,j,k),zt(i,j,k)

          enddo

          read(30,*)tot,tzs(i,j),tzp(i,j),tzd(i,j),tzt(i,j)

        write(24,*) '#kval, up, down'
      do j = 1, nband
      do i = 1, nk
          write(24,118)kval(i),feigup(i,j),feigdn(i,j)
      enddo
          write(24,*)
      enddo

      endif

      stop

      end
      
