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
      real, allocatable :: as(:,:,:), ap(:,:,:), ad(:,:,:), at(:,:,:)
      real, allocatable :: bs(:,:,:), bp(:,:,:), bd(:,:,:), bt(:,:,:)
      real, allocatable :: cs(:,:,:), cp(:,:,:), cd(:,:,:), ct(:,:,:)
      real, allocatable :: tas(:,:), tap(:,:), tad(:,:), tat(:,:)
      real, allocatable :: tbs(:,:), tbp(:,:), tbd(:,:), tbt(:,:)
      real, allocatable :: tcs(:,:), tcp(:,:), tcd(:,:), tct(:,:)
      real, allocatable :: atot(:,:), btot(:,:), ctot(:,:)
      character (len=20) :: tot
      real, allocatable :: orbitup(:,:),orbitdown(:,:)
      real :: Pi
      real, allocatable :: dotup(:,:), dotdown(:,:)
118  format(f14.7,3x,f14.7,3x,f14.7)
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

      open(30,file='PROORBMl')
      allocate (as(nk,nband,nion), ap(nk,nband,nion), ad(nk,nband,nion))
      allocate (at(nk,nband,nion), bs(nk,nband,nion), bp(nk,nband,nion))
      allocate (bd(nk,nband,nion), bt(nk,nband,nion), cs(nk,nband,nion))
      allocate (cp(nk,nband,nion), cd(nk,nband,nion), ct(nk,nband,nion))
      allocate (tas(nk,nband), tap(nk,nband), tad(nk,nband))
      allocate (tat(nk,nband), tbs(nk,nband), tbp(nk,nband))
      allocate (tbd(nk,nband), tbt(nk,nband), tcs(nk,nband))
      allocate (tcp(nk,nband), tcd(nk,nband), tct(nk,nband))

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
            read(30,*)n,as(i,j,k),ap(i,j,k),ad(i,j,k),at(i,j,k)
          enddo

          read(30,*)tot,tas(i,j),tap(i,j),tad(i,j),tat(i,j)
          
          do k = 1,nion
            read(30,*)n,bs(i,j,k),bp(i,j,k),bd(i,j,k),bt(i,j,k)
          enddo

          read(30,*)tot,tbs(i,j),tbp(i,j),tbd(i,j),tbt(i,j)

          do k = 1,nion
            read(30,*)n,cs(i,j,k),cp(i,j,k),cd(i,j,k),ct(i,j,k)
          enddo

          read(30,*)tot,tcs(i,j),tcp(i,j),tcd(i,j),tct(i,j)
        
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
            read(30,*)n,as(i,j,k),ap(i,j,k),ad(i,j,k),at(i,j,k)
          enddo

          read(30,*)tot,tas(i,j),tap(i,j),tad(i,j),tat(i,j)

          do k = 1,nion
            read(30,*)n,bs(i,j,k),bp(i,j,k),bd(i,j,k),bt(i,j,k)
          enddo

          read(30,*)tot,tbs(i,j),tbp(i,j),tbd(i,j),tbt(i,j)

          do k = 1,nion
            read(30,*)n,cs(i,j,k),cp(i,j,k),cd(i,j,k),ct(i,j,k)
          enddo
          
          read(30,*)tot,tcs(i,j),tcp(i,j),tcd(i,j),tct(i,j)

        
          read(30,*)
          read(30,*)
          read(30,*)
          read(30,*)
        enddo

      i = nk
      j = nband
          do k = 1,nion
            read(30,*)n,as(i,j,k),ap(i,j,k),ad(i,j,k),at(i,j,k)
          enddo

          read(30,*)tot,tas(i,j),tap(i,j),tad(i,j),tat(i,j)

          do k = 1,nion
            read(30,*)n,bs(i,j,k),bp(i,j,k),bd(i,j,k),bt(i,j,k)
          enddo

          read(30,*)tot,tbs(i,j),tbp(i,j),tbd(i,j),tbt(i,j)

          do k = 1,nion
            read(30,*)n,cs(i,j,k),cp(i,j,k),cd(i,j,k),ct(i,j,k)
          enddo

          read(30,*)tot,tcs(i,j),tcp(i,j),tcd(i,j),tct(i,j)

     
        allocate(atot(nk,nband),btot(nk,nband),ctot(nk,nband))

        do i = 1, nk
          do j = 1, nband
        if(tat(i,j)**2.GT.tbt(i,j)**2.and.tat(i,j)**2.GT.tct(i,j)**2)  then
              atot(i,j) = tat(i,j)
            endif
            
         if(tbt(i,j)**2.GT.tct(i,j)**2.and.tbt(i,j)**2.GT.tat(i,j)**2) then
              btot(i,j) = tbt(i,j)
            endif
            
         if(tct(i,j)**2.GT.tat(i,j)**2.and.tct(i,j)**2.GT.tbt(i,j)**2) then
              ctot(i,j) =tct(i,j)
                
            endif
          enddo
        enddo

        allocate(orbitup(nk,nband), orbitdown(nk,nband))
        allocate(dotup(nk,nband), dotdown(nk,nband))
        do i = 1, nk
          do j = 1, nband
            if(atot(i,j) .GT. 0) then
              orbitup(i,j) = feigup(i,j)
              dotup(i,j) = atot(i,j)
             endif
         
            if(atot(i,j) .LT. 0 ) then  
              orbitdown(i,j) = feigup(i,j)
              dotdown(i,j) = atot(i,j)
           endif
           
            if(btot(i,j) .GT. 0) then
              orbitup(i,j) = feigup(i,j)
              dotup(i,j) = btot(i,j)
           endif 
               
            if(btot(i,j) .LT. 0) then
              orbitdown(i,j) = feigup(i,j)
              dotdown(i,j) = btot(i,j)
            endif

            if(ctot(i,j) .GT. 0) then
              orbitup(i,j) = feigup(i,j)
              dotup(i,j) = ctot(i,j)
            endif 
               
            if(ctot(i,j) .LT. 0) then
              orbitdown(i,j) = feigup(i,j)
              dotdown(i,j) = ctot(i,j)
            endif
           
          enddo
        enddo

        do i = 1,nk
        do j = 1, nband
          if(orbitup(i,j).eq.0) then
        orbitup(i,j) = 1000000
          endif
          if(orbitdown(i,j).eq.0) then
        orbitdown(i,j) = 1000000
         endif
         if(dotup(i,j).eq.0) then
        dotup(i,j) = 1000000
         endif
         if(dotdown(i,j).eq.0) then
        dotdown(i,j) = 1000000
        endif
        enddo   
        enddo

     open(40,file='band2.txt')
 140  format(f14.7,3x,f14.7,3x,f14.7,3x,f14.7,3x,f14.7,3x,f14.7)
        write(40,*) '#kval, total'
      do j = 1, nband
      do i = 1, nk
       write(40,140)kval(i), feigup(i,j), orbitup(i,j), orbitdown(i,j),dotup(i,j),dotdown(i,j)
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

      
        open(24,file='band.txt')
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
      
