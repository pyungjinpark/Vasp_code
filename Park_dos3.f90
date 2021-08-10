!123456     
      real :: E_max, E_min
      real :: fermi
      integer :: nedos,spin, nion, soc
      real, allocatable :: s(:,:), p(:,:), d(:,:)
      real, allocatable :: s_up(:,:), s_dwn(:,:), p_up(:,:), p_dwn(:,:)
      real, allocatable :: d_up(:,:), d_dwn(:,:), s_total(:,:)
      real, allocatable :: p_total(:,:), d_total(:,:), DOS(:)
      real, allocatable :: integrated_DOS(:), DOS_up(:)
      real, allocatable :: DOS_dwn(:), integrated_DOS_up(:)
      real, allocatable :: integrated_DOS_dwn(:)
      real, allocatable :: s_x(:,:), s_y(:,:), s_z(:,:)
      real, allocatable :: p_x(:,:), p_y(:,:), p_z(:,:)
      real, allocatable :: d_x(:,:), d_y(:,:), d_z(:,:), energy(:)
      integer :: atom1, atom2

      real, allocatable :: s_sum(:,:), p_sum(:,:), d_sum(:,:)
      real, allocatable :: s_sum_up(:,:), s_sum_dwn(:,:), p_sum_up(:,:),p_sum_dwn(:,:)
      real, allocatable :: d_sum_up(:,:), d_sum_dwn(:,:)
      real, allocatable :: s_sum_x(:,:), s_sum_y(:,:), s_sum_z(:,:)
      real, allocatable :: p_sum_x(:,:), p_sum_y(:,:), p_sum_z(:,:)
      real, allocatable :: d_sum_x(:,:), d_sum_y(:,:), d_sum_z(:,:)
      real, allocatable :: s_x_p(:,:), s_x_n(:,:), s_y_p(:,:)
      real, allocatable :: s_y_n(:,:), s_z_p(:,:), s_z_n(:,:), p_x_p(:,:), p_x_n(:,:)
      real, allocatable :: p_y_p(:,:), p_y_n(:,:), p_z_p(:,:)
      real, allocatable :: p_z_n(:,:), d_x_p(:,:), d_x_n(:,:)
      real, allocatable :: d_y_p(:,:), d_y_n(:,:), d_z_p(:,:), d_z_n(:,:)
      real, allocatable :: x(:,:), y(:,:), z(:,:)

110   format(f14.10,3x,f14.10)
120   format(f14.10,3x,f14.10,3x,f14.10,3x,f14.10,3x,f14.10)
200   format(f14.10,3x,f14.10,3x,f14.10)
300   format(f14.10,3x,f14.10,3x,f14.10,3x,f14.10,3x,f14.10,3x,f14.10,3x,f14.10)
400   format(f14.10,3x,f14.10,3x,f14.10,3x,f14.10)

!123456
      OPEN(1,file='band.dat')
      read(1,*) fermi
      read(1,*) soc
      read(1,*) atom1, atom2

      open(15,file='EIGENVAL')
      read(15,*) c1, c2, c3, spin

      open(2,file='DOSCAR')
      read(2,*) nion, c1, c2, c3
      read(2,*)
      read(2,*)
      read(2,*)
      read(2,*)
      read(2,*) E_max, E_min, nedos, E_fermi, c1


!no spin no soc

      
      write(6,*) E_max
       If(spin.eq.1 .and. soc.eq.1) then
        
      allocate(energy(nedos), DOS(nedos),integrated_DOS(nedos))


       do i = 1,nedos
       
        read(2,*) energy(i), DOS(i), integrated_DOS(i)
       
         enddo

      open(4,file='result_total.out')

       write(4,*) 'energy,DOS,integrated_DOS'
      do i = 1,nedos
      write(4,200) energy(i), DOS(i), integrated_DOS(i)
      enddo
   
      

       allocate(s(nion,nedos),p(nion,nedos),d(nion,nedos))        
      
       do i = 1, atom1
       read(2,*)
       
       do j = 1, nedos
       read(2,*) energy(j), s(i,j), p(i,j), d(i,j)
       enddo
        
       enddo

   

       do i = atom1+1, nion
       read(2,*)
       
       do j = 1, nedos
       read(2,*) energy(j), s(i,j), p(i,j), d(i,j)
       enddo
        
       enddo

       


     allocate(s_sum(2,nedos), p_sum(2,nedos), d_sum(2,nedos))

       s_sum(1,1) = 0.0
       p_sum(1,1) = 0.0
       d_sum(1,1) = 0.0


       do i = 1, nedos
       do j = 1, atom1
       s_sum(1,i) = s_sum(1,i) + s(j,i)
       enddo
       enddo

       do i = 1, nedos
       do j = atom1+1, nion
       s_sum(2,i) = s_sum(2,i) + s(j,i)
       enddo
       enddo

       do i = 1, nedos
       do j = 1, atom1
       p_sum(1,i) = p_sum(1,i) + p(j,i)
       enddo
       enddo

       do i = 1, nedos
       do j = atom1+1, nion
       p_sum(2,i) = p_sum(2,i) + p(j,i)
       enddo
       enddo

       do i = 1, nedos
       do j = 1, atom1
       d_sum(1,i) = d_sum(1,i) + d(j,i)
       enddo
       enddo

       do i = 1, nedos
       do j = atom1+1, nion
       d_sum(2,i) = d_sum(2,i) + d(j,i)
       enddo
       enddo

       open(100,file='result_atom.txt')
       write(100,*) 'energy,atom1,atom2'
       do i = 1, nedos
       write(100,200)energy(i),s_sum(1,i)+p_sum(1,i)+d_sum(1,i),s_sum(2,i)+p_sum(2,i)+d_sum(2,i)
       enddo
 
       open(7,file='result_orbital')
      write(7,*) 'energy,s,p,d'
       do i = 1, nedos
        write(7,400) energy(i), s_sum(1,i)+s_sum(2,i),p_sum(1,i)+p_sum(2,i),d_sum(1,i)+d_sum(2,i)
        enddo

       open(10,file='s.out')
     write(10,*) 'energy,s_atom1,s_atom2'
      do i = 1, nedos
      write(10,200) energy(i), s_sum(1,i), s_sum(2,i)
      enddo

       open(11,file='p.out')
       write(11,*) 'energy,p_atom1,p_atom2'
      do i = 1, nedos
      write(11,200) energy(i),p_sum(1,i),p_sum(2,i)
      enddo

       open(12,file='d.out')
      write(12,*) 'energy,d_atom1,d_atom2'
      do i = 1, nedos
   
      write(12,200) energy(i),d_sum(1,i),d_sum(2,i)
      enddo





      endif  


!no spin with soc




      if(spin.eq.1 .and. soc.eq.2) then


      allocate(energy(nedos), DOS(nedos),integrated_DOS(nedos))


       do i = 1,nedos

        read(2,*) energy(i), DOS(i), integrated_DOS(i)

         enddo

      open(4,file='result_total.out')
      write(4,*) 'energy,DOS,integrated_DOS'
      do i = 1,nedos
  
      write(4,200) energy(i), DOS(i), integrated_DOS(i)
      enddo 

       allocate(s_total(nion,nedos),s_x(nion,nedos),s_y(nion,nedos),s_z(nion,nedos))
       allocate(p_total(nion,nedos),p_x(nion,nedos),p_y(nion,nedos),p_z(nion,nedos))
       allocate(d_total(nion,nedos),d_x(nion,nedos),d_y(nion,nedos),d_z(nion,nedos))

       do i = 1, atom1
       read(2,*)

       do j = 1, nedos
       read(2,*) energy(j), s_total(i,j), s_x(i,j),s_y(i,j),s_z(i,j),p_total(i,j),p_x(i,j),p_y(i,j),p_z(i,j),d_total(i,j),d_x(i,j),d_y(i,j),d_z(i,j)
       enddo

       enddo

       do i = atom1+1, nion
       read(2,*)

       do j = 1, nedos
       read(2,*) energy(j), s_total(i,j),s_x(i,j),s_y(i,j),s_z(i,j),p_total(i,j),p_x(i,j),p_y(i,j),p_z(i,j),d_total(i,j),d_x(i,j),d_y(i,j),d_z(i,j)
       enddo

       enddo



      allocate(s_sum_x(2,nedos), s_sum_y(2,nedos), s_sum_z(2,nedos))
      allocate(p_sum_x(2,nedos), p_sum_y(2,nedos), p_sum_z(2,nedos))
      allocate(d_sum_x(2,nedos), d_sum_y(2,nedos), d_sum_z(2,nedos))
      allocate(s_sum(2,nedos), p_sum(2,nedos), d_sum(2,nedos))

       do i = 1, nedos
       do j = 1, atom1
       s_sum(1,i) = s_sum(1,i) + s_total(j,i)
       enddo
       enddo

       do i = 1, nedos
       do j = atom1+1, nion
       s_sum(2,i) = s_sum(2,i) + s_total(j,i)
       enddo
       enddo

       do i = 1, nedos
       do j = 1, atom1
       p_sum(1,i) = p_sum(1,i) + p_total(j,i)
       enddo
       enddo

       do i = 1, nedos
       do j = atom1+1, nion
       p_sum(2,i) = p_sum(2,i) + p_total(j,i)
       enddo
       enddo

       do i = 1, nedos
       do j = 1, atom1
       d_sum(1,i) = d_sum(1,i) + d_total(j,i)
       enddo
       enddo

       do i = 1, nedos
       do j = atom1+1, nion
       d_sum(2,i) = d_sum(2,i) + d_total(j,i)
       enddo
       enddo


       do  i = 1, nedos
       do j = 1, atom1
       s_sum_x(1,i) = s_sum_x(1,i) + s_x(j,i)
       s_sum_y(1,i) = s_sum_y(1,i) + s_y(i,j)
       s_sum_z(1,i) = s_sum_z(1,i) + s_z(i,j)
       enddo
       enddo

       do i = 1, nedos
       do j = atom1+1, nion
       s_sum_x(2,i) = s_sum_x(2,i) + s_x(j,i)
       s_sum_y(2,i) = s_sum_y(2,i) + s_y(j,i)
       s_sum_z(2,i) = s_sum_z(2,i) + s_z(j,i)

       enddo
       enddo

       do i = 1, nedos
       do j = 1, atom1
       p_sum_x(1,i) = p_sum_x(1,i) + p_x(j,i)
       p_sum_y(1,i) = p_sum_y(1,i) + p_y(i,j)
       p_sum_z(1,i) = p_sum_z(1,i) + p_z(i,j)
       enddo
       enddo

       do i = 1, nedos
       do j = atom1+1, nion
       p_sum_x(2,i) = p_sum_x(2,i) + p_x(j,i)
       p_sum_y(2,i) = p_sum_y(2,i) + p_y(j,i)
       p_sum_z(2,i) = p_sum_z(2,i) + p_z(j,i)

       enddo
       enddo

       do i = 1, nedos
       do j = 1, atom1
       d_sum_x(1,i) = d_sum_x(1,i) + d_x(j,i)
       d_sum_y(1,i) = d_sum_y(1,i) + d_y(j,i)
       d_sum_z(1,i) = d_sum_z(1,i) + d_z(j,i)
       enddo
       enddo

       do i = 1, nedos
       do j = atom1+1, nion
       d_sum_x(2,i) = d_sum_x(2,i) + d_x(j,i)
       d_sum_y(2,i) = d_sum_y(2,i) + d_y(j,i)
       d_sum_z(2,i) = d_sum_z(2,i) + d_z(j,i)

       enddo
       enddo

       open(100,file='result_atom.txt')
       write(100,*) 'energy,atom1,atom2'
       do i = 1, nedos
       write(100,200)energy(i),s_sum(1,i)+p_sum(1,i)+d_sum(1,i),s_sum(2,i)+p_sum(2,i)+d_sum(2,i)
       enddo


       open(7,file='result_orbital')
      write(7,*) 'energy,s,p,d'
       do i = 1, nedos
        write(7,400) energy(i),s_sum(1,i)+s_sum(2,i),p_sum(1,i)+p_sum(2,i),d_sum(1,i)+d_sum(2,i)
        enddo


       open(10,file='s.out')
        write(10,*) 'E,s1_x,s1_y,s1_z,s2_x,s2_y,s2_z'
      do i = 1, nedos
      write(10,300) energy(i), s_sum_x(1,i),s_sum_y(1,i),s_sum_z(1,i),s_sum_x(2,i), s_sum_y(2,i),s_sum_z(2,i)

      enddo

       open(11,file='p.out')
      write(11,*) 'E,p1_x,p1_y,p1_z,p2_x,p2_y,p2_z'
      do i = 1, nedos

      write(11,300) energy(i),p_sum_x(1,i),p_sum_y(1,i),p_sum_z(1,i),p_sum_x(2,i),p_sum_y(2,i),p_sum_z(2,i) 
      enddo

       open(12,file='d.out')
      write(12,*) 'E,d1_x,d1_y,d1_z,d2_x,d2_y,d2_z'
      do i = 1, nedos

      write(12,300) energy(i),d_sum_x(1,i),d_sum_y(1,i),d_sum_z(1,i),d_sum_x(2,i),d_sum_y(2,i),d_sum_z(2,i)
      enddo


      allocate(x(2,nedos), y(2,nedos), z(2,nedos))
      allocate(s_x_p(2,nedos), s_x_n(2,nedos), s_y_p(2,nedos),s_y_n(2,nedos))
      allocate(s_z_p(2,nedos), s_z_n(2,nedos),p_x_p(2,nedos),p_x_n(2,nedos))
      allocate(p_y_p(2,nedos),p_y_n(2,nedos),p_z_p(2,nedos),p_z_n(2,nedos))
      allocate(d_x_p(2,nedos),d_x_n(2,nedos),d_y_p(2,nedos),d_y_n(2,nedos))
      allocate(d_z_p(2,nedos),d_z_n(2,nedos))

      do i = 1,2
      do j = 1, nedos
        x(i,j) = 1.0
        y(i,j) = 1.0
        z(i,j) = 1.0
        s_x_p(i,j) = 0.0
        s_x_n(i,j) = 0.0
        s_y_p(i,j) = 0.0
        s_y_n(i,j) = 0.0
        s_z_p(i,j) = 0.0
        s_z_n(i,j) = 0.0
        p_x_p(i,j) = 0.0
        p_x_n(i,j) = 0.0
        p_y_p(i,j) = 0.0
        p_y_n(i,j) = 0.0
        p_z_p(i,j) = 0.0
        p_z_n(i,j) = 0.0
        d_x_p(i,j) = 0.0
        d_x_n(i,j) = 0.0
        d_y_p(i,j) = 0.0
        d_y_n(i,j) = 0.0
        d_z_p(i,j) = 0.0
        d_z_n(i,j) = 0.0
      enddo
      enddo
      
      do i = 1,2
      do j = 1, nedos
       x(i,j) = x(i,j) * s_sum_x(i,j)
       y(i,j) = y(i,j) * s_sum_y(i,j)
       z(i,j) = z(i,j) * s_sum_z(i,j)
      enddo
      enddo

      do i = 1,2
      do j = 1, nedos
        if( x(i,j) .GT. 0) then
          s_x_p(i,j) = x(i,j)
        endif

        if(x(i,j) .LT. 0) then
          s_x_n(i,j) = x(i,j)
        endif

        if(y(i,j) .GT. 0) then
          s_y_p(i,j) = y(i,j)
        endif

        if(y(i,j) .LT. 0 ) then
          s_y_n(i,j) = y(i,j)
        endif

        if(z(i,j) .GT. 0 ) then
          s_z_p(i,j) = z(i,j)
        endif

        if(z(i,j) .LT. 0 ) then
          s_z_n(i,j) = z(i,j)
        endif

       enddo
       enddo


      do i = 1,2
      do j = 1, nedos
        x(i,j) = 1.0
        y(i,j) = 1.0
        z(i,j) = 1.0
      enddo
      enddo

      do i = 1,2
      do j = 1, nedos
       x(i,j) = x(i,j) * p_sum_x(i,j)
       y(i,j) = y(i,j) * p_sum_y(i,j)
       z(i,j) = z(i,j) * p_sum_z(i,j)
      enddo
      enddo

      do i = 1,2
      do j = 1, nedos
        if( x(i,j) .GT. 0) then
          p_x_p(i,j) = x(i,j)
        endif

        if(x(i,j) .LT. 0) then
          p_x_n(i,j) = x(i,j)
        endif

        if(y(i,j) .GT. 0) then
          p_y_p(i,j) = y(i,j)
        endif

        if(y(i,j) .LT. 0 ) then
          p_y_n(i,j) = y(i,j)
        endif

        if(z(i,j) .GT. 0 ) then
          p_z_p(i,j) = z(i,j)
        endif

        if(z(i,j) .LT. 0 ) then
          p_z_n(i,j) = z(i,j)
        endif

       enddo
       enddo

      do i = 1,2
      do j = 1, nedos
        x(i,j) = 1.0
        y(i,j) = 1.0
        z(i,j) = 1.0
      enddo
      enddo

      do i = 1,2
      do j = 1, nedos
       x(i,j) = x(i,j) * d_sum_x(i,j)
       y(i,j) = y(i,j) * d_sum_y(i,j)
       z(i,j) = z(i,j) * d_sum_z(i,j)
      enddo
      enddo

      do i = 1,2
      do j = 1, nedos
        if( x(i,j) .GT. 0) then
          d_x_p(i,j) = x(i,j)
        endif
       
        if(x(i,j) .LT. 0) then
          d_x_n(i,j) = x(i,j)
        endif

        if(y(i,j) .GT. 0) then
          d_y_p(i,j) = y(i,j)
        endif
        
        if(y(i,j) .LT. 0 ) then
          d_y_n(i,j) = y(i,j)
        endif

        if(z(i,j) .GT. 0 ) then
          d_z_p(i,j) = z(i,j)
        endif

        if(z(i,j) .LT. 0 ) then
          d_z_n(i,j) = z(i,j)
        endif

       enddo
       enddo

       do i = 1, 2
       do j = 1, nedos
        
        if(s_x_p(i,j) .eq. 0 ) then
          s_x_p(i,j) = 1000000000000000000000.0
        endif

        if(s_x_n(i,j) .eq. 0 ) then
          s_x_n(i,j) = 1000000000000000000000.0
        endif

        if(s_y_p(i,j) .eq. 0 ) then
          s_y_p(i,j) = 1000000000000000000000.0
        endif
        
        if(s_y_n(i,j) .eq. 0 ) then
          s_y_n(i,j) = 1000000000000000000000.0
        endif

        if(s_z_p(i,j) .eq. 0 ) then
          s_z_p(i,j) = 1000000000000000000000.0
        endif
        
        if(s_z_n(i,j) .eq. 0 ) then
          s_z_n(i,j) = 1000000000000000000000.0
        endif


        if(p_x_p(i,j) .eq. 0 ) then
          p_x_p(i,j) = 1000000000000000000000.0
        endif
        
        if(p_x_n(i,j) .eq. 0 ) then
          p_x_n(i,j) = 1000000000000000000000.0
        endif

        if(p_y_p(i,j) .eq. 0 ) then
          p_y_p(i,j) = 1000000000000000000000.0
        endif

        if(p_y_n(i,j) .eq. 0 ) then
          p_y_n(i,j) = 1000000000000000000000.0
        endif

        if(p_z_p(i,j) .eq. 0 ) then
          p_z_p(i,j) = 1000000000000000000000.0
        endif

        if(p_z_n(i,j) .eq. 0 ) then
          p_z_n(i,j) = 1000000000000000000000.0
        endif


        if(d_x_p(i,j) .eq. 0 ) then
          d_x_p(i,j) = 1000000000000000000000.0
        endif
        
        if(d_x_n(i,j) .eq. 0 ) then
          d_x_n(i,j) = 1000000000000000000000.0
        endif

        if(d_y_p(i,j) .eq. 0 ) then
          d_y_p(i,j) = 1000000000000000000000.0
        endif

        if(d_y_n(i,j) .eq. 0 ) then
          d_y_n(i,j) = 1000000000000000000000.0
        endif

        if(d_z_p(i,j) .eq. 0 ) then
          d_z_p(i,j) = 1000000000000000000000.0
        endif

        if(d_z_n(i,j) .eq. 0 ) then
          d_z_n(i,j) = 1000000000000000000000.0
        endif

        enddo 
        enddo

       open(20,file='s_x.out')
     write(20,*)'energy,s_x_p_atom1,s_x_n_atom1,s_x_p_atom2,s_x_n_atom2'
      do i = 1, nedos
      write(20,120) energy(i),s_x_p(1,i),s_x_n(1,i),s_x_p(2,i),s_x_n(2,i)

      enddo

       open(21,file='s_y.out')
     write(21,*)'energy,s_y_p_atom1,s_y_n_atom1,s_y_p_atom2,s_y_n_atom2'
      do i = 1, nedos

      write(21,120)energy(i),s_y_p(1,i),s_y_n(1,i),s_y_p(2,i),s_y_n(2,i)

      enddo

       open(22,file='s_z.out')
     write(22,*)'energy,s_z_p_atom1,s_z_n_atom1,s_z_p_atom2,s_z_n_atom2'
      do i = 1, nedos

      write(22,120)energy(i),s_z_p(1,i),s_z_n(1,i),s_z_p(2,i),s_z_n(2,i)

      enddo


       open(30,file='p_x.out')
     write(30,*)'energy,p_x_p_atom1,p_x_n_atom1,p_x_p_atom2,p_x_n_atom2'
      do i = 1, nedos

      write(30,120)energy(i),p_x_p(1,i),p_x_n(1,i),p_x_p(2,i),s_x_n(2,i)

      enddo

       open(31,file='p_y.out')
     write(31,*)'energy,p_y_p_atom1,p_y_n_atom1,p_y_p_atom2,p_y_n_atom2'
      do i = 1, nedos

      write(31,120)energy(i),p_y_p(1,i),p_y_n(1,i),p_y_p(2,i),p_y_n(2,i)

      enddo

       open(32,file='p_z.out')
     write(32,*)'energy,p_z_p_atom1,p_z_n_atom1,p_z_p_atom2,p_z_n_atom2'
      do i = 1, nedos

      write(32,120)energy(i),p_z_p(1,i),p_z_n(1,i),p_z_p(2,i),p_z_n(2,i)

      enddo



       open(40,file='d_x.out')
     write(40,*)'energy,d_x_p_atom1,d_x_n_atom1,d_x_p_atom2,d_x_n_atom2'
      do i = 1, nedos

      write(40,120)energy(i),d_x_p(1,i),d_x_n(1,i),d_x_p(2,i),d_x_n(2,i)

      enddo

       open(41,file='d_y.out')
     write(41,*)'energy,d_y_p_atom1,d_y_n_atom1,d_y_p_atom2,d_y_n_atom2'
      do i = 1, nedos

      write(41,120)energy(i),d_y_p(1,i),d_y_n(1,i),d_y_p(2,i),d_y_n(2,i)

      enddo

       open(42,file='d_z.out')
     write(42,*)'energy,d_z_p_atom1,d_z_n_atom1,d_z_p_atom2,d_z_n_atom2'
      do i = 1, nedos

      write(42,120)energy(i),d_z_p(1,i),d_z_n(1,i),d_z_p(2,i),d_z_n(2,i)

      enddo


      endif


! spin 



       If(spin.eq.2) then

      allocate(energy(nedos), DOS_up(nedos),integrated_DOS_up(nedos))
      allocate(DOS_dwn(nedos),integrated_DOS_dwn(nedos))


       do i = 1,nedos

        read(2,*) energy(i), DOS_up(i),DOS_dwn(i),integrated_DOS_up(i),integrated_DOS_dwn(i)

         enddo

      open(4,file='result_total.out')
201   format(f14.10,3x,f14.10,3x,f14.10,3x,f14.10,3x,f14.10)
     write(4,*)'E,DOS_up,DOS_dwn,integrated_DOS_up,integrated_DOS_dwn'
      do i = 1,nedos
     write(4,201) energy(i),DOS_up(i),DOS_dwn(i),integrated_DOS_up(i),integrated_DOS_dwn(i)
      enddo


       allocate(s_up(nion,nedos),p_up(nion,nedos),d_up(nion,nedos))
       allocate(s_dwn(nion,nedos),p_dwn(nion,nedos),d_dwn(nion,nedos))

       do i = 1, atom1
       read(2,*)

       do j = 1, nedos
       read(2,*) energy(j), s_up(i,j), s_dwn(i,j),p_up(i,j),p_dwn(i,j), d_up(i,j),d_dwn(i,j)
       enddo

       enddo


       do i = atom1+1,nion
       read(2,*)

       do j = 1, nedos
       read(2,*) energy(j), s_up(i,j), s_dwn(i,j),p_up(i,j),p_dwn(i,j),d_up(i,j),d_dwn(i,j)
       enddo

       enddo


      allocate(s_sum_up(2,nedos),p_sum_up(2,nedos),d_sum_up(2,nedos))
      allocate(s_sum_dwn(2,nedos),p_sum_dwn(2,nedos),d_sum_dwn(2,nedos))


       s_sum_up(1,1) = 0.0
       p_sum_up(1,1) = 0.0
       d_sum_up(1,1) = 0.0

       s_sum_dwn(1,1) = 0.0
       p_sum_dwn(1,1) = 0.0
       d_sum_dwn(1,1) = 0.0

       do i = 1, nedos
       do j = 1, atom1
       s_sum_up(1,i) = s_sum_up(1,i) + s_up(j,i)
       enddo
       enddo

       do i = 1, nedos
       do j = 1, atom1
       s_sum_dwn(1,i) = s_sum_dwn(1,i) + s_dwn(j,i)
       enddo
       enddo

       do i = 1, nedos
       do j = 1, atom1
       p_sum_up(1,i) = p_sum_up(1,i) + p_up(j,i)
       enddo
       enddo

       do i = 1, nedos
       do j = 1, atom1
       p_sum_dwn(1,i) = p_sum_dwn(1,i) + p_dwn(j,i)
       enddo
       enddo

       do i = 1, nedos
       do j = 1, atom1
       d_sum_up(1,i) = d_sum_up(1,i) + d_up(j,i)
       enddo
       enddo

       do i = 1, nedos
       do j = 1, atom1
       d_sum_dwn(1,i) = d_sum_dwn(1,i) + d_dwn(j,i)
       enddo
       enddo




       do i = 1, nedos
       do j = atom1+1, nion
       s_sum_up(2,i) = s_sum_up(2,i) + s_up(j,i)
       enddo
       enddo

       do i = 1, nedos
       do j = atom+1,nion
       s_sum_dwn(2,i) = s_sum_dwn(2,i) + s_dwn(j,i)
       enddo
       enddo

       do i = 1, nedos
       do j = atom1+1, nion
       p_sum_up(2,i) = p_sum_up(2,i) + p_up(j,i)
       enddo
       enddo

       do i = 1, nedos
       do j = atom1+1, nion
       p_sum_dwn(2,i) = p_sum_dwn(2,i) + p_dwn(j,i)
       enddo
       enddo

       do i = 1, nedos
       do j = atom1+1,nion
       d_sum_up(2,i) = d_sum_up(2,i) + d_up(j,i)
       enddo
       enddo

       do i = 1, nedos
       do j = atom1+1, nion
       d_sum_dwn(2,i) = d_sum_dwn(2,i) + d_dwn(j,i)
       enddo
       enddo

       open(100,file='result_atom.txt')
       write(100,*) 'energy,atom1_up,atom1_dwn,atom2_up,atom2_dwn'
       do i = 1, nedos
      write(100,120)energy(i),s_sum_up(1,i)+p_sum_up(1,i)+d_sum_up(1,i)&
             , s_sum_dwn(1,i)+p_sum_dwn(1,i)+d_sum_dwn(1,i),&
                    s_sum_up(2,i)+p_sum_up(2,i)+d_sum_up(2,i),&
                    s_sum_dwn(2,i)+p_sum_dwn(2,i)+d_sum_dwn(2,i)
       enddo


       open(7,file='result_orbital')
      write(7,*) 'energy,s_up,s_dwn,p_up,p_dwn,d_up,d_dwn'
       do i = 1, nedos
      write(7,300) energy(i),s_sum_up(1,i)+s_sum_up(2,i),s_sum_dwn(1,i)&
                   +s_sum_dwn(2,i),p_sum_up(1,i)+p_sum_up(2,i)&
                ,p_sum_dwn(1,i)+p_sum_dwn(2,i),&
                d_sum_up(1,i)+d_sum_up(2,i),d_sum_dwn(1,i)+d_sum_dwn(2,i)
        enddo


       open(10,file='s.out')
      write(10,*)'energy,s_up_atom1,s_dwn_aton1,s_up_atom2,s_dwn_atom2'
      do i = 1, nedos
      write(10,120) energy(i), s_sum_up(1,i), s_sum_dwn(1,i),s_sum_up(2,i), s_sum_dwn(2,i)
      enddo

       open(11,file='p.out')
      write(11,*)'energy,p_up_atom1,p_dwn_aton1,p_up_atom2,p_dwn_atom2'
      do i = 1, nedos
      write(11,120) energy(i), p_sum_up(1,i),p_sum_dwn(1,i),p_sum_up(2,i), p_sum_dwn(2,i)
      enddo

       open(12,file='d.out')                                         
      write(12,*)'energy,d_up_atom1,d_dwn_aton1,d_up_atom2,d_dwn_atom2'
      do i = 1, nedos
      write(12,120) energy(i), d_sum_up(1,i),d_sum_dwn(1,i),d_sum_up(2,i), d_sum_dwn(2,i)
      enddo


      endif


      stop

      end
