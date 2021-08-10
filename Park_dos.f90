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

!123456
      OPEN(1,file='band.dat')
      read(1,*) fermi
      read(1,*) soc

      open(15,file='EIGENVAL')
      read(15,*) c1, c2, c3, spin

      open(2,file='DOSCAR')
      read(2,*) nion, c1, c2, c3
      read(2,*)
      read(2,*)
      read(2,*)
      read(2,*)
      read(2,*) E_max, E_min, nedos, E_fermi, c1
      
      write(6,*) E_max
       If(spin.eq.1 .and. soc.eq.1) then
        
      allocate(energy(nedos), DOS(nedos),integrated_DOS(nedos))


       do i = 1,nedos
       
        read(2,*) energy(i), DOS(i), integrated_DOS(i)
       
         enddo

      open(4,file='result_total.out')
200   format(f14.10,3x,f14.10,3x,f14.10)

    
      do i = 1,nedos
      write(4,200) energy(i), DOS(i), integrated_DOS(i)
      enddo
   
      

       allocate(s(nion,nedos),p(nion,nedos),d(nion,nedos))        
      
       do i = 1, nion
       read(2,*)
       
       do j = 1, nedos
       read(2,*) energy(j), s(i,j), p(i,j), d(i,j)
       enddo
        
       enddo

      open(3,file='result.out')
100   format(f14.10,3x,f14.10,3x,f14.10,3x,f14.10)
     
      do j = 1, nion     
      do i = 1,nedos
      write(3,100) energy(i), s(j,i),p(j,i),d(j,i)
      enddo
      write(3,*)
      enddo

      endif  

      if(spin.eq.1 .and. soc.eq.2) then


      allocate(energy(nedos), DOS(nedos),integrated_DOS(nedos))


       do i = 1,nedos

        read(2,*) energy(i), DOS(i), integrated_DOS(i)

         enddo

      open(4,file='result_total.out')

      do i = 1,nedos
      write(4,200) energy(i), DOS(i), integrated_DOS(i)
      enddo 

       allocate(s_total(nion,nedos),s_x(nion,nedos),s_y(nion,nedos),s_z(nion,nedos))
       allocate(p_total(nion,nedos),p_x(nion,nedos),p_y(nion,nedos),p_z(nion,nedos))
       allocate(d_total(nion,nedos),d_x(nion,nedos),d_y(nion,nedos),d_z(nion,nedos))

       do i = 1, nion
       read(2,*)

       do j = 1, nedos
       read(2,*) energy(j), s_total(i,j), s_x(i,j),s_y(i,j),s_z(i,j),p_total(i,j),p_x(i,j),p_y(i,j),p_z(i,j),d_total(i,j),d_x(i,j),d_y(i,j),d_z(i,j)
       enddo

       enddo

      open(3,file='result.out')
101   format(f14.10,3x,f14.10,3x,f14.10,3x,f14.10,3x,f14.10,3x,f14.10,3x,f14.10,3x,f14.10,3x,f14.10,3x,f14.10,3x,f14.10,3x,f14.10,3x,f14.10)

      do i = 1, nion
      do j = 1,nedos
      write(3,101) energy(j), s_total(i,j),s_x(i,j),s_y(i,j),s_z(i,j),p_total(i,j),p_x(i,j),p_y(i,j),p_z(i,j),d_total(i,j),d_x(i,j),d_y(i,j),d_z(i,j)
      enddo
      write(3,*)
      enddo



      endif

       If(spin.eq.2) then

      allocate(energy(nedos), DOS_up(nedos),integrated_DOS_up(nedos))
      allocate(DOS_dwn(nedos),integrated_DOS_dwn(nedos))


       do i = 1,nedos

        read(2,*) energy(i), DOS_up(i),DOS_dwn(i),integrated_DOS_up(i),integrated_DOS_dwn(i)

         enddo

      open(4,file='result_total.out')
201   format(f14.10,3x,f14.10,3x,f14.10,3x,f14.10,3x,f14.10)
      do i = 1,nedos
     
     write(4,201) energy(i),DOS_up(i),DOS_dwn(i),integrated_DOS_up(i),integrated_DOS_dwn(i)
      enddo


       allocate(s_up(nion,nedos),p_up(nion,nedos),d_up(nion,nedos))
       allocate(s_dwn(nion,nedos),p_dwn(nion,nedos),d_dwn(nion,nedos))

       do i = 1, nion
       read(2,*)

       do j = 1, nedos
       read(2,*) energy(j), s_up(i,j), s_dwn(i,j),p_up(i,j),p_dwn(i,j), d_up(i,j),d_dwn(i,j)
       enddo

       enddo

      open(3,file='result.out')
102   format(f14.10,3x,f14.10,3x,f14.10,3x,f14.10,3x,f14.10,3x,f14.10,3x,f14.10)

      do i = 1, nion
      do j = 1,nedos
      write(3,102) energy(j), s_up(i,j),s_dwn(i,j),p_up(i,j),p_dwn(i,j),d_up(i,j),d_dwn(i,j)
      enddo
      write(3,*)
      enddo

      endif


      stop

      end
