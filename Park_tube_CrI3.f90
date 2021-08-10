    
        MODULE functions
                IMPLICIT NONE

          REAL, PARAMETER :: PI = 4.0 * ATAN(1.0)
          REAL, PARAMETER :: Degree180 = 180.0
          REAL, PARAMETER :: R_to_D = Degree180/PI
          REAL, PARAMETER :: D_to_R = PI/Degree180


        CONTAINS

          REAL FUNCTION RadianToDegree(Radian)
                IMPLICIT NONE
                real, intent(in) :: Radian

                RadianToDegree = Radian * R_to_D
          END FUNCTION RadianToDegree

          REAL FUNCTION DegreeToRadian(Degree)
                IMPLICIT NONE
                real, intent(in) :: Degree

                DegreeToRadian = Degree * D_to_R
          END FUNCTION DegreeToRadian


          REAL FUNCTION MySIN(x)
                IMPLICIT NONE
                REAL, intent(in) :: x

                MySIN = SIN(DegreeToRadian(x))
          END FUNCTION MySIN

          REAL FUNCTION MyCOS(x)
                IMPLICIT NONE
                REAL, INTENT(in) :: x

                MyCOS = COS(DegreeToRadian(x))
          END FUNCTION MyCOS

        END MODULE functions

          PROGRAM  main
          USE  functions        ! use a module

          IMPLICIT  NONE
         
          integer :: n1, n2, i, j, k
          real :: bind_length, total_length, lattice, total_high
          real :: radis, angle, a1, b1, c1, d1
          real :: a2, b2, c2, d2, e2,f2
          real, allocatable :: x(:,:), y(:,:), xatom(:,:,:)
          real, allocatable :: yatom(:,:,:), high(:,:,:)
          real, allocatable :: fxatom(:,:,:), fyatom(:,:,:)
          real, allocatable :: fhigh(:,:,:)
          real :: total_cell

          OPEN(1,file='tube.dat')
          
          READ(1,*)
          READ(1,*) n1, n2
          READ(1,*) bind_length

         !armchair
        
          if(n2 .NE. 0) THEN
                radis = (6*n1*((bind_length/2)*SQRT(3.0)))/(2*PI)
                lattice = radis + (radis/2)
                total_high = (bind_length/2)*6
                total_cell = lattice*2
          angle = 360.0/n1

        if(total_cell .GE. 50) THEN
                total_cell = 49
        endif
        

         allocate(x(2,12), y(2,12), xatom(2,12,n1), yatom(2,12,n1))
         allocate(high(2,12,n1))

         x(1,1) = radis
         y(1,1) = 0
       
         x(1,2) = radis*MyCOS(360.0/(n1*6)) - 0*MySIN(360.0/(n1*6))
         y(1,2) = radis*MySIN(360.0/(n1*6)) + 0*MySIN(360.0/(n1*6))

         x(1,3) = radis*MyCOS((360.0/(n1*6)*3)) - 0*MySIN((360.0/(n1*6))*3)
         y(1,3) = radis*MySIN((360.0/(n1*6)*3)) + 0*MyCOS((360.0/(n1*6))*3)

         x(1,4) = radis*MyCOS((360.0/(n1*6)*4)) -0*MySIN((360.0/(n1*6))*4)
         y(1,4) = radis*MySIN((360.0/(n1*6)*4)) +0*MyCOS((360.0/(n1*6))*4)

         a1 = radis - 1.56377985
         b1 = radis + 1.56377985

         x(2,1) = a1
         y(2,1) = 0

         x(2,2) = b1
         y(2,2) = 0

         x(2,3) = b1*MyCOS(360.0/(n1*6)) - 0*MySIN(360.0/(n1*6))
         y(2,3) = b1*MySIN(360.0/(n1*6)) + 0*MySIN(360.0/(n1*6))

         x(2,4) = a1*MyCOS((360.0/(n1*6))) -0*MySIN((360.0/(n1*6)))
         y(2,4) = a1*MySIN((360.0/(n1*6))) +0*MyCOS((360.0/(n1*6)))

         x(2,5) = b1*MyCOS((360.0/(n1*6))*2) -0*MySIN((360.0/(n1*6))*2)
         y(2,5) = b1*MySIN((360.0/(n1*6))*2) +0*MyCOS((360.0/(n1*6))*2)

         x(2,6) = a1*MyCOS((360.0/(n1*6))*2) -0*MySIN((360.0/(n1*6))*2)
         y(2,6) = a1*MySIN((360.0/(n1*6))*2) +0*MyCOS((360.0/(n1*6))*2)

         x(2,7) = b1*MyCOS((360.0/(n1*6))*3) -0*MySIN((360.0/(n1*6))*3)
         y(2,7) = b1*MySIN((360.0/(n1*6))*3) +0*MyCOS((360.0/(n1*6))*3)

         x(2,8) = a1*MyCOS((360.0/(n1*6))*3) -0*MySIN((360.0/(n1*6))*3)
         y(2,8) = a1*MySIN((360.0/(n1*6))*3) +0*MyCOS((360.0/(n1*6))*3)

         x(2,9) = a1*MyCOS((360.0/(n1*6))*4) -0*MySIN((360.0/(n1*6))*4)
         y(2,9) = a1*MySIN((360.0/(n1*6))*4) +0*MyCOS((360.0/(n1*6))*4)

         x(2,10) = b1*MyCOS((360.0/(n1*6))*4) -0*MySIN((360.0/(n1*6))*4)
         y(2,10) = b1*MySIN((360.0/(n1*6))*4) +0*MyCOS((360.0/(n1*6))*4)

         x(2,11) = b1*MyCOS((360.0/(n1*6))*5) -0*MySIN((360.0/(n1*6))*5)
         y(2,11) = b1*MySIN((360.0/(n1*6))*5) +0*MyCOS((360.0/(n1*6))*5)

         x(2,12) = a1*MyCOS((360.0/(n1*6))*5) -0*MySIN((360.0/(n1*6))*5)
         y(2,12) = a1*MySIN((360.0/(n1*6))*5) +0*MyCOS((360.0/(n1*6))*5)
        

        do i = 1,4
                xatom(1,i,1) = x(1,i)
                yatom(1,i,1) = y(1,i)
        enddo
        
        do i = 1, 12
                xatom(2,i,1) = x(2,i)
                yatom(2,i,1) = y(2,i)
        enddo

        do j = 1, 4
         do i = 2, n1
                xatom(1,j,i) = MyCOS(angle)*xatom(1,j,i-1)&
                                -  MySIN(angle)*yatom(1,j,i-1)
                yatom(1,j,i) = xatom(1,j,i-1)*MySIN(angle)&
                                + yatom(1,j,i-1)*MyCOS(angle)

        enddo
        enddo

        do j = 1, 12
         do i = 2, n1
                xatom(2,j,i) = MyCOS(angle)*xatom(2,j,i-1)&
                                -  MySIN(angle)*yatom(2,j,i-1)
                yatom(2,j,i) = xatom(2,j,i-1)*MySIN(angle)&
                                + yatom(2,j,i-1)*MyCOS(angle)



        enddo
        enddo


        do i = 1, n1

        high(1,1,i) = 4*(bind_length/2)
        high(1,2,i) = bind_length/2
        high(1,3,i) = bind_length/2
        high(1,4,i) = 4*(bind_length/2)

        enddo

        do i = 1, n1
        
        high(2,1,i) = 2*(bind_length/2)
        high(2,2,i) = 6*(bind_length/2)
        high(2,3,i) = 3*(bind_length/2)
        high(2,4,i) = 5*(bind_length/2)
        high(2,5,i) = 0
        high(2,6,i) = 2*(bind_length/2)
        high(2,7,i) = 3*(bind_length/2)
        high(2,8,i) = 5*(bind_length/2)
        high(2,9,i) = 2*(bind_length/2)
        high(2,10,i) = 6*( bind_length/2)
        high(2,11,i) = 3*(bind_length/2)
        high(2,12,i) = 5*(bind_length/2)

        enddo

        allocate(fxatom(2,12,n1), fyatom(2,12,n1),fhigh(2,12,n1))

      
        do j = 1, 4
        do k = 1, n1

          fxatom(1,j,k) = (xatom(1,j,k)+lattice)/total_cell
          fyatom(1,j,k) = (yatom(1,j,k)+lattice)/total_cell

          fhigh(1,j,k) = high(1,j,k)/total_high

        enddo
        enddo
       
        do j = 1, 12
        do k = 1, n1

          fxatom(2,j,k) = (xatom(2,j,k)+lattice)/total_cell
          fyatom(2,j,k) = (yatom(2,j,k)+lattice)/total_cell

          fhigh(2,j,k) = high(2,j,k)/total_high

        enddo
        enddo


     open(24,file='tube.out')
 100  format(f20.16,3x,f20.16,3x,f20.16)

        write(24,*) total_cell, total_high, 'atom1= ', n1*4, 'atom2 =',n1*12

      do j = 1, 4
      do k = 1, n1

         write(24,100) fxatom(1,j,k), fyatom(1,j,k), fhigh(1,j,k)
      enddo
      enddo

      do j = 1, 12
      do k = 1, n1

         write(24,100) fxatom(2,j,k), fyatom(2,j,k), fhigh(2,j,k)
      enddo
      enddo

        endif

         !zig zag

          if(n2 .EQ. 0) THEN
                radis = (bind_length*3*(n1-1))/(2*PI)
                lattice = radis + (radis/2)*3
                total_high = ((bind_length/2)*SQRT(3.0))*6
                total_cell = lattice*2
          angle = 360.0/n1

        if(total_cell .GE. 50) THEN
                total_cell = 49
        endif 
        
         allocate(x(2,12), y(2,12), xatom(2,12,n1), yatom(2,12,n1))
         allocate(high(2,12,n1))

         x(1,1) = radis
         y(1,1) = 0
        
         x(1,2) = radis
         y(1,2) = 0

         x(1,3) = radis*MyCOS((360.0/(n1*6)*3)) - 0*MySIN((360.0/(n1*6))*3)
         y(1,3) = radis*MySIN((360.0/(n1*6)*3)) + 0*MyCOS((360.0/(n1*6))*3)

         x(1,4) = radis*MyCOS((360.0/(n1*6)*3)) - 0*MySIN((360.0/(n1*6))*3)
         y(1,4) = radis*MySIN((360.0/(n1*6)*3)) + 0*MyCOS((360.0/(n1*6))*3)

         a1 = radis - 1.56377985
         b1 = radis + 1.56401799

         x(2,1) = a1*MyCOS((360.0/(n1*6)*1))&
                  -0*MySIN((360.0/(n1*6))*1)
         y(2,1) = a1*MySIN((360.0/(n1*6)*1))&
                  +0*MyCOS((360.0/(n1*6))*1)

         x(2,2) = a1*MyCOS((360.0/(n1*6)*1))&
                  -0*MySIN((360.0/(n1*6))*1)
         y(2,2) = a1*MySIN((360.0/(n1*6)*1))&
                  +0*MyCOS((360.0/(n1*6))*1)

         x(2,3) = a1*MyCOS((360.0/(n1*6)*1))&
                  -0*MySIN((360.0/(n1*6))*1)
         y(2,3) = a1*MySIN((360.0/(n1*6)*1))&
                  +0*MyCOS((360.0/(n1*6))*1)

         x(2,4) = b1*MyCOS((360.0/(n1*6)*2))&
                  -0*MySIN((360.0/(n1*6))*2)
         y(2,4) = b1*MySIN((360.0/(n1*6)*2))&
                  +0*MyCOS((360.0/(n1*6))*2)
     
         x(2,5) = b1*MyCOS((360.0/(n1*6)*2))&
                  -0*MySIN((360.0/(n1*6))*2)
         y(2,5) = b1*MySIN((360.0/(n1*6)*2))&
                  +0*MyCOS((360.0/(n1*6))*2)

         x(2,6) = b1*MyCOS((360.0/(n1*6)*2))&
                  -0*MySIN((360.0/(n1*6))*2)
         y(2,6) = b1*MySIN((360.0/(n1*6)*2))&
                  +0*MyCOS((360.0/(n1*6))*2)

         x(2,7) = a1*MyCOS((360.0/(n1*6)*4))&
                  -0*MySIN((360.0/(n1*6))*4)
         y(2,7) = a1*MySIN((360.0/(n1*6)*4))&
                  +0*MyCOS((360.0/(n1*6))*4)

         x(2,8) = a1*MyCOS((360.0/(n1*6)*4))&
                  -0*MySIN((360.0/(n1*6))*4)
         y(2,8) = a1*MySIN((360.0/(n1*6)*4))&
                  +0*MyCOS((360.0/(n1*6))*4)

         x(2,9) = a1*MyCOS((360.0/(n1*6)*4))&
                  -0*MySIN((360.0/(n1*6))*4)
         y(2,9) = a1*MySIN((360.0/(n1*6)*4))&
                  +0*MyCOS((360.0/(n1*6))*4)

         x(2,10) = b1*MyCOS((360.0/(n1*6)*5))&
                  -0*MySIN((360.0/(n1*6))*5)
         y(2,10) = b1*MySIN((360.0/(n1*6)*5))&
                  +0*MyCOS((360.0/(n1*6))*5)

         x(2,11) = b1*MyCOS((360.0/(n1*6)*5))&
                  -0*MySIN((360.0/(n1*6))*5)
         y(2,11) = b1*MySIN((360.0/(n1*6)*5))&
                  +0*MyCOS((360.0/(n1*6))*5)

         x(2,12) = b1*MyCOS((360.0/(n1*6)*5))&
                  -0*MySIN((360.0/(n1*6))*5)
         y(2,12) = b1*MySIN((360.0/(n1*6)*5))&
                  +0*MyCOS((360.0/(n1*6))*5)

        do i = 1,4
                xatom(1,i,1) = x(1,i)
                yatom(1,i,1) = y(1,i)
        enddo

        do i = 1, 12
                xatom(2,i,1) = x(2,i)
                yatom(2,i,1) = y(2,i)
        enddo

        do j = 1, 4
         do i = 2, n1
                xatom(1,j,i) = MyCOS(angle)*xatom(1,j,i-1)&
                                -  MySIN(angle)*yatom(1,j,i-1)
                yatom(1,j,i) = xatom(1,j,i-1)*MySIN(angle)&
                                + yatom(1,j,i-1)*MyCOS(angle)

        enddo
        enddo

        do j = 1, 12
         do i = 2, n1
                xatom(2,j,i) = MyCOS(angle)*xatom(2,j,i-1)&
                                -  MySIN(angle)*yatom(2,j,i-1)
                yatom(2,j,i) = xatom(2,j,i-1)*MySIN(angle)&
                                + yatom(2,j,i-1)*MyCOS(angle)



        enddo
        enddo

        c1 = (bind_length/2)*SQRT(3.0)

        do i = 1, n1

        high(1,1,i) = c1
        high(1,2,i) = c1*3
        high(1,3,i) = 0
        high(1,4,i) = c1*4

        enddo

        do i = 1, n1

        high(2,1,i) = 0
        high(2,2,i) = c1*2
        high(2,3,i) = c1*4
        high(2,4,i) = c1
        high(2,5,i) = c1*3
        high(2,6,i) = c1*5
        high(2,7,i) = c1
        high(2,8,i) = c1*3
        high(2,9,i) = c1*5
        high(2,10,i) = 0
        high(2,11,i) = c1*2
        high(2,12,i) = c1*4

        enddo

        allocate(fxatom(2,12,n1), fyatom(2,12,n1),fhigh(2,12,n1))


        do j = 1, 4
        do k = 1, n1

          fxatom(1,j,k) = (xatom(1,j,k)+lattice)/total_cell
          fyatom(1,j,k) = (yatom(1,j,k)+lattice)/total_cell

          fhigh(1,j,k) = high(1,j,k)/total_high

        enddo
        enddo

        do j = 1, 12
        do k = 1, n1

          fxatom(2,j,k) = (xatom(2,j,k)+lattice)/total_cell
          fyatom(2,j,k) = (yatom(2,j,k)+lattice)/total_cell

          fhigh(2,j,k) = high(2,j,k)/total_high

        enddo
        enddo


     open(24,file='tube.out')

        write(24,*) total_cell, total_high, 'atom1= ', n1*4, 'atom2=',n1*12

      do j = 1, 4
      do k = 1, n1

         write(24,100) fxatom(1,j,k), fyatom(1,j,k), fhigh(1,j,k)
      enddo
      enddo

      do j = 1, 12
      do k = 1, n1

         write(24,100) fxatom(2,j,k), fyatom(2,j,k), fhigh(2,j,k)
      enddo
      enddo

        endif




          END PROGRAM  main
