    
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
                radis = ((n1-1)*bind_length*3)/(2*PI)
                lattice = radis + (radis/2)*2
                total_high = (bind_length/2)*SQRT(3.0)*2
                total_cell = lattice*2
          angle = 360.0/n1

        if(total_cell .GE. 50) THEN
                total_cell = 49
        endif

         allocate(x(2,3), y(2,3), xatom(2,3,n1), yatom(2,3,n1))
         allocate(high(2,3,n1))

        
         x(1,1) = radis
         y(1,1) = 0
       
         a1 = radis - 1.3
         b1 = 0

         c1 = radis + 1.3
         d1 = 0

         x(1,2) = a1*MyCOS(360.0/(n1*3)) - b1*MySIN(360.0/(n1*3))
         y(1,2) = a1*MySIN(360.0/(n1*3)) + b1*MyCOS(360.0/(n1*3))

         x(1,3) = c1*MyCOS(360.0/(n1*3)) + d1*MySIN(360.0/(n1*3))
         y(1,3) = -c1*MySIN(360.0/(n1*3)) + d1*MyCOS(360.0/(n1*3))
        
         xatom(1,1,1) = x(1,1)
         yatom(1,1,1) = y(1,1)

         xatom(1,2,1) = x(1,2)
         yatom(1,2,1) = y(1,2)
        
         xatom(1,3,1) = x(1,3)
         yatom(1,3,1) = y(1,3)

         do i = 2, n1
                xatom(1,1,i) = MyCOS(angle)*xatom(1,1,i-1)&
                                -  MySIN(angle)*yatom(1,1,i-1)
                yatom(1,1,i) = xatom(1,1,i-1)*MySIN(angle)&
                                + yatom(1,1,i-1)*MyCOS(angle)
                
                xatom(1,2,i) = xatom(1,2,i-1)*MyCOS(angle)&
                               - yatom(1,2,i-1)*MySIN(angle)
                yatom(1,2,i) = xatom(1,2,i-1)*MySIN(angle)&
                               + yatom(1,2,i-1)*MyCOS(angle)

                xatom(1,3,i) = xatom(1,3,i-1)*MyCOS(angle) &
                               - yatom(1,3,i-1)*MySIN(angle)
                yatom(1,3,i) = xatom(1,3,i-1)*MySIN(angle) &
                               + yatom(1,3,i-1)*MyCOS(angle)
        

        enddo
        
         a2 = radis
         b2 = 0

         c2 = radis - 1.3
         d2 = 0

         e2 = radis + 1.3
         f2 = 0

         x(2,1) = a2*MyCOS((360.0/(n1*6))*3) - b2*MySIN((360.0/(n1*6))*3)
         y(2,1) = a2*MySIN((360.0/(n1*6))*3) + b2*MyCOS((360.0/(n1*6))*3)

         x(2,2) = c2*MyCOS(360.0/(n1*6)) + d2*MySIN(360.0/(n1*6))
         y(2,2) = -c2*MySIN(360.0/(n1*6)) + d2*MyCOS(360.0/(n1*6))

         x(2,3) = e2*MyCOS(360.0/(n1*6)) - f2*MySIN(360.0/(n1*6))
         y(2,3) = e2*MySIN(360.0/(n1*6)) + f2*MyCOS(360.0/(n1*6))

         xatom(2,1,1) = x(2,1)
         yatom(2,1,1) = y(2,1)

         xatom(2,2,1) = x(2,2)
         yatom(2,2,1) = y(2,2)

         xatom(2,3,1) = x(2,3)
         yatom(2,3,1) = y(2,3)

         do i = 2, n1
                xatom(2,1,i) = MyCOS(angle)*xatom(2,1,i-1)&
                                -  MySIN(angle)*yatom(2,1,i-1)
                yatom(2,1,i) = xatom(2,1,i-1)*MySIN(angle)&
                                + yatom(2,1,i-1)*MyCOS(angle)

                xatom(2,2,i) = xatom(2,2,i-1)*MyCOS(angle)&
                               - yatom(2,2,i-1)*MySIN(angle)
                yatom(2,2,i) = xatom(2,2,i-1)*MySIN(angle)&
                               + yatom(2,2,i-1)*MyCOS(angle)

                xatom(2,3,i) = xatom(2,3,i-1)*MyCOS(angle) &
                               - yatom(2,3,i-1)*MySIN(angle)
                yatom(2,3,i) = xatom(2,3,i-1)*MySIN(angle) &
                               + yatom(2,3,i-1)*MyCOS(angle)


        enddo
           
         do i = 1, n1
        high(1,1,i) = 0
        high(1,2,i) = 0
        high(1,3,i) = 0
        
        high(2,1,i) = total_high/2
        high(2,2,i) = total_high/2
        high(2,3,i) = total_high/2

         enddo

        allocate(fxatom(2,3,n1), fyatom(2,3,n1),fhigh(2,3,n1))

        do i = 1, 2
        do j = 1, 3
        do k = 1, n1
        
          fxatom(i,j,k) = (xatom(i,j,k)+lattice)/total_cell
          fyatom(i,j,k) = (yatom(i,j,k)+lattice)/total_cell

          fhigh(i,j,k) = high(i,j,k)/total_high

        enddo
        enddo
        enddo

     open(24,file='tube.out')
 100  format(f20.16,3x,f20.16,3x,f20.16)

        write(24,*) total_cell, total_high, 'atom1= ', n1*2, 'atom2 =', n1*4
     
      do j = 1, 3
      do k = 1, n1
       
         write(24,100) fxatom(1,j,k), fyatom(1,j,k), fhigh(1,j,k)
         write(24,100) fxatom(2,j,k), fyatom(2,j,k), fhigh(2,j,k)
      enddo
      enddo
        endif


      
        !zig zag

        if(n2 .eq. 0 ) THEN        
                radis = (2*(bind_length/2)*SQRT(3.0)*(n1-1))/(2*PI)
                lattice = radis + (radis/2)*3
                total_high = (bind_length/2)*6
                total_cell = lattice*2
                angle = 360.0/n1

        if(total_cell .GE. 50) THEN
                total_cell = 49
        endif

         allocate(x(2,3), y(2,3), xatom(2,3,n1), yatom(2,3,n1))
         allocate(high(2,3,n1))

        x(1,1) = radis
        y(1,1) = 0

        a1 = radis - 1.65982445
        b1 = 0

        c1 = radis + 1
        d1 = 0

        x(1,2) = a1*MyCOS((360.0/(n1*2))) - b1*MySIN(360.0/(n1*2))
        y(1,2) = a1*MySIN((360.0/(n1*2))) + b1*MyCOS(360.0/(n1*2))

        x(1,3) = c1*MyCOS((360.0/(n1*2))) - d1*MySIN(360.0/(n1*2))
        y(1,3) = c1*MySIN((360.0/(n1*2))) + d1*MyCOS(360.0/(n1*2))

         xatom(1,1,1) = x(1,1)
         yatom(1,1,1) = y(1,1)

         xatom(1,2,1) = x(1,2)
         yatom(1,2,1) = y(1,2)

         xatom(1,3,1) = x(1,3)
         yatom(1,3,1) = y(1,3)
        
         do i = 2, n1
                xatom(1,1,i) = MyCOS(angle)*xatom(1,1,i-1)&
                                -  MySIN(angle)*yatom(1,1,i-1)
                yatom(1,1,i) = xatom(1,1,i-1)*MySIN(angle)&
                                + yatom(1,1,i-1)*MyCOS(angle)

                xatom(1,2,i) = xatom(1,2,i-1)*MyCOS(angle)&
                               - yatom(1,2,i-1)*MySIN(angle)
                yatom(1,2,i) = xatom(1,2,i-1)*MySIN(angle)&
                               + yatom(1,2,i-1)*MyCOS(angle)

                xatom(1,3,i) = xatom(1,3,i-1)*MyCOS(angle) &
                               - yatom(1,3,i-1)*MySIN(angle)
                yatom(1,3,i) = xatom(1,3,i-1)*MySIN(angle) &
                               + yatom(1,3,i-1)*MyCOS(angle)


        enddo

        a2 = radis 
        b2 = 0

        x(2,2) = radis - 1.65982445
        y(2,2) = 0

        x(2,3) = radis + 1
        y(2,3) = 0

        x(2,1) = a2*MyCOS((360.0/(n1*2))) - b2*MySIN(360.0/(n1*2))
        y(2,1) = a2*MySIN((360.0/(n1*2))) + b2*MyCOS(360.0/(n1*2))

         xatom(2,1,1) = x(2,1)
         yatom(2,1,1) = y(2,1)

         xatom(2,2,1) = x(2,2)
         yatom(2,2,1) = y(2,2)

         xatom(2,3,1) = x(2,3)
         yatom(2,3,1) = y(2,3)

         do i = 2, n1
                xatom(2,1,i) = MyCOS(angle)*xatom(2,1,i-1)&
                                -  MySIN(angle)*yatom(2,1,i-1)
                yatom(2,1,i) = xatom(2,1,i-1)*MySIN(angle)&
                                + yatom(2,1,i-1)*MyCOS(angle)

                xatom(2,2,i) = xatom(2,2,i-1)*MyCOS(angle)&
                               - yatom(2,2,i-1)*MySIN(angle)
                yatom(2,2,i) = xatom(2,2,i-1)*MySIN(angle)&
                               + yatom(2,2,i-1)*MyCOS(angle)

                xatom(2,3,i) = xatom(2,3,i-1)*MyCOS(angle) &
                               - yatom(2,3,i-1)*MySIN(angle)
                yatom(2,3,i) = xatom(2,3,i-1)*MySIN(angle) &
                               + yatom(2,3,i-1)*MyCOS(angle)


        enddo
        
         do i = 1, n1
        high(1,1,i) = bind_length/2
        high(1,2,i) = 0
        high(1,3,i) = (bind_length/2)*2

        high(2,1,i) = (bind_length/2)*4
        high(2,2,i) = (bind_length/2)*3
        high(2,3,i) = (bind_length/2)*5

         enddo

        allocate(fxatom(2,3,n1), fyatom(2,3,n1),fhigh(2,3,n1))

        do i = 1, 2
        do j = 1, 3
        do k = 1, n1

          fxatom(i,j,k) = (xatom(i,j,k)+lattice)/total_cell
          fyatom(i,j,k) = (yatom(i,j,k)+lattice)/total_cell

          fhigh(i,j,k) = high(i,j,k)/total_high

        enddo
        enddo
        enddo

     open(24,file='tube.out')

        write(24,*) total_cell, total_high, 'atom1 = ', n1*2, 'atom2 = ',n1*4

      do j = 1, 3
      do k = 1, n1

         write(24,100) fxatom(1,j,k), fyatom(1,j,k), fhigh(1,j,k)
         write(24,100) fxatom(2,j,k), fyatom(2,j,k), fhigh(2,j,k)
      enddo
      enddo


        endif


        END PROGRAM  main
