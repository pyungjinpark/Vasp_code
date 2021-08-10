!123456
       integer :: atom1, atom2, totatom
       real,allocatable :: postion1(:),postion2(:),postion3(:)
       real,allocatable :: newpostion1(:),newpostion2(:)      

        OPEN(1,file='POSCAR')
        read(1,*)
        read(1,*)
        read(1,*)
        read(1,*)
        read(1,*)
        read(1,*)
        read(1,*) atom1, atom2
        read(1,*) 

        totatom = 0
        totatom = atom1 + atom2
     
        allocate(postion1(totatom), postion2(totatom),postion3(totatom))
        allocate(newpostion1(totatom), newpostion2(totatom))


        do i = 1,totatom
                
        read(1,*) postion1(i), postion2(i), postion3(i)
        
        enddo
 
        newpostion1(1) = 0
        newpostion2(1) = 0

        do i = 1,totatom
        
        newpostion1(i) = cos(30.0)*postion1(i) - sin(30.0)*postion2(i)
        newpostion2(i) = sin(30.0)*postion1(i) + cos(30.0)*postion2(i) 

        enddo
        write(6,*) postion1(1)
100     format(f14.7,3x,f14.7,3x,f14.7)
        open(2,file='postion.txt')
        
        do i = 1,totatom
        write(2,100) newpostion1(i), newpostion2(i), postion3(i)
        enddo
  
        end
