subroutine init_random_seed()
	integer :: i,n,clock
	integer,allocatable :: seed(:)
	
	call random_seed(size=n)
	allocate(seed(n))
	
	call system_clock(count=clock)
	
	seed = clock + 37*(/(i-1,i=1,n)/)
	call random_seed(put=seed)
	
	deallocate(seed)
end subroutine init_random_seed

!program test_init_randseed
!implicit none
!real :: rand1,rand2,rand
!real,parameter :: pi=3.1415926
!
!call init_random_seed() 
!!call random_seed()
!call random_number(rand1)
!call random_number(rand2)
!rand=sqrt(-2.0*log(rand1))*cos(2.0*pi*rand2)
!print *, rand
!
!end 
