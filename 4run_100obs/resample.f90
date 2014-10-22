subroutine resample(np,xk,wk,xk_new,wk_new,parent)
	use shared_data
	implicit none
	integer,intent(in) :: np
	real(kind=kind_num),dimension(np),intent(in)  :: xk,wk
	real(kind=kind_num),dimension(np),intent(out) :: xk_new,wk_new
	integer,dimension(np),intent(out) :: parent
	real(kind=kind_num) :: u1,cumsum_wk(np),u(np)
	integer :: i,j
	
	cumsum_wk(1)=wk(1)
	do i=2,np
		cumsum_wk(i)=cumsum_wk(i-1)+wk(i)
	end do
	
	i=1
	call random_seed()
	call random_number(u1)
	u1=u1/real(np)
	do j=1,np
		u(j)=u1+real(j-1)/real(np)
		do while (u(j)>cumsum_wk(i))
			i=i+1
		end do
		xk_new(j)=xk(i)
		wk_new(j)=1/real(np)
		parent(j)=i
	end do
	
	return
end subroutine resample

!program test_resample
!implicit none
!
!real :: xk(500),wk(500),xk_new(500),wk_new(500)
!integer :: parent(500)
!
!call random_number(xk)
!call random_number(wk)
!wk=wk/sum(wk)
!
!call resample(size(xk),xk,wk,xk_new,wk_new,parent)
!end program test_resample

