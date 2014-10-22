module subroutines
contains
subroutine sis_area(wk_1,obs,ObsErrStd,key,wk,neff)
use shared_data
real(kind=kind_num),dimension(:),intent(in) :: wk_1
real(kind=kind_num),dimension(:,:),intent(in) :: obs,ObsErrStd
real(kind=kind_num),dimension(:,:,:),intent(in) :: key
real(kind=kind_num),dimension(:),intent(out) :: wk
integer,intent(out) :: neff

real(kind=kind_num) :: sumw
real(kind=kind_num),allocatable,dimension(:) :: sumup
integer :: i,j,k

allocate( sumup(size(wk_1)) )

! wk(i) ~ p[yk|xk(i)]*wk_1(i)
! p[yk|xk(i)] ~ exp(-1/2 * (yk-H[xk(i)])^2/(sigma^2) )

!!=== Eq_East : 120W-115W 3S-3N
!do k=1,size(wk_1)
!	sumup(k)=0.0
!	do i=121,136
!		do j=23,29
!			sumup(k)=sumup(k)+((obs(i,j)-key(i,j,k))/ObsErrStd(i,j))**2
!		end do
!	end do
!	wk(k)= exp( -sumup(k)/2.0 )*wk_1(k)
!end do

!=== Eq_West : 180W-165W 3S-3N
do k=1,size(wk_1)
	sumup(k)=0.0
	do i=61,76
		do j=23,29
			sumup(k)=sumup(k)+((obs(i,j)-key(i,j,k))/ObsErrStd(i,j))**2
		end do
	end do
	wk(k)= exp( -sumup(k)/2.0 )*wk_1(k)
end do

sumw=sum(wk)
wk=wk/sumw

neff=1/sum(wk**2)

deallocate(sumup)
end subroutine sis_area
end module subroutines
