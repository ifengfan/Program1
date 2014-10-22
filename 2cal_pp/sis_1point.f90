subroutine sis(np,wk_1,obs,ObsErrStd,part,wk,neff)
use shared_data
implicit none
integer,intent(in) :: np
real(kind=kind_num),dimension(np),intent(in) :: wk_1,part
real(kind=kind_num),intent(in) :: obs,ObsErrStd
real(kind=kind_num),dimension(np),intent(out) :: wk
integer,intent(out) :: neff
real(kind=kind_num) :: sumw
integer :: i

! wk(i) ~ p[yk|xk(i)]*wk_1(i)
! p[yk|xk(i)] ~ exp(-1/2 * (yk-H[xk(i)])^2/(sigma^2) )

do i=1,np
	wk(i)= exp( -( (obs-part(i))/(sqrt(2.0)*ObsErrStd) )**2 )*wk_1(i)
enddo

sumw=sum(wk)
wk=wk/sumw

neff=1/sum(wk**2)

end subroutine sis
