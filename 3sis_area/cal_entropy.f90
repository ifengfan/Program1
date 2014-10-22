subroutine cal_entropy(np,vector,S)
	use shared_data
	implicit none
	integer,intent(in) :: np
	real(kind=kind_num),dimension(np),intent(in) :: vector 
	real(kind=kind_num),intent(out) :: S
	
	integer :: k,nbar
	real(kind=kind_num),allocatable :: freq(:)
	real(kind=kind_num) :: interval
	
	nbar=int(sqrt(real(np)))+1
	allocate(freq(nbar))
	interval = (maxval(vector)-minval(vector))/real(nbar)
	freq(1)=count(vector(:)<minval(vector)+interval)
	do k=2,nbar-1
		freq(k)=count( vector(:)>=minval(vector)+(k-1)*interval .AND. &
                   vector(:)< minval(vector)+ k*interval )
	enddo
	freq(nbar)=count( vector(:)>=minval(vector)+(nbar-1)*interval)
	freq=freq/real(np)
	S=0.0
	do k=1,nbar
		if (freq(k)/=0.0) then
			S=S+freq(k)*log(freq(k))
		end if
	enddo
	S=-interval*S
	
	deallocate(freq)
	return
end subroutine cal_entropy
