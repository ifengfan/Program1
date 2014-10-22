module shared_data
	integer,parameter :: kind_num=selected_real_kind(p=15,r=37)
	real(kind=kind_num),parameter :: pi=3.1415926
end module shared_data

program fengfan
use netcdf
use shared_data
implicit none 
!-----------------------NAMELIST VARs,DO NOT CHANGE----------------------------
integer :: varid_ssh,varid_sst,varid_u,varid_nino3,varid_nino4,varid_nino34,ncid
integer,dimension(4) :: start_ssh,count_ssh,start_sst,count_sst,start_u,count_u
integer,dimension(2) :: start_nino3,count_nino3,start_nino4,count_nino4,start_nino34,count_nino34
integer,dimension(20) :: EPyrs
character(len=999) :: file_name,path
namelist /mylist/ varid_ssh,start_ssh,count_ssh,&
                  varid_sst,start_sst,count_sst,&
                  varid_u,start_u,count_u,&
                  varid_nino3,start_nino3,count_nino3,&
                  varid_nino4,start_nino4,count_nino4,&
                  varid_nino34,start_nino34,count_nino34,&
                  EPyrs,path
!------------------------------------------------------------------------------

integer :: varid_key,varid_nino
integer,dimension(4) :: start_key,count_key
integer,dimension(2) :: start_nino,count_nino
real(kind=kind_num) :: fillvalue

!----Basic variables:
real(kind=kind_num),allocatable,dimension(:,:,:,:) :: key    !ssh/sst/u
real(kind=kind_num),allocatable,dimension(:,:)     :: nino   !nino3/nino4/nino3.4
logical,allocatable,dimension(:,:)     :: notnan

!----Variables used to calc ObsErrStd:
real(kind=kind_num),allocatable,dimension(:,:,:) :: sum1,sum2,key_std,ObsErrStd

!----Variables related to 'truth year':
real(kind=kind_num),allocatable,dimension(:,:,:) :: truth

!----Variables related to synthetic obs:
real(kind=kind_num) :: rand1,rand2,rand,obs

!----Variables related to weights:
real(kind=kind_num),allocatable,dimension(:) :: wk_1,wk
integer :: neff

!----Variables related to new weights after resample:
real(kind=kind_num),allocatable,dimension(:) :: xk_new,wk_new
integer,allocatable,dimension(:) :: parent

!----Variables related to new paticles after resample:
real(kind=kind_num),allocatable,dimension(:) :: key_new,nino_new

!----Variables related to Predictive Power:
real(kind=kind_num) :: Sq(12),Var_q(12)
real(kind=kind_num) :: Sp,Var_p,pp,pp_sum
real(kind=kind_num),allocatable,dimension(:,:)   :: pp_100ave
real(kind=kind_num),allocatable,dimension(:,:,:) :: pp_cases
real(kind=kind_num),allocatable,dimension(:,:)   :: pp_case_ave

integer :: i,j,k,index,cnt
integer :: m_obs,x_obs,y_obs,sele

open(unit=333,file='CCSM4_nml',status='old',action='read',delim='apostrophe')
read(unit=333,nml=mylist)
close(unit=333)

!====INPUT THE MAIN VARIABLES YOU USE IN THIS PROGRAM :
varid_key = varid_u
varid_nino= varid_nino3
start_key = start_u
count_key = count_u
start_nino= start_nino3
count_nino= count_nino3

allocate( key(count_key(1),count_key(2),count_key(3),count_key(4)),&
         nino(count_nino(1),count_nino(2)),&
       notnan(count_key(1),count_key(2)) )

allocate(     sum1(count_key(1),count_key(2),count_key(3)),&
              sum2(count_key(1),count_key(2),count_key(3)),&
           key_std(count_key(1),count_key(2),count_key(3)),&
         ObsErrStd(count_key(1),count_key(2),count_key(3)) )

allocate(    truth(count_key(1),count_key(2),count_key(3)))

allocate(   wk_1(count_key(4)) , wk(count_key(4)) )

allocate( xk_new(count_key(4)) , wk_new(count_key(4)) , parent(count_key(4)))

allocate( key_new(count_key(4)) , nino_new(count_key(4)) )

allocate(  pp_100ave(count_key(1),count_key(2)),&
            pp_cases(count_key(1),count_key(2),size(EPyrs)),&
         pp_case_ave(count_key(1),count_key(2)) )

file_name=path
write(*,*) "file is ", trim(adjustl(file_name))

call check(nf90_open(trim(adjustl(file_name)), nf90_nowrite, ncid=ncid))
call check(nf90_get_var(ncid,varid=varid_key, values=key, start=start_key, count=count_key))
call check(nf90_get_att(ncid,varid=varid_key,name="_FillValue",values=fillvalue))
call check(nf90_get_var(ncid,varid=varid_nino,values=nino,start=start_nino,count=count_nino))
call check(nf90_close(ncid))

!=======Place the NaNs=======
notnan(:,:)=.True.
where(key(:,:,1,1)==fillvalue) notnan=.False.
!print *, notnan(1,1)

!=======Calculate the Entropy======
do i=1,12
	call cal_entropy(size(nino,2),nino(i,:),Sq(i))
enddo
!open(100,file='Sq_nino3.txt')
!write(100,*) Sq
!close(100)

!======Calculate the Variance of NINO index=====
do i=1,12
	Var_q(i)= ( sum(nino(i,:)**2) - sum(nino(i,:))**2/size(nino(i,:)) )/( size(nino(i,:))-1 )
enddo
!write(*,*) Var_q

!=======Calculate the ObsErrStd======
sum1(:,:,:)=0.0
sum2(:,:,:)=0.0
key_std(:,:,:)=0.0
do i=1,size(key,4)
	sum1=sum1+key(:,:,:,i)
	sum2=sum2+key(:,:,:,i)**2
end do
key_std=sqrt((count_key(4)*sum2-sum1**2)/(count_key(4)*(count_key(4)-1)))
do i=1,count_key(3)
	where(notnan==.False.) key_std(:,:,i)=sqrt(-1.0)  !NaNs
enddo 
ObsErrStd = 0.3*key_std  !~~~~ObsErrStd now have NaNs instead of FillValue~~~~
!-----test ObsErrStd-----
!open(10,file='ObsErrStd_sst.bin',access='direct',form='unformatted',&
!    status='replace',recl=kind_num*161*51)
!do i=1,12
!	write(10,rec=i) ObsErrStd(:,:,i)
!enddo
!close(10)

!~~~~~~~~NOTICE FOR DRAWING CONTOURS~~~~~~~~
!MATLAB: USE NaNs
!NCL   : USE FillValue

!====Main====
index = 20  ! The year you select as the truth year
write(*,*) "The truth year is ", EPyrs(index)
truth = key(:,:,:,EPyrs(index))
do i=1,count_key(3)
	where(notnan==.False.) truth(:,:,i)=sqrt(-1.0)  !NaNs
enddo  !~~~~truth now have NaNs instead of FillValue~~~~

m_obs = 1  ! The first month to do the assimilation
do i=1,count_key(1)
	do j=1,count_key(2)
		if (notnan(i,j)) then
			x_obs = i
			y_obs = j
			pp_sum = 0.0
			do cnt=1,100
				call init_random_seed()
				!call randmom_seed()
				call random_number(rand1)
				call random_number(rand2)
				rand=sqrt(-2.0*log(rand1))*cos(2.0*pi*rand2)
				obs=truth(x_obs,y_obs,m_obs)+rand*ObsErrStd(x_obs,y_obs,m_obs)
				wk_1(:)=1.0/count_key(4)
				call sis(count_key(4),wk_1,obs,ObsErrStd(x_obs,y_obs,m_obs),key(x_obs,y_obs,m_obs,:),wk,neff)
				call resample(count_key(4),nino(m_obs,:),wk,xk_new,wk_new,parent)
				
				do sele=1,count_key(4)
					key_new(sele)=key(x_obs,y_obs,m_obs+1,parent(sele))
					nino_new(sele)=nino(m_obs+1,parent(sele))
				end do
				call init_random_seed()
				!call randmom_seed()
				call random_number(rand1)
				call random_number(rand2)
				rand=sqrt(-2.0*log(rand1))*cos(2.0*pi*rand2)
				obs=truth(x_obs,y_obs,m_obs+1)+rand*ObsErrStd(x_obs,y_obs,m_obs+1)
				call sis(count_key(4),wk_new,obs,ObsErrStd(x_obs,y_obs,m_obs+1),key_new,wk,neff)
				call resample(count_key(4),nino_new,wk,xk_new,wk_new,parent)
				
				do sele=1,count_key(4)
					key_new(sele)=key(x_obs,y_obs,m_obs+2,parent(sele))
					nino_new(sele)=nino(m_obs+2,parent(sele))
				end do
				call init_random_seed()
				!call randmom_seed()
				call random_number(rand1)
				call random_number(rand2)
				rand=sqrt(-2.0*log(rand1))*cos(2.0*pi*rand2)
				obs=truth(x_obs,y_obs,m_obs+2)+rand*ObsErrStd(x_obs,y_obs,m_obs+2)
				call sis(count_key(4),wk_new,obs,ObsErrStd(x_obs,y_obs,m_obs+2),key_new,wk,neff)
				call resample(count_key(4),nino_new,wk,xk_new,wk_new,parent)
				! Calculate PP using entropy:
				call cal_entropy(count_key(4),xk_new,Sp)
				pp=1.0-exp(Sp-Sq(m_obs+2))
				! OR Calculate PP using variance:
!				Var_p=( sum(xk_new**2) - sum(xk_new)**2/size(xk_new) )/( size(xk_new)-1 )
!				pp=1.0-Var_p/Var_q(m_obs+2)
				
				pp_sum=pp_sum + pp
			end do
		end if
		pp_100ave(i,j)=pp_sum/100.0
	end do
end do
where(notnan==.False.) pp_100ave=sqrt(-1.0)
print *, maxval(pp_100ave)
print *, maxloc(pp_100ave) !it only returns the fisrt element in column-major order.
!locations = maxloc(pp_100ave,mask=pp_100ave .GE. 0.13)
open(12,file='pp_u_JFM_case20_0.3*_100ave.bin',&
     access='direct',form='unformatted',&
     status='replace',recl=kind_num*161*51)
write(12,rec=1) pp_100ave
close(12)

deallocate(key,nino,notnan,sum1,sum2,key_std,ObsErrStd,truth,&
           wk_1,wk,xk_new,wk_new,parent,key_new,nino_new,&
           pp_100ave,pp_cases,pp_case_ave)

contains
subroutine check(status)
	integer,intent(in) :: status
	if(status /= nf90_noerr) then
		print *, trim(nf90_strerror(status))
		stop 2
	end if
end subroutine check
end program fengfan
