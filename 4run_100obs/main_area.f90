module shared_data
	integer,parameter :: kind_num=selected_real_kind(p=15,r=37)
	real(kind=kind_num),parameter :: pi=3.1415926
end module shared_data

program fengfan
use subroutines
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
real(kind=kind_num),allocatable,dimension(:,:,:) :: rand1,rand2,rand,obs

!----Variables related to weights:
real(kind=kind_num),allocatable,dimension(:) :: wk_1,wk
integer :: neff

!----Variables related to new weights after resample:
real(kind=kind_num),allocatable,dimension(:) :: xk_new,wk_new
integer,allocatable,dimension(:) :: parent

!----Variables related to new paticles after resample:
real(kind=kind_num),allocatable,dimension(:,:,:) :: key_new
real(kind=kind_num),allocatable,dimension(:) :: nino_new

!----Variables related to 12 months particles:
real(kind=kind_num),allocatable,dimension(:,:) :: particles_12m
real(kind=kind_num),allocatable,dimension(:,:,:) :: parts_12m_100case

integer :: i,j,k,index,cnt
integer :: m_obs,x_obs,y_obs,sele

open(unit=333,file='CCSM4_nml',status='old',action='read',delim='apostrophe')
read(unit=333,nml=mylist)
close(unit=333)

!====INPUT THE MAIN VARIABLES YOU USE IN THIS PROGRAM :
varid_key = varid_sst
varid_nino= varid_nino3
start_key = start_sst
count_key = count_sst
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

allocate( rand1(count_key(1),count_key(2),count_key(3)),&
          rand2(count_key(1),count_key(2),count_key(3)),&
           rand(count_key(1),count_key(2),count_key(3)),&
            obs(count_key(1),count_key(2),count_key(3)))

allocate(   wk_1(count_key(4)) , wk(count_key(4)) )

allocate( xk_new(count_key(4)) , wk_new(count_key(4)) , parent(count_key(4)))

allocate( key_new(count_key(1),count_key(2),count_key(4)) )

allocate( nino_new(count_key(4)) )

allocate( particles_12m(count_key(4),count_key(3)) )
allocate( parts_12m_100case(count_key(4),count_key(3),100) )

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
!do i=1,12
!	call cal_entropy(size(nino,2),nino(i,:),Sq(i))
!enddo
!open(100,file='Sq_nino3.txt')
!write(100,*) Sq
!close(100)

!======Calculate the Variance of NINO index=====
!do i=1,12
!	Var_q(i)= ( sum(nino(i,:)**2) - sum(nino(i,:))**2/size(nino(i,:)) )/( size(nino(i,:))-1 )
!enddo
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
	where(notnan==.False.) key_std(:,:,i)=sqrt(-1.0)
enddo 
!ASSIMILATE AREA OBS YOU NEED TO INCREASE ObsErrStd
!ObsErrStd = 0.4           !~~~~ObsErrStd now have NaNs instead of FillValue~~~~
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
index = 1  ! The year you select as the truth year
write(*,*) "The truth year is ", EPyrs(index)
truth = key(:,:,:,EPyrs(index))
do i=1,count_key(3)
	where(notnan==.False.) truth(:,:,i)=sqrt(-1.0)
enddo  !~~~~truth now have NaNs instead of FillValue~~~~

m_obs = 1  ! The first month to do the assimilation
do cnt=1,100
	call init_random_seed()
	!call randmom_seed()
	call random_number(rand1)
	call random_number(rand2)
	rand=sqrt(-2.0*log(rand1))*cos(2.0*pi*rand2)
	obs=truth+rand*ObsErrStd
	wk_1(:)=1.0/count_key(4)
	!WITH OR WITHOUT TRUTH PARTICLE
	!wk_1(:)=1.0/(count_key(4)-1)
	!wk_1(EPyrs(index))=0.0
	call sis_area(wk_1,obs(:,:,m_obs),ObsErrStd(:,:,m_obs),key(:,:,m_obs,:),wk,neff)
	call resample(count_key(4),nino(m_obs,:),wk,xk_new,wk_new,parent)
	particles_12m(:,1)=xk_new
!open(21,file='parent1.txt')
!write(21,*) parent
!close(21)

	do sele=1,count_key(4)
		key_new(:,:,sele)=key(:,:,m_obs+1,parent(sele))
		nino_new(sele)=nino(m_obs+1,parent(sele))
	end do
	call sis_area(wk_1,obs(:,:,m_obs+1),ObsErrStd(:,:,m_obs+1),key_new,wk,neff)
	call resample(count_key(4),nino_new,wk,xk_new,wk_new,parent)
	particles_12m(:,2)=xk_new
!open(22,file='parent2.txt')
!write(22,*) parent
!close(22)

	do sele=1,count_key(4)
		key_new(:,:,sele)=key(:,:,m_obs+2,parent(sele))
		nino_new(sele)=nino(m_obs+2,parent(sele))
	end do
	call sis_area(wk_1,obs(:,:,m_obs+2),ObsErrStd(:,:,m_obs+2),key_new,wk,neff)
	call resample(count_key(4),nino_new,wk,xk_new,wk_new,parent)
	particles_12m(:,3)=xk_new
!open(23,file='parent3.txt')
!write(23,*) parent
!close(23)

	do i=4,12
		do sele=1,count_key(4)
			particles_12m(sele,i)=nino(i,parent(sele))
		end do
	enddo
!end do
!open(24,file='particles_Apr.txt')
!write(24,*) particles_12m(:,4)
!close(24)

parts_12m_100case(:,:,cnt)=particles_12m
end do

open(12,file='parts_12m_sst_JFM_Eq_west.bin',&
    access='direct',form='unformatted',&
    status='replace',recl=kind_num*count_key(4)*12)
do i=1,100
	write(12,rec=i) parts_12m_100case(:,:,i)
enddo
close(12)

deallocate(key,nino,notnan,sum1,sum2,key_std,ObsErrStd,truth,&
           rand1,rand2,rand,obs,&
           wk_1,wk,xk_new,wk_new,parent,key_new,nino_new,&
           particles_12m,parts_12m_100case)

contains
subroutine check(status)
	integer,intent(in) :: status
	if(status /= nf90_noerr) then
		print *, trim(nf90_strerror(status))
		stop 2
	end if
end subroutine check
end program fengfan
