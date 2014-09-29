! purpose:
!         output the information about the nc file
! input :
!         the name and path of your intresting file
! output: 
!         the file(information_about_nc) containing the information about the input file _nc   
! Record of revision:
!      date          programmer    description of change
!   ==========       ==========    =====================
!   2012/12/01        lei chen         first realease
! --------------------------------------------------------------------   
! compile and link ,for instance : ifort -o readnc  -I/cvfs01/disk5/fsqueeze/netcdf.intel/include -L/cvfs01/disk5/fsqueeze/netcdf.intel/lib  -lnetcdf -lnetcdff -limf readnc.f90

program readnc
	use netcdf
	implicit none
	integer,parameter :: dbl=selected_real_kind(p=15,r=20)
	character (len = *), parameter :: file_name="CCSM4_trop_1deg_501yr.nc "
	integer:: ncid,i,j,k,l
	integer :: numdims,numvars,numglobalattrs,unlimdimid,numattrs,dimlen,xtype
	real(kind=dbl):: gattrval_nu,attrval_nu
	character (len=999):: gattrname,gattrval,varname,dimname,attrname,attrval
	character(len=99) ::errmsg
	integer,allocatable,dimension(:) :: dimids
	real(kind=dbl),dimension(17) :: lev 
	real(kind=dbl),dimension(128) :: lon 
	real(kind=dbl),dimension(60) :: lat
	open(unit=373,file='info_about_nc',form='formatted',status='replace')
	open(unit=374,file='error_about_nc',form='formatted',status='replace')
	write(*,*)file_name
	call check( nf90_open(file_name, nf90_nowrite, ncid=ncid) )
	write(*,*)"here"
	errmsg=''
	if (trim(errmsg)/='') then
		write(374,*) errmsg
		stop
	end if
	call check(nf90_inquire(ncid=ncid,ndimensions=numdims, nvariables=numvars, nattributes=numglobalattrs,unlimiteddimid=unlimdimid))
	write(373,*) "gloabl attribute:"
	if (numglobalattrs/=0) then
		do i=1,numglobalattrs
			gattrval=" " 
			call check(nf90_inq_attname(ncid=ncid,varid=nf90_global,attnum=i,name=gattrname))
			call check(nf90_get_att(ncid=ncid,varid=nf90_global, name=gattrname,values=gattrval))
			if (trim(errmsg)=='') then
				do l=1,999
					if( iachar(gattrval(l:l)) == 0) then
						gattrval(l:l)=" "
					end if
				end do
				write(373,'(1x,a,a1)') trim(gattrname),"="
				write(373,'(1x,a)') trim(gattrval)
			else 
				call check(nf90_get_att(ncid=ncid,varid=nf90_global, name=gattrname,values=gattrval_nu))
				write(373,'(1x,a,a1)') trim(gattrname),"="
				write(373,'(1x,e12.5)')gattrval_nu
				errmsg=''
			end if
		end do
	end if
	write(373,'(a82)') "----------------------------------------------------------------------------------"
	write(373,*) "variable information:"
	do i=1,numvars
		call check(nf90_inquire_variable(ncid=ncid, varid=i, name=varname, xtype=xtype, ndims=numdims,natts=numattrs))
		write(373,'(1x,a,a,1x,a,i2,1x,a,i2)') "varable=",trim(varname), "id=",i,"variable_data_type=",xtype
		allocate(dimids(numdims))
		call check(nf90_inquire_variable(ncid=ncid, varid=i,dimids=dimids))
		do j=1,numdims
			call check(nf90_inquire_dimension(ncid=ncid,dimid=dimids(j), name=dimname,len=dimlen))
			write(373,'(1x,2a,i6,1x,a)') trim(dimname),"=",dimlen, "dimesnion"
		end do
		deallocate(dimids)
		if (numattrs/=0) then
			write(373,'(1x)') 
			write(373,'(1x,a)') "attributes about the current variable:"
			do k=1,numattrs
				call check(nf90_inq_attname(ncid=ncid, varid=i, attnum=k,name=attrname))
				call check(nf90_inquire_attribute(ncid=ncid, varid=i,name=attrname,xtype=xtype))
				call check(nf90_get_att(ncid=ncid, varid=i, name=attrname, values=attrval))
				write(373,'(1x,2a)') trim(attrname),"="
				if (trim(errmsg)=='') then
					do l=1,999
						if( iachar(attrval(l:l)) == 0) then
							attrval(l:l)=" "
						end if
					end do
					write(373,'(1x,a)') trim(attrval)
					write(373,'(1x,a,i2)') "attribute_value_type=",xtype
				else 
					call check(nf90_get_att(ncid=ncid,varid=i, name=attrname,values=attrval_nu))
					write(373,'(1x,e12.5)') attrval_nu
					write(373,'(1x,a,i2)') "attribute_value_type=",xtype
					errmsg=''
				end if
			end do
		end if
		write(373,'(a82)') "----------------------------------------------------------------------------------"
	end do
	close(unit=373)
	close(unit=374)
	stop

	!-----------------------------------------------------------
	call check(nf90_get_var(ncid, varid=3, values=lev))
	call check(nf90_get_var(ncid, varid=4,values=lat)) 
	call check(nf90_get_var(ncid, varid=6,values=lon))
	open(unit=333,file='level',form='formatted',status='replace')
	write(333,*) 'lev'
	write(333,'(17(1x,e12.5))') (lev(i),i=1,17)
	close(unit=333)
	open(unit=333,file='lat',form='formatted',status='replace')
	write(333,*) 'lat'
	write(333,'(60(1x,e12.5))')(lat(j),j=1,60)
	close(unit=333)
	open(unit=333,file='lon',form='formatted',status='replace')
	write(333,*) 'lon'
	write(333,'(128(1x,e12.5))')(lon(j),j=1,128)
	close(unit=333)
	call check(nf90_close(ncid=ncid)) 
	!----------------------------------------------------------
	
	contains
	subroutine check(status)
		integer, intent ( in) :: status
		if(status /= nf90_noerr) then 
			errmsg=trim(nf90_strerror(status))
		end if
	end subroutine check  
end program readnc 

