Program1
========

Read in NetCDF files and do the SIS. 

1. rm readnc info_about_nc error_about_nc

2. input the NetCDF filename and path you need in line 18 of readnc.f90
   then compile and link: ./compile_readnc.sh

3. run ./readnc
   You'll get file "info_about_nc" and "error_about_nc"

4. Create your own NAMELIST file "ModelName_nml" according to "info_about_nc"
