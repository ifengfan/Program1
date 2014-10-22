#!/bin/bash
 ifort -o readnc readnc.f90\
 -I/cvfs01/disk5/fsqueeze/netcdf_4.3.2/include\
 -L/cvfs01/disk5/fsqueeze/netcdf_4.3.2/lib\
  -lnetcdf -lnetcdff -limf -assume byterecl
