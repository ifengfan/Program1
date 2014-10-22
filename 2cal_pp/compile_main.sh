#!/bin/bash
ifort -o main_JFM cal_entropy.f90 sis_1point.f90 resample.f90 \
 init_random_seed.f90 main_JFM.f90 \
 -I/cvfs01/disk5/fsqueeze/netcdf_4.3.2/include \
 -L/cvfs01/disk5/fsqueeze/netcdf_4.3.2/lib \
 -lnetcdf -lnetcdff -limf -assume byterecl
