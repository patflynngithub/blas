CC=gcc
FC=gfortran
CFLAGS=-lblas
FFLAGS=-lblas

all: ddot cblas_ddot ddot_f

ddot : ddot.c

	$(CC) ddot.c $(CFLAGS) -o ddot

cblas_ddot : cblas_ddot.c

	$(CC) cblas_ddot.c $(CFLAGS) -o cblas_ddot

ddot_f : ddot_f.f90

	$(FC) ddot_f.f90 $(FFLAGS) -o ddot_f


