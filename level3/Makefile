
# Makefile for blas/level3

# $@ is the name of the file to be made.
# $? is the name of the changed dependents.
# $< the name of the related file that caused the action.
# $* the prefix shared by target and dependent files.

CC=gcc
FC=gfortran

CFLAGS=-lblas
FFLAGS=-lblas

RM=/bin/rm -f

# -------------------------

all: dgemm cblas_dgemm dgemm_f

.c:

	$(CC) $@.c $(CFLAGS) -o $@

dgemm_f : dgemm_f.f90

	$(FC) $? $(FFLAGS) -o $@

# Didn't work:
#
#    make: *** No rule to make target 'dgemm_f', needed by 'all'.  Stop.
#
# .f90:
#
#	$(FC) $@.f90 $(FFLAGS) -o $@

# all: dgemm cblas_dgemm dgemm_f
# 
# dgemm : dgemm.c
# 
#     # $(CC) dgemm.c $(CFLAGS) -o dgemm
# 	$(CC) $? $(CFLAGS) -o $@
# 
# cblas_dgemm : cblas_dgemm.c
# 
#     # $(CC) cblas_dgemm.c $(CFLAGS) -o cblas_dgemm
# 	$(CC) $? $(CFLAGS) -o $@
# 
# dgemm_f : dgemm_f.f90
# 
#     # $(FC) dgemm_f.f90 $(FFLAGS) -o dgemm_f
# 	$(FC) $? $(FFLAGS) -o $@

clean:

	${RM} dgemm cblas_dgemm dgemm_f
    # find . -type f  \( ! -name "*.*" -and ! -name "Makefile" \) # -delete

