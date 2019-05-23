
# Makefile for blas/level2

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

all: dtrsv cblas_dtrsv dtrsv_f

.c:

	$(CC) $@.c $(CFLAGS) -o $@

dtrsv_f : dtrsv_f.f90

	$(FC) $? $(FFLAGS) -o $@

# Didn't work:
#
#    make: *** No rule to make target 'dtrsv_f', needed by 'all'.  Stop.
#
# .f90:
#
#	$(FC) $@.f90 $(FFLAGS) -o $@

# all: dtrsv cblas_dtrsv dtrsv_f
# 
# dtrsv : dtrsv.c
# 
#     # $(CC) dtrsv.c $(CFLAGS) -o dtrsv
# 	$(CC) $? $(CFLAGS) -o $@
# 
# cblas_dtrsv : cblas_dtrsv.c
# 
#     # $(CC) cblas_dtrsv.c $(CFLAGS) -o cblas_dtrsv
# 	$(CC) $? $(CFLAGS) -o $@
# 
# dtrsv_f : dtrsv_f.f90
# 
#     # $(FC) dtrsv_f.f90 $(FFLAGS) -o dtrsv_f
# 	$(FC) $? $(FFLAGS) -o $@

clean:

	${RM} dtrsv cblas_dtrsv dtrsv_f
    # find . -type f  \( ! -name "*.*" -and ! -name "Makefile" \) # -delete

