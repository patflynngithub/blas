
# Makefile for blas/level1


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

all: ddot cblas_ddot ddot_f

.c:

	$(CC) $@.c $(CFLAGS) -o $@

ddot_f : ddot_f.f90

	$(FC) $? $(FFLAGS) -o $@

# Didn't work:
#
#    make: *** No rule to make target 'ddot_f', needed by 'all'.  Stop.
#
# .f90:
#
#	$(FC) $@.f90 $(FFLAGS) -o $@

# all: ddot cblas_ddot ddot_f
# 
# ddot : ddot.c
# 
#     # $(CC) ddot.c $(CFLAGS) -o ddot
# 	$(CC) $? $(CFLAGS) -o $@
# 
# cblas_ddot : cblas_ddot.c
# 
#     # $(CC) cblas_ddot.c $(CFLAGS) -o cblas_ddot
# 	$(CC) $? $(CFLAGS) -o $@
# 
# ddot_f : ddot_f.f90
# 
#     # $(FC) ddot_f.f90 $(FFLAGS) -o ddot_f
# 	$(FC) $? $(FFLAGS) -o $@

clean:

	${RM} ddot cblas_ddot ddot_f
    # find . -type f  \( ! -name "*.*" -and ! -name "Makefile" \) # -delete
