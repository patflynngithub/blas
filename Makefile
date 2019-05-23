
# Makefile for BLAS programs

all: dir_level1 dir_level2 dir_level3

dir_level1:

	$(MAKE) -C level1

dir_level2:

	$(MAKE) -C level2

dir_level3:

	$(MAKE) -C level3

clean:

	$(MAKE) clean -C level1
	$(MAKE) clean -C level2
	$(MAKE) clean -C level3

