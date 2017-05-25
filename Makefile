# Generated automatically from Makefile.in by configure.
##### User configurable options #####
# This is an example Makefile.in (or Makefile configured with mpireconfig)
# for the programs cpi, pi3, and cpilog.  

DIR 	    = /usr/bin

CC          = $(DIR)/mpicc
CCC         = $(DIR)/mpicxx
F77         = $(DIR)/mpif77
CLINKER     = $(DIR)/mpicc
CCLINKER    = $(DIR)/mpicxx
FLINKER     = $(DIR)/mpif77
F90         = $(DIR)/mpif90
F90LINKER   = $(DIR)/mpif90	  
MAKE        = make --no-print-directory
SHELL       = /bin/sh
#

### End User configurable options ###

LIBS =
FLIBS =

default: mandel-parallel

mandel-parallel: mandel-parallel.o
	$(CLINKER) -o mandel-parallel mandel-parallel.o -lm
clean:
	rm -f *.o mandel-parallel
realclean: clean
	rm -f mandel.ppm
