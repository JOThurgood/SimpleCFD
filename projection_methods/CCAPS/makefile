# http://www.webalice.it/o.drofa/davide/makefile-fortran/makefile-fortran.html

#  The compiler
FC = gfortran
# flags for debugging or for maximum performance, comment as necessary
FCFLAGS = -g -fbounds-check -Wall  -Wline-truncation  -Wcharacter-truncation  -Wsurprising  -Waliasing  -Wimplicit-interface  -Wunused-parameter  -fwhole-file  -fcheck=all  -std=f2008  -pedantic      -fbacktrace -Wextra -Wconversion -Wconversion-extra 

FCFLAGS = -O2
# flags forall (e.g. look for system .mod files, required in gfortran)
FCFLAGS += -I/usr/include

# libraries needed for linking, unused in the examples
#LDFLAGS = -li_need_this_lib

# List of executables to be built within the package
PROGRAMS = ccaps 

# "make" builds all
all: $(PROGRAMS)


advance_module.o: boundary_conditions.o diagnostics.o
setup.o: control.o initial_conditions.o
ccaps.o: shared_data.o advance_module.o setup.o control.o initial_conditions.o \
  boundary_conditions.o diagnostics.o welcome.o
ccaps: shared_data.o  advance_module.o setup.o control.o initial_conditions.o \
  boundary_conditions.o diagnostics.o welcome.o


# ======================================================================
# And now the general rules, these should not require modification
# ======================================================================

# General rule for building prog from prog.o; $^ (GNU extension) is
# used in order to list additional object files on which the
# executable depends
%: %.o
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)

# General rules for building prog.o from prog.f90 or prog.F90; $< is
# used in order to list only the first prerequisite (the source file)
# and not the additional prerequisites such as module or include files
%.o: %.f90
	$(FC) $(FCFLAGS) -c $<

%.o: %.F90
	$(FC) $(FCFLAGS) -c $<

# Utility targets
.PHONY: clean cleanall 

cleanmods:
	rm -rf *.o *.mod *.MOD *.dSYM *.dat
clean:  cleanmods
	rm -f *~ $(PROGRAMS) *.png
