# http://www.webalice.it/o.drofa/davide/makefile-fortran/makefile-fortran.html

#  The compiler
FC = gfortran
# flags for debugging or for maximum performance, comment as necessary
FCFLAGS = -g -fbounds-check -Wall  -Wline-truncation  -Wcharacter-truncation  -Wsurprising  -Waliasing  -Wimplicit-interface  -Wunused-parameter  -fwhole-file  -fcheck=all  -std=f2008  -pedantic      -fbacktrace -Wextra -Wconversion -Wconversion-extra 

FCFLAGS = -O2
FCFLAGS = -O2 -fopenmp
# flags forall (e.g. look for system .mod files, required in gfortran)
FCFLAGS += -I/usr/include

# libraries needed for linking, unused in the examples
#LDFLAGS = -li_need_this_lib

# List of executables to be built within the package
PROGRAMS = test1 test2 test3 test4 test5 test6 test6_reverse \
  test6_y test6_y_reverse \
  test2_constHelmholtz

# "make" builds all
all: $(PROGRAMS)

test1: multigrid.o
test2: multigrid.o
test2_constHelmholtz: multigrid.o
test3: multigrid.o
test4: multigrid.o
test5: multigrid.o
test6: multigrid.o
test6_y: multigrid.o
test6_reverse: multigrid.o
test6_y_reverse: multigrid.o

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
cleanout: 
	rm -rf output
clean:  cleanmods 
	rm -f *~ $(PROGRAMS) *.png output/*png
