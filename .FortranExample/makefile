# Make file for the random laser (diffusion) code

SOURCES = variables.f90 prog.f90 pops.f90 matinv.f90 intens.f90 tridiag.f90

PROGRAM = prog

OBJECTS = $(SOURCES:.f90=.o)

.KEEP_STATE:
.SUFFIXES:
.SUFFIXES: .o .f90 .h

$(PROGRAM):	$(OBJECTS)
	f90 -C -o $(PROGRAM) $(OBJECTS) 

.f90.o: 
	f90 -C -c -o $@ $<

# f90 -C -c -g -o $@ $<
# f90 -C -g -o $(PROGRAM) $(OBJECTS) -multigot  bounds check, and debugging, noopt
 

