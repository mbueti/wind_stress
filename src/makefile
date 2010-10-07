
SRC = MAIN.f windsubs.f 
OBJ = ${SRC:.f=.o}
FC  = /opt/intel/Compiler/11.1/064/bin/intel64/ifort 
LDFLAGS =   
LIBS = 
CMD =  ./windonly
FFLAGS = -O3  -convert big_endian
.SUFFIXES: .o .f .f90

.f.o:
	$(FC) $(FFLAGS) -c $*.f

# Lines from here on down should not need to be changed.  They are the
# actual rules which make uses to build field.
#
all:            $(CMD)
$(CMD): $(OBJ)
	$(FC) $(FFLAGS) -o $(@) $(OBJ) $(LIBS)
	test -d ../exe || mkdir ../exe
	mv windsubs.o ../exe/
	mv MAIN.o ../exe/
	mv windonly ../exe/

clean:
	rm -f $(OBJ) 

