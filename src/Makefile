# src MAkefile

INCLUDE = -I../include -I$(TBBLOC)/include
FLAGS = -O3 -std=gnu99 -DUSING_CILK -D$(CC) $(INCLUDE)

FCFLAGS = -O3
LIB_LAP = libadap_laplace.a
LIB_YUK = libadap_yukawa.a

CSOURCES = fmm-action.c fmm-dag.c fmm-param.c
FSOURCES = ribesl.f dgamma.f
COBJECTS = $(CSOURCES:.c=.o)
FOBJECTS = $(FSOURCES:.f=.o)

ifeq ($(CC), icc)
	CCFLAGS = -pedantic $(FLAGS)
	FC      = ifort
else 
	CCFLAGS = -fcilkplus $(FLAGS)
	FC      = gfortran
endif 

####################################################

lib_lap: 
	$(CC) $(CCFLAGS) -DLAPLACE -c $(CSOURCES)
	ar -cvq $(LIB_LAP) $(COBJECTS)

lib_yuk:
	$(CC) $(CCFLAGS) -DYUKAWA -c $(CSOURCES)
	$(FC) $(FCFLAGS) -c $(FSOURCES)
	ar -cvq $(LIB_YUK) $(COBJECTS) $(FOBJECTS)

clean:
	rm -f *.o *~ $(LIB_LAP) $(LIB_YUK)


