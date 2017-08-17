# src MAkefile

INCLUDE = -I../include -I$(TBBLOC)/include 
LIBS = -L$(TBBLIB) -ltbbmalloc -lcilkrts -lm
FLAGS = -O3 -std=gnu99 -DUSING_CILK -D$(CC) $(INCLUDE) 
FCFLAGS = -O3

EXEC = fmm

ifeq ($(CC), icc)
	CCFLAGS = -Wall -pedantic $(FLAGS)
	FC = ifort 
	LINKER = $(FC) -nofor_main
else 
	CCFLAGS = -fcilkplus $(FLAGS) 
	FC = gfortran
	LINKER = $(CC)
endif 

############################################################

lap:
	$(CC) $(CCFLAGS) -DLAPLACE -c fmm-test.c -o fmm-test.o 
	$(CC) -o $(EXEC) fmm-test.o $(LIBS) -L../lib -ladap_laplace

laptest3:
	./$(EXEC) -d 1 -a 3
	./$(EXEC) -d 2 -a 3

laptest6:
	./$(EXEC) -d 1 -a 6
	./$(EXEC) -d 2 -a 6

yuk:
	$(CC) $(CCFLAGS) -DYUKAWA -c fmm-test.c -o fmm-test.o
	$(LINKER) -o $(EXEC) fmm-test.o $(LIBS) -L../lib -ladap_yukawa	

yuktest3:
	./$(EXEC) -d 1 -a 3
	./$(EXEC) -d 2 -a 3

yuktest6:
	./$(EXEC) -d 1 -a 6
	./$(EXEC) -d 2 -a 6

clean:
	rm -f $(EXEC) *.o *~
