CC      = icc
CCFLAGS = -std=gnu99 -Wall -pedantic -DLAPLACE -DUSING_CILK -g
INCLUDE = -I../include
LIBS    = -lcilkrts -ltbbmalloc -L../lib -ladap_laplace
EXEC    = example

all:
	$(CC) $(CCFLAGS) -c $(EXEC).c -o $(EXEC).o $(INCLUDE)
	$(CC) -o $(EXEC) $(EXEC).o $(LIBS)

clean:
	rm -f *.o $(EXEC) *~
