# release Makefile

#CC = icc
CC = clang
#CC = gcc

# TBBLOC = /opt/intel/composer_xe_2015/tbb
# TBBLIB = $(TBBLOC)/lib/intel64/gcc4.4
TBBLOC = /usr/local/tbb44_20151115oss
TBBLIB = $(TBBLOC)/lib

export

all:
	make clean laplace
	make clean yukawa

laplace:
	cd ./src; make lib_lap
	mkdir -p lib
	cp ./src/libadap_laplace.a ./lib/
	cd ./test; make clean lap laptest3 laptest6

yukawa:
	cd ./src; make lib_yuk
	mkdir -p lib
	cp ./src/libadap_yukawa.a ./lib/
	cd ./test; make clean yuk yuktest3 yuktest6

clean:
	(cd ./src; make clean)
	(cd ./example; make clean)
	(cd ./test; make clean)
	(cd ./lib; rm -f *.a)
	rm -rf *~ ./lib

testall:
	make all | grep "\*\*\*"
