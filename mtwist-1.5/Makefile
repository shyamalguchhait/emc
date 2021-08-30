test: mtcctest.o mtwist.o
	g++ -o test mtcctest.o mtwist.o

mtcctest.o: mtcctest.cc
	g++ -c mtcctest.cc

mtwist.o: mtwist.c
	g++ -c mtwist.c

clean:
	rm -f test m*.o