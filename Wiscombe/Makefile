# Makefile for Warren J. Wiscombe's Mie scattering code
# ftp://climate1.gsfc.nasa.gov/wiscombe/Single_Scatt/Homogen_Sphere/Exact_Mie/

objects = ErrPack.o MIEV0.o
target = test
objectst = ErrPack.o MIEV0.o testcc.o
test : $(objectst)
	g++ $(objectst) -o $@  -lgfortran
#-lgfortran

$(objects) : %.o: %.f
	gfortran -c $< -o $@

testcc.o: testcc.cpp
	g++ -c testcc.cpp

clean:
	rm -f test *.o 
	rm -f data.dat
	clear
