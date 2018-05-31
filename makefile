all: test_Prior.exe Distance_calc_wGRVS.exe

Prior.o : Prior.cc Prior.h
	g++ -c -o $@ -O3 -ffast-math $<

%.exe : %.cc Prior.o
	g++ -o $@ -O3 -ffast-math  $< Prior.o
