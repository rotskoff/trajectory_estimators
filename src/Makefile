CC = g++
DEBUG = -g
CFLAGS = -Wall $(DEBUG) -std=c++11 -O3 -ffast-math -funroll-loops
INC = -I/usr/local/include
LFLAGS = -L/usr/local/lib  -lm -larmadillo -lhdf5
#LFLAGS = -L/usr/local/lib -L/usr/lib/x86_64-linux-gnu/hdf5/serial/ -lm -larmadillo -lhdf5


quartic: particles.h dos_estimate_quartic.cc
	$(CC) $(CFLAGS) $(INC) -c particles.cc
	$(CC) $(CFLAGS) $(INC) -c dos_estimate_quartic.cc
	$(CC) $(CFLAGS) $(INC) particles.o dos_estimate_quartic.o -o run_quartic $(LFLAGS)

mixture: particles.h dos_estimate_mixture.cc
	$(CC) $(CFLAGS) $(INC) -c particles.cc
	$(CC) $(CFLAGS) $(INC) -c dos_estimate_mixture.cc
	$(CC) $(CFLAGS) $(INC) particles.o dos_estimate_mixture.o -o run_mixture $(LFLAGS)

clean:
	if [ -f run_quartic ]; then rm run_quartic; fi
	if [ -f run_mixture ]; then rm run_mixture; fi
	rm *.o
