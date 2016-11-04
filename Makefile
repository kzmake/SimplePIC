# Makefile

OPT = -std=c++14 -O3 -Wall -lboost_timer -lboost_system

all:
	make mpi

single:
	g++ $(OPT) main.cpp

mpi:
	mpic++ $(OPT) -DMPI_PIC main.cpp 

run:
	make makedir 1>/dev/null 2>/dev/null 
	mpirun -np 4 a.out

clean:
	rm -rf a.out data/

makedir:
	mkdir -p data
	mkdir -p data/f_e
	mkdir -p data/f_b
	mkdir -p data/f_j


