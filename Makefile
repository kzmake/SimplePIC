# Makefile

OPT = -std=c++14 -O3 -Wall -lboost_timer -lboost_system

all:
	mpic++ $(OPT) main.cpp

run-1:
	make makedir 1>/dev/null 2>/dev/null 
	mpirun -np 1 a.out

run-2:
	make makedir 1>/dev/null 2>/dev/null 
	mpirun -np 2 a.out

run-4:
	make makedir 1>/dev/null 2>/dev/null 
	mpirun -np 4 a.out

run-6:
	make makedir 1>/dev/null 2>/dev/null 
	mpirun -np 6 a.out

run-12:
	make makedir 1>/dev/null 2>/dev/null 
	mpirun -np 12 a.out

run:
	make run-1

clean:
	rm -rf a.out data/

makedir:
	mkdir -p data
	mkdir -p data/f_e
	mkdir -p data/f_b
	mkdir -p data/f_j


