#! /bin/bash
i=0
while ((i<$1))
    do 
	mpirun -np 4 ./build/mc-pricer ./src/perf.dat 0.1 >> result.txt
	((i=i+1))
    done


