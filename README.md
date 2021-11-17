# Shearsort
Shearsort (Parallel and Distributed Programming Uppsala University)

## Purpose
The purpose of this project is to implement a parallel Shear sort algorithm in C with the help of Message Passing Interface (MPI) and evaluate its performance through strong and weak scaling experiments. Access to a super-computer was available for the performance experiments.

## How to run
How to run the file shearsort.c :
1. Produce an executable by "make shearsort".
2. Run the program with "mpirun -np p ./shearsort input-file output-file" 
where p is the number of processes. An example input file is attached.
3. Delete the executable by "make clean".
