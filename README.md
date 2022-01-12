# Parallel Version of the Program for the Red-Black2D Task 
# && 
# Tree based token algorithm
course "Supercomputers and Parallel Data Processing" && 
course "Distributed systems"

## OpenMP version
```
v_openmp$ gcc -std=c99 -Wall -fopenmp -o run redb_2d.c  # build
v_openmp$ OMP_NUM_THREADS=8 ./run  # run
```
Time of execution
![Time of execution](./openmp_time.jpg)

## MPI version
```
v_mpi$ mpicc -std=c99 -o run_1 redb_2d.c -lm  # build
v_mpi$ mpisubmit.bg -n 8 -w 00:30:00 run # run
```
Time of execution
![Time of execution](./mpi_time.jpg)

### [Report for Supercomputers and Parallel Data Processing](./Отчет.%20Суперкомпьютеры%20и%20параллельная%20обработка%20данных.pdf)

## MPI stable version
This version parallelizes the task in the same way as "MPI version",
but if one of the processes will be killed, then the program can redistribute calculations 
among the left processes and continue to work.  
```
v_mpi_stable$ mpicc -std=c99 -o run_1 redb_2d.c -lm  # build
v_mpi_stable$ mpirun -np 6 --mca shmem posix --mca opal_event_include poll --map-by :OVERSUBSCRIBE --with-ft ulfm ./run_1
  # running on 6 processors
```

## Tree based token algorithm
The calculation synchronization algorithm 
```
tree_based_token_algorithm$ mpic++ -o run main.cpp  # build
tree_based_token_algorithm$ mpiexec -np 25 --map-by :OVERSUBSCRIBE ./run  # run
```

### [Report for Distributed systems](./Отчет.%20СРаспределенные%20Ссистемы.pdf)
