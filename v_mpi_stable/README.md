```
$ mpicc -std=c99 -o run_1 redb_2d.c -lm  # build
$ mpiexec -np 4 ./run_1  # running on 4 processors
```
