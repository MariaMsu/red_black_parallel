
```
$ mpic++ -o run main.cpp  # build
$ mpiexec -np 25 --map-by :OVERSUBSCRIBE ./run  # run
```
