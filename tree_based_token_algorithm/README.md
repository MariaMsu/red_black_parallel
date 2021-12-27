```

mpic++ -o run main.cpp
mpiexec -np 25 --map-by :OVERSUBSCRIBE ./run

```