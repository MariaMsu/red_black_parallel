ставим зависимости
```
git clone https://github.com/ICLDisco/ompi
cd ompi      
git submodule update --init --recursive 
./autogen.pl 

sudo apt-get install flex
sudo apt-get install bison
(install by hand) https://github.com/jgm/pandoc/releases/download/2.16.2/pandoc-2.16.2-1-amd64.deb

./configure --with-ft=mpi
make
sudo make install

sudo ldconfig
```

билдим и запускаем
```
$ mpicc -std=c99 -o run_1 redb_2d.c -lm  # build
$ mpirun -np 6 --mca shmem posix --mca opal_event_include poll --map-by :OVERSUBSCRIBE --with-ft ulfm ./run_1
  # running on 6 processors
```
