#!/bin/bash -x
mpisubmit.bg -n 1 -w 00:30:00 run_256
mpisubmit.bg -n 2 -w 00:30:00 run_256
mpisubmit.bg -n 4 -w 00:30:00 run_256
mpisubmit.bg -n 8 -w 00:30:00 run_256
mpisubmit.bg -n 16 -w 00:30:00 run_256
mpisubmit.bg -n 32 -w 00:30:00 run_256
mpisubmit.bg -n 64 -w 00:30:00 run_256
mpisubmit.bg -n 128 -w 00:30:00 run_256
mpisubmit.bg -n 160 -w 00:30:00 run_256

