# parallel-md-active

Author: Parswa Nath [Dept. of Physics, IIT Hyderabad]
License: GPL2

This is a massively parallel Molecular Dynamics (MD) Simulator for Active systems. I wrote this simulator while preparing my M.Sc. thesis with a hope that this might be used as a standard MD simulator for active systems 

FEATURES
* Models activity through a ABP model by Marchetti & Hagan (see http://arxiv.org/pdf/1207.1737.pdf). 
* Can be extened for other ABP models as well (Vischek model for instance).
* Parallelization is achived by decomposing the domain (so called "cell list" technique) and is handled by MPI. 
* Scalable and massively parallel

TODO
* Code cleanup
* Add some shared memory parallelism using openMP.
* Is it a really good idea to use the vector container from C++ STL for high performance codes?
* Improvement of I/O patterns
* Modification of existing data structure for particle to reduce disk usage. 

COMPILING & BUILDING 
* You will need MPICH libraries to compile this. In ubuntu, get it with sudo apt-get install mpich2
* mpic++ *.cpp -o md.x -std=c++11 -O3 -march=native

SAMPLE USAGE
* mpirun -np 48 ./md.x < some_input.txt > some_output.txt

INPUT FORMAT
```

** parallel MD input **

n_blocks_x -> 4;

start_mode -> "load"; 
input_data_path -> "/media/theo_xubuntu/old_backup/128f_eq";

output_data_path -> "try_io"; 

nx -> 32;
ny -> 32;

rho -> 1.1;
sigma -> 1.0;
rc -> 1.122;
Pe -> 99.9994;

is_confined -> 0;

max_time_steps -> 10000;
dt -> 1e-5;
config_save_interval -> 20;
wait_for_equilibration -> 10;

debug_type -> "none";
```
