# spar Benchmarks

* `reduce_bandish` - 
* `reduce_rand` - 


## Compiling

`make`


## Flags

All benchmarks

* `-d` - print header
* `-r rank`
    - rank to reduce to (-1 for allreduce)
    - default 0
* `-n order`
    - the number of rows/cols
    - default 5000

`reduce_bandish`:

`reduce_rand`:

* `-a` - Should the approx (see )
* `-p prop`
    - proportion dense
    - default 0.001



## Running the Benchmarks

All benchmarks should be launched with `mpirun` with at least 2 ranks.

Benchmarks are meant to be run across different numbers of MPI ranks as well as various other parameters (matrix size, reduce vs allreduce, etc.). These timings should take place in separate runs. Each run will produce the line for a CSV.

For example, you might do something like:

```
$ mpirun -np 2 ./reduce_rand -n 5000 -p 0.0001 -d
$ for n in `seq 3 6`; do mpirun -np $n ./reduce_rand -n 5000 -p 0.0001; done
```

Which on a desktop produces

```
benchmark,size,seed,densevec,root,n,prop_dense,bytes_index,bytes_scalar,nnz_local,len_local,time_gen,nnz,len,time_reduce
reduce_rand,2,1234,0,0,5000,0.000100,4,4,2500,2500,0.429878,2500,3264,0.006539
reduce_rand,3,1234,0,0,5000,0.000100,4,4,2500,2500,0.430399,2500,3264,0.012798 
reduce_rand,4,1234,0,0,5000,0.000100,4,4,2500,2500,0.428284,2500,3264,0.012809 
reduce_rand,5,1234,0,0,5000,0.000100,4,4,2500,2500,0.427958,2500,3264,0.018075 
reduce_rand,6,1234,0,0,5000,0.000100,4,4,2500,2500,0.503866,2500,3264,0.019085 
```

The first two lines (header and first output line) are produced by the first run. The other 4 lines are produced by the for loop.
