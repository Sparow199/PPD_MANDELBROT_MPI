# PPD TP MPI MANDELBROT

## Works

The aim of this session is to implement a
Parallel calculation of a fractal mandelbrot.

First version: split into several lines (according to the number of nodes of calculations).


## Use

```sh
$ mpirun -np [number of processors] --hostfile [hostfile] mandel-parallel [-n] [-b] [-d] [-f]
```

- hostfile: host file containing the list of machines
- number iterations [-n int].
- image terminals [-b float float float float].
- image dimensions [-d int int].
- destination file name [-f xxx.ppm]