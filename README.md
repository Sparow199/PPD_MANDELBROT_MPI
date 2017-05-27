# PPD TP MPI MANDELBROT


![MANDEL](https://cloud.githubusercontent.com/assets/22281426/26521662/19d0df30-42ee-11e7-8676-4d2cc0593fb8.png)

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
