# Black-to-White hole transition

**The Julia codes are parallelized on the available cores.** It is therefore advisable for the performance to parallelize the codes keeping into account the number of CPU cores present on the system.

A full list of the employed Julia packages can be found in `./inc/pkgs.jl`. **Before executing the source codes, all packages must be installed.**

**The Julia's Just-in-Time compiler is such that the first execution of functions is considerably slower that following ones, and it also allocates much more memory**. To avoid this, you can use the [DaemonMode package](https://github.com/dmolina/DaemonMode.jl).

## Usage

TODO

/home/frisus/Scrivania/julia-1.8.1/bin/julia bin/part_2_vertex_contraction.jl ./ /home/frisus/Scrivania/sl2cfoam-next-dev/data_sl2cfoam
