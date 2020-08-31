# q-LSKUM based Meshfree Solver for Inviscid Compressible Flows

Scientific Computing @ BITS Pilani - Hyderabad Campus.

Development of a meshfree solver for inviscid compressible fluid flows in C++. The meshfree solver is based on the Least Squares
Kinetic Upwind Method (q-LSKUM), developed by Deshpande et. al. 

For access, please contact Dr. N. Anil, BITS Pilani.

This work is done as part of my M.Sc. (Hons.) Mathematics Thesis, under the guidance of Dr. N. Anil, Asst. Prof., Department of Mathematics, BITS Pilani - Hyderabad Campus. Nischay Ram Mamidi, Kumar Prasun, Srikanth and Rupanshu have developed the Python, Julia, Fortran and Regent versions respectively.

## Dependencies:
* gcc 8.3.0 or higher
* armadillo 9.900.2

## Instructions to run:

* clone the repo
* `cd src`
* run `./batchscript`

## Directory Structure: 
```
.
+-- Meshfree_cpp
|   +-- src
|      +-- clean_main.cpp
|      +-- clean_main.hpp
|      +-- core.cpp
|      +-- core.hpp
|      +-- utils.cpp
|      +-- utils.hpp
|      +-- split_fluxes.cpp
|      +-- split_fluxes.hpp
|      +-- quadrant_fluxes.cpp
|      +-- quadrant_fluxes.hpp
|      +-- state_update.cpp
|      +-- state_update.hpp
|      +-- flux_residual.cpp
|      +-- flux_residual.hpp
|      +-- limiters.cpp
|      +-- limiters.hpp
|      +-- wall_fluxes.cpp
|      +-- wall_fluxes.hpp
|      +-- point.cpp
|      +-- point.hpp
|      +-- file.cpp
|      +-- file.hpp
|      +-- checkFileRead.txt
|      +-- Makefile   
|   +-- batchscript
|   +-- README.md
|   +-- LICENSE
```
## Roadmap:

- [x] Serial Primal Code
- [ ] Profile and Benchmark Serial Primal
- [ ] Primal CUDA Parallel Version
- [ ] Benchmark Primal CUDA Parallel on V100s
- [ ] Serial Adjoint Code with CoDiPack
- [ ] Benchmark Serial Adjoint Version
- [ ] GPU Parallel Adjoint Solver
- [ ] Benchmark Parallel Adjoint on V100s
- [ ] Begin 3D Serial Primal 


