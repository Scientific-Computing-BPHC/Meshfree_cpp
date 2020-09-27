# q-LSKUM based Meshfree Solver for Inviscid Compressible Flows

Scientific Computing @ BITS Pilani - Hyderabad Campus.

Development of a meshfree solver for inviscid compressible fluid flows in C++. The meshfree solver is based on the Least Squares
Kinetic Upwind Method (q-LSKUM), developed by Deshpande et. al. 

This work is done as part of my M.Sc. (Hons.) Mathematics Thesis, under the guidance of Dr. N. Anil, Assistant Professor, Department of Mathematics, BITS Pilani - Hyderabad Campus, and Professor S.M. Deshpande, JNCASR. 

## Dependencies:
* gcc 8.3.0 or higher
* armadillo 9.900.2
* CUDA 11

## Usage:

* Configure the parameters through the Config struct in core.hpp or cuda_core.hpp
* chmod +x batchscript.sh (within src/serial or src/CUDA)
* run `./batchscript` 

## Directory Structure: 
```
.
+-- Meshfree_cpp
|   +-- src
|       +-- serial
|           +-- clean_main.cpp
|           +-- clean_main.hpp
|           +-- core.cpp
|           +-- core.hpp
|           +-- utils.cpp
|           +-- utils.hpp
|           +-- split_fluxes.cpp
|           +-- split_fluxes.hpp
|           +-- quadrant_fluxes.cpp
|           +-- quadrant_fluxes.hpp
|           +-- state_update.cpp
|           +-- state_update.hpp
|           +-- flux_residual.cpp
|           +-- flux_residual.hpp
|           +-- limiters.cpp
|           +-- limiters.hpp
|           +-- wall_fluxes.cpp
|           +-- wall_fluxes.hpp
|           +-- point.cpp
|           +-- point.hpp
|           +-- Makefile   
|           +-- batchscript.sh
|       +-- CUDA
|           +-- main_cuda.cu
|           +-- main_cuda.hpp
|           +-- core_cuda.cu
|           +-- core_cuda.hpp
|           +-- utils.cpp
|           +-- utils.hpp
|           +-- split_fluxes_cuda.cu
|           +-- split_fluxes_cuda.hpp
|           +-- quadrant_fluxes_cuda.cu
|           +-- quadrant_fluxes_cuda.hpp
|           +-- state_update_cuda.cu
|           +-- state_update_cuda.hpp
|           +-- flux_residual_cuda.cu
|           +-- flux_residual_cuda.hpp
|           +-- limiters_cuda.cu
|           +-- limiters_cuda.hpp
|           +-- wall_fluxes.cu
|           +-- wall_fluxes.hpp
|           +-- point.cpp
|           +-- point.hpp
|           +-- Makefile   
|           +-- batchscript.sh
|   +-- README.md
|   +-- LICENSE
```
## Roadmap:

- [x] Serial Primal Code
- [x] Profile Serial Primal
- [X] Primal CUDA Parallel Version v1 
- [X] Primal CUDA Parallel Version v2 
    - [X] Reductions
    - [X] Shared Memory
    - [X] CUDA Graphs
- [ ] Serial Adjoint Code with CoDiPack
- [ ] Optimize and Benchmark Serial Primal once again
- [ ] CUDAfy CoDiPack based Adjoint Solver
- [ ] Begin 3D Serial Primal (stretch goal for thesis)

## Support:

Please contact the author for running the code, and Dr. N. Anil for access to the input grids.

## Author:

Harivallabha Rangarajan, Department of Mathematics and Computer Science, BITS Pilani - Hyderabad. 

