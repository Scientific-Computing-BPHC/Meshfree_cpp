VERSION = 0.1
CC = nvcc 
CPP = g++
JSON_PATH = /home/hari/Work/build/vcpkg/packages/jsoncpp_x64-linux/include
LDFLAGS = -lm -larmadillo -lgomp 
# CPPFLAGS  = -Wall -pg -O2 -std=c++0x -fpermissive -I $(JSON_PATH)
CUDA_PATH = /home/hari/Work/Common
CFLAGS  = -std=c++11 -Xcompiler -fopenmp -I $(CUDA_PATH)
CPPFLAGS  = -Wall -pg -O2 -fpermissive -I $(JSON_PATH) -I $(CUDA_PATH)

# ## Serial
# OBJ = clean_main.o core.o utils.o point.o state_update.o flux_residual.o wall_fluxes.o limiters.o quadrant_fluxes.o split_fluxes.o outer_fluxes.o interior_fluxes.o 
# make: $(OBJ)
# 	$(CPP) $(CPPFLAGS) -o clean_meshfree $(OBJ) $(LDFLAGS)

# Parallel
OBJ = main_cuda.o utils.o point.o state_update_cuda.o flux_residual_cuda.o wall_fluxes_cuda.o limiters_cuda.o quadrant_fluxes_cuda.o split_fluxes_cuda.o outer_fluxes_cuda.o interior_fluxes_cuda.o core_cuda.o

make: $(OBJ)
	$(CC) $(CFLAGS) -o clean_meshfree $(OBJ) $(LDFLAGS)

%.o: %.cpp
	$(CPP) $(CPPFLAGS) -c $<

%.o: %.cu
	$(CC) $(CFLAGS) -dc -c $<

clean:
	@echo "Cleaning up..."
	rm *.o clean_meshfree