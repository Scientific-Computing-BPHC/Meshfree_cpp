VERSION = 0.1
CC = nvcc -arch compute_61 
CPP = g++
JSON_PATH = /home/hari/Work/build/vcpkg/packages/jsoncpp_x64-linux/include
LDFLAGS = -lm -larmadillo -lgomp 
# CPPFLAGS  = -Wall -pg -O2 -std=c++0x -fpermissive -I $(JSON_PATH)
CFLAGS  = -std=c++11 -Xcompiler -fopenmp 
CPPFLAGS  = -Wall -pg -O2 -fpermissive -I $(JSON_PATH)


OBJ = clean_main.o core.o utils.o point.o state_update.o flux_residual.o wall_fluxes.o limiters.o quadrant_fluxes.o split_fluxes.o outer_fluxes.o interior_fluxes.o
# make: $(OBJ)
# 	$(CPP) $(CPPFLAGS) -o clean_meshfree $(OBJ) $(LDFLAGS)

make: $(OBJ)
	$(CC) $(CFLAGS) -o clean_meshfree $(OBJ) $(LDFLAGS)

%.o: %.cpp
	$(CPP) $(CPPFLAGS) -c $<

%.o: %.cu
	$(CC) $(CFLAGS) -c $<

clean:
	@echo "Cleaning up..."
	rm *.o clean_meshfree