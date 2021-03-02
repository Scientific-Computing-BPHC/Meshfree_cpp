VERSION = 0.1

CPP = g++
JSON_PATH = /home/hari/Work/build/vcpkg/packages/jsoncpp_x64-linux/include
ARMA_PATH = /opt/spack/opt/spack/linux-centos8-broadwell/gcc-10.2.0/armadillo-9.800.3-wfj72pbrwiwikrwa5end6zmwqg7qwpi3/include
LDFLAGS = -lm -lgomp 
CPPFLAGS  = -Wall -pg -O2 -std=c++0x -fpermissive -I $(JSON_PATH) -I $(ARMA_PATH)


OBJ = main.o core.o utils.o point.o state_update.o flux_residual.o wall_fluxes.o limiters.o quadrant_fluxes.o split_fluxes.o outer_fluxes.o interior_fluxes.o
make: $(OBJ)
	$(CPP) $(CPPFLAGS) -o meshfree $(OBJ) $(LDFLAGS)

%.o: %.cpp
	$(CPP) $(CPPFLAGS) -c $<

clean:
	@echo "Cleaning up..."
	rm *.o meshfree