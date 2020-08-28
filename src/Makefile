VERSION = 0.1

CPP = g++
JSON_PATH = /home/hari/Work/build/vcpkg/packages/jsoncpp_x64-linux/include
LDFLAGS = -lm -larmadillo -lgomp 
CPPFLAGS  = -Wall -O2 -std=c++0x -fpermissive -I $(JSON_PATH)


OBJ = clean_main.o core.o utils.o point.o
make: $(OBJ)
	$(CPP) $(CPPFLAGS) -o clean_meshfree $(OBJ) $(LDFLAGS)

%.o: %.cpp
	$(CPP) $(CPPFLAGS) -c $<

clean:
	@echo "Cleaning up..."
	rm *.o clean_meshfree