#!/bin/bash
# echo "Patience"
# cd src
# FILE=clean_meshfree
# if test -f "$FILE"; then
#     echo "$FILE exists."
#     rm clean_meshfree
#     echo "$FILE executable removed."
# fi
# g++ -std=c++0x clean_main.cpp -o clean_meshfree -fpermissive -lm -larmadillo -lgomp -I /home/hari/Work/build/vcpkg/packages/jsoncpp_x64-linux/include

# ./clean_meshfree /opt/grids/quadtree/part/partGrid40K 10 
# #./meshfree /opt/grids/quadtree/part/partGrid40K 10 

echo "Patience"
make clean
make 
#taskset -c 0 ./clean_meshfree /opt/grids/quadtree/part/partGrid800K 1000 &
taskset -c 1 ./clean_meshfree /opt/grids/quadtree/part/partGrid2.5M 100 &
taskset -c 2 ./clean_meshfree /opt/grids/quadtree/part/partGrid2.5M 1000

