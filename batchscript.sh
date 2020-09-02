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
cd src
rm debug_output.txt
rm debug_main_store.txt
rm debug_main_store_2.txt
rm debug_globaldata.txt
touch debug_output.txt
touch debug_main_store.txt
touch debug_main_store_2.txt
touch debug_globaldata.txt
make clean
make 
./clean_meshfree /opt/grids/quadtree/part/partGrid40K 2

