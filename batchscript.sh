#!/bin/bash
echo "Patience"
cd src
FILE=meshfree
if test -f "$FILE"; then
    echo "$FILE exists."
    rm meshfree
    echo "$FILE executable removed."
fi
g++ -std=c++0x main.cpp -o meshfree -fpermissive -lm -larmadillo -lgomp -I /home/hari/Work/build/vcpkg/packages/jsoncpp_x64-linux/include

./meshfree /opt/grids/quadtree/part/partGrid40K 10 

