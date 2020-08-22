#!/bin/bash
echo "Patience"
cd src
g++ main.cpp -o meshfree -fpermissive -lm -larmadillo -lgomp -I /home/hari/Work/build/vcpkg/packages/jsoncpp_x64-linux/include

./meshfree /opt/grids/quadtree/part/partGrid40K 10 

