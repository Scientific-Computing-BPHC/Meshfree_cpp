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
rm debug_main_store_3.txt
rm debug_globaldata.txt
rm debug_globaldata_conn.txt
rm debug_globaldata__x.txt
rm debug_globaldata__y.txt
rm debug_normals.txt
rm debug_globaldata_xpos.txt
rm debug_globaldata_ypos.txt
rm debug_globaldata_xneg.txt
rm debug_globaldata_yneg.txt
rm debug_globaldata_primal.txt
rm debug_globaldata_flux_res.txt
rm debug_globaldata_phi.txt
rm debug_globaldata_dq.txt
rm debug_Gs.txt
rm debug_Gs_again.txt
rm debug_res_sqr.txt
rm debug_state_update.txt
rm checkFileRead.txt
touch debug_output.txt
touch debug_main_store.txt
touch debug_main_store_2.txt
touch debug_main_store_3.txt
touch debug_globaldata.txt
touch debug_globaldata__x.txt
touch debug_globaldata__y.txt
touch debug_normals.txt
touch debug_globaldata_conn.txt
touch debug_globaldata_xpos.txt
touch debug_globaldata_ypos.txt
touch debug_globaldata_xneg.txt
touch debug_globaldata_yneg.txt
touch debug_globaldata_primal.txt
touch debug_globaldata_flux_res.txt
touch debug_globaldata_phi.txt
touch debug_globaldata_dq.txt
touch debug_Gs.txt
touch debug_res_sqr.txt
touch debug_state_update.txt
touch debug_Gs_again.txt
make clean
make 
./clean_meshfree /opt/grids/quadtree/part/partGrid40K 1

