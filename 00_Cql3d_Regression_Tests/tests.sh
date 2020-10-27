#~/cql3d/git/cql3d/00_Cql3d_Regression_Tests/tests

#Adjust path to the cql3d executable.  For example:
XCQL3D="../../xcql3d_gfortran64.1"


#test1: Runaway electron rates vs rho
mkdir test1
cp cqlinput_eoved_.16_multi-flux-surface_test1.0 test1/
cp cqlinput_eoved_.16_multi-flux-surface_test1.1 test1/
cd test1/
cp cqlinput_eoved_.16_multi-flux-surface_test1.0 cqlinput
time $XCQL3D > log_test1.0    # 3.3s on compx2
cp cqlinput_eoved_.16_multi-flux-surface_test1.1 cqlinput
time $XCQL3D > log_test1.1    # 1.4s on compx2
cd ..


#test2: DC electric field resistivity vs rho
mkdir test2
cp cqlinput_freidberg_full_eps.3_short_w_sxr_201024 test2/cqlinput
cd test2/
ln -s ../eqdsk_freidberg_full_eps
time $XCQL3D > log_test2   #4.7s on compx2
cd ..


#test3: D3D shot 96143 one ray ECH test case
mkdir test3
cp cqlinput_96143_one_ray.1_short_201024 test3/cqlinput
cd test3/
ln -s ../g096143.01440
ln -s ../rayech_test3 rayech
time $XCQL3D > log_test3   #6.3s on compx2
cd ..

#test4: MAST EBW OXB test case.  Ray data in companion genray test.
mkdir test4
cp cqlinput_MAST_test.0 test4/cqlinput
cp eqdsk_MASTU test4
cd test4/
#Next line is assuming genray and cql3d distribution are both under a common
#directory.
ln -s ../../../genray/00_Genray_Regression_Tests/test10/genray.nc genray_test4.nc
time $XCQL3D > log_test4   #32s on compx2
cd ..

#test5  multiURF_multiSpecies+NBI test case with D+H plasma in DIII-D
mkdir test5
cp cqlinput_H0.short_mmsv8_adjust.4 test5/cqlinput
cd test5/
ln -s ../genrayfw_18rays.nc
ln -s ../g122080.03100
time $XCQL3D > log_test5   # 33s on compx2
cd ..

