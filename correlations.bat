::correlations.exe str_path str_filename str_filename_ids b_ids_from_file b_read_flow b_periodic_boundary b_periodic_velocities real_L real_cluster_com_dist real_cluster_angle_dist real_vel_histo_start real_vel_histo_end real_vel_histo_step real_vor_histo_start real_vor_histo_end real_vor_histo_step real_rmin real_rmax real_rstep real_tstart real_tend real_t_step uint_Xvortstart uint_Xvortend uint_Xvortstep uint_times_start uint_times_end uint_times_step uint_vorticity_bins uint_FlowGridPoints
::examples for different types of input:
::flow (from opticalflow.cpp)
::bin/Release/correlations.exe "" "" 0 1 0 0 1024. 24. 0.349 -0.1 11. 0.1 -0.46 0.46 0.01 0. 1024. 16. 0 50 1 0 10 1 1 51 1 64 64
::vicsek (from vicsek.cpp)
::bin/Release/correlations.exe "" "metric.txt" "" 0 0 1 1 16. 0.75 0.349 0.0 0.02 0.001 -0.04 0.04 0.001 0. 16. 0.25 0 500 1 0 10 1 1 1999 1 64 64
::SPR (from SPR.cpp)
::bin/Release/correlations.exe "" "sprdata_N1000_l7_pf0.5.txt" "" 0 0 1 1 118. 10.5 0.349 -0.001 1. 0.001 -0.34 0.34 0.01 0. 118. 1.84 0 500 1 0 10 1 1 1001 1 64 64
::splines (experimental data from spline_generate.py)
::bin/Release/correlations.exe "" "data_exp_tracked_splines.txt" "data_exp_tracked_ids.txt" 1 0 0 0 512. 75. 0.349 -0.1 40. 0.1 -2.2 2.2 0.01 0. 512. 16. 0 200 1 0 10 1 1 1999 1 64 64


bin\Release\correlations.exe "" "sprdata.txt" "" 0 0 1 1 26. 5.0 0.349 -0.001 1. 0.001 -0.34 0.34 0.01 0. 26. 0.26 0 10 1 0 10 1 1 10 1 64 64
pause