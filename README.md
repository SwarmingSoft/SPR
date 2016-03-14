# Self-propelled rods (SPR)

tested with gcc version 4.9.2 (x86_64-win32-seh-rev2, Built by MinGW-W64 project)   
compiled with "-std=c++11", "-D_USE_MATH_DEFINES", "-DNDEBUG"  
tested with "-fexpensive-optimizations", "-O3", "-march=native"   
on Windows 7 and Linux, e.g.  
Windows:  
1) g++.exe -Wall -D_USE_MATH_DEFINES  -fexpensive-optimizations -O3 -march=native -DNDEBUG  -std=c++11 SPR.cpp -o SPR.o  
2) g++.exe -o SPR.exe SPR.o  
Linux:  
1) g++ -Wall -D_USE_MATH_DEFINES -fexpensive-optimizations -O3 -march=native -DNDEBUG -std=c++11 SPR.cpp -o SPR.o

no additional dependencies

call executable as  
Windows:  
3) SPR.exe sprdata_length sprdata U0 F f0 lambda dt rmin passive_frac t_sim_end N packing_frac seed  
Linux:  
2) /.SPR.o sprdata_length sprdata U0 F f0 lambda dt rmin passive_frac t_sim_end N packing_frac seed  

where (standard value in brackets)  
sprdata_length (sprdata_length): input file "sprdata_length.txt" with space seperated lengths, e.g. "8 8 8 8 8 8" for N=6 rods a length l=5  
sprdata (sprdata): output file "sprdata.txt"  
U0 (450): Yukawa potential amplitude  
F (1): driving force  
f0 (1): Stokesian friction  
lambda (1): width of rods  
dt (0.005): euler integrator step size  
rmin (1): cutoff of potential  
passive_frac (0.): fraction of passive particles  
t_sim_end (20000): timesteps till termination  
N (1000): number of particles  
packing_fraction (0.5): aimed for packing fraction between 0. and 1.  
seed (477): random number generator seed (integer)  
