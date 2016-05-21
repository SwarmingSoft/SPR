::example batch file for windows; if not in global path variable, copy dynamically linked librarys of compiler (e.g. libgcc_s_dw2-1.dll and libstdcc++-6.dll for gcc) in .exe folder
::"infile[lengths](string)    outputfile(string)    U0(real)    F(real)    f0(real)    lambda(real)    dt(real)    rmin(real)    passive_frac(real)    t_sim_end(real)    N(+0int)    packing_fraction(real)    seed(+0int)"
bin\Release\SPR.exe sprdata_lengths sprdata 450 1 1 1 0.005 1 0 100 25 0.2 6783
pause