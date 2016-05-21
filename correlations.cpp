//v 21.05.2016

/*
passed over arguments are:

str_path - string of path with imputfiles ("" for same folder as executable)
str_filename - string of input file (including .txt, eg. "sprdata.txt")
str_filename_ids - string of input file with id (only important for experimental data, "" otherwise)
b_ids_from_file - 1: if ids from input file (above), 0: otherwise
b_read_flow - 1: if reading in input file from opticalflow.cpp, 0: otherwise
b_periodic_boundary - 1: periodic boundaries, 0: no periodic boundaries
b_periodic_velocities - 1: periodic boundaries, 0: no periodic boundaries (kind of redundant)
real_L real_cluster_com_dist - real (= float/double, depending on declaration) of absolute center of mass distance for two particles to be considered in same cluster
real_cluster_angle_dist - real angle in radiant between two particle orientations to be considered in same cluster (0.349 is 20°)
real_vel_histo_start - lower boundary for velocities in velocity histogram
real_vel_histo_end - upper boundary for velocities in velocity histogram
real_vel_histo_step - resolution (bin width) for velocities in velocity histogram
real_vor_histo_start - lower boundary for vorticity in velocity histogram
real_vor_histo_end - upper boundary for vorticity in velocity histogram
real_vor_histo_step - resolution (bin width) for vorticity in velocity histogram
real_rmin - lower boundary for spatial distances in correlation functions
real_rmax - upper boundary for spatial distances in correlation functions
real_rstep - resolution (bin width) for spatial distances in correlation functions
real_tstart - time step where to start calculation of everything (except vorticities)
real_tend - time step where to end calculation of everything (except vorticities)
real_t_step - time step width for calculation of everything (except vorticities)
uint_Xvortstart - as above only for vorticity (choose lower, because more computation time)
uint_Xvortend - as above
uint_Xvortstep - as above
uint_times_start - what time steps to read in starting from this (the 6 parameters above are relative to this, so if you read in say 15000-17000 in steps of 1 here, you can only use 0 to 2000 above)
uint_times_end - what time steps to read in till this
uint_times_step - what time steps to read in, spacing
uint_vorticity_bins - number of bins (equally spaced) for spatial vorticity calculations
uint_FlowGridPoints  - number of flow grid points (only interesting for opticalflow.cpp output)
*/

#include "correlations.hpp"
#include <fstream>
#include <cmath>
#include <stdexcept>

typedef unsigned int ID;
typedef unsigned int TIME;
typedef real ANGLE;

//copied from SPR.cpp
//difference p1-p2 (vectors) periodic an square LxL
Vector distance_periodic_vector(Vector p1, Vector p2, real L)
{
    //from constants.cpp (only periodic ".Norm()")
    /*real absx = std::abs(p1.x-p2.x);
    real absy = std::abs(p1.y-p2.y);
    real minx = std::min(absx, L-absx);
    real miny = std::min(absy, L-absy);
    return std::sqrt(minx*minx + miny*miny);*/

    real dx = p1.x - p2.x;
    real dy = p1.y - p2.y;

    if (dx > 0.5 * L)
    {
        dx = dx - L;
    }
    else if (dx <= -0.5 * L)
    {
        dx = dx + L;
    }

    if (dy > 0.5 * L) // >= ?
    {
        dy = dy - L;
    }
    else if (dy <= -0.5 * L)
    {
        dy = dy + L;
    }

    return Vector(dx, dy);
}

real centralDifference(real a, real b, real delta)
{
    return (a - b)/delta/2.;
}

////////////////////////////////////////////////////////////////////////////////
//correlation functions as from "Collective Motion of Spherical Bacteria"
////////////////////////////////////////////////////////////////////////////////

//supply average vector of bins top (t), right (r), bottom (b), left (l) and binstep (distance between centers of bins)
real localCurl(Vector vt, Vector vr, Vector vb, Vector vl, real binstep)
{
    //(dvy/dx - dvx/dy)
    return (vl.y - vr.y - vt.x + vb.x) / binstep / 2.;
}

Vector FixBorder(Vector x, real L)
{
    real delta = 0.0001;
    if (x.x >= L)
    {
        x.x = L - delta;
    }
    else if (x.x < 0)
    {
        x.x = delta;
    }
    if (x.y >= L)
    {
        x.y = L - delta;
    }
    else if (x.y < 0)
    {
        x.y = delta;
    }
    return x;
}

//spatial velocity-velocity correlation function AND spatial orientation-orientation correlation function AND spatial direction-direction correlation function
//rmin, rmax, rstep: start, stop and width of correlation bins, tmin, tmax, tstep: start, stop und steps of time averaging
template <bool periodic>
std::vector<std::array<real, 2>> correlationPoneXVel(std::vector<std::vector<Vector>> positions, std::vector<std::vector<Vector>> velocities, real rmin, real rmax, real rstep, unsigned int tmin, unsigned int tmax, unsigned int tstep, real L)
{
    real com_dist, rleft, rright, sum, c;
    Vector v, vr, x, xr;
    int num;
    std::vector<std::array<real, 2>> ret;
    assert(positions.size() == velocities.size());
    std::vector<Vector>::iterator posit1, posit2, velit1, velit2;

    //for all bins
    for (rleft = rmin; rleft < rmax; rleft += rstep) //can get up to sqrt(2)*L?
    {
        Vector sum1(0., 0.), sum2(0., 0.);
        sum = 0.;
        num = 0;
        rright = rleft + rstep;

        std::cerr << "\rleft=" << rleft << "/" << rmax;

        //for all times
        for (TIME t = tmin; t < std::min(tmax, velocities.size()); t += tstep) //for (int t = 0; t < velocities.size(); ++t)
        {
            assert(positions[t].size() == velocities[t].size());
            //for all particles
            for (posit1 = positions[t].begin(), velit1 = velocities[t].begin();
                posit1 != positions[t].end() && velit1 != velocities[t].end();
                ++posit1, ++velit1
            )
            //for (std::pair<ID, Vector> entry1 : velocities[t])
            {
                for (posit2 = positions[t].begin(), velit2 = velocities[t].begin();
                    posit2 != positions[t].end() && velit2 != velocities[t].end();
                    ++posit2, ++velit2
                )
                //for (std::pair<ID, Vector> entry2 : velocities[t])
                {
                    v = *velit1;
                    vr = *velit2;

                    x = *posit1;
                    xr = *posit2;

                    //get distance
                    if (periodic)
                    {
                        com_dist = distance_periodic(x, xr, L);
                    }
                    else
                    {
                        com_dist = x.GetDistance(xr);
                    }

                    //test if in bin
                    if (com_dist >= rleft && com_dist < rright)
                    {
                        //calculate correlation
                        sum += v * vr;
                        sum1 += v;
                        sum2 += vr;
                        num += 1;
                    }
                }
            }
        }
        c = (num == 0 ? 0. : sum / num - (sum1 / num) * (sum2 / num));
        ret.push_back({rleft, c});
    }

    real norm = ret.front()[1];
    for (auto &cor : ret)
    {
        cor[1] /= norm;
    }

    //vector with 2-array ((rleft_1, c_1), (rleft_2, c_2), ....)
    return ret;
}


//temporal velocity-velocity correlation function AND temporal orientation-orientation correlation function AND temporal direction-direction correlation function
//tmin, tmax, tstep: time bins, vorbin: number of spatial bins
std::vector<std::array<real, 2>> correlationPoneTVel(std::vector<std::vector<Vector>> positions, std::vector<std::vector<Vector>> velocities, unsigned int tmin, unsigned int tmax, unsigned int tstep, unsigned int vorbin, real L)
{
    real c, l = L / vorbin;
    Vector v, vr, x, xr;
    unsigned int box, tminmax, tdif;
    std::vector<std::array<real, 2>> ret;

    assert(positions.size() == velocities.size());
    std::map<ID, Vector>::iterator posit1, posit2, velit1, velit2;

    //avg vel for i,j bin (vorbin x vorbin) with first index being time (from 0 to tmax-1)
    std::vector<std::vector<std::vector<Vector>>> velbins_avg(velocities.size(), std::vector<std::vector<Vector>>(vorbin, std::vector<Vector>(vorbin)));

    //for all times
    for (TIME t = 0; t < velocities.size(); ++t)
    {
        assert(velocities[t].size() == positions[t].size());
        std::vector<std::vector<Vector>> velbins(vorbin*vorbin);

        //for all particles
        for (unsigned int j = 0; j < velocities[t].size(); ++j)
        {
            v = velocities[t][j];
            x = positions[t][j];

            //order into bin
            x = FixBorder(x, L);
            assert(x.x > 0 && x.y > 0);
            assert(x.x < L && x.y < L);
            box = (static_cast<unsigned int>(x.x/l) + static_cast<unsigned int>(x.y/l)*vorbin);
            assert(box >= 0 && box < velbins.size());
            velbins[box].push_back(v);
            //std::cout << entry1.first <<" "<< x.x << " " << x.y << " " <<box << std::endl;
        }

        //average‌ bins
        for (unsigned int i = 0; i < velbins.size(); ++i)
        {
            Vector avg_vec(0, 0);

            for (unsigned int j = 0; j < velbins[i].size(); ++j)
            {
                avg_vec += velbins[i][j];
            }
            if (velbins[i].size() > 0)
            {
                avg_vec /= velbins[i].size();
            }
            assert((i-i%vorbin)/vorbin >= 0 && (i-i%vorbin)/vorbin < vorbin);
            //std::cout << t << " " << i%vorbin << " " << (i-i%vorbin)/vorbin << std::endl;
            velbins_avg[t][i%vorbin][(i-i%vorbin)/vorbin] = avg_vec;
        }
    }

    //for all time bins
    for (unsigned int dt = tmin; dt < std::min(tmax, velocities.size()); dt += tstep)
    {
        tminmax = dt + tstep;
        std::cerr << "\rtleft=" << dt << "/" << std::min(tmax, velocities.size()) << "      ";
        real sum = 0.;
        Vector sum1(0., 0.), sum2(0., 0.);
        unsigned int num = 0;
        //for all combinations of time (part 1)
        for (unsigned int t1 = 0; t1 < velocities.size(); ++t1)
        {
            //for all combinations of time (part 2)
            for (unsigned int t2 = 0; t2 < velocities.size(); ++t2)
            {
                tdif = t2 - t1;
                if (tdif >= dt && tdif < tminmax)
                {
                    //for all bins
                    for (unsigned int i = 0; i < vorbin; ++i)
                    {
                        for (unsigned int j = 0; j < vorbin; ++j)
                        {
                            //calculate correlation
                            sum += velbins_avg[t1][i][j] * velbins_avg[t2][i][j];
                            sum1 += velbins_avg[t1][i][j];
                            sum2 += velbins_avg[t2][i][j];
                            num += 1;
                        }
                    }

                }
            }
        }

        //average over space
        c = (num == 0 ? 0. : sum / num - (sum1 / num) * (sum2 / num));
        ret.push_back({dt, c});
    }

    //normalize final correlation
    real norm = ret.front()[1];
    for (auto &cor : ret)
    {
        cor[1] /= norm;
    }

    //vector with 2-array ((tleft_1, c_1), (tleft_2, c_2), ....)
    return ret;
}


//temporal vorticity-vorticity correlation function
//tmin, tmax, tstep: time bins, vorbin: number of spatial bins
std::vector<std::array<real, 2>> correlationPoneTVor(std::vector<std::vector<Vector>> positions, std::vector<std::vector<Vector>> velocities, unsigned int tmin, unsigned int tmax, unsigned int tstep, unsigned int vorbin, real L)
{
    real c, l = L / vorbin;
    Vector v, vr, x, xr;
    unsigned int box, tminmax, tdif;
    std::vector<std::array<real, 2>> ret;

    assert(positions.size() == velocities.size());
    std::map<ID, Vector>::iterator posit1, posit2, velit1, velit2;

    //avg vel for i,j bin (vorbin x vorbin) with first index being time (from 0 to tmax-1)
    //std::vector<std::vector<std::vector<Vector>>> velbins_avg(velocities.size(), std::vector<std::vector<Vector>>(vorbin, std::vector<Vector>(vorbin)));
    std::vector<std::vector<std::vector<real>>> local_curls(velocities.size(), std::vector<std::vector<real>>(vorbin-2, std::vector<real> (vorbin-2)));

    //for all times
    for (TIME t = 0; t < velocities.size(); ++t)
    {
        assert(velocities[t].size() == positions[t].size());
        std::vector<std::vector<Vector>> velbins(vorbin*vorbin);

        //for all particles
        for (unsigned int j = 0; j < velocities[t].size(); ++j)
        {
            v = velocities[t][j];
            x = positions[t][j];

            //order into bin
            x = FixBorder(x, L);
            assert(x.x > 0 && x.y > 0);
            assert(x.x < L && x.y < L);
            box = (static_cast<unsigned int>(x.x/l) + static_cast<unsigned int>(x.y/l)*vorbin);
            assert(box >= 0 && box < velbins.size());
            velbins[box].push_back(v);
            //std::cout << entry1.first <<" "<< x.x << " " << x.y << " " <<box << std::endl;
        }

        //average‌ bins
        std::vector<std::vector<Vector>> velbins_avg(vorbin, std::vector<Vector>(vorbin));
        for (unsigned int i = 0; i < velbins.size(); ++i)
        {
            Vector avg_vec(0, 0);

            for (unsigned int j = 0; j < velbins[i].size(); ++j)
            {
                avg_vec += velbins[i][j];
            }
            if (velbins[i].size() > 0)
            {
                avg_vec /= velbins[i].size();
            }
            assert((i-i%vorbin)/vorbin >= 0 && (i-i%vorbin)/vorbin < vorbin);
            //std::cout << t << " " << i%vorbin << " " << (i-i%vorbin)/vorbin << std::endl;
            velbins_avg[i%vorbin][(i-i%vorbin)/vorbin] = avg_vec;
        }

        //could be done with periodic boundaries...
        //calculate local curl for all bins (except outer layer)
        for (unsigned int i = 1; i < velbins_avg.size() - 1; ++i)
        {
            for (unsigned int j = 1; j < velbins_avg[i].size() -1; ++j)
            {
                local_curls[t][i-1][j-1] = localCurl(velbins_avg[i][j+1], velbins_avg[i+1][j], velbins_avg[i][j-1], velbins_avg[i-1][j], l);
            }
        }
    }

    //for all time bins
    for (unsigned int dt = tmin; dt < std::min(tmax, velocities.size()); dt += tstep)
    {
        tminmax = dt + tstep;
        std::cerr << "\rtleft=" << dt << "/" << std::min(tmax, velocities.size()) << "      ";
        real sum = 0., sum1 = 0., sum2 = 0.;
        unsigned int num = 0;

        //for all combinations of time (part 1)
        for (unsigned int t1 = 0; t1 < velocities.size(); ++t1)
        {
            //for all combinations of time (part 2)
            for (unsigned int t2 = 0; t2 < velocities.size(); ++t2)
            {
                tdif = t2 - t1;
                if (tdif >= dt && tdif < tminmax)
                {
                    //for all bins
                    for (unsigned int i = 0; i < local_curls[t1].size(); ++i)
                    {
                        for (unsigned int j = 0; j < local_curls[t1][i].size(); ++j)
                        {
                            //calculate correlation
                            sum += local_curls[t1][i][j] * local_curls[t2][i][j];
                            sum1 += local_curls[t1][i][j];
                            sum2 += local_curls[t2][i][j];
                            num += 1;
                        }
                    }

                }
            }
        }

        //average over space
        c = (num == 0 ? 0. : sum / num - (sum1 / num) * (sum2 / num));
        ret.push_back({dt, c});
    }

    //normalize final correlation
    real norm = ret.front()[1];
    for (auto &cor : ret)
    {
        cor[1] /= norm;
    }

    //vector with 2-array ((tleft_1, c_1), (tleft_2, c_2), ....)
    return ret;
}



//spatial vorticity-vorticity correlation function
//rmin, rmax, rstep: start, stop and width of histogram bins, tmin, tmax, tstep: start, stop und steps of time averaging, vorbin: number of spatial bins for vorticy calculation
std::vector<std::array<real, 3>> correlationPoneXVor(std::vector<std::vector<Vector>> positions, std::vector<std::vector<Vector>> velocities, real rmin, real rmax, real rstep, unsigned int tmin, unsigned int tmax, unsigned int tstep, unsigned int vorbin, real L)
{
    real rleft, rright, sum, c, l = L / vorbin;
    Vector v, vr, x, xr, avg_vec;
    int num;
    unsigned int ret2_ind = 0, box, avg_sum;
    std::vector<std::array<real, 3>> ret;
    std::vector<std::vector<real>> ret2;

    for (rleft = rmin; rleft < rmax; rleft += rstep)
    {
        ret2.push_back(std::vector<real> ());
    }
    assert(positions.size() == velocities.size());
    std::map<ID, Vector>::iterator posit1, posit2, velit1, velit2;

    //for all times
    for (unsigned int t = tmin; t < std::min(tmax, velocities.size()); t += tstep) //for (int t = 0; t < velocities.size(); ++t)
    {
        std::vector<std::vector<Vector>> velbins(vorbin*vorbin);
        //for all particles
        for (unsigned int j = 0; j < velocities[t].size(); ++j)
        {
            v = velocities[t][j];
            x = positions[t][j];

            //order into bin
            x = FixBorder(x, L);
            assert(x.x < L && x.y < L);
            box = (static_cast<unsigned int>(x.x/l) + static_cast<unsigned int>(x.y/l)*vorbin);
            velbins[box].push_back(v);
            //std::cout << entry1.first <<" "<< x.x << " " << x.y << " " <<box << std::endl;
        }


        //average‌ bins
        std::vector<std::vector<Vector>> velbins_avg (vorbin, std::vector<Vector>(vorbin));
        for (unsigned int i = 0; i < velbins.size(); ++i)
        {
            avg_vec.x = 0.;
            avg_vec.y = 0.;
            avg_sum = 0;
            for (unsigned int j = 0; j < velbins[i].size(); ++j)
            {
                avg_vec += velbins[i][j];
                avg_sum += 1;
            }
            if (avg_sum > 0)
            {
                avg_vec /= avg_sum;
            }
            velbins_avg[i%vorbin][(i-i%vorbin)/vorbin] = avg_vec;
        }

        //could be done with periodic boundaries...
        //calculate local curl for all bins (except outer layer)
        std::vector<std::vector<real>> local_curls(vorbin-2, std::vector<real> (vorbin-2));
        for (unsigned int i = 1; i < velbins_avg.size() - 1; ++i)
        {
            for (unsigned int j = 1; j < velbins_avg[i].size() -1; ++j)
            {
                local_curls[i-1][j-1] = localCurl(velbins_avg[i][j+1], velbins_avg[i+1][j], velbins_avg[i][j-1], velbins_avg[i-1][j], l);
            }
        }


        //for all bins
        ret2_ind = 0;
        for (rleft = rmin; rleft < rmax; rleft += rstep) //can get up to sqrt(2)*L?
        {
            //Vector sum1(0., 0.), sum2(0., 0.);
            real sum1 = 0., sum2 = 0.;
            sum = 0.;
            num = 0;
            rright = rleft + rstep;

            std::cerr << "\rtleft=" << t << "/" << std::min(tmax, velocities.size()) << " " << "rleft=" << rleft << "/" << rmax << "      ";


            //for all combinations of bins
            for (unsigned int i1 = 0; i1 < local_curls.size(); ++i1)
            {
                for (unsigned int j1 = 0; j1 < local_curls[i1].size(); ++j1)
                {
                    for (unsigned int i2 = 0; i2 < local_curls.size(); ++i2)
                    {
                        for (unsigned int j2 = 0; j2 < local_curls[i2].size(); ++j2)
                        {
                            real bin_dist = std::sqrt((i2-i1)*(i2-i1)+(j2-j1)*(j2-j1)) * l;
                            if (bin_dist >= rleft && bin_dist < rright)
                            {
                                //calculate correlation
                                sum += local_curls[i1][j1] * local_curls[i2][j2];
                                sum1 += local_curls[i1][j1];
                                //std::cout << local_curls[i1][j1] << std::endl;
                                sum2 += local_curls[i2][j2];
                                num += 1;
                            }
                        }
                    }
                }
            }


            c = (num == 0 ? 0. : sum / num - (sum1 / num) * (sum2 / num));
            //ret.push_back({rleft, c});
            //std::cout << ret2_ind << " " << c << std::endl;
            ret2[ret2_ind].push_back(c);
            ret2_ind += 1;
        }
    }


    //average over time
    ret2_ind = 0;
    for (rleft = rmin; rleft < rmax; rleft += rstep)
    {
        unsigned int mean_counter = 0;
        real mean_c = 0., sd_c = 0.;
        for (unsigned int i = 0; i < ret2[ret2_ind].size(); ++i)
        {
            mean_c += ret2[ret2_ind][i];
            mean_counter += 1;
            //std::cout << "rleft: " << rleft << " " << " i: " << i << " " << " c: " << ret2[ret2_ind][i] << std::endl;
        }
        if (mean_counter > 0)
        {
            mean_c /= mean_counter;
        }
        mean_counter = 0;
        for (unsigned int i = 0; i < ret2[ret2_ind].size(); ++i)
        {
            sd_c += (ret2[ret2_ind][i] - mean_c) * (ret2[ret2_ind][i] - mean_c);
            mean_counter += 1;
        }
        if (mean_counter > 1)
        {
            sd_c = std::sqrt(1./(mean_counter - 1) * sd_c);
        }
        ret.push_back({rleft, mean_c, sd_c});
        ret2_ind += 1;
    }


    //normalize final correlation
    real norm = ret.front()[1];
    for (auto &cor : ret)
    {
        cor[1] /= norm;
    }

    std::cout << std::endl;

    //vector with 3-array ((rleft_1, c_1, sd_1), (rleft_2, c_2, sd_1), ....)
    return ret;
}


//temporal velocity-velocity auto correlation function AND temporal orientation-orientation auto correlation function
//tmin, tmax, tstep: time bins
std::vector<std::array<real, 2>> correlationPoneTAutoVel(std::vector<std::vector<Vector>> velocities, std::vector<std::vector<ID>> ids, unsigned int tmin, unsigned int tmax, unsigned int tstep, real L)
{
    std::map<ID, std::vector<Vector>> per_id_velocities;

    std::cerr << "Creating trajectories map. vels: " << velocities.size() << ", ids: " << ids.size() << std::endl;

    for (TIME t = 0; t < velocities.size(); ++t)
    {
        for (unsigned int j = 0; j < velocities[t].size(); ++j)
        {
            //std::cout << t << " " << j << std::endl;
            //std::cout << ids[t][j] << std::endl;
            //std::cout << velocities[t][j].x << std::endl;
            per_id_velocities[ids[t][j]].push_back(velocities[t][j]);
        }
    }

    std::cerr << "Finding longest trajectory:";

    unsigned int maxlen = 0;
    for (const auto &vellist : per_id_velocities)
    {
        maxlen = std::max(maxlen, vellist.second.size());
    }

    std::cerr << " " << maxlen << std::endl;

    std::vector<real> correlations_ttr(maxlen, 0.); // first index is deltat, then list with correlations
    std::vector<Vector> correlations_tt(maxlen, Vector(0., 0.)); // first index is deltat, then list with correlations
    std::vector<Vector> correlations_trtr(maxlen, Vector(0., 0.)); // first index is deltat, then list with correlations
    std::vector<unsigned int> num (maxlen, 0);

    std::cerr << "Calculation time correlations" << std::endl;

    for (const auto &entry : per_id_velocities)
    {
        const auto &vellist = entry.second;
        for (TIME deltat = tmin; deltat < std::min(vellist.size(), tmax); deltat += tstep)
        {
            //TODO bin times
            for (TIME start = 0; start < vellist.size() - deltat; ++start)
            {
                assert(maxlen >= vellist.size());
                assert(start >= 0 && start <= vellist.size());
                assert(start+deltat >= 0 && start+deltat < vellist.size());
                assert(deltat >= 0 && deltat < maxlen);
                correlations_ttr[deltat] += vellist[start]*vellist[start+deltat];
                correlations_tt[deltat] += vellist[start];
                correlations_trtr[deltat] += vellist[start+deltat];
                num[deltat] += 1;
            }
        }
    }

    std::cerr << "Calculating mean and sds" << std::endl;

    std::vector<std::array<real, 2>> ret(maxlen);
    for (TIME t = 0; t < maxlen; ++t)
    {
        ret[t] = {t, (num[t] == 0 ? 0. : (correlations_ttr[t]/num[t]) - (correlations_tt[t]/num[t])*(correlations_trtr[t]/num[t]))};
    }

    real norm = ret.front()[1];
    for (auto &cor : ret)
    {
        cor[1] /= norm;
    }

    return ret;
}


//temporal vorticity auto correlation function
std::vector<std::array<real, 2>> correlationPoneTAutoVor(std::vector<std::vector<Vector>> positions, std::vector<std::vector<Vector>> velocities, std::vector<std::vector<ID>> ids, unsigned int tmin, unsigned int tmax, unsigned int tstep, unsigned int vorbin, real L)
{
    unsigned int tt = velocities.size();
    unsigned int box;
    real l = L / vorbin;

    Vector v, x;

    assert(positions.size() == velocities.size());

    std::vector<std::vector<std::vector<real>>> local_curls(tt, std::vector<std::vector<real>>(vorbin-2, std::vector<real>(vorbin-2)));
    //for all particles
    for (TIME t = 0; t < tt; ++t)
    {
        //fix
        std::vector<std::vector<Vector>> velbins(vorbin*vorbin);
        std::vector<std::vector<Vector>> velbins_avg(vorbin, std::vector<Vector>(vorbin));

        for (unsigned int j = 0; j < velocities[t].size(); ++j)
        {
            assert(velocities[t].size() == positions[t].size());
            v = velocities[t][j];
            x = positions[t][j];

            //order into bin
            x = FixBorder(x, L);
            assert(x.x < L && x.y < L);
            box = (static_cast<unsigned int>(x.x/l) + static_cast<unsigned int>(x.y/l)*vorbin);
            velbins[box].push_back(v);
            //std::cout << entry1.first <<" "<< x.x << " " << x.y << " " <<box << std::endl;
        }

        //average‌ bins

        for (unsigned int i = 0; i < velbins.size(); ++i)
        {
            Vector avg_vec(0, 0);
            for (unsigned int j = 0; j < velbins[i].size(); ++j)
            {
                avg_vec += velbins[i][j];
            }
            if (velbins[i].size() > 0)
            {
                avg_vec /= velbins[i].size();
            }
            velbins_avg[i%vorbin][(i-i%vorbin)/vorbin] = avg_vec;
        }

        //could be done with periodic boundaries...
        //calculate local curl for all bins (except outer layer)

        for (unsigned int i = 1; i < velbins_avg.size() - 1; ++i)
        {
            for (unsigned int j = 1; j < velbins_avg[i].size() -1; ++j)
            {
                local_curls[t][i-1][j-1] = localCurl(velbins_avg[i][j+1], velbins_avg[i+1][j], velbins_avg[i][j-1], velbins_avg[i-1][j], l);
            }
        }
    }

    std::map<ID, std::vector<real>> per_id_vorticities;

    std::cerr << "Creating trajectories map. vels: " << velocities.size() << ", ids: " << ids.size() << std::endl;

    unsigned int xbin, ybin;

    for (TIME t = 0; t < positions.size(); ++t)
    {
        for (unsigned int j = 0; j < positions[t].size(); ++j)
        {
            x = FixBorder(positions[t][j], L);
            xbin = static_cast<unsigned int>(x.x/l);
            ybin = static_cast<unsigned int>(x.y/l);

            //ignore particles in outer parts
            if (xbin == 0 || xbin == vorbin - 1)
            {
                continue;
            }
            if (ybin == 0 || ybin == vorbin - 1)
            {
                continue;
            }

            assert(xbin >= 1 && xbin < vorbin-1);
            assert(ybin >= 1 && ybin < vorbin-1);
            // in which bin is this particle? find corresponding local curl
            per_id_vorticities[ids[t][j]].push_back(local_curls[t][xbin-1][ybin-1]);
        }
    }

    std::cerr << "Finding longest trajectory:";

    unsigned int maxlen = 0;
    for (const auto &vellist : per_id_vorticities)
    {
        maxlen = std::max(maxlen, vellist.second.size());
    }

    std::cerr << " " << maxlen << std::endl;

    std::vector<real> correlations_ttr(maxlen, 0); // first index is deltat, then list with correlations
    std::vector<real> correlations_tt(maxlen, 0); // first index is deltat, then list with correlations
    std::vector<real> correlations_trtr(maxlen, 0); // first index is deltat, then list with correlations
    std::vector<unsigned int> num (maxlen, 0);

    std::cerr << "Calculation time correlations" << std::endl;

    for (const auto &entry : per_id_vorticities)
    {
        const auto &vellist = entry.second;
        for (TIME deltat = tmin; deltat < std::min(vellist.size(), tmax); deltat += tstep)
        {
            //TODO bin times
            for (TIME start = 0; start < vellist.size() - deltat; ++start)
            {
                assert(maxlen >= vellist.size());
                assert(start >= 0 && start <= vellist.size());
                assert(start+deltat >= 0 && start+deltat < vellist.size());
                assert(deltat >= 0 && deltat < maxlen);
                correlations_ttr[deltat] += vellist[start]*vellist[start+deltat];
                correlations_tt[deltat] += vellist[start];
                correlations_trtr[deltat] += vellist[start+deltat];
                num[deltat] += 1;
            }
        }
    }

    std::cerr << "Calculating mean and sds" << std::endl;

    std::vector<std::array<real, 2>> ret(maxlen);
    for (TIME t = 0; t < maxlen; ++t)
    {
        ret[t] = {t, (num[t] == 0 ? 0. : (correlations_ttr[t]/num[t]) - (correlations_tt[t]/num[t])*(correlations_trtr[t]/num[t]))};
    }

    real norm = ret.front()[1];
    for (auto &cor : ret)
    {
        cor[1] /= norm;
    }

    return ret;
}



//vorticity histogram
//rmin, rmax, rstep: start, stop and width of correlation bins, tmin, tmax, tstep: start, stop und steps of time averaging, vorbin: number of boxes (bins)
std::vector<std::array<real, 3>> histogramXVor(std::vector<std::vector<Vector>> positions, std::vector<std::vector<Vector>> velocities, real rmin, real rmax, real rstep, unsigned int tmin, unsigned int tmax, unsigned int tstep, unsigned int vorbin, real L)
{
    real rleft, l = L / vorbin;
    Vector v, vr, x, xr, avg_vec;
    unsigned int ret2_ind = 0, avg_sum, box;
    std::vector<std::array<real, 3>> ret;
    std::vector<unsigned int> ret2;

    for (rleft = rmin; rleft < rmax; rleft += rstep)
    {
        ret2.push_back(static_cast<unsigned int> (0));
    }
    assert(positions.size() == velocities.size());
    std::map<ID, Vector>::iterator posit1, posit2, velit1, velit2;

    //for all times
    for (unsigned int t = tmin; t < std::min(tmax, velocities.size()); t += tstep) //for (int t = 0; t < velocities.size(); ++t)
    {
        std::vector<std::vector<Vector>> velbins(vorbin*vorbin);
        //for all particles
        for (unsigned int j = 0; j < velocities[t].size(); ++j)
        {
            v = velocities[t][j];
            x = positions[t][j];

            //order into bin
            x = FixBorder(x, L);
            assert(x.x < L && x.y < L);
            box = (static_cast<unsigned int>(x.x/l) + static_cast<unsigned int>(x.y/l)*vorbin);
            velbins[box].push_back(v);
            //std::cout << entry1.first <<" "<< x.x << " " << x.y << " " <<box << std::endl;
        }


        //average‌ bins
        std::vector<std::vector<Vector>> velbins_avg (vorbin, std::vector<Vector>(vorbin));
        for (unsigned int i = 0; i < velbins.size(); ++i)
        {
            avg_vec.x = 0.;
            avg_vec.y = 0.;
            avg_sum = 0;
            for (unsigned int j = 0; j < velbins[i].size(); ++j)
            {
                avg_vec += velbins[i][j];
                avg_sum += 1;
            }
            if (avg_sum > 0)
            {
                avg_vec /= avg_sum;
            }
            velbins_avg[i%vorbin][(i-i%vorbin)/vorbin] = avg_vec;
        }


        //calculate local curl for all bins (except outer layer)
        std::vector<std::vector<real>> local_curls(vorbin-2, std::vector<real> (vorbin-2));
        for (unsigned int i = 1; i < velbins_avg.size() - 1; ++i)
        {
            for (unsigned int j = 1; j < velbins_avg[i].size() -1; ++j)
            {
                local_curls[i-1][j-1] = localCurl(velbins_avg[i][j+1], velbins_avg[i+1][j], velbins_avg[i][j-1], velbins_avg[i-1][j], l);
                box = static_cast<int> ((local_curls[i-1][j-1] - std::fmod(local_curls[i-1][j-1], rstep) - rmin) / rstep);

                //std::cout << local_curls[i-1][j-1] << " " <<  std::fmod(local_curls[i-1][j-1], rstep) << " " << box << std::endl;

                if (box > 0 && box < ret2.size())
                {
                    ret2[box] += 1;
                }
                else
                {
                    //out of all bins
                    std::cout << "warning: vorticity histogram rbins out of range " << rmin << "-" << rmax << ": " <<local_curls[i-1][j-1] << std::endl;
                }
            }
        }
    }

    //create ret
    ret2_ind = 0;
    for (rleft = rmin; rleft < rmax; rleft += rstep)
    {
        ret.push_back({rleft, ret2[ret2_ind], 0.});
        ret2_ind += 1;
    }


    //vector with 3-array ((rleft_1, c_1, sd_1), (rleft_2, c_2, sd_1), ....)
    return ret;
}


//velocity histogram
//rmin, rmax, rstep: start, stop and width of correlation bins, tmin, tmax, tstep: start, stop und steps of time averaging, vorbin: number of boxes (bins)
std::vector<std::array<real, 3>> histogramXVel(std::vector<std::vector<Vector>> velocities, real rmin, real rmax, real rstep, unsigned int tmin, unsigned int tmax, unsigned int tstep, real L)
{
    real rleft, vnorm;
    Vector v, vr, x, xr, avg_vec;
    unsigned int ret2_ind, box;
    std::vector<std::array<real, 3>> ret;
    std::vector<unsigned int> ret2;

    for (rleft = rmin; rleft < rmax; rleft += rstep)
    {
        ret2.push_back(static_cast<unsigned int> (0));
    }

    //for all times
    for (unsigned int t = tmin; t < std::min(tmax, velocities.size()); t += tstep)
    {
        //for all particles
        for (unsigned int j = 0; j < velocities[t].size(); ++j)
        {
            vnorm = velocities[t][j].Norm();
            box = static_cast<int> ((vnorm - std::fmod(vnorm, rstep) - 2* rmin) / rstep); //((vnorm - std::fmod(vnorm, rstep) - rmin) / rstep);
            if (box > 0 && box < ret2.size())
            {
                ret2[box] += 1;
            }
            else
            {
                //out of all bins
                std::cout << "warning: velocity histogram rbins out of range " << rmin << "-" << rmax << ": " << vnorm << std::endl;
            }
        }
    }

    //create ret
    ret2_ind = 0;
    for (rleft = rmin; rleft < rmax; rleft += rstep)
    {
        ret.push_back({rleft, ret2[ret2_ind], 0.});
        ret2_ind += 1;
    }


    //vector with 3-array ((rleft_1, c_1, sd_1), (rleft_2, c_2, sd_1), ....)
    return ret;
}



//per id velocity histogram
//rmin, rmax, rstep: start, stop and width of correlation bins, tmin, tmax, tstep: start, stop und steps of time averaging
std::map<ID, std::vector<std::array<real, 3>>> perIDhistogramVel(std::vector<std::vector<Vector>> velocities, std::vector<std::vector<ID>> ids, real rmin, real rmax, real rstep, unsigned int tmin, unsigned int tmax, unsigned int tstep, real L)
{
    real rleft;
    unsigned int ret2_ind, box;
    std::map<ID, std::vector<std::array<real, 3>>> ret;
    std::map<ID, std::vector<unsigned int>> ret2;
    std::map<ID, std::vector<real>> per_id_velocities;

    for (TIME t = tmin; t < std::min(tmax, velocities.size()); t += tstep)
    {
        for (unsigned int j = 0; j < velocities[t].size(); ++j)
        {
            //std::cout << t << " " << j << std::endl;
            //std::cout << ids[t][j] << std::endl;
            //std::cout << velocities[t][j].x << std::endl;
            per_id_velocities[ids[t][j]].push_back(velocities[t][j].Norm());
        }
    }

    //create ret2
    for(auto entry : per_id_velocities)
    {
        for (rleft = rmin; rleft < rmax; rleft += rstep)
        {
            ret2[entry.first].push_back(0);
        }
    }

    //for all ids
    for (std::pair<ID, std::vector<real>> entry : per_id_velocities)
    {
        //for all particles of a id
        for (real vnorm : entry.second)
        {
            box = static_cast<int> ((vnorm - std::fmod(vnorm, rstep)- 2*rmin) / rstep); //((vnorm - std::fmod(vnorm, rstep) - rmin) / rstep);
            if (box > 0 && box < ret2[entry.first].size())
            {
                ret2[entry.first][box] += 1;
            }
            else
            {
                //out of all bins
                std::cout << "warning: per id velocity histogram rbins out of range " << rmin << "-" << rmax << ": " << vnorm << std::endl;
            }
        }
    }

    //create ret
    for (auto entry : ret2)
    {
        ret2_ind = 0;
        for (rleft = rmin; rleft < rmax; rleft += rstep)
        {
            ret[entry.first].push_back({rleft, ret2[entry.first][ret2_ind], 0.});
            ret2_ind += 1;
        }
    }


    //map (1st dim: ID) with vectors (histogram) with 3-array ((rleft_1, c_1, sd_1), (rleft_2, c_2, sd_1), ....)
    return ret;
}


//per id traveled pathlength
//rmin, rmax, rstep: start, stop and width of correlation bins, tmin, tmax, tstep: start, stop und steps of time averaging
std::map<ID, std::array<real, 3>> perIDtraveledPathlength(std::vector<std::vector<Vector>> velocities, std::vector<std::vector<ID>> ids, unsigned int tmin, unsigned int tmax, unsigned int tstep)
{
    std::map<ID, std::array<real, 3>> ret;
    std::map<ID, std::vector<real>> per_id_velocities;

    for (TIME t = tmin; t < std::min(tmax, velocities.size()); t += tstep)
    {
        for (unsigned int j = 0; j < velocities[t].size(); ++j)
        {
            //std::cout << t << " " << j << std::endl;
            //std::cout << ids[t][j] << std::endl;
            //std::cout << velocities[t][j].x << std::endl;
            per_id_velocities[ids[t][j]].push_back(velocities[t][j].Norm());
        }
    }

    //for all ids
    for (std::pair<ID, std::vector<real>> entry : per_id_velocities)
    {
        ret[entry.first] = {0., 0., 0.};

        //for all particles of a id
        for (real vnorm : entry.second)
        {
            ret[entry.first][0] += tstep;
            ret[entry.first][1] += tstep * vnorm;
        }
    }

    //map IDs to 3-array ((time1, path1, 0.), (time2, path2, 0.), ....),
    //where time1 is the time of the trajectory (in tstep) and path is the
    //path length of the particle (sum over velocity * tstep)
    return ret;
}


//number fluctuations (loosely based on number_fluctuations.py)
std::vector<std::array<real, 3>> number_fluctuations(std::vector<std::vector<Vector>> positions, real L)
{
    std::vector<std::array<real, 3>> ret;
    Vector pos;
    unsigned int no_of_particles, max_no_of_particles, no_of_boxes, max_no_of_boxes, box;
    real mean_abs_dist_to_average, abs_dist_to_average, sd_abs_dist_to_average, l, avg_no_per_box;

    max_no_of_particles = positions[0].size();
    //find maximum number of particles
    for (const auto &position : positions)
    {
        no_of_particles = position.size();
        if (no_of_particles > max_no_of_particles)
        {
            max_no_of_particles = no_of_particles;
        }
    }
    max_no_of_boxes = std::sqrt(max_no_of_particles);
    //map no_of_boxes -> list of number fluctuation (abs distance to average)
    std::map<unsigned int, std::vector<real>> number_flucts;
    for (no_of_boxes = 2; no_of_boxes <= max_no_of_boxes; no_of_boxes += 1)
    {
        number_flucts[no_of_boxes] = std::vector<real>();
    }

    //for all times
    for (const auto &position : positions)
    {
        no_of_particles = position.size();
        //for all number of boxes
        for (no_of_boxes = 2; no_of_boxes <= max_no_of_boxes; no_of_boxes += 1)
        {
            assert(no_of_boxes > 0);
            avg_no_per_box = real(no_of_particles) / (no_of_boxes * no_of_boxes);
            l = L / no_of_boxes;
            std::vector<real> box_vec(no_of_boxes*no_of_boxes, 0); //save no of particles of all boxes
            //for all particles
            for (Vector pos : position)
            {
                pos = FixBorder(pos, L);
                assert(pos.x < L && pos.y < L);
                box = static_cast<unsigned int>(pos.x/l) + static_cast<unsigned int>(pos.y/l) *no_of_boxes;
                assert(box >= 0 && box < no_of_boxes*no_of_boxes);
                box_vec[box] += 1; //put particle in a box
            }
            //go over all boxes and calculate average absolute distance to average
            abs_dist_to_average = 0.;
            for (auto particles : box_vec)
            {
                abs_dist_to_average += std::abs(particles - avg_no_per_box);
            }
            //std::cout << "box vec size: " << box_vec.size() << std::endl;
            abs_dist_to_average /= box_vec.size();
            number_flucts[no_of_boxes].push_back(abs_dist_to_average);
            //std::cout << "len number_flucts " << number_flucts.size() << " abs_dist " << abs_dist_to_average << std::endl;
        }
    }
    //std::vector<real> temp_vec;
    //change keys in dictionary, such that they are actually average particle number
    /*for (no_of_boxes = 2; no_of_boxes <= max_no_of_boxes; no_of_boxes += 1)
    {
        temp_vec = number_flucts[no_of_boxes];
        number_flucts.erase(no_of_boxes);
        number_flucts[static_cast<unsigned int> (std::round(real(max_no_of_particles)/(no_of_boxes*no_of_boxes)))] = temp_vec;
    }*/

    //vector with 3-array ((<N>_L, Delta N, sd Delta N), ....)
    for (auto element = number_flucts.rbegin(); element != number_flucts.rend(); ++element)
    {
        mean_abs_dist_to_average = 0.;
        for (auto vec : element->second)
        {
            mean_abs_dist_to_average += vec;
        }
        mean_abs_dist_to_average /= element->second.size();

        sd_abs_dist_to_average = 0.;
        for (const auto vec : element->second)
        {
            sd_abs_dist_to_average += (vec - mean_abs_dist_to_average) * (vec - mean_abs_dist_to_average);
        }
        sd_abs_dist_to_average /= (element->second.size()-1);
        sd_abs_dist_to_average = std::sqrt(sd_abs_dist_to_average);

        ret.push_back({real(max_no_of_particles) / (element->first * element->first), mean_abs_dist_to_average, sd_abs_dist_to_average});
    }

    return ret;
}

//COPIED from tracking.cpp
real AngleDifference180(real a, real b)
{
    return std::min(mod(a - b, 2*PI), mod(b - a, 2*PI));
}

//copied from clustering.cpp
bool cluster_condition(Vector pos1, Vector pos2, real ang1, real ang2, real max_com_distance, real max_ang_distance, real L, bool periodic)
{
    //center of mass distance
    real com_distance;
    //PERIODIC!????????
    //com_distance = distance_periodic(positions[i], positions[j], L);
    //NON-PERIODIC??????? -> more yes than no (periodic for simulation, nonperiodic for exp?)
    if (!periodic)
    {
        //real xdiff = pos1.x-pos2.x;
        //real ydiff = pos1.y-pos2.y;
        //com_distance = std::sqrt(xdiff*xdiff + ydiff*ydiff);
        com_distance = pos1.GetDistance(pos2);
    }
    else
    {
        com_distance = distance_periodic(pos1, pos2, L);
    }

    //angle difference CORRECT????
    //real ang_distance = std::abs(ang1 - ang2);
    //ang_distance = std::min(ang_distance, 2*PI - ang_distance);
    real ang_distance = AngleDifference180(ang1, ang2); //not very very well tested....

    return (com_distance <= max_com_distance) && (ang_distance <= max_ang_distance);
}

//copied from clustering.cpp
//returns vector over time with map "size of cluster" -> "no of clusters of that size"
std::vector<std::map<unsigned int, unsigned int>> clusterSizes(std::vector<std::vector<Vector>> positions, std::vector<std::vector<real>> angles, TIME tstart, TIME tstop, TIME tstep, real max_com_distance, real max_ang_distance, real L, bool periodic_boundary_cond)
{
    std::vector<std::map<unsigned int, unsigned int>> ret;
    int delete_index, cluster_index;
    bool cluster_cond;
    for (TIME t = tstart; t < std::min(tstop, positions.size()); t += tstep)
    {
        cluster_index = 0;

        //initialize with -1s
        std::vector<int> cluster(positions[t].size(), -1);

        //later: iterate over all particles in 9 boxes
        for (unsigned int p = 0; p < positions[t].size(); ++p)
        {
            //particle not in cluster yet
            if (cluster[p] == -1)
            {
                cluster[p] = cluster_index;
                cluster_index++;
            }

            for (unsigned int p2 = 0; p2 < positions[t].size(); ++p2)
            {
                if (p == p2)
                {
                    continue;
                }
                cluster_cond = cluster_condition(positions[t][p], positions[t][p2], angles[t][p], angles[t][p2], max_com_distance, max_ang_distance, L, periodic_boundary_cond);
                //p and p2 are cluster connected
                if (cluster_cond)
                {

                    //p2 already in cluster
                    if (cluster[p2] != -1)
                    {
                        //remember id of p2
                        delete_index = cluster[p2];
                        //put p2 in cluster of p
                        cluster[p2] = cluster[p];

                        //join clusters (now the id of p2 is free?????!!!!!!)
                        //naiv: faster method??? like for example bookkeeping of all clusters
                        for (unsigned int c = 0; c < positions[t].size(); ++c)
                        {
                            if (cluster[c] == delete_index)
                            {
                                cluster[c] = cluster[p];
                            }
                        }
                    }
                    //p2 not in cluster yet
                    else
                    {
                        cluster[p2] = cluster[p];
                    }
                }
            }
        }

        //find maximum index
        int max_index = cluster.front();
        for(auto entry : cluster)
        {
            if (entry > max_index)
            {
                max_index = entry;
            }
        }


        //map index to size of cluster
        std::vector<unsigned int> count_cluster(max_index + 1, 0);
        //for all particles in cluster
        for (auto entry : cluster)
        {
            assert(entry >= 0 && entry < max_index + 1);
            count_cluster[entry] += 1;
        }


        //go over all non-0 entries in count_cluster and add 1
        //to corresponding ret-entry
        ret.push_back(std::map<unsigned int, unsigned int>());
        for (auto entry : count_cluster)
        {
            if (entry > 0)
            {
                assert(ret.back()[entry] >= 0);
                ret.back()[entry] += 1;
            }
        }
    }

    return ret;
}


//enstrophy per unit area
real Enstrophy(std::vector<std::vector<Vector>> positions, std::vector<std::vector<Vector>> velocities, unsigned int tmin, unsigned int tmax, unsigned int tstep, unsigned int vorbin, real L)
{
    real l = L / vorbin;
    Vector v, vr, x, xr, avg_vec;
    unsigned int avg_sum, box;
    real ret = 0;
    real ret_times = 0.;
    real loccurl;


    assert(positions.size() == velocities.size());

    //for all times
    for (unsigned int t = tmin; t < std::min(tmax, velocities.size()); t += tstep) //for (int t = 0; t < velocities.size(); ++t)
    {
        std::vector<std::vector<Vector>> velbins(vorbin*vorbin);
        //for all particles
        for (unsigned int j = 0; j < velocities[t].size(); ++j)
        {
            v = velocities[t][j];
            x = positions[t][j];

            //order into bin
            x = FixBorder(x, L);
            assert(x.x < L && x.y < L);
            box = (static_cast<unsigned int>(x.x/l) + static_cast<unsigned int>(x.y/l)*vorbin);
            velbins[box].push_back(v);
            //std::cout << entry1.first <<" "<< x.x << " " << x.y << " " <<box << std::endl;
        }


        //average‌ bins
        std::vector<std::vector<Vector>> velbins_avg (vorbin, std::vector<Vector>(vorbin));
        for (unsigned int i = 0; i < velbins.size(); ++i)
        {
            avg_vec.x = 0.;
            avg_vec.y = 0.;
            avg_sum = 0;
            for (unsigned int j = 0; j < velbins[i].size(); ++j)
            {
                avg_vec += velbins[i][j];
                avg_sum += 1;
            }
            if (avg_sum > 0)
            {
                avg_vec /= avg_sum;
            }
            velbins_avg[i%vorbin][(i-i%vorbin)/vorbin] = avg_vec;
        }


        //calculate local curl for all bins (except outer layer) and add to enstrophy
        for (unsigned int i = 1; i < velbins_avg.size() - 1; ++i)
        {
            for (unsigned int j = 1; j < velbins_avg[i].size() -1; ++j)
            {
                loccurl = localCurl(velbins_avg[i][j+1], velbins_avg[i+1][j], velbins_avg[i][j-1], velbins_avg[i-1][j], l);
                ret += 1./2. * loccurl * loccurl;
                ret_times += 1.;
            }
        }
    }


    //real enstrophy (averaged over time)
    return ret/ret_times/(l*l);
}



void VelocityHistogramPerParticleTrajecory()
{
    // maybe we can find histogram with two peaks for "dead" particles
}

void VelocityHistogramAveraged()
{
    // compare meso-scale appendix S4.A
}

void NormalizedEqualTimeVelocityCorrelation()
{
    // compare meso-scale appendix S4.B
}

//copied from read_data_single from constants.cpp/hpp
// times = {1,2,3} reads line 1,2,3
std::string read_ids(std::string filename, std::vector<std::vector<ID>> &ids, std::vector<unsigned int> &times)
{
    std::ifstream in;
    in.open(filename);
    assert(!in.fail() && "Could not open file");

    unsigned int id;
    unsigned int linenum = 1;
    unsigned int time = 0;

    std::string firstline, line;
    std::vector<ID> idsv; //()

    std::getline(in, firstline); //skip first line
    while (in.good())
    {
        //std::cout << time << " " << linenum << std::endl;
        std::getline(in, line);
        bool readline = linenum == times[time];
        if (readline)
        {
            std::stringstream linestream(line);
            idsv.clear();
            //if line ends with " ", we read the last value 2 times
            while (!linestream.eof())
            {
                linestream >> id;
                idsv.push_back(id);
            }
            ids.push_back(idsv);//ids[linenum] = idsv;
            ++time;
        }
        ++linenum;
    }
    assert(time == times.size());

    return firstline;
}

template <typename T>
std::ostream &PrintCorrelations(std::ostream &out, T correlations)
{
    for (auto cor : correlations)
    {
        out << cor[0] << " " << cor[1] << std::endl;
    }
    return out;
}

int main(int argc, char** argv)
{
    /*TODO

    TODO*/

    // x = yes (100%), o = yes (but still doubts), / = no
    /*
    name                        implemented        t:flow          t:vicsek         t:SPR            t: exp splined    speed(SPR,N=1000,t=100,rbins=64) in min
    spatial vel-vel                  o                o               o               o                     o                         4
    spatial dir-dir                  o                o               o               o                     o                         4
    spatial ori-ori                  o                o               o               o                     o                         4
    spatial vor-vor                  o                o               o               o                     o                         60
    temporal vel-vel                 o                o               o               o                     o                         <1
    temporal dir-dir                 o                o               o               o                     o                         <1
    temporal ori-ori                 o                o               o               o                     o                         <1
    temporal vor-vor                 o                o               o               o                     o                         <1
    temporal auto vel-vel            o                o               o               o                     o                         <1
    temporal auto dir-dir            o                o               o               o                     o                         <1
    temporal auto ori-ori            o                o               o               o                     o                         <1
    temporal auto vor-vor            o                o               o               o                     o                         <1
    temporal auto pos-pos            o                o               o               o                     o                         <1
    enstrophy                        o                o               o               o                     o                         <1
    velocity histogram               o                o               o               o                     o                         <1
    vorticity histogram              o                o               o               o                     o                         <1
    number fluctuations              o                o               o               o                     o                         <1
    per id velocity histogram        o                o               o               o                     o                         <1
    traveled pathlength              o                o               o               o                     o                         <1
    cluster size over time           o                o               o               o                     o                          1
    */

    // spatial vicsek does not crash, but has values < -1
    // temporal vorticity is not implemented yet

    ////////////////////////////////////////////////////////////////////////////
    // modify this
    ////////////////////////////////////////////////////////////////////////////

    if (argc != 31)
    {
        throw std::runtime_error("missing parameters");
    }

    std::string path = argv[1];
    std::string filename = argv[2];
    std::string filename_ids = argv[3];

    bool ids_from_file = std::stoi(argv[4]); //false for vicsek and spr, true for tracking
    bool read_flow = std::stoi(argv[5]); // true for flow
    bool periodic_boundary = std::stoi(argv[6]); //true for visek and spr, false for tracking
    bool periodic_velocities = std::stoi(argv[7]); //true for spr and vicsek

    real L = std::stof(argv[8]);

    real cluster_com_dist = std::stof(argv[9]);
    real cluster_angle_dist = std::stof(argv[10]);

    real vel_histo_start = std::stof(argv[11]);
    real vel_histo_end = std::stof(argv[12]);
    real vel_histo_step = std::stof(argv[13]);
    real vor_histo_start = std::stof(argv[14]);
    real vor_histo_end = std::stof(argv[15]);
    real vor_histo_step = std::stof(argv[16]);

    real rmin = std::stof(argv[17]);
    real rmax = std::stof(argv[18]);
    real rstep = std::stof(argv[19]);

    unsigned int tstart = std::stoi(argv[20]);
    unsigned int tend = std::stoi(argv[21]);
    unsigned int tstep = std::stoi(argv[22]);
    unsigned int Xvortstart = std::stoi(argv[23]);
    unsigned int Xvortend = std::stoi(argv[24]);
    unsigned int Xvortstep = std::stoi(argv[25]);
    unsigned int times_start = std::stoi(argv[26]);
    unsigned int times_stop = std::stoi(argv[27]);
    unsigned int times_step = std::stoi(argv[28]);
    unsigned int vorticity_bins = std::stoi(argv[29]);
    unsigned int FlowGridPoints = std::stoi(argv[30]);



    //spatial
    std::string outputfilename_xvel = path + filename + "-cor-xvel.txt";
    std::string outputfilename_xori = path + filename + "-cor-xori.txt";
    std::string outputfilename_xdir = path + filename + "-cor-xdir.txt";
    std::string outputfilename_xvor = path + filename + "-cor-xvor.txt";

    std::cout << "Saving output: " << outputfilename_xvel << std::endl;

    //temporal
    std::string outputfilename_tvel = path + filename + "-cor-tvel.txt";
    std::string outputfilename_tori = path + filename + "-cor-tori.txt";
    std::string outputfilename_tdir = path + filename + "-cor-tdir.txt";
    std::string outputfilename_tvor = path + filename + "-cor-tvor.txt";

    //temporal auto
    std::string outputfilename_tautovel = path + filename + "-cor-tautovel.txt";
    std::string outputfilename_tautoori = path + filename + "-cor-tautoori.txt";
    std::string outputfilename_tautodir = path + filename + "-cor-tautodir.txt";
    std::string outputfilename_tautovor = path + filename + "-cor-tautovor.txt";
    std::string outputfilename_tautopos = path + filename + "-cor-tautopos.txt";

    std::string outputfilename_enst = path + filename + "-cor-enstrophy.txt";
    std::string outputfilename_peridvelhisto = path + filename + "-cor-peridvelhisto.txt";
    std::string outputfilename_velhisto = path + filename + "-cor-velhisto.txt";
    std::string outputfilename_vorhisto = path + filename + "-cor-vorhisto.txt";
    std::string outputfilename_numberfluct = path + filename + "-cor-numberfluct.txt";
    std::string outputfilename_traveledpath = path + filename + "-cor-traveledpath.txt";
    std::string outputfilename_clustersizes = path + filename + "-cor-clustersizes.txt";

    //spatial
    std::ofstream out_xvel, out_xori, out_xdir, out_xvor;
    out_xvel.open(outputfilename_xvel);
    out_xori.open(outputfilename_xori);
    out_xdir.open(outputfilename_xdir);
    out_xvor.open(outputfilename_xvor);

    //temporal
    std::ofstream out_tvel, out_tori, out_tdir, out_tvor;
    out_tvel.open(outputfilename_tvel);
    out_tori.open(outputfilename_tori);
    out_tdir.open(outputfilename_tdir);
    out_tvor.open(outputfilename_tvor);

    //temporal auto
    std::ofstream out_tautovel, out_tautoori, out_tautodir, out_tautovor, out_tautopos;
    out_tautovel.open(outputfilename_tautovel);
    out_tautoori.open(outputfilename_tautoori);
    out_tautodir.open(outputfilename_tautodir);
    out_tautovor.open(outputfilename_tautovor);
    out_tautopos.open(outputfilename_tautopos);

    std::ofstream out_enst, out_peridvelhisto, out_velhisto, out_vorhisto, out_numberfluct, out_traveledpath, out_clustersizes;
    out_enst.open(outputfilename_enst);
    out_peridvelhisto.open(outputfilename_peridvelhisto);
    out_velhisto.open(outputfilename_velhisto);
    out_vorhisto.open(outputfilename_vorhisto);
    out_numberfluct.open(outputfilename_numberfluct);
    out_traveledpath.open(outputfilename_traveledpath);
    out_clustersizes.open(outputfilename_clustersizes);


    auto ALLTIMES = MakeTimes(times_start, times_stop, times_step);

    ////////////////////////////////////////////////////////////////////////////
    //read in files into data structures
    ////////////////////////////////////////////////////////////////////////////

    //read in positions and angles
    std::cout << "Read positions and angles" << std::endl;
    Posdict positions_dict;
    Angdict angles_dict;
    read_data_single(path + filename, positions_dict, angles_dict, ALLTIMES);


    std::cout << "Build positions and angles vector" << std::endl;
    std::vector<std::vector<Vector>> positions;
    std::vector<std::vector<real>> angles;
    std::vector<std::vector<Vector>> orientations;

    for (std::pair<TIME, std::vector<Vector>> element : positions_dict)
    {
        positions.push_back(element.second);
    }
    for (std::pair<TIME, std::vector<real>> element : angles_dict)
    {
        std::vector<Vector> tmp(element.second.size());
        for (unsigned int i = 0; i < tmp.size(); ++i)
        {
            tmp[i] = Vector(std::cos(element.second[i]), std::sin(element.second[i]));
        }
        orientations.push_back(tmp);
        angles.push_back(element.second);
    }

    std::cout << "Read or calculate ids" << std::endl;

    //ids: read in from file or generate automatically based on indice in
    //positions/angle/... list
    std::vector<std::vector<ID>> ids;
    int N;
    if (ids_from_file)
    {
        read_ids(path + filename_ids, ids, ALLTIMES);
        N = ids.front().size();
    }
    else
    {
        N = positions_dict.begin()->second.size();
        //ids are always 0 to N-1
        for (auto _: ALLTIMES)
        {
            ids.push_back(MakeTimes(0, N, 1));
        }
    }


    assert(positions.size() == ALLTIMES.size());
    unsigned int times = positions.size();

    std::cout << "N: " << N << ", times: " << times << std::endl; //just N of the first line (for debug)

    //need 47 byte per time and particle (400 MB @ 9000 times and 1000 particles)
    //std::vector<std::map<ID, Vector>> positions_ided(times);
    //std::vector<std::map<ID, Vector>> velocities_ided(times);
    //std::vector<std::map<ID, Vector>> directions_ided(times);
    //std::vector<std::map<ID, Vector>> orientations_ided(times);
    std::vector<std::vector<ID>> ids_vel(times-1);
    std::vector<std::vector<Vector>> velocities(times-1);
    std::vector<std::vector<Vector>> positions_vel(times-1);

    std::cout << "Merge with ids" << std::endl;

    /*
    TIME timeindex;
    ID idindex;

    // positions: make map time->(id->position)
    timeindex = 0;
    for (std::vector<Vector> element : positions) // loop over times
    {
        positions_ided[timeindex] = std::map<ID, Vector>();
        idindex = 0;
        for (Vector particle : element) // element.second is list of Vectors
        {
            positions_ided[timeindex][ids[timeindex][idindex]] = particle;
            idindex += 1;
        }
        timeindex += 1;
    }

    // orientations: make map time->(id->orientation)
    timeindex = 0;
    for (std::vector<ANGLE> element : angles) // loop over times
    {
        orientations_ided[timeindex] = std::map<ID, Vector>(); // element.first is the key, which is the time
        idindex = 0;
        for (ANGLE angle : element) // element.second is list of Vectors
        {
            orientations_ided[timeindex][ids[timeindex][idindex]] = Vector(std::cos(angle), std::sin(angle));
            idindex += 1;
        }
        timeindex += 1;
    }

    // velocities based on position differences
    Vector v1, v2;
    for (TIME t = 0; t < times - 1; ++t)
    {
        velocities_ided[t] = std::map<ID, Vector>();
        for (std::pair<ID, Vector> entry : positions_ided[t])
        {
            try
            {
                v1 = entry.second;
                v2 = positions_ided[t+1].at(entry.first);
                if (periodic_velocities)
                {
                    velocities_ided[t][entry.first] = distance_periodic_vector(v2, v1, L);
                }
                else
                {
                    velocities_ided[t][entry.first] = v2 - v1;
                }
            }
            catch (std::out_of_range)
            {
                // pass
            }
        }
    }

    if (read_flow)
    {
        velocities_ided = positions_ided;
        unsigned int pL = static_cast<unsigned int>(L);
        unsigned int FlowStep = pL/FlowGridPoints;

        // calculating positions for flow
        for (TIME t = 0; t < times; ++t)
        {
            positions_ided[t] = std::map<ID, Vector>();
            ID idindex = 0;
            for (int x = 0; x < pL; x += FlowStep)
            {
                for (int y = 0; y < pL; y += FlowStep)
                {
                    positions_ided[t][idindex] = Vector(x, y);
                    idindex += 1;
                }
            }
        }
    }

    //directions (normalized velocities)
    for (TIME t = 0; t < times - 1; ++t)
    {
        directions_ided[t] = std::map<ID, Vector>();
        for (std::pair<ID, Vector> entry : velocities_ided[t])
        {
            Vector norm = entry.second.Normal();
            real norm_real = entry.second.Norm();
            if (norm_real == 0.)
            {
                directions_ided[t][entry.first] = Vector(0, 0);
            }
            else
            {
                directions_ided[t][entry.first] = norm;
            }
        }
    }
    */


    if (read_flow)
    {
        velocities = positions;
        ids_vel = ids;
        unsigned int pL = static_cast<unsigned int>(L);
        unsigned int FlowStep = pL/FlowGridPoints;

        // calculating positions for flow
        for (TIME t = 0; t < positions.size(); ++t)
        {
            positions[t] = std::vector<Vector>(FlowGridPoints*FlowGridPoints);
            ID idindex = 0;
            for (unsigned int x = 0; x < pL; x += FlowStep)
            {
                for (unsigned int y = 0; y < pL; y += FlowStep)
                {
                    assert(idindex < FlowGridPoints*FlowGridPoints);
                    positions[t][idindex] = Vector(x, y);
                    idindex += 1;
                }
            }
        }
        positions_vel = positions;
    }
    else
    {
        // velocities based on position differences
        Vector v1, v2;
        for (TIME t = 0; t < velocities.size(); ++t)
        {
            velocities[t] = std::vector<Vector>();
            for (unsigned int j = 0; j < positions[t].size(); ++j)
            {
                v1 = positions[t][j];
                ID id = ids[t][j];
                auto idit = std::find(ids[t+1].begin(), ids[t+1].end(), id); //-> catch
                if (idit == ids[t+1].end())
                {
                    continue;
                }
                else
                {
                    v2 = positions[t+1][distance(ids[t+1].begin(), idit)];
                }

                if (periodic_velocities)
                {
                    velocities[t].push_back(distance_periodic_vector(v2, v1, L));
                }
                else
                {
                    velocities[t].push_back(v2 - v1);
                }
                positions_vel[t].push_back(v1);
                ids_vel[t].push_back(id);
            }
        }
    }

    std::vector<std::vector<Vector>> directions(velocities.size());

    //directions (normalized velocities)
    for (TIME t = 0; t < directions.size(); ++t)
    {
        directions[t] = std::vector<Vector>(velocities[t].size());
        for (unsigned int j = 0; j < velocities[t].size(); ++j)
        {
            real norm_real = velocities[t][j].Norm();
            if (norm_real == 0.)
            {
                directions[t][j] = Vector(0, 0);
            }
            else
            {
                directions[t][j] = velocities[t][j].Normal();
            }
        }
    }

    ////////////////////////////////////////////////////////////////////////////
    //calculate stuff
    ////////////////////////////////////////////////////////////////////////////


    std::cout << "Calculate spatial correlations (vel, ori, dir)" << std::endl;
    std::vector<std::array<real, 2>> cor_xvel, cor_xori, cor_xdir;

    if (periodic_boundary)
    {
        cor_xvel = correlationPoneXVel<true>(positions_vel, velocities, rmin, rmax, rstep, tstart, tend, tstep, L);
        cor_xori = correlationPoneXVel<true>(positions, orientations, rmin, rmax, rstep, tstart, tend, tstep, L);
        cor_xdir = correlationPoneXVel<true>(positions_vel, directions, rmin, rmax, rstep, tstart, tend, tstep, L);

    }
    else
    {
        cor_xvel = correlationPoneXVel<false>(positions_vel, velocities, rmin, rmax, rstep, tstart, tend, tstep, L);
        cor_xori = correlationPoneXVel<false>(positions, orientations, rmin, rmax, rstep, tstart, tend, tstep, L);
        cor_xdir = correlationPoneXVel<false>(positions_vel, directions, rmin, rmax, rstep, tstart, tend, tstep, L);
    }
    PrintCorrelations(out_xvel, cor_xvel);
    PrintCorrelations(out_xori, cor_xori);
    PrintCorrelations(out_xdir, cor_xdir);

    std::cout << "Calculate spatial correlations finished" << std::endl;


    std::cout << "Calculate spatial correlation (vor)" << std::endl;
    std::vector<std::array<real, 3>> cor_xvor;
    cor_xvor = correlationPoneXVor(positions, velocities, rmin, rmax, rstep, Xvortstart, Xvortend, Xvortstep, vorticity_bins, L);
    PrintCorrelations(out_xvor, cor_xvor);
    std::cout << "Calculate spatial correlation (vor) finished" << std::endl;

    std::cout << "Calculate temporal correlations" << std::endl;
    std::vector<std::array<real, 2>> cor_tvel, cor_tori, cor_tdir, cor_tvor;


    cor_tvel = correlationPoneTVel(positions_vel, velocities, tstart, tend, tstep, vorticity_bins, L);
    PrintCorrelations(out_tvel, cor_tvel);
    cor_tori = correlationPoneTVel(positions, orientations, tstart, tend, tstep, vorticity_bins, L);
    PrintCorrelations(out_tori, cor_tori);
    cor_tdir = correlationPoneTVel(positions_vel, directions, tstart, tend, tstep, vorticity_bins, L);
    PrintCorrelations(out_tdir, cor_tdir);
    cor_tvor = correlationPoneTVor(positions_vel, velocities, tstart, tend, tstep, vorticity_bins, L);
    PrintCorrelations(out_tvor, cor_tvor);


    std::cout << "Calculate temporal correlations finished" << std::endl;

    std::cout << "Calculate auto correlations" << std::endl;
    std::vector<std::array<real, 2>> cor_tautovel, cor_tautoori, cor_tautodir, cor_tautovor;
    std::vector<std::array<real, 2>> cor_tautopos;


    cor_tautovel = correlationPoneTAutoVel(velocities, ids_vel, tstart, tend, tstep, L);
    PrintCorrelations(out_tautovel, cor_tautovel);
    cor_tautoori = correlationPoneTAutoVel(orientations, ids, tstart, tend, tstep, L);
    PrintCorrelations(out_tautoori, cor_tautoori);
    cor_tautodir = correlationPoneTAutoVel(directions, ids_vel, tstart, tend, tstep, L);
    PrintCorrelations(out_tautodir, cor_tautodir);
    cor_tautopos = correlationPoneTAutoVel(positions, ids, tstart, tend, tstep, L);
    PrintCorrelations(out_tautopos, cor_tautopos);
    cor_tautovor = correlationPoneTAutoVor(positions_vel, velocities, ids_vel, tstart, tend, tstep, vorticity_bins, L);
    PrintCorrelations(out_tautovor, cor_tautovor);


    std::cout << "Calculate auto correlations finished" << std::endl;


    std::cout << "Calculate number fluctuations" << std::endl;
    std::vector<std::array<real, 3>> number_flucts;
    number_flucts = number_fluctuations(positions, L);
    for (auto nf : number_flucts)
    {
        out_numberfluct << nf[0] << " " << nf[1] << " " << nf[2] << std::endl;
    }
    std::cout << "Calculate number fluctuations finished" << std::endl;


    std::cout << "Calculate enstrophy" << std::endl;
    real enst;
    enst = Enstrophy(positions_vel, velocities, tstart, tend, tstep, vorticity_bins, L);
    out_enst << enst << std::endl;
    std::cout << "Calculate enstrophy finished" << std::endl;


    std::cout << "Calculate velocity histogram" << std::endl;
    std::vector<std::array<real, 3>> velhisto;
    //vicsek allways v0 * dt (=0.01)
    //flow: -0.1, 11., 0.1, vicsek: 0.0, 0.02, 0.0001, spr: -0.001, 1, 0.001, splines: -0.1, 40, 0.1
    velhisto = histogramXVel(velocities, vel_histo_start, vel_histo_end, vel_histo_step, tstart, tend, tstep, L);
    for (auto vh : velhisto)
    {
        out_velhisto << vh[0] << " " << vh[1] << " " << vh[2] << std::endl;
    }
    std::cout << "Calculate velocity histogram finished" << std::endl;


    std::cout << "Calculate vorticity histogram" << std::endl;
    std::vector<std::array<real, 3>> histo_xvor;
    //flow: -0.46, 0.46, 0.01, vicsek: -0.04, 0.04, 0.001, spr: -0.34, 0.34, 0.01, splines: -2.2, 2.2, 0.01
    histo_xvor = histogramXVor(positions_vel, velocities, vor_histo_start, vor_histo_end, vor_histo_step, tstart, tend, tstep, vorticity_bins, L);
    for (auto nf : histo_xvor)
    {
        out_vorhisto << nf[0] << " " << nf[1] << " " << nf[2] << std::endl;
    }
    std::cout << "Calculate vorticity histogram finished" << std::endl;


    std::cout << "Calculate per id vel histogram" << std::endl;
    std::map<ID, std::vector<std::array<real, 3>>> perIDvelhisto;
    //flow: -0.1, 11., 0.1, vicsek: 0.0, 0.02, 0.0001, spr: -0.001, 1, 0.001, splines: -0.1, 40, 0.1
    perIDvelhisto = perIDhistogramVel(velocities, ids_vel, vel_histo_start, vel_histo_end, vel_histo_step, tstart, tend, tstep, L);
    for (auto idhisto : perIDvelhisto)
    {
        for (auto velarray : idhisto.second)
        {
            out_peridvelhisto << velarray[0] << " " << velarray[1] << " ";
        }
        out_peridvelhisto << std::endl;
    }
    std::cout << "Calculate per id vel histogram finished" << std::endl;


    std::cout << "Calculate per id traveled path" << std::endl;
    std::map<ID, std::array<real, 3>> perIDtraveledpath;
    perIDtraveledpath = perIDtraveledPathlength(velocities, ids_vel, tstart, tend, tstep);
    for (auto path : perIDtraveledpath)
    {
        out_traveledpath << path.second[0] << " " << path.second[1] << " " << path.second[2] << std::endl;
    }
    std::cout << "Calculate per id traveled path finished" << std::endl;


    std::cout << "Calculate cluster sizes" << std::endl;
    std::vector<std::map<unsigned int, unsigned int>> clustersizes;

    //flow: 32., 15./360*2*PI, viscek: 1., 15./360*2*PI, spr: 2*7., 15./360*2*PI, splines: 94., 15./360*2*PI:
    clustersizes = clusterSizes(positions, angles, tstart, tend, tstep, cluster_com_dist, cluster_angle_dist, L, periodic_boundary);
    for (auto clusters : clustersizes)
    {
        for (auto cluster : clusters)
        {
            //size of cluster, how many
            out_clustersizes << cluster.first << " " << cluster.second << " ";
        }
        out_clustersizes << std::endl;
    }
    std::cout << "Calculate cluster sizes finished" << std::endl;

    std::cout << "Evaluation finished." << std::endl;


    return 0;
}






