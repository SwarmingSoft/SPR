//constants.cpp v16 30.10.2015

#include "constantsSPR.hpp"
#include <array>

/*
see paper: Dynamical maximum entropy approach to flocking
*/

real KroneckerDelta(int i, int j)
{
    return real(i == j);
}

void create_ni(const MATRIX &n, std::vector<real> &ni)
{
    assert(n.GetWidth() == ni.size() && n.GetHeight() == ni.size());
    for (unsigned int i = 0; i < ni.size(); ++i)
    {
        ni[i] = get_number_of_neighbours(n, i);
        //std::cout << ni[i] << "\n";
    }
}

//calculates n_i after formula directly after eq (25)
real get_number_of_neighbours(const MATRIX &n, unsigned int i)
{
    real neighbours = 0;
    for (unsigned int j = 0; j < n.GetHeight(); ++j)
    {
        if (*n.Get(i, j) != 0)
        {
            neighbours += *n.Get(i, j);
        }
    }
    return neighbours;
}

//calculates n_C after formula directly after eq (30) in paper
real get_average_number_of_neighbours(const MATRIX &n)
{
    real sum_of_neighbours = 0;
    for (unsigned int i = 0; i < n.GetWidth(); ++i)
    {
        sum_of_neighbours += get_number_of_neighbours(n, i);
    }
    return sum_of_neighbours / n.GetWidth();
}

real average(const std::vector<real> &ni)
{
    assert(ni.size() > 0);
    real sum = 0;
    for (real r : ni)
    {
        sum += r;
    }
    return sum / ni.size();
}

real distance_periodic(Vector p1, Vector p2, real L)
{
    real absx = std::abs(p1.x-p2.x);
    real absy = std::abs(p1.y-p2.y);
    real minx = std::min(absx, L-absx);
    real miny = std::min(absy, L-absy);
    return std::sqrt(minx*minx + miny*miny);
}

/*void calculate_n_voronoi(MATRIX &n, const std::vector<Vector> &positions, real L)
{
    assert(n.GetWidth() == positions.size() && n.GetHeight() == positions.size());
    voronoicell_neighbor cell;
    std::vector<int> neighbours;
	unsigned int i, j;
	int id;

	//TODO 6.4 particles per box (N/(int(x)**3 == 6.4 => solve for int(x) for arguments (int(x), int(x), 1) before true, true, false...)  (look here: http://math.lbl.gov/voro++/examples/timing_test/); already tested data.txt for same hash for different grid sizes
	container con(0, L, 0, L, 0, L, 1, 1, 1, true, true, false, positions.size());

    for (i = 0; i < positions.size(); ++i)
    {
        con.put(i, positions[i].x, positions[i].y, L/2.);
    }
	n.Zero(); //initializes matrix with n_ij = 0 for all i, j

	// Loop over all particles in the container and compute each Voronoi cell
	c_loop_all loop(con);
   	if (loop.start())
    {
        do
        {
            if (con.compute_cell(cell, loop)) //this line does all the work for us
            {
		        id = loop.pid(); //index of particle under inspection

		        // Gather information about the computed Voronoi cell
		        cell.neighbors(neighbours);

                for (j = 0; j < neighbours.size(); ++j)
                {
                    if (neighbours[j] >= 0 && id != neighbours[j]) //otherwise you sometimes get -5 and -6 as neighbouring indices; these are the boundaries in z direction; particle can be neighbour of itself in rare cases (usually when N is low)
                    {
                        //std::cout << "Set(" << id << ", " << neighbours[j] << ")\n";
                        n.Set(id, neighbours[j], 1);
                        n.Set(neighbours[j], id, 1);
                    }
                }
	        }
        }
        while (loop.inc());
    }

    //con.draw_particles("output/random_points_p.gnu");
    //con.draw_cells_gnuplot("output/random_points_v.gnu");
}*/

void calculate_n_metric(MATRIX &n, const std::vector<Vector> &positions, real L, real r0)
{
    assert(n.GetWidth() == positions.size() && n.GetHeight() == positions.size());
    real distance;
    unsigned int i, j;
    for (i = 0; i < positions.size(); ++i)
    {
        for (j = i+1; j < positions.size(); ++j)
        {
            distance = distance_periodic(positions[i], positions[j], L);
            if (distance <= r0)
            {
                n.Set(i,j,1);
                n.Set(j,i,1);
            }
            else
            {
                n.Set(i,j,0);
                n.Set(j,i,0);
            }
        }
    }
}

bool sort_func(std::tuple<real,unsigned int> i, std::tuple<real,unsigned int> j)
{
    return std::get<0>(i) < std::get<0>(j);
}

//n not symmetrical
void calculate_n_topological_unsymmetrical(MATRIX &n, const std::vector<Vector> &positions, real L, unsigned int neighbours)
{
    assert(n.GetWidth() == positions.size() && n.GetHeight() == positions.size());
    real distance;
    unsigned int i, j;
    std::vector<std::tuple<real,unsigned int>> particles(positions.size());
    n.Zero();
    unsigned int ncount;
    for (i = 0; i < positions.size(); ++i)
    {
        //particles.clear()
        for (j = 0; j < positions.size(); ++j)
        {
            distance = distance_periodic(positions[i], positions[j], L);
            particles[j] = std::make_tuple(distance, j);
        }

        std::sort(particles.begin(), particles.end(), sort_func);

        /*for (auto p : particles)
        {
            std::cout << std::get<0>(p) << " " << std::get<1>(p) << std::endl;
        }*/

        ncount = 0;
        for (auto pit = std::next(particles.begin()); pit != particles.end() && ncount < neighbours; ++pit, ++ncount)
        {
            //std::cout << std::get<0>(*pit) << " " << std::get<1>(*pit) << std::endl;
            n.Set(i,std::get<1>(*pit),1);
            //n.Set(std::get<1>(*pit),i,1);
        }

        //std::cout << "\n";
    }
}

//n symmetrical and real
void calculate_n_topological(MATRIX &n, const std::vector<Vector> &positions, real L, unsigned int neighbours)
{
    MATRIX tmp(n.GetWidth(), n.GetHeight());
    calculate_n_topological_unsymmetrical(tmp, positions, L, neighbours);

    unsigned int x, y;

    for (x = 0; x < n.GetWidth(); ++x)
    {
        for (y = 0; y < n.GetHeight(); ++y)
        {
            n.Set(x, y, (*tmp.Get(x, y) + *tmp.Get(y, x) ) / 2.);
        }
    }
}

void print_n(const MATRIX &n)
{
    unsigned int x, y;

    for (y = 0; y < n.GetHeight(); ++y)
    {
        for (x = 0; x < n.GetWidth(); ++x)
        {
            std::cout << *n.Get(y, x) << " ";
        }
        std::cout << "\n";
    }
}

real mixing_parameter(real deltat, const MATRIX &n1, const MATRIX &n2, real nc)
{
    assert(n1.GetHeight() == n1.GetWidth() && n2.GetHeight() == n2.GetWidth() && n1.GetHeight() == n2.GetHeight());
    assert(nc > 0.);
    unsigned int x, y;
    unsigned int sum = 0;

    for (y = 0; y < n1.GetHeight(); ++y)
    {
        for (x = 0; x < n1.GetWidth(); ++x)
        {
           sum += std::abs(int(*n1.Get(y, x)) - int(*n2.Get(y, x)));
        }
    }
    return real(sum) / deltat / (n1.GetWidth() * nc);
}

real order_parameter(const std::vector<Vector> &st0)
{
    Vector tmp(0, 0);

    for (Vector v : st0)
    {
        tmp += (v / v.Norm());
    }

    return tmp.Norm() / st0.size();
}

real order_parameter_angle(const std::vector<real> &angles)
{
    Vector tmp(0, 0);

    for (real a : angles)
    {
        tmp += Vector(std::cos(a), std::sin(a));
    }

    return tmp.Norm() / angles.size();
}


Pi::Pi (const Angdict &angles)
{
    unsigned int i;
    for (auto anglespair : angles)
    {
        //normalized vector of average moving direction
        Vector normn(0, 0);
        for (auto angle : anglespair.second)
        {
            normn.x += std::cos(angle);
            normn.y += std::sin(angle);
        }
        normn.Normalize();

        assert_almost_equal(normn.Norm(), real(1));

        pis[anglespair.first] = std::vector<Vector>(anglespair.second.size());
        i = 0;
        for (auto angle : anglespair.second)
        {
            //v - (v.n)*n
            Vector v(std::cos(angle), std::sin(angle));
            assert_almost_equal(v.Norm(), real(1));

            pis[anglespair.first][i] = v-(normn*(v*normn));
            //std::cout << "(" << pis[anglespair.first][i].x  << " " << pis[anglespair.first][i].y << " " << pis[anglespair.first][i].Norm() << ")  ";
            ++i;
        }
         //std::cout << "\n";
    }
    if (!pis.empty())
    {
        mysize = pis.begin()->second.size();
    }
}

Vector Pi::Get(unsigned int t, unsigned int i) const
{
    return pis.at(t)[i];
}

const std::vector<Vector> &Pi::Get(unsigned int t) const
{
    return pis.at(t);
}

unsigned int Pi::size(void) const
{
    return mysize;
}

//the following C are the correlation functions defined in TABLE 1 (page 5) of paper
real C1s(unsigned int t, Pi &pi)
{
    real sum = 0;
    Vector tmp;
    unsigned int i;
    for (i = 0; i < pi.size(); ++i)
    {
        tmp = pi.Get(t+1, i);
        sum += tmp*tmp;
    }
    return sum / pi.size();
}

real Cs(unsigned int t, Pi &pi)
{
    real sum = 0;
    Vector tmp;
    unsigned int i;
    for (i = 0; i < pi.size(); ++i)
    {
        tmp = pi.Get(t, i);
        sum += tmp*tmp;
    }
    return sum / pi.size();
}

real Gs(unsigned int t, Pi &pi)
{
    real sum = 0;
    unsigned int i;
    for (i = 0; i < pi.size(); ++i)
    {
        sum += pi.Get(t+1, i)*pi.Get(t, i);
    }
    return sum / pi.size();
}

real CTs(unsigned int t, Pi &pi, real nc, const std::vector<real> &ni)
{
    assert(nc > 0.);
    real sum = 0;
    Vector tmp;
    unsigned int i;
    for (i = 0; i < pi.size(); ++i)
    {
        tmp = pi.Get(t, i);
        sum += tmp*tmp*ni[i];
    }
    return sum / (pi.size() * nc);
}

real GTs(unsigned int t, Pi &pi, real nc, const std::vector<real> &ni)
{
    assert(nc > 0.);
    real sum = 0;
    unsigned int i;
    for (i = 0; i < pi.size(); ++i)
    {
        sum += pi.Get(t+1, i)*pi.Get(t, i)*ni[i];
    }
    return sum / (pi.size() * nc);
}

real Cint(unsigned int t, Pi &pi, real nc, MATRIX &n)
{
    assert(nc > 0.);
    real sum = 0.;
    unsigned int i, j;
    //std::cout << "cint summands:\n";
    for (i = 0; i < pi.size(); ++i)
    {
        for (j = 0; j < pi.size(); ++j)
        {
            sum += pi.Get(t, i)*pi.Get(t, j)*real(*n.Get(i, j));
            //std::cout << pi.Get(t, i)*pi.Get(t, j)*real(*n.Get(i, j)) << "\n";
        }
    }
    //std::cout << "sum: " << sum << "\n";
    //std::cout << "\n";
    return sum / (pi.size() * nc);
}

real CPint(unsigned int t, Pi &pi, real nc, const MATRIX &n)
{
    assert(nc > 0.);
    real sum = 0;
    unsigned int i, j, k;
    for (i = 0; i < pi.size(); ++i)
    {
        for (j = 0; j < pi.size(); ++j)
        {
            for (k = 0; k < pi.size(); ++k)
            {
                sum += pi.Get(t, j)*pi.Get(t, k)*real(*n.Get(i, j))*real(*n.Get(i, k));
            }
        }
    }
    return sum / (pi.size() * nc * nc);
}

real Gint(unsigned int t, Pi &pi, real nc, const MATRIX &n)
{
    assert(nc > 0.);
    real sum = 0;
    unsigned int i, j;
    for (i = 0; i < pi.size(); ++i)
    {
        for (j = 0; j < pi.size(); ++j)
        {
            sum += pi.Get(t+1, i)*pi.Get(t, j)*real(*n.Get(i, j));
        }
    }
    return sum / (pi.size() * nc);
}

real CHs(unsigned int t, Pi &pi, real nc, const std::vector<real> &ni)
{
    assert(nc > 0.);
    real sum = 0;
    Vector tmp;
    unsigned int i;
    for (i = 0; i < pi.size(); ++i)
    {
        tmp = pi.Get(t, i)*ni[i]; //NO j USED !!! ERROR IN PAPER
        sum += tmp*tmp;
    }
    return sum / (pi.size() * nc * nc);
}

real CTint(unsigned int t, Pi &pi, real nc, const std::vector<real> &ni, const MATRIX &n)
{
    assert(nc > 0.);
    real sum = 0;
    unsigned int i, j;
    for (i = 0; i < pi.size(); ++i)
    {
        for (j = 0; j < pi.size(); ++j)
        {
            sum += pi.Get(t, i)*pi.Get(t, j)*ni[i]*real(*n.Get(i,j));
        }
    }
    return sum / (pi.size() * nc * nc);
}

//order of parameter
//t, pi, nc, ni, n
//C1s, Cs, Gs, CTs, GTs, Cint, CPint, Gint, CHs, CTint

//calculates J after eq (31) in paper
real J(real nc, real T0, real Omega, real CPint, real CHs, real CTint)
{
    assert(nc > 0.);
    return 1./nc * (Omega + T0) / (CPint + CHs - 2.*CTint);
}

//calculates T after eq (32) in paper
/*
real T(real nc, real J, real T0, real Omega, real C1s, real Cs, real CTs, real GTs)
{
    return T0 + (C1s - Cs) / (2.*deltat) - (J*nc*deltat) / 2. * ((CTs - GTs)/deltat + Omega);
}*/

//is equal to eq. 32 in paper, see: Appendix-B.nb
real T(real deltat, real nc, real J, real T0, real Omega, real C1s, real Cs, real Gs, real CTs, real GTs, real Cpint, real Chs, real Ctint)
{
    return T0/2. + (C1s - Cs + Cs - Gs) / (2.*deltat) - (J*nc*deltat) / 2. * (2.*(CTs - GTs)/deltat + 2*Omega - J*nc*(Chs + Cpint - 2*Ctint));
}

real T0(real deltat, real Cs, real Gs)
{
    return (Cs - Gs) / deltat;
}

real T0t(real deltat, real Cts, real Gts)
{
    return (Cts - Gts) / deltat;
}

//calculates Omega after eq (33) in paper; Omega is related to dynamics network
real Omega(real deltat, real Cint, real Gint)
{
    assert(deltat >= 0.);
    return (Gint - Cint) / deltat;
}

real eta_of_t(real temperature)
{
    assert_log(temperature >= 0.);
    return std::sqrt(temperature*6./(PI*PI));
}

real JV_of_J(real deltat, real interaction, real nv)
{
    return interaction / (1. - deltat * interaction * nv);
}

//Dyn. max ent. B1-B2
real DynamicLogLikelihood(real deltat, real N, real nc, real C1s, real Cs, real Gs, real CTs, real GTs, real Cint, real CPint, real Gint, real CHs, real CTint, real J, real T)
{
    real alpha = J*nc*deltat;
    real lh = C1s+Cs-2*alpha*CTs+alpha*alpha*CHs+2*alpha*(Cint-alpha*CTint)+alpha*alpha*CPint-2*alpha*Gint-2*(Gs-alpha*GTs);
    return -(N/2.) * std::log(2*PI*2*T*deltat) - N*lh/(4*T*deltat);
}

//Dyn. max ent. B6
real DynamicLogLikelihoodMaximized(real N, real C1s, real Cs, real Gs, real CTs, real GTs, real Cint, real CPint, real Gint, real CHs, real CTint)
{
    real lt = (Cint-CTs+GTs-Gint)*(Cint-CTs+GTs-Gint)/(2*CTint-CPint-CHs);
    real lh = C1s + Cs - 2*Gs + lt;
    return -(N/2.)*(std::log(2*PI*lh) + 1);
}

void output(Posdict &positions, Angdict &angles, std::vector<unsigned int> &times)
{
    std::cout << "sizeof(times) = " << times.size() << "\n";
    std::array<unsigned int, 2> step = {0, 1};

    std::vector<Vector> positionsv;

    for (auto time : times)
    {
        for (auto s : step)
        {
            positionsv = positions[time+s];
            //std::cout << "times = " << time+s << ", positionsv.size() = " << positionsv.size() << "\n";
            for (auto p : positionsv)
            {
                std::cout << p.x << " " << p.y << " ";
            }
            std::cout << "\n";
            for (auto a : angles[time+s])
            {
                std::cout << a << " ";
            }
            std::cout << "\n";
        }
    }
}

// times = {1,2} reads line {1,2,3}
// times = {1,3} reads line {1,2,3,4}
// times = {1,4} reads line {1,2,4,5}
std::string read_data_consecutive(std::string filename, Posdict &positions, Angdict &angles, std::vector<size_t> &times)
{
    std::ifstream in;
    in.open(filename);
    assert(!in.fail() && "Could not open file");

    real pos1, pos2, angle;
    size_t linenum = 1;
    size_t time = 0;

    std::string firstline, line;
    std::vector<Vector> positionsv; //()
    std::vector<real> anglesv; //()

    std::getline(in, firstline); //skip first line
    while (in.good())
    {
        std::getline(in, line);
        bool readline = linenum == times[time] || linenum == times[time] + 1;
        #ifndef NDEBUG
        std::cerr << "linenum: " << linenum << ", readline: " << readline << ", time: " << time << ", times[time]: " << times[time] << std::endl;
        #endif // NDEBUG
        if (readline)
        {
            std::stringstream linestream(line);
            positionsv.clear();
            //if line ends with " ", we read the last value 2 times
            while (!linestream.eof())
            {
                linestream >> pos1;
                linestream >> pos2;
                positionsv.push_back(Vector(pos1, pos2));
            }
            positions[linenum] = positionsv;
        }
        std::getline(in, line);
        //assert(in.good()); failes why?
        if (readline)
        {
            std::stringstream linestream(line);
            anglesv.clear();
            //if line ends with " ", we read the last value 2 times
            while (!linestream.eof())
            {
                linestream >> angle;
                anglesv.push_back(angle);
            }
            angles[linenum] = anglesv;
        }
        if (readline)
        {
            if (linenum == times[time] + 1)
            {
                ++time;
                if (time == times.size())
                {
                    break;
                }
            }
        }
        ++linenum;
    }
    return firstline;
}

// times = {1,2,3} reads line 1,2,3
std::string read_data_single(std::string filename, Posdict &positions, Angdict &angles, std::vector<size_t> &times)
{
    std::ifstream in;
    in.open(filename);
    assert(!in.fail() && "Could not open file");

    real pos1, pos2, angle;
    size_t linenum = 1;
    size_t time = 0;

    std::string firstline, line;
    std::vector<Vector> positionsv; //()
    std::vector<real> anglesv; //()

    std::getline(in, firstline); //skip first line
    while (in.good())
    {
        //std::cout << time << " " << linenum << std::endl;
        std::getline(in, line);
        bool readline = linenum == times[time];
        if (readline)
        {
            std::stringstream linestream(line);
            positionsv.clear();
            //if line ends with " ", we read the last value 2 times
            while (!linestream.eof())
            {
                linestream >> pos1;
                linestream >> pos2;
                positionsv.push_back(Vector(pos1, pos2));
            }
            positions[linenum] = positionsv;
            ++time;
        }
        //assert(in.good()); failes why?
        std::getline(in, line);
        if (readline)
        {
            std::stringstream linestream(line);
            anglesv.clear();
            //if line ends with " ", we read the last value 2 times
            while (!linestream.eof())
            {
                linestream >> angle;
                anglesv.push_back(angle);
            }
            angles[linenum] = anglesv;
        }
        ++linenum;
    }
    assert(time == times.size());
    return firstline;
}

void set_diagonal_true(MATRIX &n)
{
    assert(n.GetWidth() == n.GetHeight());
    for (unsigned int i = 0; i < n.GetHeight(); ++i)
    {
        n.Set(i, i, true);
    }
}
