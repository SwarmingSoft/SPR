//constants.hpp v16 30.10.2015

#ifndef CONSTANTS_HEADER
#define CONSTANTS_HEADER

#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>
#include <tuple>
#include <cassert>
#include <map>
#include <sstream>


#include "Buffer.hpp"
#include "Vector2D.hpp"

typedef double real; //float, double, ...

const real PI = 3.141592653589793238463;

/* constants in the paper "Dynamical maximum entropy approach to flocking":
2-dimensional torus
L = 32
N = 1024
Voronoi adjacency matrix of mean nonuniform degree n_V = 6
deltat = 0.01
mixing parameter mu = 0.18, 0.35, 0.76
eta = 0.3, 0.2, 0.12
v0 = 0.5, 1.0, 2.0
JV = 1.0, 1.0, 0.1
all display polarization (order parameter) of approx 0.97
*/

//TODO
/*
- r0 -> nc in entropy.cpp
- marginal likelihood?
*/

/*//const unsigned int d = 2; //possible overrides with temporary d; "DIM"?
const unsigned int SEED = 205743;
const unsigned int N = 1024; // 256
const unsigned int L = 32; // 16; //only quadratic or l_1, l_2?

const real deltat = 0.01;
const real ceta = 0.12;
const real v0 = 1.0;
const real JV = 1.0; //better int; dependent on particle i?
//mu = ? (mixing parameter under eq (36) has to be defined a prior, but HOW???)
const int ttl = 110000; //110000 //deltat steps till termination*/

enum class Geometry {voronoi, metric, topological};

//need to use heap memory for matrices larger than 719 x 719. if N > 19853 => error code "0xC00000FD" might occur
typedef GenericBuffer_dyn<real, MemoryConst::Dynamic, MemoryConst::Dynamic> MATRIX;

typedef Vector2DT<real> Vector;

template <typename T, typename C> void _assert_almost_equal(T a, T b, C delta, const char *file, int line)
{
    if (std::abs(a - b) <= std::numeric_limits<T>::epsilon()*delta)
    {
        //good
    }
    else
    {
        std::cout << "Assertion failed!\nFile: " << file << ":" << line << "\n"
        << a << " not almost equal to "<< b << std::endl;
        abort();
    }
}

template <typename T> void _assert_log(bool cond, T s, const char *file, int line)
{
    if (cond)
    {
        //good
    }
    else
    {
        std::cout << s << " (" << file << ":" << line << ")" << std::endl;
    }
}

#ifdef NDEBUG
#define assert_almost_equal(a, b) ((void) 0)
#define assert_log(s) ((void) 0)
#else  /* NDEBUG */
#define assert_almost_equal(a, b) _assert_almost_equal((a), (b), 100, __FILE__, __LINE__)
#define assert_log(s) _assert_log((s), (#s), __FILE__, __LINE__)
#endif  /* NDEBUG */

real KroneckerDelta(int i, int j);
void create_ni(const MATRIX &n, std::vector<real> &ni);
real get_number_of_neighbours(const MATRIX &n, unsigned int i);
real average(const std::vector<real> &ni);
real distance_periodic(Vector p1, Vector p2, real L);
void calculate_n_voronoi(MATRIX &n, const std::vector<Vector> &positions, real L);
void calculate_n_metric(MATRIX &n, const std::vector<Vector> &positions, real L, real r0);
bool sort_func(std::tuple<real,unsigned int> i, std::tuple<real,unsigned int> j);
void calculate_n_topological(MATRIX &n, const std::vector<Vector> &positions, real L, unsigned int neighbours);
void calculate_n_topological_unsymmetrical(MATRIX &n, const std::vector<Vector> &positions, real L, unsigned int neighbours);
void print_n(const MATRIX &n);

real mixing_parameter(real deltat, const MATRIX &n1, const MATRIX &n2, real nc);
real order_parameter(const std::vector<Vector> &st0);
real order_parameter_angle(const std::vector<real> &angles);

real get_average_number_of_neighbours(const MATRIX &n);

typedef std::map<unsigned int, std::vector<Vector>> Posdict;
typedef std::map<unsigned int, std::vector<real>> Angdict;

struct OutOpt {
    enum OutOptEnum : unsigned int {
        Main = 1,
        Extra = 2,
        Neighbours = 4
    };
};

//stores all \pi_{i} of all times specified (and the timestep after)
class Pi
{
    Posdict pis;
    unsigned int mysize = 0;

public:

    Pi (const Angdict &angles);
    Vector Get(unsigned int t, unsigned int i) const;
    const std::vector<Vector> &Get(unsigned int t) const;
    unsigned int size(void) const;
};

//the following C are the correlation functions defined in TABLE 1 (page 5) of paper
real C1s(unsigned int t, Pi &pi);
real Cs(unsigned int t, Pi &pi);
real Gs(unsigned int t, Pi &pi);
real CTs(unsigned int t, Pi &pi, real nc, const std::vector<real> &ni);
real GTs(unsigned int t, Pi &pi, real nc, const std::vector<real> &ni);
real Cint(unsigned int t, Pi &pi, real nc, MATRIX &n);
real CPint(unsigned int t, Pi &pi, real nc, const MATRIX &n);
real Gint(unsigned int t, Pi &pi, real nc, const MATRIX &n);
real CHs(unsigned int t, Pi &pi, real nc, const std::vector<real> &ni);
real CTint(unsigned int t, Pi &pi, real nc, const std::vector<real> &ni, const MATRIX &n);

real J(real nc, real T0, real Omega, real CPint, real CHs, real CTint);
real T(real deltat, real nc, real J, real T0, real Omega, real C1s, real Cs, real Gs, real CTs, real GTs, real Cpint, real Chs, real Ctint);
real T0(real deltat, real Cs, real Gs);
real T0t(real deltat, real Cts, real Gts);
real Omega(real deltat, real Cint, real Gint);
real eta_of_t(real temperature);
real JV_of_J(real deltat, real interaction, real nv);
real DynamicLogLikelihood(real deltat, real N, real nc, real C1s, real Cs, real Gs, real CTs, real GTs, real Cint, real CPint, real Gint, real CHs, real CTint, real J, real T);
real DynamicLogLikelihoodMaximized(real N, real C1s, real Cs, real Gs, real CTs, real GTs, real Cint, real CPint, real Gint, real CHs, real CTint);

void output(Posdict &positions, Angdict &angles, std::vector<unsigned int> &times);

std::string read_data_consecutive(std::string filename, Posdict &positions, Angdict &angles, std::vector<unsigned int> &times);
std::string read_data_single(std::string filename, Posdict &positions, Angdict &angles, std::vector<unsigned int> &times);

std::vector<unsigned int> MakeTimes(unsigned int start, unsigned int stop, unsigned int step);

#endif //CONSTANTS_HEADER
