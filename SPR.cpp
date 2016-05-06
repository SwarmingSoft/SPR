//v 06.05.2016

#include <cmath>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <cassert>
#include <vector>
#include <utility>
#include <random>
#include <tuple>
#include <sstream>

#include "Vector2D.hpp"
#include "Buffer.hpp"

typedef float real;
typedef Vector2DT<real> Vector;

//TODO
/*
-potential for potential dependent friction: shift by U(rmin) for r<rmin (masoud abkenar phd page 36)
-F passive_fraction adaptive
*/

/*KNOWN OR SUSPECTED ISSUES
-since the small box optimization it is possible that pointers for runge-kutta und co can become invalid (due to copying)
*/


//Notes and conventions
/*
Nomenclature and most details as in "Meso-scale turbulence in living fluids"
(Wensink et al, 2012). "eq. [x]" means [x] in this paper, "SI" supportive info.
*/

/* passive mode:
0 : passive static, fraction of particles
1 : passive adaptive, fraction of particles
2 : passive adaptive, fraction of maximum potential
3 : passive activate, fraction of time
*/
const unsigned int passive_mode = 0;


//OLD Constants (mainly for big box optimization)
/*
real maxlength = 0; // (Maximum) Length of particle //OBSOLETE: from reading in sprdata_leghts.txt
//now calculated by GetNumberOfParts()
//unsigned int n = 9; // Parts per particle // calculate as: 1 if length == 1, 3 if 1 < length <= 3, [9/8 * length] if 3 < length where [x] denotes the nearest int
real boxlength = maxlength - lambda + rmin; //length of boxes used for optimization //OBSOLETE: from reading in sprdata_leghts.txt
unsigned int number_of_boxes = 40; // number of boxes in one dimension; change this to change packing fraction
//real L = real(number_of_boxes) * boxlength; // *length //length of simulation square //recalculate after maxlength from sprdata_lengths.txt is set correctly
//4 cases: 0.0,false,false: no passive particles; x,false,false: static passive fraction; x,true,false: passive fraction constant, but dynamically the x with lowest potential; x,true,true: passive fraction non constant and is every rod that has potential below x times the maximum potential
*/


//return maximum element of std::vector
template <typename T> T MaxElement(const std::vector<T> &elements)
{
    T max = elements.front();

    // TODO first loop not needed
    for (auto elm : elements)
    {
        if (elm > max)
        {
            max = elm;
        }
    }
    return max;
}


//calculates the needed number of boxes for uniform grid optimisation (and therefore simulation box lengths) from aimed for packing fraction
//in: aimed for packing fraction, lengths, lambda, rmin
//out: tuple with (number of boxes, real packing fraction)
std::tuple<unsigned int, real> CalcNumberOfBoxes(real packing_fraction, real rmin, const std::vector<real> &lengths, real lambda)
{
    real sum_areas = 0;
    unsigned int num_boxes = 0;
    for (unsigned int i = 0; i < lengths.size(); ++i)
    {
        sum_areas += lengths[i];
    }
    sum_areas *= lambda;
    num_boxes = (unsigned int)(std::round(std::sqrt(sum_areas/packing_fraction)/rmin));
    return std::make_tuple(num_boxes, sum_areas/(rmin*rmin*num_boxes*num_boxes));
}


//average aspect ratio from lengths and lambda
real CalcAvgAspectRatio(const std::vector<real> &lengths, real lambda)
{
    real avg_ar = 0;
    for (unsigned int i = 0; i < lengths.size(); ++i)
    {
        avg_ar += lengths[i];
    }
    avg_ar /= lambda;
    return avg_ar/lengths.size();
}


//class storing all parameters passed by command line arguments
class CMDParameterSPR
{
public:

    std::string in_file, out_file;
    real U0; //potential amplitude //250,455,500,555,625,1250
    real F; //self-motility force // 1
    real f0; // Stokesian friction coefficient // 1
    real lambda; // Screening length // 1
    real dt; //real(0.005) / std::sqrt(N * lambda*lambda / (length * length)); //Euler update steps 0.002*density^(-1/2)
    real rmin; //cutoff potential // 1 * lambda
    real passive_frac; //fraction of particles with no self-motility force F (=passive particles) // 0.5

    real t_sim_end; //time after which the simulation terminates; at the moment cannot be < 20 // 100
    unsigned int N; // Number of particles // 100
    real aimed_for_packing_fraction;
    unsigned int seed; //RNGenerator seed // 477

    //----

    unsigned int number_of_boxes;
    real packing_fraction;
    real L; //*length //length of simulation square //recalculate after maxlength from sprdata_lengths.txt is set correctly
    real maxlength;
    real avg_aspect_ratio;

    CMDParameterSPR(int argc, char** argv)
    {
        if (argc == 1)
        {
            std::cerr << "infile[lengths](string)\toutputfile(string)\tU0(real)\tF(real)\tf0(real)\tlambda(real)\tdt(real)\trmin(real)\tpassive_frac(real)\tt_sim_end(real)\tN(+0int)\tpacking_fraction(real)\tseed(+-int)" << std::endl;
            throw std::runtime_error("Parameters missing");
        }
        assert(argc == 1 + 13);

        in_file = std::string(argv[1]) + ".txt";
        out_file = std::string(argv[2]) + ".txt";

        U0 = std::atof(argv[3]);
        F = std::atof(argv[4]);
        f0 = std::atof(argv[5]);
        lambda = std::atof(argv[6]);
        dt = std::atof(argv[7]);
        rmin = std::atof(argv[8]);
        passive_frac = std::atof(argv[9]);

        t_sim_end = std::atof(argv[10]);
        N = std::atoi(argv[11]);
        aimed_for_packing_fraction = std::atof(argv[12]);

        seed = std::atoi(argv[13]);

        assert(U0 > 0 && F >= 0 && f0 > 0 && lambda > 0 && dt > 0 && N > 20);
    }

    void SetRemainingParameters(const std::vector<real> &lengths)
    {
        assert(lengths.size() == N); // ensure number of lengths read equals number of rods
        auto ret = CalcNumberOfBoxes(aimed_for_packing_fraction, rmin, lengths, lambda);
        number_of_boxes = std::get<0>(ret);
        packing_fraction = std::get<1>(ret);
        L = number_of_boxes * rmin;
        maxlength = MaxElement(lengths);
        avg_aspect_ratio = CalcAvgAspectRatio(lengths, lambda);

        std::cout << "Aimed for packing fraction: " << aimed_for_packing_fraction << ", Real packing fraction: " << packing_fraction << std::endl;
    }

    CMDParameterSPR(std::string in_file, std::string out_file, real U0, real F, real f0, real lambda, real dt, real rmin, real passive_frac, real t_sim_end, unsigned int N, real aimed_for_packing_fraction, unsigned int seed) : in_file(in_file), out_file(out_file), U0(U0), F(F), f0(f0), lambda(lambda), dt(dt), rmin(rmin),
        passive_frac(passive_frac), t_sim_end(t_sim_end), N(N), aimed_for_packing_fraction(aimed_for_packing_fraction), seed(seed)
    {

    }

    std::string header(void)
    {
        std::stringstream ss;

        ss << "In-file: " << in_file << " N: " << N << " L: " << L << " lambda: " << lambda << " maxlength: " << maxlength << " avg_aspect_ratio: " << avg_aspect_ratio << " number_of_boxes: " << number_of_boxes << " rmin: " << rmin << " aimed_for_packing_fraction: " << aimed_for_packing_fraction << " packing_fraction: " << packing_fraction << " U0: " << U0 << " F: " << F << " f0: " << f0 << " passive_frac: " << passive_frac << " passive_mode: " << passive_mode << " dt: " << dt << " SEED " << seed;
        return ss.str();
    }
};


//n (number of segments per rod) from aspect ratio
unsigned int GetNumberOfParts(real aspect_ratio)
{
    assert(aspect_ratio >= real(1));
    if (aspect_ratio == real(1)) //real vs int!
    {
        return 1;
    }
    else if (aspect_ratio <= real(3)) //care: a < 1 falls wrongly into this, should be exception/error
    {
        return 3;
    }
    else
    {
        return std::nearbyint(9*aspect_ratio/real(8));
    }
}

// START COPIED FROM CONSTANTS.HPP
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

// END COPIED


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


//simple 2x2 matrix with multiplication
template <typename T> class Matrix2DT
{
public:
    T a, b, c, d;

    Matrix2DT(void) noexcept
    {

    }

    Matrix2DT(T a, T b, T c, T d) noexcept : a(a), b(b), c(c), d(d)
    {

    }

    /*
    (a b) (x)
    (c d) (y)
    */
    Vector2DT<T> operator*(const Vector2DT<T> &v) noexcept
    {
        return Vector2DT<T>(a*v.x + b*v.y, c*v.x + d*v.y);
    }
};

typedef Matrix2DT<real> Matrix;

class Segment;

#define INIT_NONE 0
#define INIT_BASIC 3475831
#define INIT_FULL 93848659

//single rod, stores all important information (position r, heading u, ...)
//override + and * operator: adds/multiplicates r and u,
//asserts all other field variables are equal
class Rod
{
//private:
public:
    //global
    real lambda;

#ifndef NDEBUG
    unsigned int init = INIT_NONE;
#endif

public:
    //local
    Vector r, u; //center of mass and unit vector of direction
    real length;
    real F, L;
    unsigned int n; // number of segments
    std::vector<real> ls;
    std::vector<Segment> segments;

    unsigned int ticker = 0;
    bool active;

    real a; // aspect ratio
    real fpara, fperp, fr;
    real finvtmp, ftmp0, ftmp2, fdiff;

    int boxx, boxy; // position in optimization box

    /*Rod(const Rod &);
    Rod &operator=(const Rod &);
    Rod(const Rod &&);*/

    Rod(void) noexcept;
    Rod(Vector r, Vector u) noexcept;
    Rod(real lambda, Vector r, Vector u, real length, real F, real f0, real L);

    void PrecomputeLs(void);
    void PrecomputeSegments(void);
    void UpdateSegments(void);
    Matrix Getfinv(void);
    real l(unsigned int i) const;
    Rod operator+(const Rod &rhs) const;
    Rod &operator+=(const Rod &rhs);
    Rod operator*(real dt) const;
    Rod operator/(real dt) const;
    Rod &operator*=(real dt);
    Rod &operator/=(real dt);
};


//single yukawa potential (or similar) point particle -> disk
//store segment identity in rod, identity in box (for optimisation) grid and position
class Segment
{
public:
    Rod *rod; // rod the segment belongs to
    unsigned int i;
    real l; // position relativ to center of rod along its axis

    int boxx, boxy; // position in optimization box

    Vector pos; // absolut position of the segment

    Segment(void);
    Segment(Rod *rod, unsigned int i);
    void UpdatePos(void);
};

/*
//copy constructor
Rod::Rod(const Rod &rhs)
{
    std::cerr << "Rod::Rod() copy constructor" << std::endl;
    lambda = rhs.lambda;
    r = rhs.r; u = rhs.u;
    length = rhs.length;
    F = rhs.F;
    L = rhs.L;
    n = rhs.n;
    ls = rhs.ls;
    segments = rhs.segments;
    //fix pointers

    ticker = rhs.ticker;
    active = rhs.active;

    a = rhs.a;
    fpara = rhs.fpara; fperp = rhs.fperp; fr = rhs.fr;
    finvtmp = rhs.finvtmp; ftmp0 = rhs.ftmp0; ftmp2 = rhs.ftmp2; fdiff = rhs.fdiff;

    boxx = rhs.boxx;
    boxy = rhs.boxy;

    #ifndef NDEBUG
    init = rhs.init;
    #endif
}

//assignment constructor
Rod &Rod::operator=(const Rod &rhs)
{
    std::cerr << "Rod::operator=" << std::endl;
    lambda = rhs.lambda;
    r = rhs.r; u = rhs.u;
    length = rhs.length;
    F = rhs.F;
    L = rhs.L;
    n = rhs.n;
    ls = rhs.ls;
    segments = rhs.segments;
    //fix pointers

    ticker = rhs.ticker;
    active = rhs.active;

    a = rhs.a;
    fpara = rhs.fpara; fperp = rhs.fperp; fr = rhs.fr;
    finvtmp = rhs.finvtmp; ftmp0 = rhs.ftmp0; ftmp2 = rhs.ftmp2; fdiff = rhs.fdiff;

    boxx = rhs.boxx;
    boxy = rhs.boxy;

    #ifndef NDEBUG
    init = rhs.init;
    #endif

    return *this;
}

//move constructor
Rod::Rod(const Rod &&rhs)
{
    std::cerr << "Rod::Rod() move constructor" << std::endl;
}
*/

Rod::Rod(void) noexcept
{

}

Rod::Rod(Vector r, Vector u) noexcept : r(r), u(u)
{
    //uninitialized Rod, except r and u
    #ifndef NDEBUG
    init = INIT_BASIC;
    #endif
}

Rod::Rod(real lambda, Vector r, Vector u, real length, real F, real f0, real L) : lambda(lambda), r(r), u(u), length(length), F(F), L(L), n(GetNumberOfParts(length/lambda)), ls(n), segments(n), active(F > 0.001) /*std::nextafter(0)*/ //n(GetNumberOfParts(length))
{
    assert(n > 0 && length > 0 && lambda > 0 && f0 != 0);

    a = length/lambda; //aspect ratio
    assert(a >= 1);

    //standard expressions for rod-like macromolecules: "Comparison of theories for the translational and rotational diffusion coefficients of rod-like macromolecules"
    fpara = 2*M_PI/(std::log(a) - 0.207 + 0.980/a - 0.133/(a*a));
    fperp = 4*M_PI/(std::log(a) + 0.839 + 0.185/a + 0.233/(a*a));
    fr = M_PI*a*a/3/(std::log(a) - 0.662 + 0.917/a - 0.050/(a*a));

    // precompute some parameters
    ftmp0 = f0*fpara;
    finvtmp = ftmp0*fperp;
    ftmp2 = 1/(f0*fr);
    fdiff = fperp - fpara;

    PrecomputeLs();
    #ifndef NDEBUG
    init = INIT_FULL;
    #endif
    PrecomputeSegments();
}

// save some memory by deleting this and only use PrecomputeSegments
void Rod::PrecomputeLs(void)
{
    assert(ls.size() == n);
    if (n != 1)
    {
        for (unsigned int i = 0; i < n; ++i)
        {
            ls[i] = -(length - lambda)/2. + real(i) * (length-lambda)/(n-1);
        }
    }
    else
    {
        ls[0] = 0;
    }
}

void Rod::PrecomputeSegments(void)
{
    assert(init == INIT_FULL);
    assert(segments.size() == n);
    for (unsigned int i = 0; i < n; ++i)
    {
        segments[i] = Segment(this, i);
    }
}

// Updates positions of all segments within the rod
void Rod::UpdateSegments(void)
{
    assert(init == INIT_FULL);
    for (Segment &seg : segments)
    {
        assert(seg.rod == this);
        seg.UpdatePos();
    }
}


//inverse translational friction tensor f_T^{-1}, eq. [6] SI
inline Matrix Rod::Getfinv(void)
{
    assert(init == INIT_FULL);
    assert(finvtmp != 0);
    real tmp = fdiff/finvtmp * u.x*u.y;
    real tmp2 = fdiff * u.y * u.y;

    return Matrix(
        (fperp - tmp2)/finvtmp,
        tmp,
        tmp,
        (fpara + tmp2)/finvtmp
    );
}

// getter for l
real Rod::l(unsigned int i) const
{
    assert(init == INIT_FULL);
    assert(ls.size() == n);
    return ls[i];
}

// add two rods while preserving periodic bounaries
// should only be used to add the same two rods at different times
Rod Rod::operator+(const Rod &rhs) const
{
    assert(init == INIT_FULL);
    assert(rhs.init == INIT_BASIC);

    Rod ret(*this); //copy self
    ret.r.x = mod(r.x + rhs.r.x, L);
    ret.r.y = mod(r.y + rhs.r.y, L);
    assert(ret.r.x < L && ret.r.y < L);
    ret.u += rhs.u;
    return ret;

    //return Rod(lambda, Vector(mod(r.x + rhs.r.x, L), mod(r.y + rhs.r.y, L)), u+rhs.u, length, F);
    //return Rod(lambda, Vector(r.x + rhs.r.x, r.y + rhs.r.y), u+rhs.u, length, F); //non-periodic
}

// appends a rod while preserving periodic boundaries
// should only be used to append the same two rods at a different time
Rod &Rod::operator+=(const Rod &rhs)
{
    assert(init == INIT_FULL);
    assert(rhs.init == INIT_BASIC);

    r += rhs.r;

    //periodic boundary conditions
    r.x = mod(r.x, L);
    r.y = mod(r.y, L);

    assert(r.x < L && r.y < L);

    u += rhs.u;
    assert(u.Norm() != 0.); //assures length of u is always 1
    u.Normalize();
    return *this;
}

Rod Rod::operator*(real dt) const
{
    assert(init == INIT_BASIC);
    Rod ret(*this);
    ret.r *= dt;
    ret.u *= dt;
    return ret;
    //return Rod(lambda, r*dt, u*dt, length, F);
}

Rod Rod::operator/(real dt) const
{
    assert(init == INIT_BASIC);
    Rod ret(*this);
    ret.r /= dt;
    ret.u /= dt;
    return ret;
    //return Rod(lambda, r/dt, u/dt, length, F);
}

Rod &Rod::operator*=(real dt)
{
    assert(init == INIT_FULL);
    r *= dt;
    u *= dt;
    return *this;
}

Rod &Rod::operator/=(real dt)
{
    assert(init == INIT_FULL);
    r /= dt;
    u /= dt;
    return *this;
}

Segment::Segment(void) { }

// Construct segment with parent rod and index within rod
Segment::Segment(Rod *rod, unsigned int i) : rod(rod), i(i), l(rod->l(i))
{

}

// update position based on rod
void Segment::UpdatePos(void)
{
    assert(rod != nullptr);
    assert(rod->init == INIT_FULL);
    assert(l == rod->l(i));
    //pos = rod->r + rod->l(i)*rod->u;
    pos = rod->r + l*rod->u;
}

typedef std::vector<Rod> Rods;
typedef GenericBuffer_dyn<std::vector<Rod*>, MemoryConst::Dynamic, MemoryConst::Dynamic> RodsMatrix;
typedef GenericBuffer_dyn<std::vector<Segment*>, MemoryConst::Dynamic, MemoryConst::Dynamic> SegmentsMatrix;


//override operators for multiplication and addition for vector rods of Rod
//objects, so that all the responsibility is now in Rods class
Rods &operator+=(Rods &lhs, const Rods &rhs)
{
    assert(lhs.size() == rhs.size());
    for (unsigned int i = 0; i < lhs.size(); ++i)
    {
        assert(lhs[i].init == INIT_FULL);
        assert(rhs[i].init == INIT_BASIC);
        lhs[i] += rhs[i]; //just adds r and u, takes all other parameters from lhs
    }
    return lhs;
}

/*Rods operator+(const Rods &lhs, const Rods &rhs)
{
    assert(lhs.size() == rhs.size());

    Rods ret(lhs.size());
    for (unsigned int i = 0; i < lhs.size(); ++i)
    {
        assert(lhs[i].init == INIT_FULL);
        assert(rhs[i].init == INIT_BASIC);
        ret[i] = lhs[i] + rhs[i]; //just adds r and u, takes all other parameters from lhs
    }
    return ret;
}*/

/*Rods operatorEuler*(real dt, Rods rods) // should be the same as normal *operator
{
    for (unsigned int i = 0; i < rods.size(); ++i)
    {
        rods[i] *= dt;
    }
    return rods;
}*/

Rods operator*(real dt, const Rods &rods)
{
    Rods ret(rods.size());
    for (unsigned int i = 0; i < rods.size(); ++i)
    {
        assert(rods[i].init == INIT_BASIC);
        ret[i] = rods[i] * dt;
    }
    return ret;
}

Rods operator*(const Rods &rods, real dt)
{
    return dt*rods;
}

Rods operator/(const Rods &rods, real c)
{
    Rods ret(rods.size());
    for (unsigned int i = 0; i < rods.size(); ++i)
    {
        assert(rods[i].init == INIT_BASIC);
        ret[i] = rods[i] / c;
    }
    return ret;
}

// normalize all rods in list
void Normalize(Rods &rods)
{
    for (auto &rod : rods)
    {
        rod.u.Normalize();
    }
}

//Yukawa Potential non-periodic
real YukawaPotential(Rod a, Rod b, real U0, real lambda)
{
    assert(lambda > 0 && U0 > 0);
    real sum = 0, rij;

    for (unsigned int i = 0; i < a.n; ++i)
    {
        for (unsigned int j = 0; j < b.n; ++j)
        {
            rij = ((a.r-b.r)+(a.l(i)*a.u+b.l(j)*b.u)).Norm();
            sum += std::exp(-rij/lambda)/rij;
        }
    }

    return U0/a.n/b.n * sum;
}

//Lennard-Jones Potential non-periodic
real LennardJonesPotential(Rod a, Rod b, real epsilon, real sigma)
{
    assert(sigma > 0 && epsilon > 0);
    real sum = 0, rij;

    for (unsigned int i = 0; i < a.n; ++i)
    {
        for (unsigned int j = 0; j < b.n; ++j)
        {
            rij = ((a.r-b.r)+(a.l(i)*a.u+b.l(j)*b.u)).Norm();
            sum += (std::pow(sigma/rij, 12) - std::pow(sigma/rij, 6));
        }
    }

    return 4*epsilon/a.n/b.n * sum;
}

//sum of Yukawa potentials of all particles i to all other particles j
//excluding self-self potentials
real U(const Rods &rods, real U0, real lambda)
{
    assert(lambda > 0.);
    real sum = 0.;
    for (unsigned int i = 0; i < rods.size(); ++i)
    {
        for (unsigned int j = 0; j < rods.size(); ++j)
        {
            if (i != j) // exclude self-self potential
            {
                sum += YukawaPotential(rods[i], rods[j], U0, lambda);
            }

        }
    }
    return sum/real(2.);
}

//partial derivatives of U (sum of Yukawa potentials) wrt r and u as a tuple, periodic
//use get<0>, get<1> on tuple
//TODO optimise using symmetry in i and j ??? -> should not be possible???
std::tuple<Vector,Vector> DerivativeRaUaYukawaPotential(Rod a, Rod b, real U0, real lambda, real L, real rmin)
{
    assert(lambda > 0 && U0 > 0);
    Vector r(0, 0), u(0, 0), t;

    real norm, tmp;

    for (unsigned int i = 0; i < a.n; ++i)
    {
        for (unsigned int j = 0; j < b.n; ++j)
        {
            //assert(a.l(i) == b.l(i) && a.l(j) == b.l(j));
            //t = a.r - b.r + a.l(i)*a.u - b.l(j)*b.u;
            //t = distance_periodic_vector(a.r, b.r, L) + a.l(i)*a.u - b.l(j)*b.u;
            t = distance_periodic_vector(a.r + a.l(i)*a.u, b.r + b.l(j)*b.u, L); // (a.r + a.l(i)*a.u, b.r + a.l(j)*b.u, L)

            norm = t.Norm();
            assert(norm != 0);
            if (norm <= rmin)
            {
                tmp = std::exp(-norm/lambda)*(norm + lambda)/(norm*norm*norm*lambda);

                r -= t*tmp;
                u -= t*(a.l(i)*tmp);
            }
            //else {}
        }
    }

    return std::make_tuple(U0/a.n/b.n * r, U0/a.n/b.n * u);
}

//partial derivatives of sum of Lennard-Jones Potential wrt r and u as tuples, periodic
std::tuple<Vector,Vector> DerivativeRaUaLennardJonesPotential(Rod a, Rod b, real epsilon, real sigma, real L, real rmin)
{
    assert(epsilon > 0 && sigma > 0);
    Vector r(0, 0), u(0, 0), t;

    real norm, tmp;
    real sigma6 = std::pow(sigma, 6);

    for (unsigned int i = 0; i < a.n; ++i)
    {
        for (unsigned int j = 0; j < b.n; ++j)
        {
            t = distance_periodic_vector(a.r + a.l(i)*a.u, b.r + b.l(j)*b.u, L); // (a.r + a.l(i)*a.u, b.r + a.l(j)*b.u, L)

            norm = t.Norm();
            assert(norm != 0);

            if (norm <= rmin)
            {
                tmp = (std::pow(norm, 6) - 2*sigma6)/std::pow(norm, 14);

                r += t*tmp;
                u += t*(a.l(i)*tmp);
            }
            //else {}
        }
    }

    return std::make_tuple(24*epsilon*sigma6*a.n/b.n * r, 24*epsilon*sigma6*a.n/b.n * u);
}


//partial derivative of U_alpha,beta_i,j (Yukawa potential for rod alpha and beta only for segment i (of alpha) and j (of beta)) wrt r and u as a tuple
//rmin is cutoff, L width of simulation box - needed for periodicity
//use get<0>, get<1> on tuple
std::tuple<Vector,Vector> DerivativeRaUaYukawaPotentialSegment(Segment a, Segment b, real lambda, real L, real rmin)
{
    assert(lambda > 0);
    Vector r(0, 0), u(0, 0), t;

    real norm, tmp;

    t = distance_periodic_vector(a.pos, b.pos, L);

    norm = t.Norm();
    assert(norm != 0);
    if (norm <= rmin)
    {
        tmp = std::exp(-norm/lambda)*(norm + lambda)/(norm*norm*norm*lambda);

        r -= t*tmp;
        u -= t*(a.l*tmp);
    }

    return std::make_tuple(r, u);
}

//partial derivatives of single segment Lennard-Jones Potential wrt r and u, periodic
//rmin is cutoff, L width of simulation box - needed for periodicity
std::tuple<Vector,Vector> DerivativeRaUaLennardJonesPotentialSegment(Segment a, Segment b, real sigma6, real L, real rmin)
{
    assert(sigma6 > 0);
    Vector r(0, 0), u(0, 0), t;

    real norm, tmp;

    t = distance_periodic_vector(a.pos, b.pos, L);

    norm = t.Norm();
    assert(norm != 0);
    if (norm <= rmin)
    {
        tmp = (std::pow(norm, 6) - 2*sigma6)/std::pow(norm, 14);

        r += t*tmp;
        u += t*(a.l*tmp);
    }

    return std::make_tuple(r, u);
}

std::tuple<Vector,Vector,real> DerivativeRaUaYukawaPotentialSegmentAdaptive(Segment a, Segment b, real lambda, real L, real rmin)
{
    assert(lambda > 0);
    Vector r(0, 0), u(0, 0), t;

    real norm, tmp, pot = 0;

    t = distance_periodic_vector(a.pos, b.pos, L);

    norm = t.Norm();
    assert(norm != 0);
    if (norm <= rmin)
    {
        //tmp = std::exp(-norm/lambda)*(norm + lambda)/(norm*norm*norm*lambda);
        pot = std::exp(-norm/lambda)/norm;
        tmp = pot/(norm*norm*lambda)*(norm + lambda);

        r -= t*tmp;
        u -= t*(a.l*tmp);
    }

    return std::make_tuple(r, u, pot);
}


//gradients of U (sum of Yukawa potentials) wrt r and u as a tuple
//access .first and .second (now get<0> get<1>)
std::tuple<Vector,Vector> GradientRaUa(const Rods &rods, unsigned int a, real U0, real lambda, real L, real rmin)
{
    Vector rsum(0, 0);
    Vector usum(0, 0);
    real searchlen;
    for (unsigned int b = 0; b < rods.size(); ++b)
    {
        searchlen = rods[a].length+2*rods[a].lambda; //finetune
        //if (a != b)
        if (a != b && rods[a].r.GetDistance(rods[b].r) < searchlen) //test cutoff parameter
        {
            auto derivative = DerivativeRaUaYukawaPotential(rods[a], rods[b], U0, lambda, L, rmin);
            rsum += std::get<0>(derivative);
            usum += std::get<1>(derivative);
        }
    }
    return std::make_tuple(rsum/2, usum/2);
}

//same as GradientRaUa but with 9 cells instead of global (rod boxes)
std::tuple<Vector,Vector> GradientRaUaBoxed(const RodsMatrix &matrix, const Rod &a, real U0, real lambda, real L, real rmin)
{
    Vector rsum(0, 0);
    Vector usum(0, 0);

    //int counter = 0;
    for (int x = a.boxx - 1; x < a.boxx + 2; ++x)
    {
        for (int y = a.boxy - 1; y < a.boxy + 2; ++y)
        {
            const std::vector<Rod*> *rods = matrix.PeriodicGet(x, y);
            //const std::vector<Rod*> *rods = matrix.PeriodicGet(x, y); //additionally returns x,y
            for (unsigned int p = 0; p < rods->size(); ++p)
            {
                const Rod &b = *rods->at(p);
                //if (x == a.boxx && y == a.boxy && a.boxpos == b.boxpos) {}
                //if (mod(x, matric.GetWidth()) == a.boxx && mod(y,yboxes) == a.boxy && a.boxpos == b.boxpos) {}
                if (&a != &b)
                {
                    auto derivative = DerivativeRaUaYukawaPotential(a, b, U0, lambda, L, rmin);
                    rsum += std::get<0>(derivative);
                    usum += std::get<1>(derivative);
                }
                //else {std::cout << x << " " << y << std::endl; counter += 1;}
            }
        }
    }
    //std::cout << "counter: " << counter << std::endl;
    return std::make_tuple(rsum/2, usum/2);
}

//same as GradientRaUa but with 9 cells instead of global and for segment boxes
std::tuple<Vector,Vector> GradientRaUaSegment(const SegmentsMatrix &matrix, const Segment &a, real lambda, real L, real rmin)
{
    Vector rsum(0, 0);
    Vector usum(0, 0);

    //int counter = 0;
    for (int x = a.boxx - 1; x < a.boxx + 2; ++x)
    {
        for (int y = a.boxy - 1; y < a.boxy + 2; ++y)
        {
            const std::vector<Segment*> *segs = matrix.PeriodicGet(x, y);
            for (unsigned int p = 0; p < segs->size(); ++p)
            {
                const Segment &b = *segs->at(p);
                //if (x == a.boxx && y == a.boxy && a.boxpos == b.boxpos) {}
                //if (mod(x, matric.GetWidth()) == a.boxx && mod(y,yboxes) == a.boxy && a.boxpos == b.boxpos) {}
                if (a.rod != b.rod) // same segment: &a != &b
                {
                    auto derivative = DerivativeRaUaYukawaPotentialSegment(a, b, lambda, L, rmin); // * U0 / a.n / b.n outsourced
                    rsum += std::get<0>(derivative) / b.rod->n; // /b.n from before; U0 / a.n in SPR operator()
                    usum += std::get<1>(derivative) / b.rod->n; // /b.n from before; U0 / a.n in SPR operator()
                }
                //else {std::cout << x << " " << y << std::endl; counter += 1;}
            }
        }
    }
    //std::cout << "counter: " << counter << std::endl;
    return std::make_tuple(rsum/2, usum/2);
}

//same as GradientRaUa but with 9 cells instead of global and for segment boxes; adaptive returns tuple instead of tuple (with the potential as 3rd)
std::tuple<Vector,Vector,real> GradientRaUaSegmentAdaptive(const SegmentsMatrix &matrix, const Segment &a, real lambda, real L, real rmin)
{
    Vector rsum(0, 0);
    Vector usum(0, 0);
    real potsum = 0.;

    //int counter = 0;
    for (int x = a.boxx - 1; x < a.boxx + 2; ++x)
    {
        for (int y = a.boxy - 1; y < a.boxy + 2; ++y)
        {
            const std::vector<Segment*> *segs = matrix.PeriodicGet(x, y);
            for (unsigned int p = 0; p < segs->size(); ++p)
            {
                const Segment &b = *segs->at(p);
                //if (x == a.boxx && y == a.boxy && a.boxpos == b.boxpos) {}
                //if (mod(x, matric.GetWidth()) == a.boxx && mod(y,yboxes) == a.boxy && a.boxpos == b.boxpos) {}
                if (a.rod != b.rod) // same segment: &a != &b
                {
                    auto derivative = DerivativeRaUaYukawaPotentialSegmentAdaptive(a, b, lambda, L, rmin); // * U0 / a.n / b.n outsourced
                    rsum += std::get<0>(derivative) / b.rod->n; // /b.n from before; U0 / a.n in SPR operator()
                    usum += std::get<1>(derivative) / b.rod->n; // /b.n from before; U0 / a.n in SPR operator()
                    potsum += std::get<2>(derivative) / b.rod->n;
                }
                //else {std::cout << x << " " << y << std::endl; counter += 1;}
            }
        }
    }
    //std::cout << "counter: " << counter << std::endl;
    return std::make_tuple(rsum/2, usum/2, potsum/2); // *2 OR /2 OR *1 ?????
}


void operator+=(std::tuple<Vector, Vector> &a, const std::tuple<Vector, Vector> &b)
{
    std::get<0>(a) += std::get<0>(b);
    std::get<1>(a) += std::get<1>(b);
}

void operator*=(std::tuple<Vector, Vector> &a, real s)
{
    std::get<0>(a) *= s;
    std::get<1>(a) *= s;
}

void operator+=(std::tuple<Vector, Vector, real> &a, const std::tuple<Vector, Vector, real> &b)
{
    std::get<0>(a) += std::get<0>(b);
    std::get<1>(a) += std::get<1>(b);
    std::get<2>(a) += std::get<2>(b);
}

void operator*=(std::tuple<Vector, Vector, real> &a, real s)
{
    std::get<0>(a) *= s;
    std::get<1>(a) *= s;
    std::get<2>(a) *= s;
}

//self propelled rod: stores global constants (such as friction coefficients)
//does all the calculation, eg division into boxes and gradients in operator()
class SPR
{
    unsigned int N;
    real U0, F, L, rmin, lambda, passive_frac, dt;
    Rods ret;

    Matrix finv;
    unsigned int potentials_frac;

    public:

    /*SPR(unsigned int N, real U0, real lambda) : N(N), U0(U0), lambda(lambda), ret(N)
    {
        unsigned int numboxes = std::lround(L / boxlength); //what if rounding down and particle outside of all boxes

        std::cout << "Boxes: " << numboxes << "**2" << std::endl;

        boxes = RodsMatrix(numboxes, numboxes);

        for (unsigned int i = 0; i < boxes.GetWidth(); ++i)
        {
            for (unsigned int j = 0; j < boxes.GetHeight(); ++j)
            {
                boxes.Set(i, j, std::vector<Rod*>());
            }
        }
    }*/

    SPR(const CMDParameterSPR &args) : N(args.N), U0(args.U0), F(args.F), L(args.L), rmin(args.rmin), lambda(args.lambda), passive_frac(args.passive_frac), dt(args.dt), ret(args.N)
    {
        //unsigned int numboxes = std::lround(args.L / args.rmin); //what if rounding down and particle outside of all boxes
        unsigned int numboxes = args.number_of_boxes;

        std::cout << "Boxes: " << numboxes << "**2" << std::endl;

        boxes = SegmentsMatrix(numboxes, numboxes);

        for (unsigned int i = 0; i < boxes.GetWidth(); ++i)
        {
            for (unsigned int j = 0; j < boxes.GetHeight(); ++j)
            {
                boxes.Set(i, j, std::vector<Segment*>());
            }
        }

        potentials_frac = std::lround(args.passive_frac*N);
    }

    //obsolete: unoptimised (without boxes)
    /*Rods UNOPT_operator()(const Rods &rods, real)
    {
        assert(N == rods.size());
        Vector r, u;

        //gradients wrt r and u
        for (unsigned int a = 0; a < rods.size(); ++a)
        {
            finv = Buildfinv(rods[a].u);
            auto gradient = GradientRaUa(rods, a, U0, lambda);
            r = rods[a].F/ftmp0*rods[a].u - finv*std::get<0>(gradient);
            u = -std::get<1>(gradient)*ftmp2;
            u.Normalize();

            //std::cout << "dr: " << std::get<0>(gradient).x << " " << std::get<0>(gradient).y << " du: " << std::get<1>(gradient).x << " " << std::get<1>(gradient).y << " ";
            ret[a] = Rod(rods[a].lambda, r, u, rods[a].length);
        }
        return ret;
    }*/

    //does all the work, calculates gradients with respect to r and u and returns
    //them stored in a "Rod" object (as r and u)
    //confer eq. [4] and [5] in SI
    //MUSTN'T change data stored in Rods rods other than the index in box ints
    //Rod Boxes
    /*RodsMatrix boxes;
    Rods operator()(Rods &rods, real)
    {
        assert(N == rods.size());
        unsigned int i, j;

        //reset particles per box
        for (i = 0; i < boxes.GetWidth(); ++i)
        {
            for (j = 0; j < boxes.GetHeight(); ++j)
            {
                boxes.Get(i, j)->clear();
            }
        }

        int x, y; //same type as boxx, boxy
        unsigned int ri;

        //save particle's index (position in vector in boxes.Get(i, j))
        //in rods field variables boxx, boxy, boxpos
        for (ri = 0; ri < rods.size(); ++ri)
        {
            assert(rods[ri].r.x >= 0 && rods[ri].r.y >= 0);
            assert(rods[ri].r.x < L && rods[ri].r.y < L);
            assert_almost_equal(rods[ri].u.Norm(), real(1)); // needed?

            x = std::floor(rods[ri].r.x / boxlength);
            y = std::floor(rods[ri].r.y / boxlength);
            //std::cout << "395 "  << a << " x: " << x << " y: " << y << " " << rods[a].r.x << " " << rods[a].r.y << " " << length << std::endl;
            rods[ri].boxx = x;
            rods[ri].boxy = y;
            //rods[ri].boxpos = boxes.Get(x, y)->size();
            boxes.Get(x, y)->push_back(&rods[ri]);
        }

        Vector r, u;

        //gradients for all rods wrt r and u
        for (ri = 0; ri < rods.size(); ++ri)
        {
            finv = rods[ri].Getfinv(); //Buildfinv(rods[a]);
            auto gradient = GradientRaUaBoxed(boxes, rods[ri], U0, lambda);
            r = rods[ri].F/rods[ri].ftmp0*rods[ri].u - finv*std::get<0>(gradient);
            assert(!std::isnan(r.x) && !std::isnan(r.y));
            u = -std::get<1>(gradient)*rods[ri].ftmp2;
            assert(!std::isnan(u.x) && !std::isnan(u.y));

            //std::cout << "dr: " << std::get<0>(gradient).x << " " << std::get<0>(gradient).y << " du: " << std::get<1>(gradient).x << " " << std::get<1>(gradient).y << " ";

            ret[ri] = Rod(r, u); // return only partially initialized rod, will only be used on rhs of addition
        }
        return ret;
    }*/

    //Segment Boxes
    SegmentsMatrix boxes;
    Rods operator()(Rods &rods, real t)
    {
        assert(N == rods.size());

        unsigned int i, j;
        //reset segments in boxes
        for (i = 0; i < boxes.GetWidth(); ++i)
        {
            for (j = 0; j < boxes.GetHeight(); ++j)
            {
                boxes.Get(i, j)->clear();
            }
        }

        int x, y; //same type as boxx, boxy
        unsigned int ri, pi;

        //add segment ptr to boxes and save box coords to segment
        for (ri = 0; ri < rods.size(); ++ri)
        {
            // ensure that all rods lie in simulation box
            assert(rods[ri].r.x >= 0 && rods[ri].r.y >= 0);
            assert(rods[ri].r.x < L && rods[ri].r.y < L);
            assert_almost_equal(rods[ri].u.Norm(), real(1)); // needed?
            rods[ri].UpdateSegments(); // Update positions of all segments
            for (Segment &seg : rods[ri].segments)
            {
                x = std::floor(seg.pos.x / rmin);
                y = std::floor(seg.pos.y / rmin);
                seg.boxx = x;
                seg.boxy = y;
                //std::cout << "box: "<< boxes.GetWidth() << "x" << boxes.GetHeight() << " x " << x << " y " << y  << std::endl;
                //boxes.Get(x, y)->push_back(&seg);
                boxes.PeriodicGet(x, y)->push_back(&seg);
            }
        }

        Vector r, u;

        assert(N == ret.size());

        //PASSIVE PARTICLES NON ADAPTIVE
        if (passive_mode == 0)
        {
            for (ri = 0; ri < rods.size(); ++ri)
            {
                std::tuple<Vector, Vector> gradient(Vector(0, 0), Vector(0, 0));
                for (Segment &seg : rods[ri].segments)
                {
                    gradient += GradientRaUaSegment(boxes, seg, lambda, L, rmin);
                }
                gradient *= U0 / rods[ri].n; // /b.n is in inner loop
                finv = rods[ri].Getfinv(); //Buildfinv(rods[a]);
                r = rods[ri].F/rods[ri].ftmp0*rods[ri].u - finv*std::get<0>(gradient);
                assert(!std::isnan(r.x) && !std::isnan(r.y));
                u = -std::get<1>(gradient)*rods[ri].ftmp2;
                assert(!std::isnan(u.x) && !std::isnan(u.y));
                ret[ri] = Rod(r, u);
            }
        }

        //PASSIVE PARTICLES ADAPTIVE relative to number of particles
        else if (passive_mode == 1)
        {
            std::vector<std::tuple<unsigned int, real>> potentials(N);
            std::vector<std::tuple<Vector, Vector, real>> gradients(N);

            for (ri = 0; ri < rods.size(); ++ri)
            {
                std::tuple<Vector, Vector, real> gradient(Vector(0, 0), Vector(0, 0), 0);
                for (Segment &seg : rods[ri].segments)
                {
                    gradient += GradientRaUaSegmentAdaptive(boxes, seg, lambda, L, rmin);
                }
                gradient *= U0 / rods[ri].n; // /b.n is in inner loop
                gradients[ri] = gradient; //TODO third element of gradient not needed, could be left out
                potentials[ri] = std::make_tuple(ri, std::get<2>(gradient)); // (index, potential)
            }

            std::sort(potentials.begin(), potentials.end(), [](const std::tuple<unsigned int, real> &a, const std::tuple<unsigned int, real> &b) { return std::get<1>(a) < std::get<1>(b); });

            for (pi = 0; pi < potentials.size(); ++pi)
            {
                ri = std::get<0>(potentials[pi]);

                if (pi < potentials_frac)
                {
                    rods[ri].F = 0; //should be 0
                    rods[ri].active = false;
                    //if (t > 1) {std::cout << pi << " Index: " << ri << " n: " << rods[ri].n << " Potential (low): " << std::get<1>(potentials[pi]) << std::endl;}
                }
                else
                {
                    rods[ri].F = F; //should be F
                    rods[ri].active = true;
                    //if (t > 1) {std::cout << pi << " Index: " << ri << " n: " << rods[ri].n << " Potential (hig): " << std::get<1>(potentials[pi]) << std::endl;}
                }
                finv = rods[ri].Getfinv(); //Buildfinv(rods[a]);
                r = rods[ri].F/rods[ri].ftmp0*rods[ri].u - finv*std::get<0>(gradients[ri]);
                assert(!std::isnan(r.x) && !std::isnan(r.y));
                u = -std::get<1>(gradients[ri])*rods[ri].ftmp2;
                assert(!std::isnan(u.x) && !std::isnan(u.y));
                ret[ri] = Rod(r, u);
            }
        }

        //PASSIVE PARTICLES ADAPTIVE relative to maximum potential
        else if (passive_mode == 2)
        {
            std::vector<std::tuple<Vector, Vector, real>> gradients(N);

            for (ri = 0; ri < rods.size(); ++ri)
            {
                std::tuple<Vector, Vector, real> gradient(Vector(0, 0), Vector(0, 0), 0);
                for (Segment &seg : rods[ri].segments)
                {
                    gradient += GradientRaUaSegmentAdaptive(boxes, seg, lambda, L, rmin);
                }
                gradient *= U0 / rods[ri].n; // /b.n is in inner loop
                gradients[ri] = gradient;
            }

            real max_potential = std::get<2>(gradients.front());

            // TODO first loop iteration not needed
            for (auto elm : gradients)
            {
                if (std::get<2>(elm) > max_potential)
                {
                    max_potential = std::get<2>(elm);
                }
            }

            real potential_frac = passive_frac*max_potential;

            for (ri = 0; ri < rods.size(); ++ri)
            {
                //below threshold
                if (std::get<2>(gradients[ri]) < potential_frac)
                {
                    rods[ri].F = 0;
                    rods[ri].active = false;
                }
                else
                {
                    rods[ri].F = F;
                    rods[ri].active = true;
                }
                finv = rods[ri].Getfinv(); //Buildfinv(rods[a]);
                r = rods[ri].F/rods[ri].ftmp0*rods[ri].u - finv*std::get<0>(gradients[ri]);
                assert(!std::isnan(r.x) && !std::isnan(r.y));
                u = -std::get<1>(gradients[ri])*rods[ri].ftmp2;
                assert(!std::isnan(u.x) && !std::isnan(u.y));
                ret[ri] = Rod(r, u);
            }
        }

        //PASSIVE PARTICLES ADAPTIVE relative to number of particles
        else if (passive_mode == 3)
        {
            for (ri = 0; ri < rods.size(); ++ri)
            {
                std::tuple<Vector, Vector, real> gradient(Vector(0, 0), Vector(0, 0), 0);
                for (Segment &seg : rods[ri].segments)
                {
                    gradient += GradientRaUaSegmentAdaptive(boxes, seg, lambda, L, rmin);
                }
                gradient *= U0 / rods[ri].n; // /b.n is in inner loop
                finv = rods[ri].Getfinv(); //Buildfinv(rods[a]);
                r = rods[ri].F/rods[ri].ftmp0*rods[ri].u - finv*std::get<0>(gradient);
                assert(!std::isnan(r.x) && !std::isnan(r.y));
                u = -std::get<1>(gradient)*rods[ri].ftmp2;
                assert(!std::isnan(u.x) && !std::isnan(u.y));

                bool has_potential = (std::get<2>(gradient) > 0); //0.001
                if (rods[ri].active != has_potential and rods[ri].ticker < passive_frac/dt) // passive_frac == 100
                {
                    rods[ri].ticker += 1;
                }
                else
                {
                    rods[ri].active = has_potential;
                    rods[ri].ticker = 0;
                }
                if (rods[ri].active)
                {
                    rods[ri].F = F;
                }
                else
                {
                    rods[ri].F = 0;
                }

                ret[ri] = Rod(r, u);
            }
        }

        //PASSIVE PARTICLES: not a valid mode
        else
        {
            static_assert(passive_mode >= 0 && passive_mode <= 3, "Unknown passive mode");
        }

        return ret;
    }
};


//fix pointers for the small box segment optimization (Segment.*rod of rod[ri]  has to be &rods[ri])
void FixRodPointers(Rods &rods)
{
    for (unsigned int ri = 0; ri < rods.size(); ++ri)
    {
        assert(rods[ri].n == rods[ri].segments.size());
        for (unsigned int j = 0; j < rods[ri].n; ++j)
        {
            rods[ri].segments[j].rod = &rods[ri];
        }
    }
}

//initialise rods at start of simulation
void InitRods(Rods &rods, const std::vector<real> &lengths, const CMDParameterSPR &args)
{
    //initialise RNGenerator
    std::default_random_engine gen(args.seed);
    std::uniform_real_distribution<real> realdist(0, args.L);
    std::discrete_distribution<int> booldist {0.5,0.5};
    std::discrete_distribution<int> Fdist {args.passive_frac,real(1)-args.passive_frac};

    real x_tmp, y_tmp, customF;
    int failed_placements = 0;

    //initialise rods at uniformly random positions in box with orientation upwards or downwards
    for (unsigned int i = 0; i < rods.size(); ++i)
    {
        if (passive_mode == 3) //passive_frac must be time fraction (therefore > 1 etc)
        {
            customF = args.F;
        }
        else
        {
            customF = Fdist(gen)*args.F;
        }

        tryagain:
        x_tmp = realdist(gen);
        y_tmp = realdist(gen);
        Vector rod_pos_tmp(x_tmp, y_tmp);

        if (x_tmp < args.L/2. - lengths[i]/2. && x_tmp > lengths[i]/2.) // left side
        {
            for (unsigned int j = 0; j < i; ++j)
            {
                auto diff = distance_periodic_vector(rod_pos_tmp, rods[j].r, args.L);
                if (std::abs(diff.x) <= ((lengths[i] + rods[j].length)/2) && std::abs(diff.y) < args.lambda/2)
                {
                    ++failed_placements;
                    goto tryagain;
                }
            }
            rods[i] = Rod(args.lambda, rod_pos_tmp, Vector(booldist(gen)*2-1, 0), lengths[i], customF, args.f0, args.L); //Vector(booldist(gen)*2-1, 0)
        }
        else if (x_tmp > args.L/2. && x_tmp < args.L) // right side
        {
            for (unsigned int j = 0; j < i; ++j)
            {
                auto diff = distance_periodic_vector(rod_pos_tmp, rods[j].r, args.L);
                if (std::abs(diff.y) <= ((lengths[i] + rods[j].length)/2) && std::abs(diff.x) < args.lambda/2)
                {
                    ++failed_placements;
                    goto tryagain;
                }
            }
            rods[i] = Rod(args.lambda, rod_pos_tmp, Vector(0, booldist(gen)*2-1), lengths[i], customF, args.f0, args.L); //Vector(booldist(gen)*2-1, 0)
        }
        else //nothing; clear zone; assures non-crossing rods
        {
            goto tryagain;
        }
        //rods[i].UpdateSegments(); //test
    }
    std::cout << "Failed placements in random InitRods: " << failed_placements << std::endl;
    FixRodPointers(rods);
    /*rods[0] = Rod(lambda, Vector(0, 3*length), Vector(1, 0), n, length);
    rods[1] = Rod(lambda, Vector(3*length-4*lambda, 3*length+3*length*lambda), Vector(0, -1), length);
    rods[2] = Rod(lambda, Vector(3*length-3*lambda, 3*length+3*length*lambda), Vector(0, -1), length);
    rods[3] = Rod(lambda, Vector(3*length-2*lambda, 3*length+3*length*lambda), Vector(0, -1), length);
    rods[4] = Rod(lambda, Vector(3*length-lambda, 3*length+3*length*lambda), Vector(0, -1), length);
    rods[5] = Rod(lambda, Vector(3*length, 3*length+3*length*lambda+0.9*length), Vector(0, -1), length);
    rods[6] = Rod(lambda, Vector(3*length+lambda, 3*length+3*length*lambda+0.9*length), Vector(0, -1), length);
    rods[7] = Rod(lambda, Vector(3*length+2*lambda, 3*length+3*length*lambda), Vector(0, -1), length);*/
}

//output center of mass positions x and orientations (in rad wrt to x axis)
//into stream (eg file)
//makes sure no space appears at eol
void output(std::ostream &out, std::ostream &outextras, const Rods &x)
{
    for (unsigned int i = 0; i < x.size() - 1; ++i)
    {
        out << x[i].r.x << " " << (x[i].L - x[i].r.y) << " ";
    }
    out << x.back().r.x << " " << (x.back().L - x.back().r.y) << std::endl;
    for (unsigned int i = 0; i < x.size() - 1; ++i)
    {
        out << std::atan2(x[i].u.y, x[i].u.x) << " ";
    }
    out << std::atan2(x.back().u.y, x.back().u.x) << std::endl;

    int active = 0;
    for (const auto &rod : x)
    {
        active += int(rod.active);
    }
    outextras << real(active) / x.size() << std::endl;
}

/// START FUNCTION TESTS

//test the PeriodicGet from Buffer.hpp
void main_test_periodicget(void) //periodic_test
{
    int w = 3, h = 3;
    GenericBuffer_dyn<int, MemoryConst::Dynamic, MemoryConst::Dynamic> t(w, h);

    for (unsigned int x = 0; x < t.GetWidth(); ++x)
    {
        for (unsigned int y = 0; y < t.GetHeight(); ++y)
        {
            t.Set(x, y, (y*t.GetWidth()) + x);
        }
    }

    //all pairs should be equal
    std::cout << *t.PeriodicGet(-1, -1) << " " << *t.Get(w-1, h-1) << std::endl;
    std::cout << *t.PeriodicGet(w, h) << " " << *t.Get(0, 0) << std::endl;
    std::cout << *t.PeriodicGet(0, -1) << " " << *t.Get(0, h-1) << std::endl;
    std::cout << *t.PeriodicGet(-1, 0) << " " << *t.Get(w-1, 0) << std::endl;
    std::cout << *t.PeriodicGet(-1, h) << " " << *t.Get(w-1, 0) << std::endl;
    std::cout << *t.PeriodicGet(w, -1) << " " << *t.Get(0, h-1) << std::endl;

    std::cout << *t.PeriodicGet(-2, -2) << " " << *t.Get(w-2, h-2) << std::endl;
    std::cout << *t.PeriodicGet(0, -2) << " " << *t.Get(0, h-2) << std::endl;
    std::cout << *t.PeriodicGet(-2, 0) << " " << *t.Get(w-2, 0) << std::endl;
    std::cout << *t.PeriodicGet(-2, h) << " " << *t.Get(w-2, 0) << std::endl;
    std::cout << *t.PeriodicGet(w, -2) << " " << *t.Get(0, h-2) << std::endl;

    std::cout << *t.PeriodicGet(-1, h-2) << " " << *t.Get(w-1, h-2) << std::endl;
}


//test the distance_periodic_vector function
void main_test_distance_periodic_vector(void)
{
    assert(Vector(2, 0) == distance_periodic_vector(Vector(5, 3), Vector(3, 3), 6));
    assert(Vector(0, -2) == distance_periodic_vector(Vector(1, 2), Vector(1, 4), 6));

    assert(Vector(2, 0) == distance_periodic_vector(Vector(1, 3), Vector(5, 3), 6));
    assert(Vector(-2, 0) == distance_periodic_vector(Vector(5, 3), Vector(1, 3), 6));
    assert(Vector(0, 2) == distance_periodic_vector(Vector(3, 1), Vector(3, 5), 6));
    assert(Vector(0, -2) == distance_periodic_vector(Vector(3, 5), Vector(3, 1), 6));

    assert(Vector(-2, -2) == distance_periodic_vector(Vector(5, 5), Vector(1, 1), 6));
    assert(Vector(2, -2) == distance_periodic_vector(Vector(1, 1), Vector(5, 3), 6));
    assert(Vector(-2, 2) == distance_periodic_vector(Vector(5, 5), Vector(1, 3), 6));
    assert(Vector(2, 2) == distance_periodic_vector(Vector(1, 1), Vector(5, 5), 6));
}

//test the mod function (Vector2D.hpp)
void main_test_mod(void)
{
    assert(1 == mod(1, 10));
    assert(1 == mod(11, 10));
    assert(1 == mod(21, 10));
    assert(9 == mod(-1, 10));
    assert(9 == mod(-11, 10));
    assert(9 == mod(-21, 10));
}

int main_() //_test_normalization(void)
{
    real U0 = 100;
    real F = 0;
    real f0 = 1;
    real L = 10;
    real lambda = 1;
    real dt = 0.002;
    real rmin = 10;
    real passive_frac = 0;
    real t_sim_end = 100;
    unsigned int N = 2;
    unsigned int number_of_boxes = 12123123;
    unsigned int seed = 1337;

    CMDParameterSPR args("", "", U0, F, f0, lambda, dt, rmin, passive_frac, t_sim_end, N, number_of_boxes, seed);

    SPR f1(args);
    Rods rods1(2);
    rods1[0] = Rod(lambda, Vector(3, 3), Vector(1, 0), 4, F, f0, L);
    rods1[1] = Rod(lambda, Vector(6, 6), Vector(1, 0), 8, F, f0, L);

    FixRodPointers(rods1);

    f1(rods1, 2);


    SPR f2(args);
    Rods rods2(2);
    rods2[0] = Rod(lambda, Vector(3, 3), Vector(1, 0), 4, F, f0, L);
    rods2[1] = Rod(lambda, Vector(6, 6), Vector(1, 0), 4, F, f0, L);

    FixRodPointers(rods2);

    f2(rods2, 2);


    SPR f3(args);
    Rods rods3(2);
    rods3[0] = Rod(lambda, Vector(3, 3), Vector(1, 0), 8, F, f0, L);
    rods3[1] = Rod(lambda, Vector(6, 6), Vector(1, 0), 8, F, f0, L);

    FixRodPointers(rods3);

    f3(rods3, 2);


    SPR f4(args);
    Rods rods4(2);
    rods4[0] = Rod(lambda, Vector(5, 5), Vector(1, 0), 4, F, f0, L);
    rods4[1] = Rod(lambda, Vector(5, 6), Vector(1, 0), 4, F, f0, L);

    FixRodPointers(rods4);

    f4(rods4, 2);


    SPR f5(args);
    Rods rods5(2);
    rods5[0] = Rod(lambda, Vector(5, 5), Vector(1, 0), 8, F, f0, L);
    rods5[1] = Rod(lambda, Vector(5, 6), Vector(1, 0), 8, F, f0, L);

    FixRodPointers(rods5);

    f5(rods5, 2);

    return 0;
}

/// END FUNCTION TESTS

// reads lengths/aspect ratios from filename
// file is expected to contain one line with lengths separated by one space
std::vector<real> ReadLengths(std::string filename)
{
    std::vector<real> ret;

    std::ifstream in;
    in.open(filename);
    if (!in)
    {
        throw std::runtime_error("Could not open lengths file!");
    }

    std::string item;

    while (std::getline(in, item, ' '))
    {
        ret.push_back(std::atof(item.c_str()));
    }

    return ret;
}

//MAIN
int main(int argc, char** argv)
{
    CMDParameterSPR args(argc, argv);

    std::vector<real> lengths = ReadLengths(args.in_file);

    args.SetRemainingParameters(lengths);

    //SetMaxlength(MaxElement(lengths));
    //assert(boxlength >= maxlength);

    Rods x(args.N);
    SPR f(args);
    InitRods(x, lengths, args);

    //std::ostream &out = std::cout;
    std::ofstream out;
    out.open(args.out_file);
    if (!out)
    {
        throw std::runtime_error("Could not open output file!");
    }

    //header in sprdata.txt
    out << args.header() << std::endl;

    std::ofstream outextras;
    outextras.open(args.out_file.substr(0, args.out_file.size() - 4) + ".extras.txt"); //TODO kill ".txt.extras.txt"
    if (!outextras)
    {
        throw std::runtime_error("Could not open extras output file!");
    }

    outextras << args.header() << std::endl;
    outextras << "active" << std::endl;

    std::cout << "Initialisation complete." << std::endl;

    std::cout << "t: 0" << std::endl;
    output(out, outextras, x);

    int t_threshold = 1; //when to output

    real t = real(0);

    //initilisation bug: if too close -> integration step yields entangled rods
    //-> fix: make steps even smaller for first dt
    //Euler update steps
    while (t <= 0.1) //0.005
    {
        t += 0.005/100.;

        //Euler update steps
        x += 0.005/100.*f(x, t);
    }

    //initilisation bug: if too close -> integration step yields entangled rods
    //-> fix: make steps even smaller for first dt
    //Euler update steps
    while (t <= 20.)
    {
        t += 0.005;

        if (t >= t_threshold)
        {
            std::cout << "t: " << std::round(t) << std::endl;
        }

        //Euler update steps
        x += 0.005*f(x, t);

        //output
        if (t >= t_threshold)
        {
            output(out, outextras, x);
            t_threshold += 1;
        }
    }

    for (; t < args.t_sim_end; t += args.dt) //real t = 0
    {
        if (t >= t_threshold)
        {
            std::cout << "t: " << std::round(t) << std::endl;
        }

        //Euler update steps //+= operator normalizes u of all rods, but + not
        x += args.dt*f(x, t);

        //Runge Kutta update steps
        /*Rods k1 = dt*f(x, t);
        Rods tmp1 = x + k1/real(2);
        Normalize(tmp1); // todo: and translate for nonperiodic+
        Rods k2 = dt*f(tmp1, t + dt/real(2));
        Rods tmp2 = x + k2/real(2);
        Normalize(tmp2);
        Rods k3 = dt*f(tmp2, t + dt/real(2));
        Rods tmp3 = x + k3;
        Normalize(tmp3);
        Rods k4 = dt*f(tmp3, t + dt);

        //+= operator normalizes u of all rods, but + not
        //x += (k1 + k2*real(2) + k3*real(2) + k4) / real(6); //averaging does not work because of periodic boundary condition
        x += (k1/6 + k2/3 + k3/3 + k4/6); //averaging should still not work much better...*/

        // explicit midpoint method
        //x += k2;

        // Heun method
        /*Rods k1 = dt*f(x, t);
        Rods tmp = x + k1;
        Normalize(tmp);
        Rods k2 = dt*f(tmp, t + dt);
        x += k1/2 + k2/2;*/

        // Ralston method
        /*Rods k1 = dt*f(x, t);
        Rods tmp = x + 2*k1/3;
        Normalize(tmp);
        Rods k2 = dt*f(tmp, t + 2*dt/3);

        //x += (k1 + 3*k2)/4; // DOES NOT WORK !?!?
        x += k1/4 + 3*k2/4;*/


        //output
        if (t >= t_threshold)
        {
            output(out, outextras, x);
            t_threshold += 1;
        }
    }

    std::cout << "Evaluation complete." << std::endl;

    out.close();
    outextras.close();
}
