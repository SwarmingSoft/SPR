#ifndef Vector2D_HEADER
#define Vector2D_HEADER

#include <algorithm>
#include <cmath>
#include <limits>

template <typename X, typename Y> auto mod(X x, Y y) -> decltype(std::fmod(std::fmod(x, y) + y, y))
{
    return std::fmod(std::fmod(x, y) + y, y);
}

/*
absolute vector: implicit start is at zero, gives position: SquarePeriodicVector2DT
relative vector: gives relative translation: Vector2DT

abs - abs: cycled rel
abs - rel: cycled abs
abs + abs: ?
abs + rel: cycled abs

abs * scalar: cycle
rel * scalar: none

abs / scalar: cycle
rel / scalar: none

abs * rel: not possible

rel*rel: none
abs.norm(): none
rel.norm(): none


*/

template <class T> class Vector2DT
{
public:

    union
    {
        struct
        {
            T x, y;
        };
        T v[2];
    };

    Vector2DT(void) noexcept;
    Vector2DT(T x, T y) noexcept;
    T Norm(void) const;
    T atan2(void) const;
    T GetDistance(const Vector2DT<T> &rhs) const;
    T GetAngle(const Vector2DT<T> &rhs) const;
    void Normalize(void);
    Vector2DT<T> Normal(void) const;
    Vector2DT<T> &Periodify(T edge);
    Vector2DT<T> Periodified(T edge) const;

    bool operator==(const Vector2DT<T> &rhs) const;
    bool operator!=(const Vector2DT<T> &rhs) const;
    Vector2DT<T> operator+(const Vector2DT<T> &rhs) const;
    Vector2DT<T> operator-(void) const;
    Vector2DT<T> operator-(const Vector2DT<T> &rhs) const;
    Vector2DT<T> operator*(T scalar) const;
    Vector2DT<T> operator/(T scalar) const;
    T operator*(const Vector2DT<T> &rhs) const;
    Vector2DT<T> &operator+=(const Vector2DT<T> &rhs);
    Vector2DT<T> &operator-=(const Vector2DT<T> &rhs);
    Vector2DT<T> &operator*=(T scalar);
    Vector2DT<T> &operator/=(T scalar);
    bool Compare(const Vector2DT<T> &rhs, T delta = 10) const;
};

//typedef Vector2DT<double> Vector2D;

template <class T> Vector2DT<T>::Vector2DT(void) noexcept {} //uninitialized vector

template <class T> Vector2DT<T>::Vector2DT(T x, T y) noexcept : x(x), y(y) {}

template <class T> T Vector2DT<T>::Norm(void) const
{
    return std::sqrt(x * x + y * y);
}

template <class T> T Vector2DT<T>::atan2(void) const
{
    return std::atan2(y, x);
}

template <class T> T Vector2DT<T>::GetDistance(const Vector2DT<T> &rhs) const
{
    return (rhs - *this).Norm();
}

template <class T> T Vector2DT<T>::GetAngle(const Vector2DT<T> &rhs) const
{
    return std::acos(*this * rhs) / (Norm() * rhs.Norm());
}

template <class T> void Vector2DT<T>::Normalize(void)
{
    *this /= Norm();
}

template <class T>  Vector2DT<T> Vector2DT<T>::Normal(void) const
{
    return *this / Norm();
}

template <class T> bool Vector2DT<T>::operator==(const Vector2DT<T> &rhs) const
{
    return x == rhs.x && y == rhs.y;
}

template <class T> bool Vector2DT<T>::operator!=(const Vector2DT<T> &rhs) const
{
    return x != rhs.x || y != rhs.y;
}

template <class T> Vector2DT<T> Vector2DT<T>::operator+(const Vector2DT<T> &rhs) const
{
    return Vector2DT<T>(x + rhs.x, y + rhs.y);
}

template <class T> Vector2DT<T> Vector2DT<T>::operator-(const Vector2DT<T> &rhs) const
{
    return Vector2DT<T>(x - rhs.x, y - rhs.y);
}

template <class T> Vector2DT<T> Vector2DT<T>::operator-(void) const
{
    return Vector2DT<T>(-x, -y);
}

template <class T> Vector2DT<T> Vector2DT<T>::operator*(T scalar) const
{
    return Vector2DT<T>(x * scalar, y * scalar);
}

template <class T> Vector2DT<T> Vector2DT<T>::operator/(T scalar) const
{
    return Vector2DT<T>(x / scalar, y / scalar);
}

template <class T> T Vector2DT<T>::operator*(const Vector2DT<T> &rhs) const
{
    return x * rhs.x + y * rhs.y;
}

template <class T> Vector2DT<T> &Vector2DT<T>::operator+=(const Vector2DT<T> &rhs)
{
    x += rhs.x;
    y += rhs.y;
    return *this;
}

template <class T> Vector2DT<T> &Vector2DT<T>::operator-=(const Vector2DT<T> &rhs)
{
    x -= rhs.x;
    y -= rhs.y;
    return *this;
}

template <class T> Vector2DT<T> &Vector2DT<T>::operator*=(T scalar)
{
    x *= scalar;
    y *= scalar;
    return *this;
}

template <class T> Vector2DT<T> &Vector2DT<T>::operator/=(T scalar)
{
    x /= scalar;
    y /= scalar;
    return *this;
}

template <class T> bool Vector2DT<T>::Compare(const Vector2DT<T> &rhs, T delta) const
{
    return std::abs(x - rhs.x) < std::numeric_limits<T>::epsilon()*delta
        && std::abs(y - rhs.y) < std::numeric_limits<T>::epsilon()*delta;
}

template <class T> Vector2DT<T> operator*(T scalar, const Vector2DT<T> &rhs)
{
    return rhs * scalar;
}

template <class T> Vector2DT<T> &Vector2DT<T>::Periodify(T edge)
{
    *this = Periodified(edge);
    return *this;
}

template <class T> Vector2DT<T> Vector2DT<T>::Periodified(T edge) const
{
    return Vector2DT<T>(mod(x, edge), mod(y, edge));
}

template <class T, T edge> class SquarePeriodicVector2DT : public Vector2DT<T>
{
public:

    SquarePeriodicVector2DT(T x, T y);
    T GetDistance(const SquarePeriodicVector2DT<T,edge> &rhs) const;
    SquarePeriodicVector2DT<T,edge> &Periodify(void);
    SquarePeriodicVector2DT<T,edge> Periodified(void) const;

    //Vector2DT<T> operator+(const SquarePeriodicVector2DT<T,edge> &rhs) const;
    SquarePeriodicVector2DT<T,edge> operator+(const Vector2DT<T> &rhs) const;
    Vector2DT<T> operator-(const SquarePeriodicVector2DT<T,edge> &rhs) const;
    SquarePeriodicVector2DT<T,edge> operator-(const Vector2DT<T> &rhs) const;
    SquarePeriodicVector2DT<T,edge> operator-(void) const;
    SquarePeriodicVector2DT<T,edge> operator*(T scalar) const;
    SquarePeriodicVector2DT<T,edge> operator/(T scalar) const;
    SquarePeriodicVector2DT<T,edge> &operator+=(const Vector2DT<T> &rhs);
    SquarePeriodicVector2DT<T,edge> &operator-=(const Vector2DT<T> &rhs);
    SquarePeriodicVector2DT<T,edge> &operator*=(T scalar);
    SquarePeriodicVector2DT<T,edge> &operator/=(T scalar);
    bool Compare(const SquarePeriodicVector2DT<T,edge> &rhs, T delta = 10) const;
};

//double is sadly not valid as template param
//template <double edge> using SquarePeriodicVector2D = SquarePeriodicVector2DT<double, edge>;

template <class T, T edge> SquarePeriodicVector2DT<T,edge>::SquarePeriodicVector2DT(T x, T y) : Vector2DT<T>(x, y)
{
    assert(x >= 0 && y >= 0 && x < edge && y < edge);
}

// faster than (a - b).Normalized()
template <class T, T edge> T SquarePeriodicVector2DT<T,edge>::GetDistance(const SquarePeriodicVector2DT<T,edge> &rhs) const
{
    T absx = std::abs(this->x-rhs.x);
    T absy = std::abs(this->y-rhs.y);
    T minx = std::min(absx, edge-absx);
    T miny = std::min(absy, edge-absy);
    return std::sqrt(minx*minx + miny*miny);
}

template <class T, T edge> SquarePeriodicVector2DT<T,edge> SquarePeriodicVector2DT<T,edge>::operator+(const Vector2DT<T> &rhs) const
{
    return SquarePeriodicVector2DT<T,edge>(mod(this->x + rhs.x, edge), mod(this->y + rhs.y, edge));
}

template <class T, T edge> Vector2DT<T> SquarePeriodicVector2DT<T,edge>::operator-(const SquarePeriodicVector2DT<T,edge> &rhs) const
{
    T dx = this->x - rhs.x;
    T dy = this->y - rhs.y;

    if (dx > 0.5 * edge)
    {
        dx = dx - edge;
    }
    else if (dx <= -0.5 * edge)
    {
        dx = dx + edge;
    }

    if (dy > 0.5 * edge) // >= ?
    {
        dy = dy - edge;
    }
    else if (dy <= -0.5 * edge)
    {
        dy = dy + edge;
    }

    return Vector2DT<T>(dx, dy);
}

template <class T, T edge> SquarePeriodicVector2DT<T,edge> SquarePeriodicVector2DT<T,edge>::operator-(const Vector2DT<T> &rhs) const
{
    return SquarePeriodicVector2DT<T,edge>(mod(this->x - rhs.x, edge), mod(this->y - rhs.y, edge));
}

template <class T, T edge> SquarePeriodicVector2DT<T,edge> SquarePeriodicVector2DT<T,edge>::operator-(void) const
{
    return SquarePeriodicVector2DT<T,edge>(mod(-this->x, edge), mod(-this->y, edge));
}

template <class T, T edge> SquarePeriodicVector2DT<T,edge> SquarePeriodicVector2DT<T,edge>::operator*(T scalar) const
{
    return SquarePeriodicVector2DT<T,edge>(mod(this->x * scalar, edge), mod(this->y * scalar, edge));
}

template <class T, T edge> SquarePeriodicVector2DT<T,edge> SquarePeriodicVector2DT<T,edge>::operator/(T scalar) const
{
    return SquarePeriodicVector2DT<T,edge>(mod(this->x / scalar, edge), mod(this->y / scalar, edge));
}

template <class T, T edge> SquarePeriodicVector2DT<T,edge> &SquarePeriodicVector2DT<T,edge>::operator+=(const Vector2DT<T> &rhs)
{
    Vector2DT<T>::operator+=(rhs);

    Periodify();

    return *this;
}

template <class T, T edge> SquarePeriodicVector2DT<T,edge> &SquarePeriodicVector2DT<T,edge>::operator-=(const Vector2DT<T> &rhs)
{
    Vector2DT<T>::operator-=(rhs);

    Periodify();

    return *this;
}

template <class T, T edge> SquarePeriodicVector2DT<T,edge> &SquarePeriodicVector2DT<T,edge>::operator*=(T scalar)
{
    Vector2DT<T>::operator*=(scalar);

    Periodify();

    return *this;
}

template <class T, T edge> SquarePeriodicVector2DT<T,edge> &SquarePeriodicVector2DT<T,edge>::operator/=(T scalar)
{
    Vector2DT<T>::operator/=(scalar);

    Periodify();

    return *this;
}

template <class T, T edge> bool SquarePeriodicVector2DT<T,edge>::Compare(const SquarePeriodicVector2DT<T,edge> &rhs, T delta) const
{
    return std::abs(mod(this->x - rhs.x, edge)) < std::numeric_limits<T>::epsilon()*delta
        && std::abs(mod(this->y - rhs.y, edge)) < std::numeric_limits<T>::epsilon()*delta;
}

template <class T, T edge> SquarePeriodicVector2DT<T,edge> &SquarePeriodicVector2DT<T,edge>::Periodify(void)
{
    this->x = mod(this->x, edge);
    this->y = mod(this->y, edge);
    return *this;
}

template <class T, T edge> SquarePeriodicVector2DT<T,edge> SquarePeriodicVector2DT<T,edge>::Periodified(void) const
{
    return SquarePeriodicVector2DT<T,edge>(mod(this->x, edge), mod(this->y, edge));
}

namespace std
{
    /*template <class T> T abs(Vector2DT<T> v)
    {
        return v.Norm();
    }*/

    /*template <class T> Vector2DT<T> abs(Vector2DT<T> v)
    {
        return Vector2DT<T>(std::abs(x), std::abs(y));
    }*/
}

#endif //Vector2D_HEADER
