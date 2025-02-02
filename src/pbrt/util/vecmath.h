// Mostly copied from PBRT Code
#pragma once
#include <string>
#include <iostream>
#include <cmath>
#include <array>
#include "check.h"
#include "float.h"
#include "../pbrt.h"
#include "math.h"

namespace {
// TupleLength Definition
// TODO: Figure out what TupleLength is for
template <typename T>
struct TupleLength {
    using type = Float;
};

template <>
struct TupleLength<double> {
    using type = double;
};

template <>
struct TupleLength<long double> {
    using type = long double;
};
}

namespace pbrt{

template <typename T>
inline bool IsNaN(T x) {
    // TODO: Protect against NaN and Infinity
    return false;
}

// 2D
// Tuple2 Definition
template <template <typename> class Child, typename T>
class Tuple2 {
public:
    // Tuple2 Public Members
    T x{}, y{};

    // Tuple2 Public Methods
    //static const int nDimensions = 2;

    Tuple2() = default;
    Tuple2(T x, T y) : x(x), y(y) {}
    bool HasNaN() const { return IsNaN(x) || IsNaN(y); }

    template <typename U>
    auto operator+(Child<U> c) const -> Child<decltype(T{} + U{})> {
        DCHECK(!c.HasNaN());
        return {x + c.x, y + c.y};
    }

    template <typename U>
    Child<T> &operator+=(Child<U> c) {
        DCHECK(!c.HasNaN());
        x += c.x;
        y += c.y;
        return static_cast<Child<T> &>(*this);
    }

    template <typename U>
    auto operator-(Child<U> c) const -> Child<decltype(T{} - U{})> {
        DCHECK(!c.HasNaN());
        return {x - c.x, y - c.y};
    }

    template <typename U>
    Child<T> &operator-=(Child<U> c) {
        DCHECK(!c.HasNaN());
        x -= c.x;
        y -= c.y;
        return static_cast<Child<T> &>(*this);
    }

    bool operator==(Child<T> c) const { return x == c.x && y == c.y; }
    bool operator!=(Child<T> c) const { return x != c.x || y != c.y; }

    template <typename U>
    auto operator*(U s) const -> Child<decltype(T{} * U{})> {
        return {s * x, s * y};
    }

    template <typename U>
    Child<T> &operator*=(U s) {
        DCHECK(!IsNaN(s));
        x *= s;
        y *= s;
        return static_cast<Child<T> &>(*this);
    }

    template <typename U>
    auto operator/(U d) const -> Child<decltype(T{} / U{})> {
        DCHECK(d != 0 && !IsNaN(d));
        return {x / d, y / d};
    }

    template <typename U>
    Child<T> &operator/=(U d) {
        DCHECK_NE(d, 0);
        DCHECK(!IsNaN(d));
        x /= d;
        y /= d;
        return static_cast<Child<T> &>(*this);
    }

    Child<T> operator-() const { return {-x, -y}; }

    T operator[](int i) const {
        DCHECK(i >= 0 && i <= 1);
        return (i == 0) ? x : y;
    }

    T &operator[](int i) {
        DCHECK(i >= 0 && i <= 1);
        return (i == 0) ? x : y;
    }

    friend std::ostream& operator<<(std::ostream& os, const Child<T>& obj) {
        os << "[ " << obj.x << ", " << obj.y << " ]";
        return os;
    }
};

// Tuple2 Inline Functions
template <template <class> class C, typename T, typename U>
inline auto operator*(U s, Tuple2<C, T> t) -> C<decltype(T{} * U{})> {
    DCHECK(!t.HasNaN());
    return t * s;
}

template <template <class> class C, typename T>
inline C<T> Abs(Tuple2<C, T> t) {
    // "argument-dependent lookup..." (here and elsewhere)
    using std::abs;
    return {abs(t.x), abs(t.y)};
}

template <template <class> class C, typename T>
inline C<T> Ceil(Tuple2<C, T> t) {
    using std::ceil;
    return {ceil(t.x), ceil(t.y)};
}

template <template <class> class C, typename T>
inline C<T> Floor(Tuple2<C, T> t) {
    using std::floor;
    return {floor(t.x), floor(t.y)};
}

template <template <class> class C, typename T>
inline auto Lerp(Float t, Tuple2<C, T> t0, Tuple2<C, T> t1) {
    return (1 - t) * t0 + t * t1;
}

template <template <class> class C, typename T>
inline C<T> FMA(Float a, Tuple2<C, T> b, Tuple2<C, T> c) {
    return {FMA(a, b.x, c.x), FMA(a, b.y, c.y)};
}

template <template <class> class C, typename T>
inline C<T> FMA(Tuple2<C, T> a, Float b, Tuple2<C, T> c) {
    return FMA(b, a, c);
}

template <template <class> class C, typename T>
inline C<T> Min(Tuple2<C, T> t0, Tuple2<C, T> t1) {
    using std::min;
    return {min(t0.x, t1.x), min(t0.y, t1.y)};
}

template <template <class> class C, typename T>
inline T MinComponentValue(Tuple2<C, T> t) {
    using std::min;
    return min({t.x, t.y});
}

template <template <class> class C, typename T>
inline int MinComponentIndex(Tuple2<C, T> t) {
    return (t.x < t.y) ? 0 : 1;
}

template <template <class> class C, typename T>
inline C<T> Max(Tuple2<C, T> t0, Tuple2<C, T> t1) {
    using std::max;
    return {max(t0.x, t1.x), max(t0.y, t1.y)};
}

template <template <class> class C, typename T>
inline T MaxComponentValue(Tuple2<C, T> t) {
    using std::max;
    return max({t.x, t.y});
}

template <template <class> class C, typename T>
inline int MaxComponentIndex(Tuple2<C, T> t) {
    return (t.x > t.y) ? 0 : 1;
}

template <template <class> class C, typename T>
inline C<T> Permute(Tuple2<C, T> t, std::array<int, 2> p) {
    return {t[p[0]], t[p[1]]};
}

template <template <class> class C, typename T>
inline T HProd(Tuple2<C, T> t) {
    return t.x * t.y;
}

// Vector2 Definition
template <typename T>
class Vector2 : public Tuple2<Vector2, T> {
  public:
    // Vector2 Public Methods
    using Tuple2<Vector2, T>::x;
    using Tuple2<Vector2, T>::y;

    Vector2() = default;
    Vector2(T x, T y) : Tuple2<pbrt::Vector2, T>(x, y) {}

    template <typename U>
    explicit Vector2(Vector2<U> v)
        : Tuple2<pbrt::Vector2, T>(T(v.x), T(v.y)) {}

    // TODO: Fix this conversion from Point2 to Vector2
    template <typename U>
    explicit Vector2(Point2<U> p);
};

// Vector2* Definitions
using Vector2f = Vector2<Float>;
using Vector2i = Vector2<int>;

// Vector2 Inline Functions
template <typename T>
template <typename U>
Vector2<T>::Vector2(Point2<U> p) : Tuple2<pbrt::Vector2, T>(T(p.x), T(p.y)) {}

template <typename T>
inline auto Dot(Vector2<T> v1, Vector2<T> v2) ->
    typename TupleLength<T>::type {
    DCHECK(!v1.HasNaN() && !v2.HasNaN());
    return SumOfProducts(v1.x, v2.x, v1.y, v2.y);
}

template <typename T>
inline auto AbsDot(Vector2<T> v1, Vector2<T> v2) ->
    typename TupleLength<T>::type {
    DCHECK(!v1.HasNaN() && !v2.HasNaN());
    return std::abs(Dot(v1, v2));
}

template <typename T>
inline auto LengthSquared(Vector2<T> v) -> typename TupleLength<T>::type {
    return Sqr(v.x) + Sqr(v.y);
}

template <typename T>
inline auto Length(Vector2<T> v) -> typename TupleLength<T>::type {
    using std::sqrt;
    return sqrt(LengthSquared(v));
}

template <typename T>
inline auto Normalize(Vector2<T> v) {
    return v / Length(v);
}

template <typename T>
inline auto Distance(Point2<T> p1, Point2<T> p2) ->
    typename TupleLength<T>::type {
    return Length(p1 - p2);
}

template <typename T>
inline auto DistanceSquared(Point2<T> p1, Point2<T> p2) ->
    typename TupleLength<T>::type {
    return LengthSquared(p1 - p2);
}

// Point2 Definition
template <typename T>
class Point2 : public Tuple2<Point2, T> {
  public:
    // Point2 Public Methods
    using Tuple2<Point2, T>::x;
    using Tuple2<Point2, T>::y;
    using Tuple2<Point2, T>::HasNaN;
    using Tuple2<Point2, T>::operator+;
    using Tuple2<Point2, T>::operator+=;
    using Tuple2<Point2, T>::operator*;
    using Tuple2<Point2, T>::operator*=;

    Point2() { x = y = 0; }
    Point2(T x, T y) : Tuple2<pbrt::Point2, T>(x, y) {}
    template <typename U>
    explicit Point2(Point2<U> v) : Tuple2<pbrt::Point2, T>(T(v.x), T(v.y)) {}
    template <typename U>
    explicit Point2(Vector2<U> v)
        : Tuple2<pbrt::Point2, T>(T(v.x), T(v.y)) {}

    template <typename U>
    auto operator+(Vector2<U> v) const -> Point2<decltype(T{} + U{})> {
        DCHECK(!v.HasNaN());
        return {x + v.x, y + v.y};
    }
    template <typename U>
    Point2<T> &operator+=(Vector2<U> v) {
        DCHECK(!v.HasNaN());
        x += v.x;
        y += v.y;
        return *this;
    }
    // There is no Point -= Point because Point-Point = Vector
    // We can use Point += -Point instead
    
    Point2<T> operator-() const { return {-x, -y}; }

    template <typename U>
    auto operator-(Point2<U> p) const -> Vector2<decltype(T{} - U{})> {
        DCHECK(!p.HasNaN());
        return {x - p.x, y - p.y};
    }
    template <typename U>
    auto operator-(Vector2<U> v) const -> Point2<decltype(T{} - U{})> {
        DCHECK(!v.HasNaN());
        return {x - v.x, y - v.y};
    }
    template <typename U>
    Point2<T> &operator-=(Vector2<U> v) {
        DCHECK(!v.HasNaN());
        x -= v.x;
        y -= v.y;
        return *this;
    }
};

// TODO: Point2 Inline Functions

// Point2* Definitions
using Point2f = Point2<Float>;
using Point2i = Point2<int>;

// 3D
// Tuple3 Definition
template <template <typename> class Child, typename T>
class Tuple3 {
  public:
    // Tuple3 Public Methods
    Tuple3() = default;
    Tuple3(T x): x(x), y(x), z(x) { DCHECK(!HasNaN());}    
    Tuple3(T x, T y, T z) : x(x), y(y), z(z) { DCHECK(!HasNaN()); }

    bool HasNaN() const { return IsNaN(x) || IsNaN(y) || IsNaN(z); }

    T operator[](int i) const {
        DCHECK(i >= 0 && i <= 2);
        if (i == 0)
            return x;
        if (i == 1)
            return y;
        return z;
    }

    
    T &operator[](int i) {
        DCHECK(i >= 0 && i <= 2);
        if (i == 0)
            return x;
        if (i == 1)
            return y;
        return z;
    }

    template <typename U>
     auto operator+(Child<U> c) const -> Child<decltype(T{} + U{})> {
        DCHECK(!c.HasNaN());
        return {x + c.x, y + c.y, z + c.z};
    }

    static const int nDimensions = 3;

    Tuple3(Child<T> c) {
        DCHECK(!c.HasNaN());
        x = c.x;
        y = c.y;
        z = c.z;
    }

    Child<T> &operator=(Child<T> c) {
        DCHECK(!c.HasNaN());
        x = c.x;
        y = c.y;
        z = c.z;
        return static_cast<Child<T> &>(*this);
    }

    template <typename U>
     Child<T> &operator+=(Child<U> c) {
        DCHECK(!c.HasNaN());
        x += c.x;
        y += c.y;
        z += c.z;
        return static_cast<Child<T> &>(*this);
    }

    template <typename U>
     auto operator-(Child<U> c) const -> Child<decltype(T{} - U{})> {
        DCHECK(!c.HasNaN());
        return {x - c.x, y - c.y, z - c.z};
    }
    template <typename U>
     Child<T> &operator-=(Child<U> c) {
        DCHECK(!c.HasNaN());
        x -= c.x;
        y -= c.y;
        z -= c.z;
        return static_cast<Child<T> &>(*this);
    }

    
    bool operator==(Child<T> c) const { return x == c.x && y == c.y && z == c.z; }
    
    bool operator!=(Child<T> c) const { return x != c.x || y != c.y || z != c.z; }

    template <typename U>
     auto operator*(U s) const -> Child<decltype(T{} * U{})> {
        return {s * x, s * y, s * z};
    }
    template <typename U>
     Child<T> &operator*=(U s) {
        DCHECK(!IsNaN(s));
        x *= s;
        y *= s;
        z *= s;
        return static_cast<Child<T> &>(*this);
    }

    template <typename U>
     auto operator/(U d) const -> Child<decltype(T{} / U{})> {
        DCHECK_NE(d, 0);
        return {x / d, y / d, z / d};
    }
    template <typename U>
     Child<T> &operator/=(U d) {
        DCHECK_NE(d, 0);
        x /= d;
        y /= d;
        z /= d;
        return static_cast<Child<T> &>(*this);
    }
    
    Child<T> operator-() const { return {-x, -y, -z}; }

    friend std::ostream& operator<<(std::ostream& os, const Child<T>& obj) {
        os << "[ " << obj.x << ", " << obj.y << ", " << obj.z << " ]";
        return os;
    }

    // Tuple3 Public Members
    T x{}, y{}, z{};
};

// Tuple3 Inline Functions
template <template <class> class C, typename T, typename U>
 inline auto operator*(U s, Tuple3<C, T> t) -> C<decltype(T{} * U{})> {
    return t * s;
}

template <template <class> class C, typename T>
 inline C<T> Abs(Tuple3<C, T> t) {
    using std::abs;
    return {abs(t.x), abs(t.y), abs(t.z)};
}

template <template <class> class C, typename T>
 inline C<T> Ceil(Tuple3<C, T> t) {
    using std::ceil;
    return {ceil(t.x), ceil(t.y), ceil(t.z)};
}

template <template <class> class C, typename T>
 inline C<T> Floor(Tuple3<C, T> t) {
    using std::floor;
    return {floor(t.x), floor(t.y), floor(t.z)};
}

template <template <class> class C, typename T>
 inline auto Lerp(Float t, Tuple3<C, T> t0, Tuple3<C, T> t1) {
    return (1 - t) * t0 + t * t1;
}

template <template <class> class C, typename T>
 inline C<T> FMA(Float a, Tuple3<C, T> b, Tuple3<C, T> c) {
    return {FMA(a, b.x, c.x), FMA(a, b.y, c.y), FMA(a, b.z, c.z)};
}

template <template <class> class C, typename T>
 inline C<T> FMA(Tuple3<C, T> a, Float b, Tuple3<C, T> c) {
    return FMA(b, a, c);
}

template <template <class> class C, typename T>
 inline C<T> Min(Tuple3<C, T> t1, Tuple3<C, T> t2) {
    using std::min;
    return {min(t1.x, t2.x), min(t1.y, t2.y), min(t1.z, t2.z)};
}

template <template <class> class C, typename T>
 inline T MinComponentValue(Tuple3<C, T> t) {
    using std::min;
    return min({t.x, t.y, t.z});
}

template <template <class> class C, typename T>
 inline int MinComponentIndex(Tuple3<C, T> t) {
    return (t.x < t.y) ? ((t.x < t.z) ? 0 : 2) : ((t.y < t.z) ? 1 : 2);
}

template <template <class> class C, typename T>
 inline C<T> Max(Tuple3<C, T> t1, Tuple3<C, T> t2) {
    using std::max;
    return {max(t1.x, t2.x), max(t1.y, t2.y), max(t1.z, t2.z)};
}

template <template <class> class C, typename T>
 inline T MaxComponentValue(Tuple3<C, T> t) {
    using std::max;
    return max({t.x, t.y, t.z});
}

template <template <class> class C, typename T>
 inline int MaxComponentIndex(Tuple3<C, T> t) {
    return (t.x > t.y) ? ((t.x > t.z) ? 0 : 2) : ((t.y > t.z) ? 1 : 2);
}

template <template <class> class C, typename T>
 inline C<T> Permute(Tuple3<C, T> t, std::array<int, 3> p) {
    return {t[p[0]], t[p[1]], t[p[2]]};
}

template <template <class> class C, typename T>
 inline T HProd(Tuple3<C, T> t) {
    return t.x * t.y * t.z;
}

// Vector3 Definition
template <typename T>
class Vector3 : public Tuple3<Vector3, T> {
  public:
    // Vector3 Public Methods
    using Tuple3<Vector3, T>::x;
    using Tuple3<Vector3, T>::y;
    using Tuple3<Vector3, T>::z;

    Vector3() = default;
    Vector3(T x) : Tuple3<pbrt::Vector3, T>(x) {};    
    Vector3(T x, T y, T z) : Tuple3<pbrt::Vector3, T>(x, y, z) {}

    template <typename U>
     explicit Vector3(Vector3<U> v)
        : Tuple3<pbrt::Vector3, T>(T(v.x), T(v.y), T(v.z)) {}

    template <typename U>
    explicit Vector3(Point3<U> p);
    template <typename U>
    explicit Vector3(Normal3<U> n);
};

// Vector3* Definitions
using Vector3f = Vector3<Float>;
using Vector3i = Vector3<int>;

// Vector3 Inline Functions
template <typename T>
template <typename U>
Vector3<T>::Vector3(Point3<U> p) : Tuple3<pbrt::Vector3, T>(T(p.x), T(p.y), T(p.z)) {}

template <typename T>
inline Vector3<T> HorizontalProduct(Vector3<T> v1, Vector3<T> v2) {
    return {v1.x*v2.x, v1.y*v2.y, v1.z*v2.z};
}

template <typename T>
inline Vector3<T> Cross(Vector3<T> v1, Normal3<T> v2) {
    DCHECK(!v1.HasNaN() && !v2.HasNaN());
    return {DifferenceOfProducts(v1.y, v2.z, v1.z, v2.y),
            DifferenceOfProducts(v1.z, v2.x, v1.x, v2.z),
            DifferenceOfProducts(v1.x, v2.y, v1.y, v2.x)};
}

template <typename T>
inline Vector3<T> Cross(Normal3<T> v1, Vector3<T> v2) {
    DCHECK(!v1.HasNaN() && !v2.HasNaN());
    return {DifferenceOfProducts(v1.y, v2.z, v1.z, v2.y),
            DifferenceOfProducts(v1.z, v2.x, v1.x, v2.z),
            DifferenceOfProducts(v1.x, v2.y, v1.y, v2.x)};
}

template <typename T>
inline T LengthSquared(Vector3<T> v) {
    return Sqr(v.x) + Sqr(v.y) + Sqr(v.z);
}

template <typename T>
inline auto Length(Vector3<T> v) -> typename TupleLength<T>::type {
    using std::sqrt;
    return sqrt(LengthSquared(v));
}

template <typename T>
inline auto Normalize(Vector3<T> v) {
    return v / Length(v);
}

template <typename T>
inline T Dot(Vector3<T> v, Vector3<T> w) {
    DCHECK(!v.HasNaN() && !w.HasNaN());
    return v.x * w.x + v.y * w.y + v.z * w.z;
}

// Equivalent to std::acos(Dot(a, b)), but more numerically stable.
// via http://www.plunk.org/~hatch/rightway.html
template <typename T>
inline Float AngleBetween(Vector3<T> v1, Vector3<T> v2) {
    if (Dot(v1, v2) < 0)
        return Pi - 2 * SafeASin(Length(v1 + v2) / 2);
    else
        return 2 * SafeASin(Length(v2 - v1) / 2);
}

template <typename T>
inline T AbsDot(Vector3<T> v1, Vector3<T> v2) {
    DCHECK(!v1.HasNaN() && !v2.HasNaN());
    return std::abs(Dot(v1, v2));
}

template <typename T>
inline Float AngleBetween(Normal3<T> a, Normal3<T> b) {
    if (Dot(a, b) < 0)
        return Pi - 2 * SafeASin(Length(a + b) / 2);
    else
        return 2 * SafeASin(Length(b - a) / 2);
}

template <typename T>
inline Vector3<T> GramSchmidt(Vector3<T> v, Vector3<T> w) {
    return v - Dot(v, w) * w;
}

template <typename T>
inline Vector3<T> Cross(Vector3<T> v, Vector3<T> w) {
    DCHECK(!v.HasNaN() && !w.HasNaN());
    return {DifferenceOfProducts(v.y, w.z, v.z, w.y),
            DifferenceOfProducts(v.z, w.x, v.x, w.z),
            DifferenceOfProducts(v.x, w.y, v.y, w.x)};
}

template <typename T>
inline void CoordinateSystem(Vector3<T> v1, Vector3<T> *v2, Vector3<T> *v3) {
    Float sign = std::copysign(Float(1), v1.z);
    Float a = -1 / (sign + v1.z);
    Float b = v1.x * v1.y * a;
    *v2 = Vector3<T>(1 + sign * Sqr(v1.x) * a, sign * b, -sign * v1.x);
    *v3 = Vector3<T>(b, sign + Sqr(v1.y) * a, -v1.y);
}

template <typename T>
inline void CoordinateSystem(Normal3<T> v1, Vector3<T> *v2, Vector3<T> *v3) {
    Float sign = std::copysign(Float(1), v1.z);
    Float a = -1 / (sign + v1.z);
    Float b = v1.x * v1.y * a;
    *v2 = Vector3<T>(1 + sign * Sqr(v1.x) * a, sign * b, -sign * v1.x);
    *v3 = Vector3<T>(b, sign + Sqr(v1.y) * a, -v1.y);
}

template <typename T>
template <typename U>
Vector3<T>::Vector3(Normal3<U> n) : Tuple3<pbrt::Vector3, T>(T(n.x), T(n.y), T(n.z)) {}

template <typename T>
inline Vector3<T> Reflect(const Vector3<T>& v, const Normal3<T>& n) {
    return v - Vector3<T>(2*Dot(v,n)*n);
}
// TODO: Vector3fi

// Point3 Definition
template <typename T>
class Point3 : public Tuple3<Point3, T> {
  public:
    // Point3 Public Methods
    using Tuple3<Point3, T>::x;
    using Tuple3<Point3, T>::y;
    using Tuple3<Point3, T>::z;
    using Tuple3<Point3, T>::HasNaN;
    using Tuple3<Point3, T>::operator+;
    using Tuple3<Point3, T>::operator+=;
    using Tuple3<Point3, T>::operator*;
    using Tuple3<Point3, T>::operator*=;

    Point3() = default;
    Point3(T x) : Tuple3<pbrt::Point3, T>(x) {}    
    Point3(T x, T y, T z) : Tuple3<pbrt::Point3, T>(x, y, z) {}

    // We can't do using operator- above, since we don't want to pull in
    // the Point-Point -> Point one so that we can return a vector
    // instead...
    
    Point3<T> operator-() const { return {-x, -y, -z}; }

    template <typename U>
     explicit Point3(Point3<U> p)
        : Tuple3<pbrt::Point3, T>(T(p.x), T(p.y), T(p.z)) {}
    template <typename U>
     explicit Point3(Vector3<U> v)
        : Tuple3<pbrt::Point3, T>(T(v.x), T(v.y), T(v.z)) {}

    template <typename U>
     auto operator+(Vector3<U> v) const -> Point3<decltype(T{} + U{})> {
        DCHECK(!v.HasNaN());
        return {x + v.x, y + v.y, z + v.z};
    }
    template <typename U>
     Point3<T> &operator+=(Vector3<U> v) {
        DCHECK(!v.HasNaN());
        x += v.x;
        y += v.y;
        z += v.z;
        return *this;
    }

    template <typename U>
     auto operator-(Vector3<U> v) const -> Point3<decltype(T{} - U{})> {
        DCHECK(!v.HasNaN());
        return {x - v.x, y - v.y, z - v.z};
    }
    template <typename U>
     Point3<T> &operator-=(Vector3<U> v) {
        DCHECK(!v.HasNaN());
        x -= v.x;
        y -= v.y;
        z -= v.z;
        return *this;
    }

    template <typename U>
     auto operator-(Point3<U> p) const -> Vector3<decltype(T{} - U{})> {
        DCHECK(!p.HasNaN());
        return {x - p.x, y - p.y, z - p.z};
    }
};

// Point3* Definitions
using Point3f = Point3<Float>;
using Point3i = Point3<int>;
 
// Point3 Inline Functions
template <typename T>
inline auto Distance(Point3<T> p1, Point3<T> p2) {
    return Length(p1 - p2);
}

template <typename T>
inline auto DistanceSquared(Point3<T> p1, Point3<T> p2) {
    return LengthSquared(p1 - p2);
}

// TODO: Point3fi

// Normal3 Definition
template <typename T>
class Normal3 : public Tuple3<Normal3, T> {
  public:
    // Normal3 Public Methods
    using Tuple3<Normal3, T>::x;
    using Tuple3<Normal3, T>::y;
    using Tuple3<Normal3, T>::z;
    using Tuple3<Normal3, T>::HasNaN;
    using Tuple3<Normal3, T>::operator+;
    using Tuple3<Normal3, T>::operator*;
    using Tuple3<Normal3, T>::operator*=;

    Normal3() = default;
    Normal3(T x) : Tuple3<pbrt::Normal3, T>(x) {}    
    Normal3(T x, T y, T z) : Tuple3<pbrt::Normal3, T>(x, y, z) {}
    template <typename U>
    explicit Normal3<T>(Normal3<U> v)
        : Tuple3<pbrt::Normal3, T>(T(v.x), T(v.y), T(v.z)) {}

    template <typename U>
    explicit Normal3<T>(Vector3<U> v)
        : Tuple3<pbrt::Normal3, T>(T(v.x), T(v.y), T(v.z)) {}
};

using Normal3f = Normal3<Float>;

// Normal3 Inline Functions
template <typename T>
inline auto LengthSquared(Normal3<T> n) -> typename TupleLength<T>::type {
    return Sqr(n.x) + Sqr(n.y) + Sqr(n.z);
}

template <typename T>
inline auto Length(Normal3<T> n) -> typename TupleLength<T>::type {
    using std::sqrt;
    return sqrt(LengthSquared(n));
}

template <typename T>
inline auto Normalize(Normal3<T> n) {
    return n / Length(n);
}

template <typename T>
inline auto Dot(Normal3<T> n, Vector3<T> v) ->
    typename TupleLength<T>::type {
    DCHECK(!n.HasNaN() && !v.HasNaN());
    //return FMA(n.x, v.x, SumOfProducts(n.y, v.y, n.z, v.z));
    return n.x*v.x + n.y*v.y + n.z*v.z;
}

template <typename T>
inline auto Dot(Vector3<T> v, Normal3<T> n) ->
    typename TupleLength<T>::type {
    DCHECK(!v.HasNaN() && !n.HasNaN());
    //return FMA(n.x, v.x, SumOfProducts(n.y, v.y, n.z, v.z));
    return n.x*v.x + n.y*v.y + n.z*v.z;
}

template <typename T>
inline auto Dot(Normal3<T> n1, Normal3<T> n2) ->
    typename TupleLength<T>::type {
    DCHECK(!n1.HasNaN() && !n2.HasNaN());
    //return FMA(n1.x, n2.x, SumOfProducts(n1.y, n2.y, n1.z, n2.z));
    return n1.x*n2.x + n1.y*n2.y + n1.z*n2.z;
}

template <typename T>
inline auto AbsDot(Normal3<T> n, Vector3<T> v) ->
    typename TupleLength<T>::type {
    DCHECK(!n.HasNaN() && !v.HasNaN());
    return std::abs(Dot(n, v));
}

template <typename T>
inline auto AbsDot(Vector3<T> v, Normal3<T> n) ->
    typename TupleLength<T>::type {
    using std::abs;
    DCHECK(!v.HasNaN() && !n.HasNaN());
    return abs(Dot(v, n));
}

template <typename T>
inline auto AbsDot(Normal3<T> n1, Normal3<T> n2) ->
    typename TupleLength<T>::type {
    using std::abs;
    DCHECK(!n1.HasNaN() && !n2.HasNaN());
    return abs(Dot(n1, n2));
}

template <typename T>
inline Normal3<T> FaceForward(Normal3<T> n, Vector3<T> v) {
    return (Dot(n, v) < 0.f) ? -n : n;
}

template <typename T>
inline Normal3<T> FaceForward(Normal3<T> n, Normal3<T> n2) {
    return (Dot(n, n2) < 0.f) ? -n : n;
}

template <typename T>
inline Vector3<T> FaceForward(Vector3<T> v, Vector3<T> v2) {
    return (Dot(v, v2) < 0.f) ? -v : v;
}

template <typename T>
inline Vector3<T> FaceForward(Vector3<T> v, Normal3<T> n2) {
    return (Dot(v, n2) < 0.f) ? -v : v;
}
}


