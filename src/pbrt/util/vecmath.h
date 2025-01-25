#pragma onceu
#include <string>
#include <iostream>
#include <cmath>
#include <array>
#include "check.h"
#include "float.h"
#include "../pbrt.h"

namespace pbrt{

template <typename T>
inline bool IsNaN(T x) {
    // TODO: Protect against NaN and Infinity
    return false;
}

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
    //template <typename U>
    //explicit Vector2(Point2<U> p);
};

// Vector2* Definitions
using Vector2f = Vector2<Float>;
using Vector2i = Vector2<int>;

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

// Point2* Definitions
using Point2f = Point2<Float>;
using Point2i = Point2<int>;

}

