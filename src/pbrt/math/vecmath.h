// Mostly copied from PBRT Code
#pragma once
#include <array>
#include <cmath>
#include <iostream>
#include <string>
#include "pbrt/math/math.h"
#include "pbrt/pbrt.h"
#include "pbrt/test/check.h"

/**
Required by:
 * ray.h
 * interaction.h
Classes
 * Tuple2
 * Tuple3
 * Vector2
 * Vector3
 * Normal3
 * Point2
 * Point3
**/

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

namespace pbrt {

template <typename T>
inline bool IsNaN(T x)
{
  // TODO: Protect against NaN and Infinity
  return false;
}

// 2D
// Tuple2 Definition
template <template <typename> class Child, typename T>
class Tuple2 {
 public:
  // Tuple2 Public Members
  T x {}, y {};

  // Tuple2 Public Methods
  // static const int nDimensions = 2;

  Tuple2() = default;
  Tuple2(T x, T y): x(x), y(y) {}
  bool HasNaN() const { return IsNaN(x) || IsNaN(y); }

  template <typename U>
  auto operator+(Child<U> c) const
      -> Child<decltype(T {} + U {})>
  {
    DCHECK(!c.HasNaN());
    return {x + c.x, y + c.y};
  }

  template <typename U>
  Child<T> &operator+=(Child<U> c)
  {
    DCHECK(!c.HasNaN());
    x += c.x;
    y += c.y;
    return static_cast<Child<T> &>(*this);
  }

  template <typename U>
  auto operator-(Child<U> c) const
      -> Child<decltype(T {} - U {})>
  {
    DCHECK(!c.HasNaN());
    return {x - c.x, y - c.y};
  }

  template <typename U>
  Child<T> &operator-=(Child<U> c)
  {
    DCHECK(!c.HasNaN());
    x -= c.x;
    y -= c.y;
    return static_cast<Child<T> &>(*this);
  }

  bool operator==(Child<T> c) const
  {
    return x == c.x && y == c.y;
  }
  bool operator!=(Child<T> c) const
  {
    return x != c.x || y != c.y;
  }

  template <typename U>
  auto operator*(U s) const -> Child<decltype(T {} * U {})>
  {
    return {s * x, s * y};
  }

  template <typename U>
  Child<T> &operator*=(U s)
  {
    DCHECK(!IsNaN(s));
    x *= s;
    y *= s;
    return static_cast<Child<T> &>(*this);
  }

  template <typename U>
  auto operator/(U d) const -> Child<decltype(T {} / U {})>
  {
    DCHECK(d != 0 && !IsNaN(d));
    return {x / d, y / d};
  }

  template <typename U>
  Child<T> &operator/=(U d)
  {
    DCHECK_NE(d, 0);
    DCHECK(!IsNaN(d));
    x /= d;
    y /= d;
    return static_cast<Child<T> &>(*this);
  }

  Child<T> operator-() const { return {-x, -y}; }

  T operator[](int i) const
  {
    DCHECK(i >= 0 && i <= 1);
    return (i == 0) ? x : y;
  }

  T &operator[](int i)
  {
    DCHECK(i >= 0 && i <= 1);
    return (i == 0) ? x : y;
  }

  friend std::ostream &operator<<(
      std::ostream &os, const Child<T> &obj
  )
  {
    os << "[ " << obj.x << ", " << obj.y << " ]";
    return os;
  }
};

// Tuple2 Inline Functions
template <template <class> class C, typename T, typename U>
inline auto operator*(U s, Tuple2<C, T> t)
    -> C<decltype(T {} * U {})>
{
  DCHECK(!t.HasNaN());
  return t * s;
}

template <template <class> class C, typename T>
inline C<T> Abs(Tuple2<C, T> t)
{
  // "argument-dependent lookup..." (here and elsewhere)
  using std::abs;
  return {abs(t.x), abs(t.y)};
}

template <template <class> class C, typename T>
inline C<T> Ceil(Tuple2<C, T> t)
{
  using std::ceil;
  return {ceil(t.x), ceil(t.y)};
}

template <template <class> class C, typename T>
inline C<T> Floor(Tuple2<C, T> t)
{
  using std::floor;
  return {floor(t.x), floor(t.y)};
}

template <template <class> class C, typename T>
inline auto Lerp(Float t, Tuple2<C, T> t0, Tuple2<C, T> t1)
{
  return (1 - t) * t0 + t * t1;
}

template <template <class> class C, typename T>
inline C<T> FMA(Float a, Tuple2<C, T> b, Tuple2<C, T> c)
{
  return {FMA(a, b.x, c.x), FMA(a, b.y, c.y)};
}

template <template <class> class C, typename T>
inline C<T> FMA(Tuple2<C, T> a, Float b, Tuple2<C, T> c)
{
  return FMA(b, a, c);
}

template <template <class> class C, typename T>
inline C<T> Min(Tuple2<C, T> t0, Tuple2<C, T> t1)
{
  using std::min;
  return {min(t0.x, t1.x), min(t0.y, t1.y)};
}

template <template <class> class C, typename T>
inline T MinComponentValue(Tuple2<C, T> t)
{
  using std::min;
  return min({t.x, t.y});
}

template <template <class> class C, typename T>
inline int MinComponentIndex(Tuple2<C, T> t)
{
  return (t.x < t.y) ? 0 : 1;
}

template <template <class> class C, typename T>
inline C<T> Max(Tuple2<C, T> t0, Tuple2<C, T> t1)
{
  using std::max;
  return {max(t0.x, t1.x), max(t0.y, t1.y)};
}

template <template <class> class C, typename T>
inline T MaxComponentValue(Tuple2<C, T> t)
{
  using std::max;
  return max({t.x, t.y});
}

template <template <class> class C, typename T>
inline int MaxComponentIndex(Tuple2<C, T> t)
{
  return (t.x > t.y) ? 0 : 1;
}

template <template <class> class C, typename T>
inline C<T> Permute(Tuple2<C, T> t, std::array<int, 2> p)
{
  return {t[p[0]], t[p[1]]};
}

template <template <class> class C, typename T>
inline T HProd(Tuple2<C, T> t)
{
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
  Vector2(T x, T y): Tuple2<pbrt::Vector2, T>(x, y) {}

  template <typename U>
  explicit Vector2(Vector2<U> v)
      : Tuple2<pbrt::Vector2, T>(T(v.x), T(v.y))
  {}

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
Vector2<T>::Vector2(Point2<U> p)
    : Tuple2<pbrt::Vector2, T>(T(p.x), T(p.y))
{}

template <typename T>
inline auto Dot(Vector2<T> v1, Vector2<T> v2) ->
    typename TupleLength<T>::type
{
  DCHECK(!v1.HasNaN() && !v2.HasNaN());
  return SumOfProducts(v1.x, v2.x, v1.y, v2.y);
}

template <typename T>
inline auto AbsDot(Vector2<T> v1, Vector2<T> v2) ->
    typename TupleLength<T>::type
{
  DCHECK(!v1.HasNaN() && !v2.HasNaN());
  return std::abs(Dot(v1, v2));
}

template <typename T>
inline auto LengthSquared(Vector2<T> v) ->
    typename TupleLength<T>::type
{
  return Sqr(v.x) + Sqr(v.y);
}

template <typename T>
inline auto Length(Vector2<T> v) ->
    typename TupleLength<T>::type
{
  using std::sqrt;
  return sqrt(LengthSquared(v));
}

template <typename T>
inline auto Normalize(Vector2<T> v)
{
  return v / Length(v);
}

template <typename T>
inline auto Distance(Point2<T> p1, Point2<T> p2) ->
    typename TupleLength<T>::type
{
  return Length(p1 - p2);
}

template <typename T>
inline auto DistanceSquared(Point2<T> p1, Point2<T> p2) ->
    typename TupleLength<T>::type
{
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
  Point2(T x, T y): Tuple2<pbrt::Point2, T>(x, y) {}
  template <typename U>
  explicit Point2(Point2<U> v)
      : Tuple2<pbrt::Point2, T>(T(v.x), T(v.y))
  {}
  template <typename U>
  explicit Point2(Vector2<U> v)
      : Tuple2<pbrt::Point2, T>(T(v.x), T(v.y))
  {}

  template <typename U>
  auto operator+(Vector2<U> v) const
      -> Point2<decltype(T {} + U {})>
  {
    DCHECK(!v.HasNaN());
    return {x + v.x, y + v.y};
  }
  template <typename U>
  Point2<T> &operator+=(Vector2<U> v)
  {
    DCHECK(!v.HasNaN());
    x += v.x;
    y += v.y;
    return *this;
  }
  // There is no Point -= Point because Point-Point = Vector
  // We can use Point += -Point instead

  Point2<T> operator-() const { return {-x, -y}; }

  template <typename U>
  auto operator-(Point2<U> p) const
      -> Vector2<decltype(T {} - U {})>
  {
    DCHECK(!p.HasNaN());
    return {x - p.x, y - p.y};
  }
  template <typename U>
  auto operator-(Vector2<U> v) const
      -> Point2<decltype(T {} - U {})>
  {
    DCHECK(!v.HasNaN());
    return {x - v.x, y - v.y};
  }
  template <typename U>
  Point2<T> &operator-=(Vector2<U> v)
  {
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

// Bounds2 Definition
template <typename T>
class Bounds2 {
 public:
  // Bounds2 Public Methods
  Bounds2()
  {
    T minNum = std::numeric_limits<T>::lowest();
    T maxNum = std::numeric_limits<T>::max();
    pMin = Point2<T>(maxNum, maxNum);
    pMax = Point2<T>(minNum, minNum);
  }
  explicit Bounds2(Point2<T> p): pMin(p), pMax(p) {}
  Bounds2(Point2<T> p1, Point2<T> p2)
      : pMin(Min(p1, p2)), pMax(Max(p1, p2))
  {}
  template <typename U>
  explicit Bounds2(const Bounds2<U> &b)
  {
    if (b.IsEmpty())
      // Be careful about overflowing float->int conversions and
      // the like.
      *this = Bounds2<T>();
    else {
      pMin = Point2<T>(b.pMin);
      pMax = Point2<T>(b.pMax);
    }
  }

  Vector2<T> Diagonal() const { return pMax - pMin; }

  T Area() const
  {
    Vector2<T> d = pMax - pMin;
    return d.x * d.y;
  }

  bool IsEmpty() const
  {
    return pMin.x >= pMax.x || pMin.y >= pMax.y;
  }

  bool IsDegenerate() const
  {
    return pMin.x > pMax.x || pMin.y > pMax.y;
  }

  int MaxDimension() const
  {
    Vector2<T> diag = Diagonal();
    if (diag.x > diag.y)
      return 0;
    else
      return 1;
  }
  Point2<T> operator[](int i) const
  {
    DCHECK(i == 0 || i == 1);
    return (i == 0) ? pMin : pMax;
  }
  Point2<T> &operator[](int i)
  {
    DCHECK(i == 0 || i == 1);
    return (i == 0) ? pMin : pMax;
  }
  bool operator==(const Bounds2<T> &b) const
  {
    return b.pMin == pMin && b.pMax == pMax;
  }
  bool operator!=(const Bounds2<T> &b) const
  {
    return b.pMin != pMin || b.pMax != pMax;
  }
  Point2<T> Corner(int corner) const
  {
    DCHECK(corner >= 0 && corner < 4);
    return Point2<T>(
        (*this)[(corner & 1)].x, (*this)[(corner & 2) ? 1 : 0].y
    );
  }
  Point2<T> Lerp(Point2f t) const
  {
    return Point2<T>(
        pbrt::Lerp(t.x, pMin.x, pMax.x),
        pbrt::Lerp(t.y, pMin.y, pMax.y)
    );
  }
  Vector2<T> Offset(Point2<T> p) const
  {
    Vector2<T> o = p - pMin;
    if (pMax.x > pMin.x) o.x /= pMax.x - pMin.x;
    if (pMax.y > pMin.y) o.y /= pMax.y - pMin.y;
    return o;
  }
  void BoundingSphere(Point2<T> *c, Float *rad) const
  {
    *c = (pMin + pMax) / 2;
    *rad = Inside(*c, *this) ? Distance(*c, pMax) : 0;
  }

  std::string ToString() const
  {
    return StringPrintf("[ %s - %s ]", pMin, pMax);
  }

  // Bounds2 Public Members
  Point2<T> pMin, pMax;
};

// Bounds[23][fi] Definitions
using Bounds2f = Bounds2<Float>;
using Bounds2i = Bounds2<int>;

class Bounds2iIterator : public std::forward_iterator_tag {
 public:
  Bounds2iIterator(const Bounds2i &b, const Point2i &pt)
      : p(pt), bounds(&b)
  {}
  Bounds2iIterator operator++()
  {
    advance();
    return *this;
  }
  Bounds2iIterator operator++(int)
  {
    Bounds2iIterator old = *this;
    advance();
    return old;
  }
  bool operator==(const Bounds2iIterator &bi) const
  {
    return p == bi.p && bounds == bi.bounds;
  }
  bool operator!=(const Bounds2iIterator &bi) const
  {
    return p != bi.p || bounds != bi.bounds;
  }

  Point2i operator*() const { return p; }

 private:
  void advance()
  {
    ++p.x;
    if (p.x == bounds->pMax.x) {
      p.x = bounds->pMin.x;
      ++p.y;
    }
  }
  Point2i p;
  const Bounds2i *bounds;
};

// Bounds2 Inline Functions
template <typename T>
inline Bounds2<T> Union(
    const Bounds2<T> &b1, const Bounds2<T> &b2
)
{
  // Be careful to not run the two-point Bounds constructor.
  Bounds2<T> ret;
  ret.pMin = Min(b1.pMin, b2.pMin);
  ret.pMax = Max(b1.pMax, b2.pMax);
  return ret;
}

template <typename T>
inline Bounds2<T> Intersect(
    const Bounds2<T> &b1, const Bounds2<T> &b2
)
{
  // Be careful to not run the two-point Bounds constructor.
  Bounds2<T> b;
  b.pMin = Max(b1.pMin, b2.pMin);
  b.pMax = Min(b1.pMax, b2.pMax);
  return b;
}

template <typename T>
inline bool Overlaps(const Bounds2<T> &ba, const Bounds2<T> &bb)
{
  bool x = (ba.pMax.x >= bb.pMin.x) && (ba.pMin.x <= bb.pMax.x);
  bool y = (ba.pMax.y >= bb.pMin.y) && (ba.pMin.y <= bb.pMax.y);
  return (x && y);
}

template <typename T>
inline bool Inside(Point2<T> pt, const Bounds2<T> &b)
{
  return (
      pt.x >= b.pMin.x && pt.x <= b.pMax.x &&
      pt.y >= b.pMin.y && pt.y <= b.pMax.y
  );
}

template <typename T>
inline bool Inside(const Bounds2<T> &ba, const Bounds2<T> &bb)
{
  return (
      ba.pMin.x >= bb.pMin.x && ba.pMax.x <= bb.pMax.x &&
      ba.pMin.y >= bb.pMin.y && ba.pMax.y <= bb.pMax.y
  );
}

template <typename T>
inline bool InsideExclusive(Point2<T> pt, const Bounds2<T> &b)
{
  return (
      pt.x >= b.pMin.x && pt.x < b.pMax.x && pt.y >= b.pMin.y &&
      pt.y < b.pMax.y
  );
}

template <typename T, typename U>
inline Bounds2<T> Expand(const Bounds2<T> &b, U delta)
{
  Bounds2<T> ret;
  ret.pMin = b.pMin - Vector2<T>(delta, delta);
  ret.pMax = b.pMax + Vector2<T>(delta, delta);
  return ret;
}

// 3D
// Tuple3 Definition
template <template <typename> class Child, typename T>
class Tuple3 {
 public:
  // Tuple3 Public Methods
  Tuple3() = default;
  Tuple3(T x): x(x), y(x), z(x) { DCHECK(!HasNaN()); }
  Tuple3(T x, T y, T z): x(x), y(y), z(z) { DCHECK(!HasNaN()); }

  bool HasNaN() const
  {
    return IsNaN(x) || IsNaN(y) || IsNaN(z);
  }

  T operator[](int i) const
  {
    DCHECK(i >= 0 && i <= 2);
    if (i == 0) return x;
    if (i == 1) return y;
    return z;
  }


  T &operator[](int i)
  {
    DCHECK(i >= 0 && i <= 2);
    if (i == 0) return x;
    if (i == 1) return y;
    return z;
  }

  template <typename U>
  auto operator+(Child<U> c) const
      -> Child<decltype(T {} + U {})>
  {
    DCHECK(!c.HasNaN());
    return {x + c.x, y + c.y, z + c.z};
  }

  static const int nDimensions = 3;

  Tuple3(Child<T> c)
  {
    DCHECK(!c.HasNaN());
    x = c.x;
    y = c.y;
    z = c.z;
  }

  Child<T> &operator=(Child<T> c)
  {
    DCHECK(!c.HasNaN());
    x = c.x;
    y = c.y;
    z = c.z;
    return static_cast<Child<T> &>(*this);
  }

  template <typename U>
  Child<T> &operator+=(Child<U> c)
  {
    DCHECK(!c.HasNaN());
    x += c.x;
    y += c.y;
    z += c.z;
    return static_cast<Child<T> &>(*this);
  }

  template <typename U>
  auto operator-(Child<U> c) const
      -> Child<decltype(T {} - U {})>
  {
    DCHECK(!c.HasNaN());
    return {x - c.x, y - c.y, z - c.z};
  }
  template <typename U>
  Child<T> &operator-=(Child<U> c)
  {
    DCHECK(!c.HasNaN());
    x -= c.x;
    y -= c.y;
    z -= c.z;
    return static_cast<Child<T> &>(*this);
  }


  bool operator==(Child<T> c) const
  {
    return x == c.x && y == c.y && z == c.z;
  }

  bool operator!=(Child<T> c) const
  {
    return x != c.x || y != c.y || z != c.z;
  }

  template <typename U>
  auto operator*(U s) const -> Child<decltype(T {} * U {})>
  {
    return {s * x, s * y, s * z};
  }
  template <typename U>
  Child<T> &operator*=(U s)
  {
    DCHECK(!IsNaN(s));
    x *= s;
    y *= s;
    z *= s;
    return static_cast<Child<T> &>(*this);
  }

  template <typename U>
  auto operator/(U d) const -> Child<decltype(T {} / U {})>
  {
    DCHECK_NE(d, 0);
    return {x / d, y / d, z / d};
  }
  template <typename U>
  Child<T> &operator/=(U d)
  {
    DCHECK_NE(d, 0);
    x /= d;
    y /= d;
    z /= d;
    return static_cast<Child<T> &>(*this);
  }

  Child<T> operator-() const { return {-x, -y, -z}; }

  friend std::ostream &operator<<(
      std::ostream &os, const Child<T> &obj
  )
  {
    os << "[ " << obj.x << ", " << obj.y << ", " << obj.z
       << " ]";
    return os;
  }

  // Tuple3 Public Members
  T x {}, y {}, z {};
};

// Tuple3 Inline Functions
template <template <class> class C, typename T, typename U>
inline auto operator*(U s, Tuple3<C, T> t)
    -> C<decltype(T {} * U {})>
{
  return t * s;
}

template <template <class> class C, typename T>
inline C<T> Abs(Tuple3<C, T> t)
{
  using std::abs;
  return {abs(t.x), abs(t.y), abs(t.z)};
}

template <template <class> class C, typename T>
inline C<T> Ceil(Tuple3<C, T> t)
{
  using std::ceil;
  return {ceil(t.x), ceil(t.y), ceil(t.z)};
}

template <template <class> class C, typename T>
inline C<T> Floor(Tuple3<C, T> t)
{
  using std::floor;
  return {floor(t.x), floor(t.y), floor(t.z)};
}

template <template <class> class C, typename T>
inline auto Lerp(Float t, Tuple3<C, T> t0, Tuple3<C, T> t1)
{
  return (1 - t) * t0 + t * t1;
}

template <template <class> class C, typename T>
inline C<T> FMA(Float a, Tuple3<C, T> b, Tuple3<C, T> c)
{
  return {FMA(a, b.x, c.x), FMA(a, b.y, c.y), FMA(a, b.z, c.z)};
}

template <template <class> class C, typename T>
inline C<T> FMA(Tuple3<C, T> a, Float b, Tuple3<C, T> c)
{
  return FMA(b, a, c);
}

template <template <class> class C, typename T>
inline C<T> Min(Tuple3<C, T> t1, Tuple3<C, T> t2)
{
  using std::min;
  return {min(t1.x, t2.x), min(t1.y, t2.y), min(t1.z, t2.z)};
}

template <template <class> class C, typename T>
inline T MinComponentValue(Tuple3<C, T> t)
{
  using std::min;
  return min({t.x, t.y, t.z});
}

template <template <class> class C, typename T>
inline int MinComponentIndex(Tuple3<C, T> t)
{
  return (t.x < t.y) ? ((t.x < t.z) ? 0 : 2)
                     : ((t.y < t.z) ? 1 : 2);
}

template <template <class> class C, typename T>
inline C<T> Max(Tuple3<C, T> t1, Tuple3<C, T> t2)
{
  using std::max;
  return {max(t1.x, t2.x), max(t1.y, t2.y), max(t1.z, t2.z)};
}

template <template <class> class C, typename T>
inline C<T> ClampMax(Tuple3<C, T> t, T val)
{
  using std::min;
  return {min(t.x, val), min(t.y, val), min(t.z, val)};
}

template <template <class> class C, typename T>
inline T MaxComponentValue(Tuple3<C, T> t)
{
  using std::max;
  return max({t.x, t.y, t.z});
}

template <template <class> class C, typename T>
inline int MaxComponentIndex(Tuple3<C, T> t)
{
  return (t.x > t.y) ? ((t.x > t.z) ? 0 : 2)
                     : ((t.y > t.z) ? 1 : 2);
}

template <template <class> class C, typename T>
inline C<T> Permute(Tuple3<C, T> t, std::array<int, 3> p)
{
  return {t[p[0]], t[p[1]], t[p[2]]};
}

template <template <class> class C, typename T>
inline T HProd(Tuple3<C, T> t)
{
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

  Vector3() { CREATED("Vector3"); }
  Vector3(T x): Tuple3<pbrt::Vector3, T>(x)
  {
    CREATED("Vector3");
  }
  Vector3(T x, T y, T z): Tuple3<pbrt::Vector3, T>(x, y, z)
  {
    CREATED("Vector3");
  }

  template <typename U>
  explicit Vector3(Vector3<U> v)
      : Tuple3<pbrt::Vector3, T>(T(v.x), T(v.y), T(v.z))
  {
    CREATED("Vector3");
  }

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
Vector3<T>::Vector3(Point3<U> p)
    : Tuple3<pbrt::Vector3, T>(T(p.x), T(p.y), T(p.z))
{}

template <typename T>
inline Vector3<T> HorizontalProduct(
    Vector3<T> v1, Vector3<T> v2
)
{
  return {v1.x * v2.x, v1.y * v2.y, v1.z * v2.z};
}

template <typename T>
inline Vector3<T> Cross(Vector3<T> v1, Normal3<T> v2)
{
  DCHECK(!v1.HasNaN() && !v2.HasNaN());
  return {
      DifferenceOfProducts(v1.y, v2.z, v1.z, v2.y),
      DifferenceOfProducts(v1.z, v2.x, v1.x, v2.z),
      DifferenceOfProducts(v1.x, v2.y, v1.y, v2.x)
  };
}

template <typename T>
inline Vector3<T> Cross(Normal3<T> v1, Vector3<T> v2)
{
  DCHECK(!v1.HasNaN() && !v2.HasNaN());
  return {
      DifferenceOfProducts(v1.y, v2.z, v1.z, v2.y),
      DifferenceOfProducts(v1.z, v2.x, v1.x, v2.z),
      DifferenceOfProducts(v1.x, v2.y, v1.y, v2.x)
  };
}

template <typename T>
inline T LengthSquared(Vector3<T> v)
{
  return Sqr(v.x) + Sqr(v.y) + Sqr(v.z);
}

template <typename T>
inline auto Length(Vector3<T> v) ->
    typename TupleLength<T>::type
{
  using std::sqrt;
  return sqrt(LengthSquared(v));
}

template <typename T>
inline auto Normalize(Vector3<T> v)
{
  return v / Length(v);
}

template <typename T>
inline T Dot(Vector3<T> v, Vector3<T> w)
{
  DCHECK(!v.HasNaN() && !w.HasNaN());
  return v.x * w.x + v.y * w.y + v.z * w.z;
}

// Equivalent to std::acos(Dot(a, b)), but more numerically
// stable. via http://www.plunk.org/~hatch/rightway.html
template <typename T>
inline Float AngleBetween(Vector3<T> v1, Vector3<T> v2)
{
  if (Dot(v1, v2) < 0)
    return Pi - 2 * SafeASin(Length(v1 + v2) / 2);
  else
    return 2 * SafeASin(Length(v2 - v1) / 2);
}

template <typename T>
inline T AbsDot(Vector3<T> v1, Vector3<T> v2)
{
  DCHECK(!v1.HasNaN() && !v2.HasNaN());
  return std::abs(Dot(v1, v2));
}

template <typename T>
inline Float AngleBetween(Normal3<T> a, Normal3<T> b)
{
  if (Dot(a, b) < 0)
    return Pi - 2 * SafeASin(Length(a + b) / 2);
  else
    return 2 * SafeASin(Length(b - a) / 2);
}

template <typename T>
inline Vector3<T> GramSchmidt(Vector3<T> v, Vector3<T> w)
{
  return v - Dot(v, w) * w;
}

template <typename T>
inline Vector3<T> Cross(Vector3<T> v, Vector3<T> w)
{
  DCHECK(!v.HasNaN() && !w.HasNaN());
  return {
      DifferenceOfProducts(v.y, w.z, v.z, w.y),
      DifferenceOfProducts(v.z, w.x, v.x, w.z),
      DifferenceOfProducts(v.x, w.y, v.y, w.x)
  };
}

template <typename T>
inline void CoordinateSystem(
    Vector3<T> v1, Vector3<T> *v2, Vector3<T> *v3
)
{
  Float sign = std::copysign(Float(1), v1.z);
  Float a = -1 / (sign + v1.z);
  Float b = v1.x * v1.y * a;
  *v2 = Vector3<T>(
      1 + sign * Sqr(v1.x) * a, sign * b, -sign * v1.x
  );
  *v3 = Vector3<T>(b, sign + Sqr(v1.y) * a, -v1.y);
}

template <typename T>
inline void CoordinateSystem(
    Normal3<T> v1, Vector3<T> *v2, Vector3<T> *v3
)
{
  Float sign = std::copysign(Float(1), v1.z);
  Float a = -1 / (sign + v1.z);
  Float b = v1.x * v1.y * a;
  *v2 = Vector3<T>(
      1 + sign * Sqr(v1.x) * a, sign * b, -sign * v1.x
  );
  *v3 = Vector3<T>(b, sign + Sqr(v1.y) * a, -v1.y);
}

template <typename T>
template <typename U>
Vector3<T>::Vector3(Normal3<U> n)
    : Tuple3<pbrt::Vector3, T>(T(n.x), T(n.y), T(n.z))
{}

template <typename T>
inline Vector3<T> Reflect(
    const Vector3<T> &v, const Normal3<T> &n
)
{
  return v - Vector3<T>(2 * Dot(v, n) * n);
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
  Point3(T x): Tuple3<pbrt::Point3, T>(x) {}
  Point3(T x, T y, T z): Tuple3<pbrt::Point3, T>(x, y, z) {}

  // We can't do using operator- above, since we don't want to
  // pull in the Point-Point -> Point one so that we can return
  // a vector instead...

  Point3<T> operator-() const { return {-x, -y, -z}; }

  template <typename U>
  explicit Point3(Point3<U> p)
      : Tuple3<pbrt::Point3, T>(T(p.x), T(p.y), T(p.z))
  {}
  template <typename U>
  explicit Point3(Vector3<U> v)
      : Tuple3<pbrt::Point3, T>(T(v.x), T(v.y), T(v.z))
  {}

  template <typename U>
  auto operator+(Vector3<U> v) const
      -> Point3<decltype(T {} + U {})>
  {
    DCHECK(!v.HasNaN());
    return {x + v.x, y + v.y, z + v.z};
  }
  template <typename U>
  Point3<T> &operator+=(Vector3<U> v)
  {
    DCHECK(!v.HasNaN());
    x += v.x;
    y += v.y;
    z += v.z;
    return *this;
  }

  template <typename U>
  auto operator-(Vector3<U> v) const
      -> Point3<decltype(T {} - U {})>
  {
    DCHECK(!v.HasNaN());
    return {x - v.x, y - v.y, z - v.z};
  }
  template <typename U>
  Point3<T> &operator-=(Vector3<U> v)
  {
    DCHECK(!v.HasNaN());
    x -= v.x;
    y -= v.y;
    z -= v.z;
    return *this;
  }

  template <typename U>
  auto operator-(Point3<U> p) const
      -> Vector3<decltype(T {} - U {})>
  {
    DCHECK(!p.HasNaN());
    return {x - p.x, y - p.y, z - p.z};
  }
};

// Point3* Definitions
using Point3f = Point3<Float>;
using Point3i = Point3<int>;

// Point3 Inline Functions
template <typename T>
inline auto Distance(Point3<T> p1, Point3<T> p2)
{
  return Length(p1 - p2);
}

template <typename T>
inline auto DistanceSquared(Point3<T> p1, Point3<T> p2)
{
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
  Normal3(T x): Tuple3<pbrt::Normal3, T>(x) {}
  Normal3(T x, T y, T z): Tuple3<pbrt::Normal3, T>(x, y, z) {}
  template <typename U>
  explicit Normal3(Normal3<U> v)
      : Tuple3<pbrt::Normal3, T>(T(v.x), T(v.y), T(v.z))
  {}

  template <typename U>
  explicit Normal3(Vector3<U> v)
      : Tuple3<pbrt::Normal3, T>(T(v.x), T(v.y), T(v.z))
  {}
};

using Normal3f = Normal3<Float>;

// Normal3 Inline Functions
template <typename T>
inline auto LengthSquared(Normal3<T> n) ->
    typename TupleLength<T>::type
{
  return Sqr(n.x) + Sqr(n.y) + Sqr(n.z);
}

template <typename T>
inline auto Length(Normal3<T> n) ->
    typename TupleLength<T>::type
{
  using std::sqrt;
  return sqrt(LengthSquared(n));
}

template <typename T>
inline auto Normalize(Normal3<T> n)
{
  return n / Length(n);
}

template <typename T>
inline auto Dot(Normal3<T> n, Vector3<T> v) ->
    typename TupleLength<T>::type
{
  DCHECK(!n.HasNaN() && !v.HasNaN());
  // return FMA(n.x, v.x, SumOfProducts(n.y, v.y, n.z, v.z));
  return n.x * v.x + n.y * v.y + n.z * v.z;
}

template <typename T>
inline auto Dot(Vector3<T> v, Normal3<T> n) ->
    typename TupleLength<T>::type
{
  DCHECK(!v.HasNaN() && !n.HasNaN());
  // return FMA(n.x, v.x, SumOfProducts(n.y, v.y, n.z, v.z));
  return n.x * v.x + n.y * v.y + n.z * v.z;
}

template <typename T>
inline auto Dot(Normal3<T> n1, Normal3<T> n2) ->
    typename TupleLength<T>::type
{
  DCHECK(!n1.HasNaN() && !n2.HasNaN());
  // return FMA(n1.x, n2.x, SumOfProducts(n1.y, n2.y, n1.z, n2.z));
  return n1.x * n2.x + n1.y * n2.y + n1.z * n2.z;
}

template <typename T>
inline auto AbsDot(Normal3<T> n, Vector3<T> v) ->
    typename TupleLength<T>::type
{
  DCHECK(!n.HasNaN() && !v.HasNaN());
  return std::abs(Dot(n, v));
}

template <typename T>
inline auto AbsDot(Vector3<T> v, Normal3<T> n) ->
    typename TupleLength<T>::type
{
  using std::abs;
  DCHECK(!v.HasNaN() && !n.HasNaN());
  return abs(Dot(v, n));
}

template <typename T>
inline auto AbsDot(Normal3<T> n1, Normal3<T> n2) ->
    typename TupleLength<T>::type
{
  using std::abs;
  DCHECK(!n1.HasNaN() && !n2.HasNaN());
  return abs(Dot(n1, n2));
}

template <typename T>
inline Normal3<T> FaceForward(Normal3<T> n, Vector3<T> v)
{
  return (Dot(n, v) < 0.f) ? -n : n;
}

template <typename T>
inline Normal3<T> FaceForward(Normal3<T> n, Normal3<T> n2)
{
  return (Dot(n, n2) < 0.f) ? -n : n;
}

template <typename T>
inline Vector3<T> FaceForward(Vector3<T> v, Vector3<T> v2)
{
  return (Dot(v, v2) < 0.f) ? -v : v;
}

template <typename T>
inline Vector3<T> FaceForward(Vector3<T> v, Normal3<T> n2)
{
  return (Dot(v, n2) < 0.f) ? -v : v;
}

// Bounds3 Definition
template <typename T>
class Bounds3 {
 public:
  // Bounds3 Public Methods
  Bounds3()
  {
    T minNum = std::numeric_limits<T>::lowest();
    T maxNum = std::numeric_limits<T>::max();
    pMin = Point3<T>(maxNum, maxNum, maxNum);
    pMax = Point3<T>(minNum, minNum, minNum);
  }

  explicit Bounds3(Point3<T> p): pMin(p), pMax(p) {}

  Bounds3(Point3<T> p1, Point3<T> p2)
      : pMin(Min(p1, p2)), pMax(Max(p1, p2))
  {}

  Point3<T> operator[](int i) const
  {
    DCHECK(i == 0 || i == 1);
    return (i == 0) ? pMin : pMax;
  }
  Point3<T> &operator[](int i)
  {
    DCHECK(i == 0 || i == 1);
    return (i == 0) ? pMin : pMax;
  }

  Point3<T> Corner(int corner) const
  {
    DCHECK(corner >= 0 && corner < 8);
    return Point3<T>(
        (*this)[(corner & 1)].x,
        (*this)[(corner & 2) ? 1 : 0].y,
        (*this)[(corner & 4) ? 1 : 0].z
    );
  }

  Vector3<T> Diagonal() const { return pMax - pMin; }

  T SurfaceArea() const
  {
    Vector3<T> d = Diagonal();
    return 2 * (d.x * d.y + d.x * d.z + d.y * d.z);
  }

  T Volume() const
  {
    Vector3<T> d = Diagonal();
    return d.x * d.y * d.z;
  }

  int MaxDimension() const
  {
    Vector3<T> d = Diagonal();
    if (d.x > d.y && d.x > d.z)
      return 0;
    else if (d.y > d.z)
      return 1;
    else
      return 2;
  }

  Point3f Lerp(Point3f t) const
  {
    return Point3f(
        pbrt::Lerp(t.x, pMin.x, pMax.x),
        pbrt::Lerp(t.y, pMin.y, pMax.y),
        pbrt::Lerp(t.z, pMin.z, pMax.z)
    );
  }

  Vector3f Offset(Point3f p) const
  {
    Vector3f o = p - pMin;
    if (pMax.x > pMin.x) o.x /= pMax.x - pMin.x;
    if (pMax.y > pMin.y) o.y /= pMax.y - pMin.y;
    if (pMax.z > pMin.z) o.z /= pMax.z - pMin.z;
    return o;
  }

  void BoundingSphere(Point3<T> *center, Float *radius) const
  {
    *center = (pMin + pMax) / 2;
    *radius =
        Inside(*center, *this) ? Distance(*center, pMax) : 0;
  }

  bool IsEmpty() const
  {
    return pMin.x >= pMax.x || pMin.y >= pMax.y ||
           pMin.z >= pMax.z;
  }
  bool IsDegenerate() const
  {
    return pMin.x > pMax.x || pMin.y > pMax.y ||
           pMin.z > pMax.z;
  }

  template <typename U>
  explicit Bounds3(const Bounds3<U> &b)
  {
    if (b.IsEmpty())
      // Be careful about overflowing float->int conversions and
      // the like.
      *this = Bounds3<T>();
    else {
      pMin = Point3<T>(b.pMin);
      pMax = Point3<T>(b.pMax);
    }
  }
  bool operator==(const Bounds3<T> &b) const
  {
    return b.pMin == pMin && b.pMax == pMax;
  }
  bool operator!=(const Bounds3<T> &b) const
  {
    return b.pMin != pMin || b.pMax != pMax;
  }
  bool IntersectP(
      Point3f o, Vector3f d, Float tMax = Infinity,
      Float *hitt0 = nullptr, Float *hitt1 = nullptr
  ) const;
  bool IntersectP(
      Point3f o, Vector3f d, Float tMax, Vector3f invDir,
      const int dirIsNeg[3]
  ) const;

  std::string ToString() const
  {
    return StringPrintf("[ %s - %s ]", pMin, pMax);
  }

  // Bounds3 Public Members
  Point3<T> pMin, pMax;
};


using Bounds3f = Bounds3<Float>;
using Bounds3i = Bounds3<int>;

// Bounds3 Inline Functions
template <typename T>
inline Bounds3<T> Union(const Bounds3<T> &b, Point3<T> p)
{
  Bounds3<T> ret;
  ret.pMin = Min(b.pMin, p);
  ret.pMax = Max(b.pMax, p);
  return ret;
}

template <typename T>
inline Bounds3<T> Union(
    const Bounds3<T> &b1, const Bounds3<T> &b2
)
{
  Bounds3<T> ret;
  ret.pMin = Min(b1.pMin, b2.pMin);
  ret.pMax = Max(b1.pMax, b2.pMax);
  return ret;
}

template <typename T>
inline Bounds3<T> Intersect(
    const Bounds3<T> &b1, const Bounds3<T> &b2
)
{
  Bounds3<T> b;
  b.pMin = Max(b1.pMin, b2.pMin);
  b.pMax = Min(b1.pMax, b2.pMax);
  return b;
}

template <typename T>
inline bool Overlaps(const Bounds3<T> &b1, const Bounds3<T> &b2)
{
  bool x = (b1.pMax.x >= b2.pMin.x) && (b1.pMin.x <= b2.pMax.x);
  bool y = (b1.pMax.y >= b2.pMin.y) && (b1.pMin.y <= b2.pMax.y);
  bool z = (b1.pMax.z >= b2.pMin.z) && (b1.pMin.z <= b2.pMax.z);
  return (x && y && z);
}

template <typename T>
inline bool Inside(Point3<T> p, const Bounds3<T> &b)
{
  return (
      p.x >= b.pMin.x && p.x <= b.pMax.x && p.y >= b.pMin.y &&
      p.y <= b.pMax.y && p.z >= b.pMin.z && p.z <= b.pMax.z
  );
}

template <typename T>
inline bool InsideExclusive(Point3<T> p, const Bounds3<T> &b)
{
  return (
      p.x >= b.pMin.x && p.x < b.pMax.x && p.y >= b.pMin.y &&
      p.y < b.pMax.y && p.z >= b.pMin.z && p.z < b.pMax.z
  );
}

template <typename T, typename U>
inline auto DistanceSquared(Point3<T> p, const Bounds3<U> &b)
{
  using TDist = decltype(T {} - U {});
  TDist dx =
      std::max<TDist>({0, b.pMin.x - p.x, p.x - b.pMax.x});
  TDist dy =
      std::max<TDist>({0, b.pMin.y - p.y, p.y - b.pMax.y});
  TDist dz =
      std::max<TDist>({0, b.pMin.z - p.z, p.z - b.pMax.z});
  return Sqr(dx) + Sqr(dy) + Sqr(dz);
}

template <typename T, typename U>
inline auto Distance(Point3<T> p, const Bounds3<U> &b)
{
  auto dist2 = DistanceSquared(p, b);
  using TDist = typename TupleLength<decltype(dist2)>::type;
  return std::sqrt(TDist(dist2));
}

template <typename T, typename U>
inline Bounds3<T> Expand(const Bounds3<T> &b, U delta)
{
  Bounds3<T> ret;
  ret.pMin = b.pMin - Vector3<T>(delta, delta, delta);
  ret.pMax = b.pMax + Vector3<T>(delta, delta, delta);
  return ret;
}

template <typename T>
inline bool Bounds3<T>::IntersectP(
    Point3f o, Vector3f d, Float tMax, Float *hitt0,
    Float *hitt1
) const
{
  Float t0 = 0, t1 = tMax;
  for (int i = 0; i < 3; ++i) {
    // Update interval for _i_th bounding box slab
    Float invRayDir = 1 / d[i];
    Float tNear = (pMin[i] - o[i]) * invRayDir;
    Float tFar = (pMax[i] - o[i]) * invRayDir;
    // Update parametric interval from slab intersection $t$ values
    if (tNear > tFar) std::swap(tNear, tFar);
    // Update _tFar_ to ensure robust ray--bounds intersection
    // tFar *= 1 + 2 * gamma(3); // TODO: What is this for?  it
    // creates a bug

    t0 = tNear > t0 ? tNear : t0;
    t1 = tFar < t1 ? tFar : t1;
    if (t0 > t1) return false;
  }
  if (hitt0) *hitt0 = t0;
  if (hitt1) *hitt1 = t1;
  return true;
}

template <typename T>
inline bool Bounds3<T>::IntersectP(
    Point3f o, Vector3f d, Float raytMax, Vector3f invDir,
    const int dirIsNeg[3]
) const
{
  const Bounds3f &bounds = *this;
  // Check for ray intersection against $x$ and $y$ slabs
  Float tMin = (bounds[dirIsNeg[0]].x - o.x) * invDir.x;
  Float tMax = (bounds[1 - dirIsNeg[0]].x - o.x) * invDir.x;
  Float tyMin = (bounds[dirIsNeg[1]].y - o.y) * invDir.y;
  Float tyMax = (bounds[1 - dirIsNeg[1]].y - o.y) * invDir.y;
  // Update _tMax_ and _tyMax_ to ensure robust bounds intersection
  tMax *= 1 + 2 * gamma(3);
  tyMax *= 1 + 2 * gamma(3);

  if (tMin > tyMax || tyMin > tMax) return false;
  if (tyMin > tMin) tMin = tyMin;
  if (tyMax < tMax) tMax = tyMax;

  // Check for ray intersection against $z$ slab
  Float tzMin = (bounds[dirIsNeg[2]].z - o.z) * invDir.z;
  Float tzMax = (bounds[1 - dirIsNeg[2]].z - o.z) * invDir.z;
  // Update _tzMax_ to ensure robust bounds intersection
  tzMax *= 1 + 2 * gamma(3);

  if (tMin > tzMax || tzMin > tMax) return false;
  if (tzMin > tMin) tMin = tzMin;
  if (tzMax < tMax) tMax = tzMax;

  return (tMin < raytMax) && (tMax > 0);
}

inline Bounds2iIterator begin(const Bounds2i &b)
{
  return Bounds2iIterator(b, b.pMin);
}

inline Bounds2iIterator end(const Bounds2i &b)
{
  // Normally, the ending point is at the minimum x value and
  // one past the last valid y value.
  Point2i pEnd(b.pMin.x, b.pMax.y);
  // However, if the bounds are degenerate, override the end
  // point to equal the start point so that any attempt to
  // iterate over the bounds exits out immediately.
  if (b.pMin.x >= b.pMax.x || b.pMin.y >= b.pMax.y)
    pEnd = b.pMin;
  return Bounds2iIterator(b, pEnd);
}

template <typename T>
inline Bounds2<T> Union(const Bounds2<T> &b, Point2<T> p)
{
  // Be careful to not run the two-point Bounds constructor.
  Bounds2<T> ret;
  ret.pMin = Min(b.pMin, p);
  ret.pMax = Max(b.pMax, p);
  return ret;
}


// Quaternion Definition
class Quaternion {
 public:
  // Quaternion Public Methods
  Quaternion() = default;
  Quaternion(Vector3f v, Float w): v(v), w(w) {}


  Quaternion &operator+=(Quaternion q)
  {
    v += q.v;
    w += q.w;
    return *this;
  }


  Quaternion operator+(Quaternion q) const
  {
    return {v + q.v, w + q.w};
  }

  Quaternion &operator-=(Quaternion q)
  {
    v -= q.v;
    w -= q.w;
    return *this;
  }

  Quaternion operator-() const { return {-v, -w}; }

  Quaternion operator-(Quaternion q) const
  {
    return {v - q.v, w - q.w};
  }

  Quaternion &operator*=(Float f)
  {
    v *= f;
    w *= f;
    return *this;
  }

  Quaternion operator*(Float f) const { return {v * f, w * f}; }

  Quaternion &operator/=(Float f)
  {
    DCHECK_NE(0, f);
    v /= f;
    w /= f;
    return *this;
  }

  Quaternion operator/(Float f) const
  {
    DCHECK_NE(0, f);
    return {v / f, w / f};
  }

  std::string ToString() const;

  // Quaternion Public Members
  Vector3f v;
  Float w = 1;
};

// Quaternion Inline Functions

inline Quaternion operator*(Float f, Quaternion q)
{
  return q * f;
}

inline Float Dot(Quaternion q1, Quaternion q2)
{
  return Dot(q1.v, q2.v) + q1.w * q2.w;
}

inline Float Length(Quaternion q)
{
  return std::sqrt(Dot(q, q));
}
inline Quaternion Normalize(Quaternion q)
{
  DCHECK_GT(Length(q), 0);
  return q / Length(q);
}

inline Float AngleBetween(Quaternion q1, Quaternion q2)
{
  if (Dot(q1, q2) < 0)
    return Pi - 2 * SafeASin(Length(q1 + q2) / 2);
  else
    return 2 * SafeASin(Length(q2 - q1) / 2);
}

// http://www.plunk.org/~hatch/rightway.html
inline Quaternion Slerp(Float t, Quaternion q1, Quaternion q2)
{
  Float theta = AngleBetween(q1, q2);
  Float sinThetaOverTheta = SinXOverX(theta);
  return q1 * (1 - t) * SinXOverX((1 - t) * theta) /
             sinThetaOverTheta +
         q2 * t * SinXOverX(t * theta) / sinThetaOverTheta;
}

// Frame Definition
class Frame {
 public:
  // Frame Public Methods

  Frame(): x(1, 0, 0), y(0, 1, 0), z(0, 0, 1) {}

  Frame(Vector3f x, Vector3f y, Vector3f z);


  static Frame FromXZ(Vector3f x, Vector3f z)
  {
    return Frame(x, Cross(z, x), z);
  }

  static Frame FromXY(Vector3f x, Vector3f y)
  {
    return Frame(x, y, Cross(x, y));
  }


  static Frame FromZ(Vector3f z)
  {
    Vector3f x, y;
    CoordinateSystem(z, &x, &y);
    return Frame(x, y, z);
  }


  static Frame FromX(Vector3f x)
  {
    Vector3f y, z;
    CoordinateSystem(x, &y, &z);
    return Frame(x, y, z);
  }


  static Frame FromY(Vector3f y)
  {
    Vector3f x, z;
    CoordinateSystem(y, &z, &x);
    return Frame(x, y, z);
  }


  static Frame FromX(Normal3f x)
  {
    Vector3f y, z;
    CoordinateSystem(x, &y, &z);
    return Frame(Vector3f(x), y, z);
  }


  static Frame FromY(Normal3f y)
  {
    Vector3f x, z;
    CoordinateSystem(y, &z, &x);
    return Frame(x, Vector3f(y), z);
  }


  static Frame FromZ(Normal3f z) { return FromZ(Vector3f(z)); }


  Vector3f ToLocal(Vector3f v) const
  {
    return Vector3f(Dot(v, x), Dot(v, y), Dot(v, z));
  }


  Normal3f ToLocal(Normal3f n) const
  {
    return Normal3f(Dot(n, x), Dot(n, y), Dot(n, z));
  }


  Vector3f FromLocal(Vector3f v) const
  {
    return v.x * x + v.y * y + v.z * z;
  }


  Normal3f FromLocal(Normal3f n) const
  {
    return Normal3f(n.x * x + n.y * y + n.z * z);
  }

  // Frame Public Members
  Vector3f x, y, z;
};

// Frame Inline Functions
inline Frame::Frame(Vector3f x, Vector3f y, Vector3f z)
    : x(x), y(y), z(z)
{
  DCHECK_LT(std::abs(LengthSquared(x) - 1), 1e-4);
  DCHECK_LT(std::abs(LengthSquared(y) - 1), 1e-4);
  DCHECK_LT(std::abs(LengthSquared(z) - 1), 1e-4);
  DCHECK_LT(std::abs(Dot(x, y)), 1e-4);
  DCHECK_LT(std::abs(Dot(y, z)), 1e-4);
  DCHECK_LT(std::abs(Dot(z, x)), 1e-4);
}
}  // namespace pbrt
