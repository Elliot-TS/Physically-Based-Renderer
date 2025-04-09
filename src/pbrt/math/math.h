// Mostly copied from PBRT
#pragma once
#include <cmath>
#include <limits>
#include <optional>
#include <span>
#include <string>
#include "pbrt/test/check.h"


namespace pbrt {
/** Float **/
// FIXME: Things don't work if this is set to float
#define Float double

// Floating-point Constants
static constexpr Float Infinity =
    std::numeric_limits<Float>::infinity();

/** Mathematical Constants **/
constexpr Float ShadowEpsilon = 0.0001f;

constexpr Float Pi = 3.14159265358979323846;
constexpr Float InvPi = 0.31830988618379067154;
constexpr Float Inv2Pi = 0.15915494309189533577;
constexpr Float Inv4Pi = 0.07957747154594766788;
constexpr Float PiOver2 = 1.57079632679489661923;
constexpr Float PiOver4 = 0.78539816339744830961;
constexpr Float Sqrt2 = 1.41421356237309504880;

/** Math Functions **/

template <typename T>
inline T FMA(T a, T b, T c)
{
  return a * b + c;
}

template <typename T>
inline constexpr T Sqr(T x)
{
  return x * x;
}

inline float SafeSqrt(float x)
{
  DCHECK_GE(x, -1e-3f);  // not too negative
  return std::sqrt(std::max(0.f, x));
}

inline Float Radians(Float deg)
{
  return (Pi / 180) * deg;
}

// TODO: PBRT used intervals for error correction
template <typename T>
inline T DifferenceOfProducts(T a, T b, T c, T d)
{
  return a * b - c * d;
}

// InnerProduct Helper Functions
inline Float InnerProduct(Float a, Float b)
{
  return a * b;
}

template <typename... T>
inline Float InnerProduct(Float a, Float b, T... terms)
{
  Float ab = a * b;
  Float tp = InnerProduct(terms...);
  Float sum = ab + tp;
  return sum;
}

// SquareMatrix
template <int N>
class SquareMatrix;

template <int N>
std::optional<SquareMatrix<N>> LinearLeastSquares(
    const Float A[][N], const Float B[][N], int rows
);

template <int N>
std::optional<SquareMatrix<N>> LinearLeastSquares(
    const Float A[][N], const Float B[][N], int rows
)
{
  SquareMatrix<N> AtA = SquareMatrix<N>::Zero();
  SquareMatrix<N> AtB = SquareMatrix<N>::Zero();

  for (int i = 0; i < N; ++i)
    for (int j = 0; j < N; ++j)
      for (int r = 0; r < rows; ++r) {
        AtA[i][j] += A[r][i] * A[r][j];
        AtB[i][j] += A[r][i] * B[r][j];
      }

  auto AtAi = Inverse(AtA);
  if (!AtAi) return {};
  return Transpose(*AtAi * AtB);
}

// SquareMatrix Definition
template <int N>
class SquareMatrix {
 public:
  // SquareMatrix Public Methods
  static SquareMatrix Zero()
  {
    SquareMatrix m;
    for (int i = 0; i < N; ++i)
      for (int j = 0; j < N; ++j) m.m[i][j] = 0;
    return m;
  }

  SquareMatrix()
  {
    for (int i = 0; i < N; ++i)
      for (int j = 0; j < N; ++j) m[i][j] = (i == j) ? 1 : 0;
  }
  SquareMatrix(const Float mat[N][N])
  {
    for (int i = 0; i < N; ++i)
      for (int j = 0; j < N; ++j) m[i][j] = mat[i][j];
  }

  SquareMatrix(std::span<const Float> t);
  template <typename... Args>
  SquareMatrix(Float v, Args... args)
  {
    static_assert(
        1 + sizeof...(Args) == N * N,
        "Incorrect number of values provided to SquareMatrix constructor"
    );
    init<N>(m, 0, 0, v, args...);
  }
  template <typename... Args>
  static SquareMatrix Diag(Float v, Args... args)
  {
    static_assert(
        1 + sizeof...(Args) == N,
        "Incorrect number of values provided to SquareMatrix::Diag"
    );
    SquareMatrix m;
    initDiag<N>(m.m, 0, v, args...);
    return m;
  }

  SquareMatrix operator+(const SquareMatrix &m) const
  {
    SquareMatrix r = *this;
    for (int i = 0; i < N; ++i)
      for (int j = 0; j < N; ++j) r.m[i][j] += m.m[i][j];
    return r;
  }

  SquareMatrix operator*(Float s) const
  {
    SquareMatrix r = *this;
    for (int i = 0; i < N; ++i)
      for (int j = 0; j < N; ++j) r.m[i][j] *= s;
    return r;
  }
  SquareMatrix operator/(Float s) const
  {
    DCHECK_NE(s, 0);
    SquareMatrix r = *this;
    for (int i = 0; i < N; ++i)
      for (int j = 0; j < N; ++j) r.m[i][j] /= s;
    return r;
  }

  bool operator==(const SquareMatrix<N> &m2) const
  {
    for (int i = 0; i < N; ++i)
      for (int j = 0; j < N; ++j)
        if (m[i][j] != m2.m[i][j]) return false;
    return true;
  }


  bool operator!=(const SquareMatrix<N> &m2) const
  {
    for (int i = 0; i < N; ++i)
      for (int j = 0; j < N; ++j)
        if (m[i][j] != m2.m[i][j]) return true;
    return false;
  }


  bool operator<(const SquareMatrix<N> &m2) const
  {
    for (int i = 0; i < N; ++i)
      for (int j = 0; j < N; ++j) {
        if (m[i][j] < m2.m[i][j]) return true;
        if (m[i][j] > m2.m[i][j]) return false;
      }
    return false;
  }


  bool IsIdentity() const;

  std::string ToString() const;


  std::span<const Float> operator[](int i) const
  {
    return m[i];
  }

  std::span<Float> operator[](int i)
  {
    return std::span<Float>(m[i]);
  }

 private:
  Float m[N][N];
};

// SquareMatrix Inline Methods
template <int N>
inline bool SquareMatrix<N>::IsIdentity() const
{
  for (int i = 0; i < N; ++i)
    for (int j = 0; j < N; ++j) {
      if (i == j) {
        if (m[i][j] != 1) return false;
      }
      else if (m[i][j] != 0)
        return false;
    }
  return true;
}

// SquareMatrix Inline Functions
template <int N>
inline SquareMatrix<N> operator*(
    Float s, const SquareMatrix<N> &m
)
{
  return m * s;
}

template <typename Tresult, int N, typename T>
inline Tresult Mul(const SquareMatrix<N> &m, const T &v)
{
  Tresult result;
  for (int i = 0; i < N; ++i) {
    result[i] = 0;
    for (int j = 0; j < N; ++j) result[i] += m[i][j] * v[j];
  }
  return result;
}

template <int N>
Float Determinant(const SquareMatrix<N> &m);

template <>
inline Float Determinant(const SquareMatrix<3> &m)
{
  Float minor12 =
      DifferenceOfProducts(m[1][1], m[2][2], m[1][2], m[2][1]);
  Float minor02 =
      DifferenceOfProducts(m[1][0], m[2][2], m[1][2], m[2][0]);
  Float minor01 =
      DifferenceOfProducts(m[1][0], m[2][1], m[1][1], m[2][0]);
  return m[0][2] * minor01 +
         DifferenceOfProducts(
             m[0][0], minor12, m[0][1], minor02
         );
}

template <int N>
inline SquareMatrix<N> Transpose(const SquareMatrix<N> &m);
template <int N>
std::optional<SquareMatrix<N>> Inverse(const SquareMatrix<N> &);

template <int N>
SquareMatrix<N> InvertOrExit(const SquareMatrix<N> &m)
{
  std::optional<SquareMatrix<N>> inv = Inverse(m);
  CHECK(inv.has_value());
  return *inv;
}

template <int N>
inline SquareMatrix<N> Transpose(const SquareMatrix<N> &m)
{
  SquareMatrix<N> r;
  for (int i = 0; i < N; ++i)
    for (int j = 0; j < N; ++j) r[i][j] = m[j][i];
  return r;
}

template <>
inline std::optional<SquareMatrix<3>> Inverse(
    const SquareMatrix<3> &m
)
{
  Float det = Determinant(m);
  if (det == 0) return {};
  Float invDet = 1 / det;

  SquareMatrix<3> r;

  r[0][0] =
      invDet *
      DifferenceOfProducts(m[1][1], m[2][2], m[1][2], m[2][1]);
  r[1][0] =
      invDet *
      DifferenceOfProducts(m[1][2], m[2][0], m[1][0], m[2][2]);
  r[2][0] =
      invDet *
      DifferenceOfProducts(m[1][0], m[2][1], m[1][1], m[2][0]);
  r[0][1] =
      invDet *
      DifferenceOfProducts(m[0][2], m[2][1], m[0][1], m[2][2]);
  r[1][1] =
      invDet *
      DifferenceOfProducts(m[0][0], m[2][2], m[0][2], m[2][0]);
  r[2][1] =
      invDet *
      DifferenceOfProducts(m[0][1], m[2][0], m[0][0], m[2][1]);
  r[0][2] =
      invDet *
      DifferenceOfProducts(m[0][1], m[1][2], m[0][2], m[1][1]);
  r[1][2] =
      invDet *
      DifferenceOfProducts(m[0][2], m[1][0], m[0][0], m[1][2]);
  r[2][2] =
      invDet *
      DifferenceOfProducts(m[0][0], m[1][1], m[0][1], m[1][0]);

  return r;
}

template <int N, typename T>
inline T operator*(const SquareMatrix<N> &m, const T &v)
{
  return Mul<T>(m, v);
}

template <>
inline SquareMatrix<4> operator*(
    const SquareMatrix<4> &m1, const SquareMatrix<4> &m2
)
{
  SquareMatrix<4> r;
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j)
      r[i][j] = InnerProduct(
          m1[i][0], m2[0][j], m1[i][1], m2[1][j], m1[i][2],
          m2[2][j], m1[i][3], m2[3][j]
      );
  return r;
}

template <>
inline SquareMatrix<3> operator*(
    const SquareMatrix<3> &m1, const SquareMatrix<3> &m2
)
{
  SquareMatrix<3> r;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      r[i][j] = InnerProduct(
          m1[i][0], m2[0][j], m1[i][1], m2[1][j], m1[i][2],
          m2[2][j]
      );
  return r;
}

template <int N>
inline SquareMatrix<N> operator*(
    const SquareMatrix<N> &m1, const SquareMatrix<N> &m2
)
{
  SquareMatrix<N> r;
  for (int i = 0; i < N; ++i)
    for (int j = 0; j < N; ++j) {
      r[i][j] = 0;
      for (int k = 0; k < N; ++k)
        r[i][j] = FMA(m1[i][k], m2[k][j], r[i][j]);
    }
  return r;
}

template <int N>
inline SquareMatrix<N>::SquareMatrix(std::span<const Float> t)
{
  CHECK_EQ(N * N, t.size());
  for (int i = 0; i < N * N; ++i) m[i / N][i % N] = t[i];
}

template <int N>
SquareMatrix<N> operator*(
    const SquareMatrix<N> &m1, const SquareMatrix<N> &m2
);

template <>
inline Float Determinant(const SquareMatrix<1> &m)
{
  return m[0][0];
}

template <>
inline Float Determinant(const SquareMatrix<2> &m)
{
  return DifferenceOfProducts(
      m[0][0], m[1][1], m[0][1], m[1][0]
  );
}

template <>
inline Float Determinant(const SquareMatrix<4> &m)
{
  Float s0 =
      DifferenceOfProducts(m[0][0], m[1][1], m[1][0], m[0][1]);
  Float s1 =
      DifferenceOfProducts(m[0][0], m[1][2], m[1][0], m[0][2]);
  Float s2 =
      DifferenceOfProducts(m[0][0], m[1][3], m[1][0], m[0][3]);

  Float s3 =
      DifferenceOfProducts(m[0][1], m[1][2], m[1][1], m[0][2]);
  Float s4 =
      DifferenceOfProducts(m[0][1], m[1][3], m[1][1], m[0][3]);
  Float s5 =
      DifferenceOfProducts(m[0][2], m[1][3], m[1][2], m[0][3]);

  Float c0 =
      DifferenceOfProducts(m[2][0], m[3][1], m[3][0], m[2][1]);
  Float c1 =
      DifferenceOfProducts(m[2][0], m[3][2], m[3][0], m[2][2]);
  Float c2 =
      DifferenceOfProducts(m[2][0], m[3][3], m[3][0], m[2][3]);

  Float c3 =
      DifferenceOfProducts(m[2][1], m[3][2], m[3][1], m[2][2]);
  Float c4 =
      DifferenceOfProducts(m[2][1], m[3][3], m[3][1], m[2][3]);
  Float c5 =
      DifferenceOfProducts(m[2][2], m[3][3], m[3][2], m[2][3]);

  return (
      DifferenceOfProducts(s0, c5, s1, c4) +
      DifferenceOfProducts(s2, c3, -s3, c2) +
      DifferenceOfProducts(s5, c0, s4, c1)
  );
}

template <int N>
inline Float Determinant(const SquareMatrix<N> &m)
{
  SquareMatrix<N - 1> sub;
  Float det = 0;
  // Inefficient, but we don't currently use N>4 anyway..
  for (int i = 0; i < N; ++i) {
    // Sub-matrix without row 0 and column i
    for (int j = 0; j < N - 1; ++j)
      for (int k = 0; k < N - 1; ++k)
        sub[j][k] = m[j + 1][k < i ? k : k + 1];

    Float sign = (i & 1) ? -1 : 1;
    det += sign * m[0][i] * Determinant(sub);
  }
  return det;
}

template <>
inline std::optional<SquareMatrix<4>> Inverse(
    const SquareMatrix<4> &m
)
{
  // Via:
  // https://github.com/google/ion/blob/master/ion/math/matrixutils.cc,
  // (c) Google, Apache license.

  // For 4x4 do not compute the adjugate as the transpose of the
  // cofactor matrix, because this results in extra work.
  // Several calculations can be shared across the
  // sub-determinants.
  //
  // This approach is explained in David Eberly's Geometric
  // Tools book, excerpted here:
  //   http://www.geometrictools.com/Documentation/LaplaceExpansionTheorem.pdf
  Float s0 =
      DifferenceOfProducts(m[0][0], m[1][1], m[1][0], m[0][1]);
  Float s1 =
      DifferenceOfProducts(m[0][0], m[1][2], m[1][0], m[0][2]);
  Float s2 =
      DifferenceOfProducts(m[0][0], m[1][3], m[1][0], m[0][3]);

  Float s3 =
      DifferenceOfProducts(m[0][1], m[1][2], m[1][1], m[0][2]);
  Float s4 =
      DifferenceOfProducts(m[0][1], m[1][3], m[1][1], m[0][3]);
  Float s5 =
      DifferenceOfProducts(m[0][2], m[1][3], m[1][2], m[0][3]);

  Float c0 =
      DifferenceOfProducts(m[2][0], m[3][1], m[3][0], m[2][1]);
  Float c1 =
      DifferenceOfProducts(m[2][0], m[3][2], m[3][0], m[2][2]);
  Float c2 =
      DifferenceOfProducts(m[2][0], m[3][3], m[3][0], m[2][3]);

  Float c3 =
      DifferenceOfProducts(m[2][1], m[3][2], m[3][1], m[2][2]);
  Float c4 =
      DifferenceOfProducts(m[2][1], m[3][3], m[3][1], m[2][3]);
  Float c5 =
      DifferenceOfProducts(m[2][2], m[3][3], m[3][2], m[2][3]);

  Float determinant = InnerProduct(
      s0, c5, -s1, c4, s2, c3, s3, c2, s5, c0, -s4, c1
  );
  if (determinant == 0) return {};
  Float s = 1 / determinant;

  Float inv[4][4] = {
      {s * InnerProduct(m[1][1], c5, m[1][3], c3, -m[1][2],         c4),
       s * InnerProduct(
       -m[0][1], c5, m[0][2], c4, -m[0][3], c3
 ),s * InnerProduct(m[3][1], s5, m[3][3], s3, -m[3][2],         s4),
       s * InnerProduct(-m[2][1], s5, m[2][2], s4, -m[2][3], s3)                                       },

      { s * InnerProduct(
 -m[1][0], c5, m[1][2], c2, -m[1][3], c1
 ), s * InnerProduct(m[0][0], c5, m[0][3], c1, -m[0][2],         c2),
       s * InnerProduct(
       -m[3][0], s5, m[3][2], s2, -m[3][3], s1
 ), s * InnerProduct(m[2][0], s5, m[2][3], s1, -m[2][2], s2)},

      {s * InnerProduct(m[1][0], c4, m[1][3], c0, -m[1][1],         c2),
       s * InnerProduct(
       -m[0][0], c4, m[0][1], c2, -m[0][3], c0
 ), s * InnerProduct(m[3][0], s4, m[3][3], s0, -m[3][1],         s2),
       s * InnerProduct(-m[2][0], s4, m[2][1], s2, -m[2][3], s0)                                       },

      { s * InnerProduct(
 -m[1][0], c3, m[1][1], c1, -m[1][2], c0
 ), s * InnerProduct(m[0][0], c3, m[0][2], c0, -m[0][1],         c1),
       s * InnerProduct(
       -m[3][0], s3, m[3][1], s1, -m[3][2], s0
 ), s * InnerProduct(m[2][0], s3, m[2][2], s0, -m[2][1], s1)}
  };

  return SquareMatrix<4>(inv);
}

extern template class SquareMatrix<2>;
extern template class SquareMatrix<3>;
extern template class SquareMatrix<4>;


template <typename T, typename U, typename V>
inline constexpr T Clamp(T val, U low, V high)
{
  if (val < low)
    return T(low);
  else if (val > high)
    return T(high);
  else
    return val;
}

// http://www.plunk.org/~hatch/rightway.html
inline Float SinXOverX(Float x)
{
  if (1 - x * x == 1) return 1;
  return std::sin(x) / x;
}

inline float SafeASin(float x)
{
  DCHECK(x >= -1.0001 && x <= 1.0001);
  return std::asin(Clamp(x, -1, 1));
}
inline float SafeACos(float x)
{
  DCHECK(x >= -1.0001 && x <= 1.0001);
  return std::acos(Clamp(x, -1, 1));
}


inline double SafeASin(double x)
{
  DCHECK(x >= -1.0001 && x <= 1.0001);
  return std::asin(Clamp(x, -1, 1));
}

inline double SafeACos(double x)
{
  DCHECK(x >= -1.0001 && x <= 1.0001);
  return std::acos(Clamp(x, -1, 1));
}
}  // namespace pbrt
