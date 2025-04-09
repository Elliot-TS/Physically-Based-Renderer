// pbrt is Copyright(c) 1998-2020 Matt Pharr, Wenzel Jakob, and
// Greg Humphreys. The pbrt source code is licensed under the
// Apache License, Version 2.0. SPDX: Apache-2.0

#ifndef PBRT_UTIL_TRANSFORM_H
#define PBRT_UTIL_TRANSFORM_H

#include <pbrt/interaction.h>
#include <pbrt/math/math.h>
#include <pbrt/math/vecmath.h>
#include <pbrt/pbrt.h>
#include <pbrt/ray.h>
#include <pbrt/util/hash.h>
#include <cmath>
#include <functional>
#include <limits>

namespace pbrt {

// Transform Definition
class Transform {
 public:
  // Transform Public Methods

  inline Ray ApplyInverse(const Ray &r, Float *tMax = nullptr)
      const;

  inline RayDifferential ApplyInverse(
      const RayDifferential &r, Float *tMax = nullptr
  ) const;
  template <typename T>
  inline Vector3<T> ApplyInverse(Vector3<T> v) const;
  template <typename T>
  inline Normal3<T> ApplyInverse(Normal3<T>) const;

  std::string ToString() const;

  Transform() = default;


  Transform(const SquareMatrix<4> &m): m(m)
  {
    std::optional<SquareMatrix<4>> inv = Inverse(m);
    if (inv)
      mInv = *inv;
    else {
      // Initialize _mInv_ with not-a-number values
      Float NaN =
          std::numeric_limits<Float>::has_signaling_NaN
              ? std::numeric_limits<Float>::signaling_NaN()
              : std::numeric_limits<Float>::quiet_NaN();
      for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j) mInv[i][j] = NaN;
    }
  }


  Transform(const Float mat[4][4])
      : Transform(SquareMatrix<4>(mat))
  {}


  Transform(
      const SquareMatrix<4> &m, const SquareMatrix<4> &mInv
  )
      : m(m), mInv(mInv)
  {}


  const SquareMatrix<4> &GetMatrix() const { return m; }

  const SquareMatrix<4> &GetInverseMatrix() const
  {
    return mInv;
  }


  bool operator==(const Transform &t) const { return t.m == m; }

  bool operator!=(const Transform &t) const { return t.m != m; }

  bool IsIdentity() const { return m.IsIdentity(); }


  bool HasScale(Float tolerance = 1e-3f) const
  {
    Float la2 = LengthSquared((*this)(Vector3f(1, 0, 0)));
    Float lb2 = LengthSquared((*this)(Vector3f(0, 1, 0)));
    Float lc2 = LengthSquared((*this)(Vector3f(0, 0, 1)));
    return (
        std::abs(la2 - 1) > tolerance ||
        std::abs(lb2 - 1) > tolerance ||
        std::abs(lc2 - 1) > tolerance
    );
  }

  template <typename T>
  Point3<T> operator()(Point3<T> p) const;

  template <typename T>
  inline Point3<T> ApplyInverse(Point3<T> p) const;

  template <typename T>
  Vector3<T> operator()(Vector3<T> v) const;

  template <typename T>
  Normal3<T> operator()(Normal3<T>) const;


  Ray operator()(const Ray &r, Float *tMax = nullptr) const;

  RayDifferential operator()(
      const RayDifferential &r, Float *tMax = nullptr
  ) const;


  Bounds3f operator()(const Bounds3f &b) const;


  Transform operator*(const Transform &t2) const;


  bool SwapsHandedness() const;


  explicit Transform(const Frame &frame);


  explicit Transform(Quaternion q);


  explicit operator Quaternion() const;

  void Decompose(
      Vector3f *T, SquareMatrix<4> *R, SquareMatrix<4> *S
  ) const;


  Interaction operator()(const Interaction &in) const;

  Interaction ApplyInverse(const Interaction &in) const;

  // SurfaceInteraction operator()(const SurfaceInteraction &si
  //) const;

  // SurfaceInteraction ApplyInverse(const SurfaceInteraction
  // &in ) const;


  /*
  Point3fi operator()(const Point3fi &p) const
  {
    Float x = Float(p.x), y = Float(p.y), z = Float(p.z);
    // Compute transformed coordinates from point _(x, y, z)_
    Float xp =
        (m[0][0] * x + m[0][1] * y) + (m[0][2] * z + m[0][3]);
    Float yp =
        (m[1][0] * x + m[1][1] * y) + (m[1][2] * z + m[1][3]);
    Float zp =
        (m[2][0] * x + m[2][1] * y) + (m[2][2] * z + m[2][3]);
    Float wp =
        (m[3][0] * x + m[3][1] * y) + (m[3][2] * z + m[3][3]);

    // Compute absolute error for transformed point, _pError_
    Vector3f pError;
    if (p.IsExact()) {
      // Compute error for transformed exact _p_
      pError.x =
          gamma(3) *
          (std::abs(m[0][0] * x) + std::abs(m[0][1] * y) +
           std::abs(m[0][2] * z) + std::abs(m[0][3]));
      pError.y =
          gamma(3) *
          (std::abs(m[1][0] * x) + std::abs(m[1][1] * y) +
           std::abs(m[1][2] * z) + std::abs(m[1][3]));
      pError.z =
          gamma(3) *
          (std::abs(m[2][0] * x) + std::abs(m[2][1] * y) +
           std::abs(m[2][2] * z) + std::abs(m[2][3]));
    }
    else {
      // Compute error for transformed approximate _p_
      Vector3f pInError = p.Error();
      pError.x =
          (gamma(3) + 1) * (std::abs(m[0][0]) * pInError.x +
                            std::abs(m[0][1]) * pInError.y +
                            std::abs(m[0][2]) * pInError.z) +
          gamma(3) *
              (std::abs(m[0][0] * x) + std::abs(m[0][1] * y) +
               std::abs(m[0][2] * z) + std::abs(m[0][3]));
      pError.y =
          (gamma(3) + 1) * (std::abs(m[1][0]) * pInError.x +
                            std::abs(m[1][1]) * pInError.y +
                            std::abs(m[1][2]) * pInError.z) +
          gamma(3) *
              (std::abs(m[1][0] * x) + std::abs(m[1][1] * y) +
               std::abs(m[1][2] * z) + std::abs(m[1][3]));
      pError.z =
          (gamma(3) + 1) * (std::abs(m[2][0]) * pInError.x +
                            std::abs(m[2][1]) * pInError.y +
                            std::abs(m[2][2]) * pInError.z) +
          gamma(3) *
              (std::abs(m[2][0] * x) + std::abs(m[2][1] * y) +
               std::abs(m[2][2] * z) + std::abs(m[2][3]));
    }

    if (wp == 1)
      return Point3fi(Point3f(xp, yp, zp), pError);
    else
      return Point3fi(Point3f(xp, yp, zp), pError) / wp;
  }


  Vector3fi operator()(const Vector3fi &v) const;

  Point3fi ApplyInverse(const Point3fi &p) const;
  */

 private:
  // Transform Private Members
  SquareMatrix<4> m, mInv;
};

// Transform Function Declarations

Transform Translate(Vector3f delta);


Transform Scale(Float x, Float y, Float z);


Transform RotateX(Float theta);

Transform RotateY(Float theta);

Transform RotateZ(Float theta);


Transform LookAt(Point3f pos, Point3f look, Vector3f up);


Transform Orthographic(Float znear, Float zfar);


Transform Perspective(Float fov, Float znear, Float zfar);

// Transform Inline Functions
inline Transform Inverse(const Transform &t)
{
  return Transform(t.GetInverseMatrix(), t.GetMatrix());
}

inline Transform Transpose(const Transform &t)
{
  return Transform(
      Transpose(t.GetMatrix()), Transpose(t.GetInverseMatrix())
  );
}

inline Transform Rotate(
    Float sinTheta, Float cosTheta, Vector3f axis
)
{
  Vector3f a = Normalize(axis);
  SquareMatrix<4> m;
  // Compute rotation of first basis vector
  m[0][0] = a.x * a.x + (1 - a.x * a.x) * cosTheta;
  m[0][1] = a.x * a.y * (1 - cosTheta) - a.z * sinTheta;
  m[0][2] = a.x * a.z * (1 - cosTheta) + a.y * sinTheta;
  m[0][3] = 0;

  // Compute rotations of second and third basis vectors
  m[1][0] = a.x * a.y * (1 - cosTheta) + a.z * sinTheta;
  m[1][1] = a.y * a.y + (1 - a.y * a.y) * cosTheta;
  m[1][2] = a.y * a.z * (1 - cosTheta) - a.x * sinTheta;
  m[1][3] = 0;

  m[2][0] = a.x * a.z * (1 - cosTheta) - a.y * sinTheta;
  m[2][1] = a.y * a.z * (1 - cosTheta) + a.x * sinTheta;
  m[2][2] = a.z * a.z + (1 - a.z * a.z) * cosTheta;
  m[2][3] = 0;

  return Transform(m, Transpose(m));
}

inline Transform Rotate(Float theta, Vector3f axis)
{
  Float sinTheta = std::sin(Radians(theta));
  Float cosTheta = std::cos(Radians(theta));
  return Rotate(sinTheta, cosTheta, axis);
}

inline Transform RotateFromTo(Vector3f from, Vector3f to)
{
  // Compute intermediate vector for vector reflection
  Vector3f refl;
  if (std::abs(from.x) < 0.72f && std::abs(to.x) < 0.72f)
    refl = Vector3f(1, 0, 0);
  else if (std::abs(from.y) < 0.72f && std::abs(to.y) < 0.72f)
    refl = Vector3f(0, 1, 0);
  else
    refl = Vector3f(0, 0, 1);

  // Initialize matrix _r_ for rotation
  Vector3f u = refl - from, v = refl - to;
  SquareMatrix<4> r;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      // Initialize matrix element _r[i][j]_
      r[i][j] =
          ((i == j) ? 1 : 0) - 2 / Dot(u, u) * u[i] * u[j] -
          2 / Dot(v, v) * v[i] * v[j] +
          4 * Dot(u, v) / (Dot(u, u) * Dot(v, v)) * v[i] * u[j];

  return Transform(r, Transpose(r));
}

/*
inline Vector3fi Transform::operator()(const Vector3fi &v) const
{
  Float x = Float(v.x), y = Float(v.y), z = Float(v.z);
  Vector3f vOutError;
  if (v.IsExact()) {
    vOutError.x = gamma(3) * (std::abs(m[0][0] * x) +
                              std::abs(m[0][1] * y) +
                              std::abs(m[0][2] * z));
    vOutError.y = gamma(3) * (std::abs(m[1][0] * x) +
                              std::abs(m[1][1] * y) +
                              std::abs(m[1][2] * z));
    vOutError.z = gamma(3) * (std::abs(m[2][0] * x) +
                              std::abs(m[2][1] * y) +
                              std::abs(m[2][2] * z));
  }
  else {
    Vector3f vInError = v.Error();
    vOutError.x =
        (gamma(3) + 1) * (std::abs(m[0][0]) * vInError.x +
                          std::abs(m[0][1]) * vInError.y +
                          std::abs(m[0][2]) * vInError.z) +
        gamma(3) *
            (std::abs(m[0][0] * x) + std::abs(m[0][1] * y) +
             std::abs(m[0][2] * z));
    vOutError.y =
        (gamma(3) + 1) * (std::abs(m[1][0]) * vInError.x +
                          std::abs(m[1][1]) * vInError.y +
                          std::abs(m[1][2]) * vInError.z) +
        gamma(3) *
            (std::abs(m[1][0] * x) + std::abs(m[1][1] * y) +
             std::abs(m[1][2] * z));
    vOutError.z =
        (gamma(3) + 1) * (std::abs(m[2][0]) * vInError.x +
                          std::abs(m[2][1]) * vInError.y +
                          std::abs(m[2][2]) * vInError.z) +
        gamma(3) *
            (std::abs(m[2][0] * x) + std::abs(m[2][1] * y) +
             std::abs(m[2][2] * z));
  }

  Float xp = m[0][0] * x + m[0][1] * y + m[0][2] * z;
  Float yp = m[1][0] * x + m[1][1] * y + m[1][2] * z;
  Float zp = m[2][0] * x + m[2][1] * y + m[2][2] * z;

  return Vector3fi(Vector3f(xp, yp, zp), vOutError);
}

*/
// Transform Inline Methods
template <typename T>
inline Point3<T> Transform::operator()(Point3<T> p) const
{
  T xp =
      m[0][0] * p.x + m[0][1] * p.y + m[0][2] * p.z + m[0][3];
  T yp =
      m[1][0] * p.x + m[1][1] * p.y + m[1][2] * p.z + m[1][3];
  T zp =
      m[2][0] * p.x + m[2][1] * p.y + m[2][2] * p.z + m[2][3];
  T wp =
      m[3][0] * p.x + m[3][1] * p.y + m[3][2] * p.z + m[3][3];
  if (wp == 1)
    return Point3<T>(xp, yp, zp);
  else
    return Point3<T>(xp, yp, zp) / wp;
}

template <typename T>
inline Vector3<T> Transform::operator()(Vector3<T> v) const
{
  return Vector3<T>(
      m[0][0] * v.x + m[0][1] * v.y + m[0][2] * v.z,
      m[1][0] * v.x + m[1][1] * v.y + m[1][2] * v.z,
      m[2][0] * v.x + m[2][1] * v.y + m[2][2] * v.z
  );
}

template <typename T>
inline Normal3<T> Transform::operator()(Normal3<T> n) const
{
  T x = n.x, y = n.y, z = n.z;
  return Normal3<T>(
      mInv[0][0] * x + mInv[1][0] * y + mInv[2][0] * z,
      mInv[0][1] * x + mInv[1][1] * y + mInv[2][1] * z,
      mInv[0][2] * x + mInv[1][2] * y + mInv[2][2] * z
  );
}

inline Ray Transform::operator()(const Ray &r, Float *tMax)
    const
{
  Point3f o = (*this)(Point3f(r.origin));
  Vector3f d = (*this)(r.direction);
  // Offset ray origin to edge of error bounds and compute _tMax_
  /*
  if (Float lengthSquared = LengthSquared(d); lengthSquared > 0)
  {
    Float dt = Dot(Abs(d), o.Error()) / lengthSquared;
    o += d * dt;
    if (tMax) *tMax -= dt;
  }
  */

  return Ray(Point3f(o), d, r.time /*, r.medium*/);
}

inline RayDifferential Transform::operator()(
    const RayDifferential &r, Float *tMax
) const
{
  Ray tr = (*this)(Ray(r), tMax);
  RayDifferential
      ret(tr.origin, tr.direction, tr.time /*, tr.medium*/);
  ret.hasDifferentials = r.hasDifferentials;
  ret.rxOrigin = (*this)(r.rxOrigin);
  ret.ryOrigin = (*this)(r.ryOrigin);
  ret.rxDirection = (*this)(r.rxDirection);
  ret.ryDirection = (*this)(r.ryDirection);
  return ret;
}

inline Transform::Transform(const Frame &frame)
    : Transform(SquareMatrix<4>(
          frame.x.x, frame.x.y, frame.x.z, 0, frame.y.x,
          frame.y.y, frame.y.z, 0, frame.z.x, frame.z.y,
          frame.z.z, 0, 0, 0, 0, 1
      ))
{}

inline Transform::Transform(Quaternion q)
{
  Float xx = q.v.x * q.v.x, yy = q.v.y * q.v.y,
        zz = q.v.z * q.v.z;
  Float xy = q.v.x * q.v.y, xz = q.v.x * q.v.z,
        yz = q.v.y * q.v.z;
  Float wx = q.v.x * q.w, wy = q.v.y * q.w, wz = q.v.z * q.w;

  mInv[0][0] = 1 - 2 * (yy + zz);
  mInv[0][1] = 2 * (xy + wz);
  mInv[0][2] = 2 * (xz - wy);
  mInv[1][0] = 2 * (xy - wz);
  mInv[1][1] = 1 - 2 * (xx + zz);
  mInv[1][2] = 2 * (yz + wx);
  mInv[2][0] = 2 * (xz + wy);
  mInv[2][1] = 2 * (yz - wx);
  mInv[2][2] = 1 - 2 * (xx + yy);

  // Transpose since we are left-handed.  Ugh.
  m = Transpose(mInv);
}

template <typename T>
inline Point3<T> Transform::ApplyInverse(Point3<T> p) const
{
  T x = p.x, y = p.y, z = p.z;
  T xp = (mInv[0][0] * x + mInv[0][1] * y) +
         (mInv[0][2] * z + mInv[0][3]);
  T yp = (mInv[1][0] * x + mInv[1][1] * y) +
         (mInv[1][2] * z + mInv[1][3]);
  T zp = (mInv[2][0] * x + mInv[2][1] * y) +
         (mInv[2][2] * z + mInv[2][3]);
  T wp = (mInv[3][0] * x + mInv[3][1] * y) +
         (mInv[3][2] * z + mInv[3][3]);
  CHECK_NE(wp, 0);
  if (wp == 1)
    return Point3<T>(xp, yp, zp);
  else
    return Point3<T>(xp, yp, zp) / wp;
}

template <typename T>
inline Vector3<T> Transform::ApplyInverse(Vector3<T> v) const
{
  T x = v.x, y = v.y, z = v.z;
  return Vector3<T>(
      mInv[0][0] * x + mInv[0][1] * y + mInv[0][2] * z,
      mInv[1][0] * x + mInv[1][1] * y + mInv[1][2] * z,
      mInv[2][0] * x + mInv[2][1] * y + mInv[2][2] * z
  );
}

template <typename T>
inline Normal3<T> Transform::ApplyInverse(Normal3<T> n) const
{
  T x = n.x, y = n.y, z = n.z;
  return Normal3<T>(
      m[0][0] * x + m[1][0] * y + m[2][0] * z,
      m[0][1] * x + m[1][1] * y + m[2][1] * z,
      m[0][2] * x + m[1][2] * y + m[2][2] * z
  );
}

inline Ray Transform::ApplyInverse(const Ray &r, Float *tMax)
    const
{
  Point3f o = ApplyInverse(Point3f(r.origin));
  Vector3f d = ApplyInverse(r.direction);
  // Offset ray origin to edge of error bounds and compute _tMax_
  /*
  Float lengthSquared = LengthSquared(d);
  if (lengthSquared > 0) {
    Vector3f oError(
        o.x.Width() / 2, o.y.Width() / 2, o.z.Width() / 2
    );
    Float dt = Dot(Abs(d), oError) / lengthSquared;
    o += d * dt;
    if (tMax) *tMax -= dt;
  }
  */
  return Ray(Point3f(o), d, r.time /*, r.medium*/);
}

inline RayDifferential Transform::ApplyInverse(
    const RayDifferential &r, Float *tMax
) const
{
  Ray tr = ApplyInverse(Ray(r), tMax);
  RayDifferential
      ret(tr.origin, tr.direction, tr.time /*, tr.medium*/);
  ret.hasDifferentials = r.hasDifferentials;
  ret.rxOrigin = ApplyInverse(r.rxOrigin);
  ret.ryOrigin = ApplyInverse(r.ryOrigin);
  ret.rxDirection = ApplyInverse(r.rxDirection);
  ret.ryDirection = ApplyInverse(r.ryDirection);
  return ret;
}

}  // namespace pbrt

namespace std {

template <>
struct hash<pbrt::Transform> {
  size_t operator()(const pbrt::Transform &t) const
  {
    pbrt::SquareMatrix<4> m = t.GetMatrix();
    return pbrt::Hash(m);
  }
};

}  // namespace std

#endif  // PBRT_UTIL_TRANSFORM_H
