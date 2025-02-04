// Mostly copied from PBRT
#pragma once
#include <limits>

/**
Required by
 * vecmath.h
Sections:
 * Float
 * Mathematical Constants
 * Math Functions
**/


/** Float **/
// FIXME: Things don't work if this is set to float
#define Float double

// Floating-point Constants
static constexpr Float Infinity = std::numeric_limits<Float>::infinity();

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
template<typename T>
inline constexpr T Sqr(T x) {
    return x*x;
}

// TODO: PBRT used intervals for error correction
template<typename T>
inline T DifferenceOfProducts(T a, T b, T c, T d) {
    return a*b - c*d;
}
