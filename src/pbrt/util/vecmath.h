#pragma once
#include <string>
#include <iostream>
//#include "./check.h"


namespace pbrt{
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
    //bool HasNaN() const { return IsNaN(x) || IsNaN(y); }

    //template <typename U>
    //auto operator+(Child<U> c) const -> Child<decltype(T{} + U{})> {
        //DCHECK(!c.HasNaN());
        //return {x + c.x, y + c.y};
    //}

    //template <typename U>
    //Child<T> &operator+=(Child<U> c) {
        //DCHECK(!c.HasNaN());
        //x += c.x;
        //y += c.y;
        //return static_cast<Child<T> &>(*this);
    //}

    //template <typename U>
    //auto operator-(Child<U> c) const -> Child<decltype(T{} - U{})> {
        //DCHECK(!c.HasNaN());
        //return {x - c.x, y - c.y};
    //}

    //template <typename U>
    //Child<T> &operator-=(Child<U> c) {
        //DCHECK(!c.HasNaN());
        //x -= c.x;
        //y -= c.y;
        //return static_cast<Child<T> &>(*this);
    //}

    //bool operator==(Child<T> c) const { return x == c.x && y == c.y; }
    //bool operator!=(Child<T> c) const { return x != c.x || y != c.y; }

    //template <typename U>
    //auto operator*(U s) const -> Child<decltype(T{} * U{})> {
        //return {s * x, s * y};
    //}

    //template <typename U>
    //Child<T> &operator*=(U s) {
        //DCHECK(!IsNaN(s));
        //x *= s;
        //y *= s;
        //return static_cast<Child<T> &>(*this);
    //}

    //template <typename U>
    //auto operator/(U d) const -> Child<decltype(T{} / U{})> {
        //DCHECK(d != 0 && !IsNaN(d));
        //return {x / d, y / d};
    //}

    //template <typename U>
    //Child<T> &operator/=(U d) {
        //DCHECK_NE(d, 0);
        //DCHECK(!IsNaN(d));
        //x /= d;
        //y /= d;
        //return static_cast<Child<T> &>(*this);
    //}

    //Child<T> operator-() const { return {-x, -y}; }

    //T operator[](int i) const {
        //DCHECK(i >= 0 && i <= 1);
        //return (i == 0) ? x : y;
    //}

    //T &operator[](int i) {
        //DCHECK(i >= 0 && i <= 1);
        //return (i == 0) ? x : y;
    //}

    friend std::ostream& operator<<(std::ostream& os, const Child<T>& obj) {
        os << "[ " << obj.x << ", " << obj.y << " ]";
        return os;
    }
};


template <typename T>
class Vector2 : public Tuple2<Vector2, T> {
  public:
    // Vector2 Public Methods
    using Tuple2<Vector2, T>::x;
    using Tuple2<Vector2, T>::y;

    Vector2() = default;
    Vector2(T x, T y) : Tuple2<pbrt::Vector2, T>(x, y) {}
    //template <typename U>
    //explicit Vector2(Point2<U> p);
    //template <typename U>
    //explicit Vector2(Vector2<U> v)
        //: Tuple2<pbrt::Vector2, T>(T(v.x), T(v.y)) {}
};


}

