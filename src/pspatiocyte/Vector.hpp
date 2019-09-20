//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
//        This file is part of Motocyte
//
//        Copyright (C) 2019 Satya N.V. Arjunan
//
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
//
// Motocyte is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public
// License as published by the Free Software Foundation; either
// version 2 of the License, or (at your option) any later version.
// 
// Motocyte is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// See the GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public
// License along with Motocyte -- see the file COPYING.
// If not, write to the Free Software Foundation, Inc.,
// 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
// 
//END_HEADER
//
// written by Satya N. V. Arjunan <satya.arjunan@gmail.com>
//


#ifndef __Vector_hpp
#define __Vector_hpp

#include <cmath>

template<typename T>
struct Vector
{
  Vector(const T a=0, const T b=0, const T c=0);
  Vector<T> operator + (const T val) const;
  Vector<T> operator - (const T val) const;
  Vector<T> operator * (const T val) const;
  Vector<T> operator / (const T val) const;
  const Vector<T>& operator += (const T val);
  const Vector<T>& operator -= (const T val);
  const Vector<T>& operator *= (const T val);
  const Vector<T>& operator /= (const T val);
  Vector<T> operator + (const Vector<T>& vector) const;
  Vector<T> operator - (const Vector<T>& vector) const;
  Vector<T> operator * (const Vector<T>& vector) const;
  Vector<T> operator / (const Vector<T>& vector) const;
  const Vector<T>& operator += (const Vector<T>& vector);
  const Vector<T>& operator -= (const Vector<T>& vector);
  const Vector<T>& operator *= (const Vector<T>& vector);
  const Vector<T>& operator /= (const Vector<T>& vector);
  const bool operator == (const Vector<T>& vector) const;
  const bool operator != (const Vector<T>& vector) const;
  T distance(const Vector<T>& vector) const;
  Vector<T> cross_product(const Vector<T>& a) const;
  T angle_between_vectors(const Vector<T>& a) const;
  T dot_product(const Vector<T>& a) const;
  void mod(const Vector<T>& vector);
  Vector<T> norm() const;
  T magnitude() const;
  T magnitude_squared() const;
  T x;
  T y;
  T z;
};

template<class T>
T Vector<T>::angle_between_vectors(const Vector<T>& a) const {
  T dp(this->dot_product( a));
  T cosAng(dp/(this->magnitude()*a.magnitude()));
  if (abs(cosAng) > 1.0) { 
    cosAng /= cosAng;
  }
  return acos(cosAng);	// return value is in radians. 
}

template<class T>
Vector<T> Vector<T>::cross_product(const Vector<T>& a) const {
  return Vector<T>(a.y*z - a.z*y, a.z*x - a.x*z, a.x*y - a.y*x);
}

template<class T>
Vector<T> Vector<T>::norm() const {
  T denom(sqrt(x*x + y*y + z*z));
  return Vector<T>(x/denom, y/denom, z/denom);
}

template<class T>
T Vector<T>::magnitude() const {
  return sqrt(x*x + y*y + z*z);
}

template<class T>
T Vector<T>::dot_product(const Vector<T>& a) const {
  return (a.x*x + a.y*y + a.z*z);
}

template<class T>
T Vector<T>::distance(const Vector<T>& vector) const {
  return sqrt(pow(vector.x-x, 2) + pow(vector.y-y, 2) + pow(vector.z-z, 2));
}

template<class T>
T Vector<T>::magnitude_squared() const {
  return x*x + y*y + z*z;
}

template<class T>
Vector<T>::Vector(const T a, const T b, const T c)
  : x(a),
    y(b),
    z(c) {}

template<class T>
Vector<T> Vector<T>::operator + (const T val) const {
  return Vector<T>(x+val, y+val, z+val);
}

template<class T>
Vector<T> Vector<T>::operator - (const T val) const {
  return Vector<T>(x-val, y-val, z-val);
}

template<class T>
Vector<T> Vector<T>::operator * (const T val) const {
  return Vector<T>(x*val, y*val, z*val);
}

template<class T>
Vector<T> Vector<T>::operator / (const T val) const {
  return Vector<T>(x/val, y/val, z/val);
}

template<class T>
const Vector<T>& Vector<T>::operator += (const T val) {
  x += val;
  y += val;
  z += val;
  return *this;
}

template<class T>
const Vector<T>& Vector<T>::operator -= (const T val) {
  x -= val;
  y -= val;
  z -= val;
  return *this;
}

template<class T>
const Vector<T>& Vector<T>::operator *= (const T val) {
  x *= val;
  y *= val;
  z *= val;
  return *this;
}

template<class T>
const Vector<T>& Vector<T>::operator /= (const T val) {
  x /= val;
  y /= val;
  z /= val;
  return *this;
}

template<class T>
Vector<T> Vector<T>::operator + (const Vector<T>& vector) const {
  return Vector<T>(x+vector.x, y+vector.y, z+vector.z);
}

template<class T>
Vector<T> Vector<T>::operator - (const Vector<T>& vector) const {
  return Vector<T>(x-vector.x, y-vector.y, z-vector.z);
}

template<class T>
Vector<T> Vector<T>::operator * (const Vector<T>& vector) const {
  return Vector<T>(x*vector.x, y*vector.y, z*vector.z);
}

template<class T>
Vector<T> Vector<T>::operator / (const Vector<T>& vector) const {
  return Vector<T>(x/vector.x, y/vector.y, z/vector.z);
}

template<class T>
const Vector<T>& Vector<T>::operator += (const Vector<T>& vector) {
  x += vector.x;
  y += vector.y;
  z += vector.z;
  return *this;
}

template<class T>
const Vector<T>& Vector<T>::operator -= (const Vector<T>& vector) {
  x -= vector.x;
  y -= vector.y;
  z -= vector.z;
  return *this;
}

template<class T>
const Vector<T>& Vector<T>::operator *= (const Vector<T>& vector) {
  x *= vector.x;
  y *= vector.y;
  z *= vector.z;
  return *this;
}

template<class T>
const Vector<T>& Vector<T>::operator /= (const Vector<T>& vector) {
  x /= vector.x;
  y /= vector.y;
  z /= vector.z;
  return *this;
}

template<class T>
const bool Vector<T>::operator == (const Vector<T>& vector) const {
  return (x == vector.x && y == vector.y && z == vector.z);
}

template<class T>
const bool Vector<T>::operator != (const Vector<T>& vector) const {
  return (x != vector.x || y != vector.y || z != vector.z);
}

template<class T>
void Vector<T>::mod(const Vector<T>& vector) {
  x -= vector.x*floor(x/vector.x);
  y -= vector.y*floor(y/vector.y);
  z -= vector.z*floor(z/vector.z);
  /*
  x = fmod(x, vector.x);
  x = x >= 0 ? x : x + vector.x;
  y = fmod(y, vector.y);
  y = y >= 0 ? y : y + vector.y;
  z = fmod(z, vector.z);
  z = z >= 0 ? z : z + vector.z;
  */
}


#endif /* __Vector_hpp */
