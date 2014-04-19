#pragma once

#include <iostream>

class Vector3D {
 public:
  double x, y, z;

  explicit Vector3D(double x = 0.0, double y = 0.0, double z = 0.0);
  Vector3D(const Vector3D &a);

  static double dot(const Vector3D &a, const Vector3D &b);
  static Vector3D cross(const Vector3D &a, const Vector3D &b);
  static double angle(const Vector3D &a, const Vector3D &b);
  static Vector3D normal(const Vector3D &a, const Vector3D &b);

  double len() const;
  Vector3D normalized() const;
  Vector3D operator -() const;
  Vector3D operator +() const;
  double operator [](int i) const;

  friend std::ostream & operator <<(std::ostream &stream, const Vector3D &t);
};

bool operator ==(const Vector3D &a, const Vector3D &b);
bool operator !=(const Vector3D &a, const Vector3D &b);

Vector3D operator +(const Vector3D &a, const Vector3D &b);
Vector3D operator -(const Vector3D &a, const Vector3D &b);

Vector3D operator *(double a, const Vector3D &b);
Vector3D operator *(const Vector3D &a, double b);
Vector3D operator *(const Vector3D &a, const Vector3D &b);
Vector3D operator /(const Vector3D &a, double b);
Vector3D operator /(const Vector3D &a, const Vector3D &b);

std::ostream & operator <<(std::ostream &stream, const Vector3D &t);
