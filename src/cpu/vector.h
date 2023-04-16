#ifndef VECTOR_H
#define VECTOR_H
#include <iostream>

class Vector
{
public:
    double x;
    double y;

    Vector() : Vector(0, 0) {}
    Vector(double x, double y) : x(x), y(y) {}
    Vector(const Vector &vt);
    Vector(Vector &&vt) noexcept;
    double getDistance(Vector to);
    double mod();
    Vector operator+(const Vector &rhs);
    Vector operator-(const Vector &rhs);
    Vector operator*(double v);
    Vector operator*(Vector v);
    Vector operator/(double v);
    Vector &operator+=(const Vector &rhs);
    bool operator==(const Vector &rhs);
    bool operator!=(const Vector &rhs);
    Vector &operator=(const Vector &other);
    Vector &operator=(Vector &&other) noexcept;
    friend std::ostream &operator<<(std::ostream &os, const Vector &vt);
};

#endif