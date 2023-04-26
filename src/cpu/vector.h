/*
   Copyright 2023 Hsin-Hung Wu

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
*/

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
    double dot(Vector v);
    Vector operator+(const Vector &rhs);
    Vector operator-(const Vector &rhs);
    Vector operator*(double v);
    Vector operator*(Vector v);
    Vector operator/(double v);
    Vector &operator+=(const Vector &rhs);
    Vector &operator-=(const Vector &rhs);
    bool operator==(const Vector &rhs);
    bool operator!=(const Vector &rhs);
    Vector &operator=(const Vector &other);
    Vector &operator=(Vector &&other) noexcept;
    friend std::ostream &operator<<(std::ostream &os, const Vector &vt);
};

#endif