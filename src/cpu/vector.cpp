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

#include <math.h>
#include <iomanip>
#include "vector.h"

Vector::Vector(const Vector &vt) : x(vt.x), y(vt.y)
{
}

Vector::Vector(Vector &&vt) noexcept : x(vt.x), y(vt.y)
{
    vt.x = 0;
    vt.y = 0;
}

double Vector::getDistance(Vector to)
{

    return sqrt(pow(to.x - x, 2) + pow(to.y - y, 2));
}

double Vector::mod()
{

    return sqrt(x * x + y * y);
}

double Vector::dot(Vector v)
{
    return x * v.x + y * v.y;
}

Vector &Vector::operator=(const Vector &vt)
{
    if (this != &vt)
    {
        x = vt.x;
        y = vt.y;
    }

    return *this;
}

Vector &Vector::operator=(Vector &&vt) noexcept
{
    if (this != &vt)
    {

        x = vt.x;
        y = vt.y;
        vt.x = 0;
        vt.y = 0;
    }
    return *this;
}

Vector Vector::operator+(const Vector &rhs)
{

    return Vector(x + rhs.x, y + rhs.y);
}

Vector Vector::operator-(const Vector &rhs)
{

    return Vector(x - rhs.x, y - rhs.y);
}

Vector Vector::operator*(double v)
{
    return Vector(x * v, y * v);
}

Vector Vector::operator*(Vector v)
{

    return Vector(x * v.x, y * v.y);
}

Vector Vector::operator/(double v)
{

    return Vector(x / v, y / v);
}

Vector &Vector::operator+=(const Vector &rhs)
{
    x += rhs.x;
    y += rhs.y;
    return *this;
}

Vector &Vector::operator-=(const Vector &rhs)
{
    x -= rhs.x;
    y -= rhs.y;
    return *this;
}

bool Vector::operator==(const Vector &rhs)
{
    return x == rhs.x && y == rhs.y;
}

bool Vector::operator!=(const Vector &rhs)
{
    return x != rhs.x || y != rhs.y;
}

std::ostream &operator<<(std::ostream &os, const Vector &vt)
{

    os << "(" << std::fixed << std::setprecision(6) << std::setw(9) << vt.x << ",";
    os << std::fixed << std::setprecision(6) << std::setw(9) << vt.y << ")";

    return os;
}