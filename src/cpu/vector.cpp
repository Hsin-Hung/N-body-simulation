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
double Vector::mod()
{

    return sqrt(x * x + y * y);
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