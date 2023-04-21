#include "body.h"

std::ostream &operator<<(std::ostream &os, const Body &b)
{

    os << "Body(" << b.mass << "," << b.radius << "," << b.position << "," << b.velocity << "," << b.acceleration << "," << b.isDynamic;
    return os;
}