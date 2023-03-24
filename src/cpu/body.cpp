#include "body.h"

std::ostream &operator<<(std::ostream &os, const Body &b)
{

    os << "Body(" << b.mass << ",(" << b.position << "))";
    return os;
}