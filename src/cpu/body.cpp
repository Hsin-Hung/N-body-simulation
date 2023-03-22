#include "body.h"

void Body::UpdatePosition(Vector p) {}
void Body::UpdateVelocity(Vector v) {}
std::ostream &operator<<(std::ostream &os, const Body &b)
{

    os << "Body(" << b.mass << ",(" << b.position << "))";
    return os;
}