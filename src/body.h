#ifndef BODY_H
#define BODY_H

#include "vector.h"
#include <iostream>

class Body
{
public:
    double mass;
    double radius;
    Vector position;
    Vector velocity;
    Vector acceleration;

    Body(double m, double r, Vector p, Vector v, Vector a) : mass(m), radius(r), position(p), velocity(v), acceleration(a) {}

    void UpdatePosition(Vector p);
    void UpdateVelocity(Vector v);
    friend std::ostream &operator<<(std::ostream &os, const Body &b);
};

#endif