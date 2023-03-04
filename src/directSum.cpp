#include "directSum.h"
#include <algorithm>
#include <math.h>

DirectSum::DirectSum(std::vector<std::shared_ptr<Body>> &bodies) : bodies(bodies), n(bodies.size()) {}

void DirectSum::resolveCollision()
{

    for (int i = 0; i < n; ++i)
    {
        for (int j = i + 1; j < n; ++j)
        {
            if (bodies[i]->position == bodies[j]->position)
            {
                std::swap(bodies[i]->velocity, bodies[j]->velocity);
            }
        }
    }
}

void DirectSum::calculateAcceleration()
{

    for (int i = 0; i < n; ++i)
    {
        Body &bi = *bodies[i];
        bi.acceleration = Vector(0, 0);
        for (int j = 0; j < n; ++j)
        {
            if (i != j)
            {
                Body &bj = *bodies[j];
                Vector rij = bj.position - bi.position;
                Vector force = rij * ((GRAVITY * bi.mass * bj.mass) / pow((bi.position - bj.position).mod(), 3));
                bi.acceleration = bi.acceleration + force;
            }
        }
    }
}

void DirectSum::calculateVelocity()
{
    for (auto &body : bodies)
    {
        body->velocity = body->velocity + body->acceleration;
    }
}

void DirectSum::calculatePosition()
{

    for (auto &body : bodies)
    {
        body->position = body->position + body->velocity + body->acceleration * 0.5;
    }
}
void DirectSum::update()
{
    calculateAcceleration();
    calculatePosition();
    calculateVelocity();
    resolveCollision();
}