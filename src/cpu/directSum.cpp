#include "directSum.h"
#include "constants.h"
#include <algorithm>
#include <math.h>

DirectSum::DirectSum(std::vector<std::shared_ptr<Body>> &bs) : Algorithm(bs, bs.size()) {}

void DirectSum::calculateAcceleration()
{

    for (int i = 0; i < nBodies; ++i)
    {
        Body &bi = *bodies[i];
        bi.acceleration = Vector(0, 0);
        for (int j = 0; j < nBodies; ++j)
        {
            Body &bj = *bodies[j];
            Vector rij = bj.position - bi.position;
            double inv_r3 = pow(rij.x * rij.x + rij.y * rij.y + epsilon * epsilon, -1.5);
            if (i != j && bi.isDynamic && !isCollide(bi, bj))
            {
                Vector force = rij * ((GRAVITY * bj.mass) / inv_r3);
                bi.acceleration += (force / bi.mass);
            }
        }
    }
}

void DirectSum::calculateVelocity()
{

    for (auto &body : bodies)
    {
        body->velocity += body->acceleration * dt / 2.0;
    }
}

void DirectSum::calculatePosition()
{
    int boundaryWidth = WINDOW_WIDTH, boundaryHeight = WINDOW_HEIGHT;

    // check if body is at boundary
    for (auto &body : bodies)
    {
        body->position += body->velocity * dt;
        if (body->position.x < 0)
        {
            body->position.x = 0;
            body->velocity.x *= -1.0f;
        }
        if (body->position.x > boundaryWidth)
        {
            body->position.x = boundaryWidth;
            body->velocity.x *= -1.0f;
        }
        if (body->position.y < 0)
        {
            body->position.y = 0;
            body->velocity.y *= -1.0f;
        }
        if (body->position.y > boundaryHeight)
        {
            body->position.y = boundaryHeight;
            body->velocity.y *= -1.0f;
        }
    }
}

bool DirectSum::isCollide(Body b1, Body b2)
{

    return b1.radius + b2.radius > b1.position.getDistance(b2.position);
}

void DirectSum::update()
{
    calculateVelocity();
    calculatePosition();
    calculateAcceleration();
    calculateVelocity();
}