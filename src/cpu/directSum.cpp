#include "directSum.h"
#include <algorithm>
#include <math.h>

DirectSum::DirectSum(std::vector<std::shared_ptr<Body>> &bodies) : bodies(bodies), n(bodies.size()) {}

void DirectSum::calculateCollision(Body &b1, Body &b2)
{
    double vx = (b1.mass * b1.position.x) + (b2.mass * b2.position.x) / (b1.mass + b2.mass);
    double vy = (b1.mass * b1.position.y) + (b2.mass * b2.position.y) / (b1.mass + b2.mass);
    b1.velocity = Vector(vx, vy);
    b1.mass = b1.mass + b2.mass;
    b1.radius = std::max(b1.radius, b2.radius);
}

void DirectSum::calculateAcceleration()
{

    for (int i = 0; i < n; ++i)
    {
        Body &bi = *bodies[i];
        bi.acceleration = Vector(0, 0);
        for (int j = 0; j < n; ++j)
        {
            Body &bj = *bodies[j];
            Vector rij = bj.position - bi.position;
            double inv_r3 = pow(rij.x * rij.x + rij.y * rij.y + epsilon * epsilon, -1.5);
            if (i != j)
            {
                if (isCollide(bi, bj))
                {
                    // calculateCollision(bi, bj);
                }
                else
                {
                    Vector force = rij * ((GRAVITY * bj.mass) / inv_r3);
                    bi.acceleration += (force / bi.mass);
                }
            }
        }
    }
}

void DirectSum::calculateVelocity()
{
    for (auto &body : bodies)
    {
        body->velocity += body->acceleration * dt/2.0;
    }
}

void DirectSum::calculatePosition()
{
    int boundaryWidth = 1600, boundaryHeight = 1600;
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

void DirectSum::update()
{
    calculateVelocity();
    calculatePosition();
    calculateAcceleration();
    calculateVelocity();
}

bool DirectSum::isCollide(Body b1, Body b2)
{

    return b1.radius + b2.radius > b1.position.getDistance(b2.position);
}
