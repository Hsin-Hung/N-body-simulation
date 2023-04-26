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

#include <algorithm>
#include <math.h>
#include "directSum.h"
#include "constants.h"

DirectSum::DirectSum(std::vector<std::shared_ptr<Body>> &bs) : Algorithm(bs, bs.size()) {}

void DirectSum::calculateAcceleration()
{

    for (int i = 0; i < nBodies; ++i)
    {
        Body &bi = *bodies[i];
        bi.acceleration = Vector(0, 0);
        Vector force(0, 0);
        for (int j = 0; j < nBodies; ++j)
        {
            Body &bj = *bodies[j];
            if (i != j && bi.isDynamic && !isCollide(bi, bj))
            {

                Vector rij = bj.position - bi.position;
                double r = sqrt((rij.x * rij.x) + (rij.y * rij.y) + (epsilon * epsilon));
                double f = (GRAVITY * bi.mass * bj.mass) / (r * r * r + (epsilon * epsilon));
                force += rij * f;
            }
        }

        bi.acceleration += (force / bi.mass);
    }
}

void DirectSum::calculateVelocity()
{

    for (auto &body : bodies)
    {
        body->velocity += (body->acceleration * dt);
    }
}

void DirectSum::calculatePosition()
{
    double boundaryWidth = NBODY_WIDTH, boundaryHeight = NBODY_HEIGHT;

    // check if body is at boundary
    for (auto &body : bodies)
    {
        body->position += body->velocity * dt;
    }
}

bool DirectSum::isCollide(Body b1, Body b2)
{
    return b1.radius + b2.radius + COLLISION_TH >= b1.position.getDistance(b2.position);
}

void DirectSum::update()
{
    calculateAcceleration();
    calculateVelocity();
    calculatePosition();
}