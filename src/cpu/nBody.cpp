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

#include <iostream>
#include <algorithm>
#include <stdlib.h>
#include <math.h>
#include "nBody.h"
#include "directSum.h"
#include "barnesHut.h"
#include "constants.h"

NBody::NBody(int n, int a, int s) : nBodies(n)
{

    if (s == 0)
    {
        initSpiralBodies();
    }
    else if (s == 1)
    {
        initRandomBodies();
    }
    else
    {
        nBodies = 5;
        initSolarSystem();
    }

    if (a == 0)
    {
        alg = std::make_unique<DirectSum>(bodies);
    }
    else
    {
        alg = std::make_unique<BarnesHut>(bodies);
    }
}

void NBody::initRandomBodies()
{
    srand(time(NULL));

    bodies.clear();
    double maxDistance = MAX_DIST;
    double minDistance = MIN_DIST;
    Vector centerPos(CENTERX, CENTERY);
    for (int i = 0; i < nBodies - 1; ++i)
    {
        double angle = 2 * M_PI * (rand() / (double)RAND_MAX);
        // Generate random distance from center within the given max distance
        double radius = (maxDistance - minDistance) * (rand() / (double)RAND_MAX) + minDistance;

        // Calculate coordinates of the point
        double x = centerPos.x + radius * std::cos(angle);
        double y = centerPos.y + radius * std::sin(angle);
        Vector position(x, y);
        bodies.push_back(std::make_shared<Body>(EARTH_MASS, EARTH_DIA, position, Vector(0, 0), Vector(0, 0)));
    }
    bodies.push_back(std::make_shared<Body>(SUN_MASS, SUN_DIA, centerPos, Vector(0, 0), Vector(0, 0), false));
}

void NBody::initSpiralBodies()
{
    srand(time(NULL));

    bodies.clear();
    double maxDistance = MAX_DIST;
    double minDistance = MIN_DIST;
    Vector centerPos(CENTERX, CENTERY);
    for (int i = 0; i < nBodies - 1; ++i)
    {

        double angle = 2 * M_PI * (rand() / (double)RAND_MAX);
        // Generate random distance from center within the given max distance
        double radius = (maxDistance - minDistance) * (rand() / (double)RAND_MAX) + minDistance;

        // Calculate coordinates of the point
        double x = centerPos.x + radius * std::cos(angle);
        double y = centerPos.y + radius * std::sin(angle);

        Vector position(x, y);
        double distance = position.getDistance(centerPos);
        Vector r = position - centerPos;
        Vector a = r / distance;

        // Calculate velocity vector components
        double esc = sqrt((GRAVITY * SUN_MASS) / (distance));
        Vector velocity(-a.y * esc, a.x * esc);

        bodies.push_back(std::make_shared<Body>(EARTH_MASS, EARTH_DIA, position, velocity, Vector(0, 0)));
    }

    bodies.push_back(std::make_shared<Body>(SUN_MASS, SUN_DIA, centerPos, Vector(0, 0), Vector(0, 0), false));
}

void NBody::initSolarSystem()
{
    bodies.clear();

    bodies.push_back(std::make_shared<Body>(5.9740e24, 1.3927e6, Vector(1.4960e11, 0), Vector(0, 2.9800e4), Vector(0, 0)));
    bodies.push_back(std::make_shared<Body>(6.4190e23, 1.3927e6, Vector(2.2790e11, 0), Vector(0, 2.4100e4), Vector(0, 0)));
    bodies.push_back(std::make_shared<Body>(3.3020e23, 1.3927e6, Vector(5.7900e10, 0), Vector(0, 4.7900e4), Vector(0, 0)));
    bodies.push_back(std::make_shared<Body>(4.8690e24, 1.3927e6, Vector(1.0820e11, 0), Vector(0, 3.5000e4), Vector(0, 0)));
    bodies.push_back(std::make_shared<Body>(1.9890e30, 1.3927e6, Vector(CENTERX, CENTERY), Vector(0, 0), Vector(0, 0), false));
}

void NBody::update()
{
    alg->update();
}