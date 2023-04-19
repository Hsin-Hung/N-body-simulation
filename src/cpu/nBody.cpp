#include <iostream>
#include <algorithm>
#include <stdlib.h>
#include <math.h>
#include "nBody.h"
#include "directSum.h"
#include "barnesHut.h"
#include "constants.h"

NBody::NBody(int n, int i) : nBodies(n)
{
    // initBodies();
    initSpiralBodies();
    if (i == 0)
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
    double maxDistance = 2.2790e11;
    double minDistance = 1.4960e11;
    for (int i = 0; i < nBodies - 1; ++i)
    {
        double angle = 2 * M_PI * (rand() / (double)RAND_MAX);
        // Generate random distance from center within the given max distance
        double radius = (maxDistance - minDistance) * (rand() / (double)RAND_MAX) + minDistance;

        // Calculate coordinates of the point
        double x = CENTERX + radius * std::cos(angle);
        double y = CENTERY + radius * std::sin(angle);
        Vector position(x, y);
        bodies.push_back(std::make_shared<Body>(5.974e24, 1.3927e6, position, Vector(0, 0), Vector(0, 0)));
    }
    bodies.push_back(std::make_shared<Body>(1.9890e30, 1.3927e6, Vector(CENTERX, CENTERY), Vector(0, 0), Vector(0, 0), false));
}

void NBody::initSpiralBodies()
{
    srand(time(NULL));
    bodies.clear();

    double maxDistance = 2.2790e11;
    double minDistance = 1.4960e11;
    Vector centerPos(CENTERX, CENTERY);
    for (int i = 0; i < nBodies - 1; ++i)
    {

        double angle = 2 * M_PI * (rand() / (double)RAND_MAX);
        // Generate random distance from center within the given max distance
        double radius = (maxDistance - minDistance) * (rand() / (double)RAND_MAX) + minDistance;

        // Calculate coordinates of the point
        double x = CENTERX + radius * std::cos(angle);
        double y = CENTERY + radius * std::sin(angle);

        Vector position(x, y);

        double distance = position.getDistance(centerPos);
        Vector r = position - centerPos;
        Vector a = r / distance;

        // Calculate velocity vector components
        double esc = sqrt((GRAVITY * 1.9891e30) / (distance));
        Vector velocity(-a.y * esc, a.x * esc);

        bodies.push_back(std::make_shared<Body>(5.974e24, 1.3927e6, position, velocity, Vector(0, 0)));
    }

    bodies.push_back(std::make_shared<Body>(1.9891e30, 1.3927e6, centerPos, Vector(0, 0), Vector(0, 0), false));
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
    // for (auto &b : bodies)
    // {

    //     std::cout << b->position << " " << b->velocity << " " << b->mass << std::endl;
    // }
}