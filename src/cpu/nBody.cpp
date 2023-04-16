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

    for (int i = 0; i < nBodies; ++i)
    {
        int randPx = rand() % (1000 - 1 + 1) + 500;
        int randPy = rand() % (1000 - 1 + 1) + 500;
        Vector position = Vector(randPx, randPy);
        bodies.push_back(std::make_shared<Body>(1.0 / nBodies, 50, position, Vector(0, 0), Vector(0, 0)));
    }
}

void NBody::initSpiralBodies()
{

    bodies.clear();

    int maxDistance = 200;
    for (int i = 0; i < nBodies; ++i)
    {

        double angle = 2 * M_PI * (rand() / (double)RAND_MAX);
        // Generate random distance from center within the given max distance
        double distance = maxDistance * (rand() / (double)RAND_MAX);

        // Calculate coordinates of the point
        double x = CENTERX + distance * std::cos(angle);
        double y = CENTERY + distance * std::sin(angle);

        Vector position(x, y);
        Vector r(x - CENTERX, y - CENTERY);
        Vector velocity(r.y, -r.x);
        bodies.push_back(std::make_shared<Body>(1.0 / (double)nBodies, 1, position, velocity, Vector(0, 0)));
    }
}

void NBody::update()
{

    alg->update();
}