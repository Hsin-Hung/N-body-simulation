#include <iostream>
#include <algorithm>
#include <stdlib.h>
#include <math.h>
#include "nBody.h"
#include "constants.h"

NBody::NBody(const int n) : n(n)
{
    initBodies();
}

void NBody::initBodies()
{
    srand(time(NULL));

    bodies.push_back(std::make_shared<Body>(200.0 / n, 50, Vector(CENTERX, CENTERY), Vector(0, 0), Vector(0, 0), false));

    for (int i = 0; i < n; ++i)
    {
        int randPx = rand() % (1000 - 1 + 1) + 500;
        int randPy = rand() % (1000 - 1 + 1) + 500;
        // int randVx = rand() % (500 - 1 + 1) + 1;
        // int randVy = rand() % (500 - 1 + 1) + 1;
        Vector position = Vector(randPx, randPy);
        Vector r = Vector(randPx - CENTERX, randPy - CENTERY);
        Vector velocity = Vector(r.y, -r.x);
        bodies.push_back(std::make_shared<Body>(1.0 / n, 50, position, velocity, Vector(0, 0)));
    }
}

void NBody::display()
{
    for (auto &body : bodies)
    {
        std::cout << body->velocity << " ";
    }
    std::cout << std::endl;
}