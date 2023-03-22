#include <iostream>
#include <algorithm>
#include <stdlib.h>
#include <math.h>
#include "nBody.h"

NBody::NBody()
{
    initBodies();
}

void NBody::initBodies()
{   
    srand(time(NULL));
    for (int i = 0; i < n; ++i)
    {
        int randX = rand() % (1000 - 1 + 1) + 500;
        int randY = rand() % (1000 - 1 + 1) + 500;
        bodies.push_back(std::make_shared<Body>(1.899e12, 10, Vector(randX, randY), Vector(0, 0), Vector(0, 0)));
    }

}

void NBody::simulate()
{
}
void NBody::display()
{
    for (auto &body : bodies)
    {
        std::cout << body->position << " ";
    }
    std::cout << std::endl;
}