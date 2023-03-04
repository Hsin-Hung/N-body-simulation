#include <iostream>
#include <algorithm>
#include <math.h>
#include "nBody.h"

NBody::NBody()
{
    bodies.push_back(std::make_shared<Body>(1, 1, Vector(1, 8), Vector(0, 0), Vector(0, 0)));
    bodies.push_back(std::make_shared<Body>(1, 1, Vector(6, 9), Vector(0, 0), Vector(0, 0)));
    bodies.push_back(std::make_shared<Body>(1, 1, Vector(6, 6), Vector(0, 0), Vector(0, 0)));
    bodies.push_back(std::make_shared<Body>(1, 1, Vector(9, 6), Vector(0, 0), Vector(0, 0)));
    bodies.push_back(std::make_shared<Body>(1, 1, Vector(9, 1), Vector(0, 0), Vector(0, 0)));

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