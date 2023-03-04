#ifndef NBODY_H
#define NBODY_H
#include <vector>
#include <memory>
#include "vector.h"
#include "body.h"

class NBody
{
public:
    int n;
    int timeSteps;
    std::vector<std::shared_ptr<Body>> bodies;

    NBody();
    void simulate();
    void display();
};

#endif