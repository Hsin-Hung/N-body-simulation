#ifndef NBODY_H
#define NBODY_H
#include <vector>
#include <memory>
#include "vector.h"
#include "body.h"

class NBody
{
    const int n;

public:
    int timeSteps;
    std::vector<std::shared_ptr<Body>> bodies;

    NBody(const int n);
    void display();
    void initBodies();
};

#endif