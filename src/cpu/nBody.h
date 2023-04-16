#ifndef NBODY_H
#define NBODY_H
#include <vector>
#include <memory>
#include "vector.h"
#include "body.h"
#include "algorithm.h"

class NBody
{
    std::unique_ptr<Algorithm> alg;
    int nBodies;
    void initRandomBodies();
    void initSpiralBodies();

public:
    std::vector<std::shared_ptr<Body>> bodies;
    NBody(int n, int i);
    void update();
};

#endif