#ifndef DIRECT_SUM_H
#define DIRECT_SUM_H
#include <vector>
#include <memory>
#include "body.h"
#include "constants.h"

class DirectSum
{
    int n;
    std::vector<std::shared_ptr<Body>> &bodies;

    void resolveCollision();
    void calculateAcceleration();
    void calculateVelocity();
    void calculatePosition();

public:
    DirectSum(std::vector<std::shared_ptr<Body>> &bodies);
    void update();
};

#endif