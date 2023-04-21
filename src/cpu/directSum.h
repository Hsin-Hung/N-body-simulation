#ifndef DIRECT_SUM_H
#define DIRECT_SUM_H
#include <vector>
#include <memory>
#include "body.h"
#include "constants.h"
#include "algorithm.h"

class DirectSum : public Algorithm
{
    const double epsilon = 0.5;
    const double dt = 25000.0;

    void calculateAcceleration();
    void calculateVelocity();
    void calculatePosition();
    bool isCollide(Body b1, Body b2);

public:
    DirectSum(std::vector<std::shared_ptr<Body>> &bs);
    void update() override;
};

#endif