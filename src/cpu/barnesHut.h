#ifndef BARNES_HUT_H
#define BARNES_HUT_H
#include <vector>
#include <memory>
#include "body.h"
#include "quadtree.h"
#include "constants.h"
#include "algorithm.h"

#define THETA 0.5

class BarnesHut : public Algorithm
{
    const double epsilon = 0.5;
    const double dt = 0.005;

    std::unique_ptr<QuadTree> quadTree;
    void constructQuadTree();
    void computeCenterMass();
    void calculateForceHelper(std::unique_ptr<QuadTree> &root, std::shared_ptr<Body> body);
    void computeBoundingBox();
    void calculateForce(std::shared_ptr<Body> b);
    void calculateAcceleration();
    void calculateVelocity();
    void calculatePosition();
    bool isCollide(Body b1, Body b2);

public:
    BarnesHut(std::vector<std::shared_ptr<Body>> &bs);
    void update() override;
};

#endif