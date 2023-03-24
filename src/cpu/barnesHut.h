#ifndef BARNES_HUT_H
#define BARNES_HUT_H
#include <vector>
#include <memory>
#include "body.h"
#include "quadtree.h"
#include "constants.h"
#define THETA 0.3
class BarnesHut
{
    const double epsilon = 0.5;
    const double dt = 0.005;
    const int n;

    std::vector<std::shared_ptr<Body>> &bodies;
    void constructQuadTree();
    void calculateForceHelper(std::unique_ptr<QuadTree> &root, std::shared_ptr<Body> body);
    void calculateForce(std::shared_ptr<Body> b);
    void calculateAcceleration();
    void calculateVelocity();
    void calculatePosition();
    bool isCollide(Body b1, Body b2);

public:
    std::unique_ptr<QuadTree> quadTree;
    BarnesHut(std::vector<std::shared_ptr<Body>> &bs);
    void update();
};

#endif