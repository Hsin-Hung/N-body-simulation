#include <math.h>
#include "barnesHut.h"

BarnesHut::BarnesHut(std::vector<std::shared_ptr<Body>> &bodies) : bodies(bodies), n(bodies.size()) {}

void BarnesHut::constructQuadTree()
{
    quadTree = std::make_unique<QuadTree>(Vector(0, 10), Vector(10, 0));
    for (auto &body : bodies)
    {
        quadTree->insert(body);
    }
    updateCM(quadTree);

}
void BarnesHut::calculateForceHelper(std::unique_ptr<QuadTree> &root, std::shared_ptr<Body> body)
{
    if (!root)
        return;

    if (root->b && root->b != body)
    {
        Body &bi = *body, &bj = *root->b;
        Vector rij = bj.position - bi.position;
        Vector force = rij * ((GRAVITY * bi.mass * bj.mass) / pow((bi.position - bj.position).mod(), 3));
        bi.acceleration = bi.acceleration + force;
        return;
    }

    double sd = root->getWidth() / body->position.getDistance(root->centerMass);
    if (sd < THETA)
    {
        Body &bi = *body;
        Vector rij = root->centerMass - bi.position;
        Vector force = rij * ((GRAVITY * bi.mass * root->totalMass) / pow((bi.position - root->centerMass).mod(), 3));
        bi.acceleration = bi.acceleration + force;
        return;
    }

    calculateForceHelper(root->topLeftTree, body);
    calculateForceHelper(root->topRightTree, body);
    calculateForceHelper(root->botLeftTree, body);
    calculateForceHelper(root->botRightTree, body);
}

void BarnesHut::calculateForce(std::shared_ptr<Body> b)
{
    calculateForceHelper(quadTree, b);
}
void BarnesHut::calculateAcceleration()
{
    for (auto &body : bodies)
    {
        calculateForce(body);
    }
}
void BarnesHut::calculateVelocity()
{
    for (auto &body : bodies)
    {
        body->velocity = body->velocity + body->acceleration;
    }
}

void BarnesHut::calculatePosition()
{

    for (auto &body : bodies)
    {
        body->position = body->position + body->velocity + body->acceleration * 0.5;
    }
}
void BarnesHut::update()
{
    constructQuadTree();
    calculateAcceleration();
    calculatePosition();
    calculateVelocity();
}
