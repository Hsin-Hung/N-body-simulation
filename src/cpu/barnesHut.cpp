#include <math.h>
#include "barnesHut.h"

BarnesHut::BarnesHut(std::vector<std::shared_ptr<Body>> &bodies) : bodies(bodies), n(bodies.size()) {}

void BarnesHut::constructQuadTree()
{
    quadTree = std::make_unique<QuadTree>(Vector(0, 1650), Vector(1650, 0));
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
        double inv_r3 = pow(rij.x * rij.x + rij.y * rij.y + epsilon * epsilon, -1.5);
        Vector force = rij * ((GRAVITY * bj.mass) / inv_r3);
        bi.acceleration += (force / bi.mass);
        return;
    }

    double sd = root->getWidth() / body->position.getDistance(root->centerMass);
    if (sd < THETA)
    {
        Body &bi = *body;
        Vector rij = root->centerMass - bi.position;
        double inv_r3 = pow(rij.x * rij.x + rij.y * rij.y + epsilon * epsilon, -1.5);
        Vector force = rij * ((GRAVITY * root->totalMass) / inv_r3);
        bi.acceleration += (force / bi.mass);
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
        body->velocity += body->acceleration * dt/2.0;
    }
}

void BarnesHut::calculatePosition()
{

    int boundaryWidth = 1600, boundaryHeight = 1600;
    for (auto &body : bodies)
    {
        body->position += body->velocity * dt;
        if (body->position.x < 0)
        {
            body->position.x = 0;
            body->velocity.x *= -1.0f;
        }
        if (body->position.x > boundaryWidth)
        {
            body->position.x = boundaryWidth;
            body->velocity.x *= -1.0f;
        }
        if (body->position.y < 0)
        {
            body->position.y = 0;
            body->velocity.y *= -1.0f;
        }
        if (body->position.y > boundaryHeight)
        {
            body->position.y = boundaryHeight;
            body->velocity.y *= -1.0f;
        }
    }
}
void BarnesHut::update()
{
    calculateVelocity();
    calculatePosition();
    constructQuadTree();
    calculateAcceleration();
    calculateVelocity();
}
