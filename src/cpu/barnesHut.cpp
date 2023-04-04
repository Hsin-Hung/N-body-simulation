#include <math.h>
#include "barnesHut.h"
#include "constants.h"

BarnesHut::BarnesHut(std::vector<std::shared_ptr<Body>> &bs) : bodies(bs), n(bs.size()) {}

void BarnesHut::computeBoundingBox()
{
    Vector topLeft = Vector(INFINITY, -INFINITY), botRight = Vector(-INFINITY, INFINITY);
    for (auto &body : bodies)
    {
        topLeft.x = fminf(topLeft.x, body->position.x);
        topLeft.y = fmaxf(topLeft.y, body->position.y);
        botRight.x = fmaxf(botRight.x, body->position.x);
        botRight.y = fminf(botRight.y, body->position.y);
    }

    quadTree = std::make_unique<QuadTree>(topLeft, botRight);
}

void BarnesHut::constructQuadTree()
{
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

    if (root->b)
    {

        Body &bi = *body, &bj = *root->b;
        if (isCollide(bi, bj) || root->b == body)
            return;

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
        if (bi.radius * 2 > bi.position.getDistance(root->centerMass))
            return;
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
        if (body->isDynamic)
        {
            body->acceleration = Vector(0, 0);
            calculateForce(body);
        }
    }
}
void BarnesHut::calculateVelocity()
{
    double mag = 0.0;
    for (auto &body : bodies)
    {
        body->velocity += body->acceleration * dt / 2.0;
        mag = std::max(mag, sqrt(body->velocity.x * body->velocity.x + body->velocity.y * body->velocity.y));
    }
    std::cout << mag << std::endl;
}

void BarnesHut::calculatePosition()
{

    int boundaryWidth = WINDOW_WIDTH, boundaryHeight = WINDOW_HEIGHT;

    // check if body is at boundary
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
    computeBoundingBox();
    constructQuadTree();
    calculateAcceleration();
    calculateVelocity();
}
bool BarnesHut::isCollide(Body b1, Body b2)
{

    return b1.radius + b2.radius > b1.position.getDistance(b2.position);
}