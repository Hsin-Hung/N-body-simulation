/*
   Copyright 2023 Hsin-Hung Wu

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
*/

#include <math.h>
#include "barnesHut.h"
#include "constants.h"

BarnesHut::BarnesHut(std::vector<std::shared_ptr<Body>> &bs) : Algorithm(bs, bs.size()) {}

void BarnesHut::computeBoundingBox()
{
    Vector topLeft = Vector(INFINITY, -INFINITY), botRight = Vector(-INFINITY, INFINITY);
    for (auto &body : bodies)
    {
        topLeft.x = fminf(topLeft.x, body->position.x - 1.0e10);
        topLeft.y = fmaxf(topLeft.y, body->position.y + 1.0e10);
        botRight.x = fmaxf(botRight.x, body->position.x + 1.0e10);
        botRight.y = fminf(botRight.y, body->position.y - 1.0e10);
    }

    quadTree = std::make_unique<QuadTree>(topLeft, botRight);
}

void BarnesHut::constructQuadTree()
{
    for (auto &body : bodies)
    {
        quadTree->insert(body);
    }
}

void BarnesHut::computeCenterMass()
{
    updateCenterMass(quadTree);
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

        double r = sqrt((rij.x * rij.x) + (rij.y * rij.y) + (epsilon * epsilon));
        Vector force = rij * ((GRAVITY * bi.mass * bj.mass) / (r * r * r + (epsilon * epsilon)));
        bi.acceleration += (force / bi.mass);
        return;
    }

    double sd = root->getWidth() / body->position.getDistance(root->centerMass);
    if (sd < theta)
    {
        Body &bi = *body;
        Vector rij = root->centerMass - bi.position;
        if (!isCollide(bi, root->centerMass))
        {
            double r = sqrt((rij.x * rij.x) + (rij.y * rij.y) + (epsilon * epsilon));
            Vector force = rij * ((GRAVITY * bi.mass * root->totalMass) / (r * r * r + (epsilon * epsilon)));
            bi.acceleration += (force / bi.mass);
        }

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
    for (auto &body : bodies)
    {
        body->velocity += body->acceleration * dt;
    }
}

void BarnesHut::calculatePosition()
{
    // check if body is at boundary
    for (auto &body : bodies)
    {
        body->position += body->velocity * dt;
    }
}

bool BarnesHut::isCollide(Body b1, Body b2)
{

    return b1.radius + b2.radius + COLLISION_TH > b1.position.getDistance(b2.position);
}

bool BarnesHut::isCollide(Body b, Vector cm)
{
    return b.radius * 2 + COLLISION_TH > b.position.getDistance(cm);
}

void BarnesHut::update()
{
    computeBoundingBox();
    constructQuadTree();
    computeCenterMass();
    calculateAcceleration();
    calculateVelocity();
    calculatePosition();
}
