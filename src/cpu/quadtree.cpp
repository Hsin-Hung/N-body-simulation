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
#include <iostream>
#include "quadtree.h"

QuadTree::QuadTree() : QuadTree(Vector(-1, -1), Vector(-1, -1))
{
}
QuadTree::QuadTree(Vector topL, Vector botR)
{
    topLeft = topL;
    botRight = botR;
    centerMass = Vector(-1, -1);
    totalMass = 0.0;
    isLeaf = true;
    b = nullptr;
    topLeftTree = nullptr;
    topRightTree = nullptr;
    botLeftTree = nullptr;
    botRightTree = nullptr;
}
QuadTree::~QuadTree()
{
}

int QuadTree::getQuadrant(Vector pos)
{

    if ((topLeft.x + botRight.x) / 2.0 >= pos.x)
    {
        // Indicates topLeftTree
        if ((topLeft.y + botRight.y) / 2.0 <= pos.y)
        {
            return 2;
        }

        // Indicates botLeftTree
        else
        {
            return 3;
        }
    }
    else
    {
        // Indicates topRightTree
        if ((topLeft.y + botRight.y) / 2.0 <= pos.y)
        {
            return 1;
        }

        // Indicates botRightTree
        else
        {
            return 4;
        }
    }
}
void QuadTree::insert(std::shared_ptr<Body> body)
{
    if (body == nullptr)
        return;

    if (!inBoundary(body->position))
    {
        std::cout << "ERROR: body out of bound" << std::endl;
        return;
    }

    // If node x does not contain a body, put the new body here.
    if (b == nullptr && isLeaf)
    {
        b = body;
        return;
    }

    if (b != nullptr)
    {
        if (b == body)
        {
            b = nullptr;
        }
        else
        {
            this->insert(b);
        }
    }

    isLeaf = false;

    int q = getQuadrant(body->position);

    if (q == 1)
    {
        if (topRightTree == nullptr)
            topRightTree = std::make_unique<QuadTree>(Vector((topLeft.x + botRight.x) / 2,
                                                             topLeft.y),
                                                      Vector(botRight.x,
                                                             (topLeft.y + botRight.y) / 2));

        topRightTree->insert(body);
    }
    else if (q == 2)
    {
        if (topLeftTree == nullptr)
            topLeftTree = std::make_unique<QuadTree>(Vector(topLeft.x, topLeft.y), Vector((topLeft.x + botRight.x) / 2,
                                                                                          (topLeft.y + botRight.y) / 2));

        topLeftTree->insert(body);
    }
    else if (q == 3)
    {
        if (botLeftTree == nullptr)
            botLeftTree = std::make_unique<QuadTree>(Vector(topLeft.x,
                                                            (topLeft.y + botRight.y) / 2),
                                                     Vector((topLeft.x + botRight.x) / 2,
                                                            botRight.y));

        botLeftTree->insert(body);
    }
    else
    {
        if (botRightTree == nullptr)
            botRightTree = std::make_unique<QuadTree>(Vector((topLeft.x + botRight.x) / 2,
                                                             (topLeft.y + botRight.y) / 2),
                                                      Vector(botRight.x, botRight.y));

        botRightTree->insert(body);
    }
}
std::shared_ptr<Body> QuadTree::search(Vector p)
{

    if (!inBoundary(p))
        return nullptr;

    if (b != nullptr)
        return b;

    int q = getQuadrant(p);

    if (q == 1)
    {
        if (topRightTree == nullptr)
            return nullptr;
        return topRightTree->search(p);
    }
    else if (q == 2)
    {
        if (topLeftTree == nullptr)
            return nullptr;
        return topLeftTree->search(p);
    }
    else if (q == 3)
    {
        if (botLeftTree == nullptr)
            return nullptr;
        return botLeftTree->search(p);
    }
    else
    {
        if (botRightTree == nullptr)
            return nullptr;
        return botRightTree->search(p);
    }
}
double getTotalMass(std::unique_ptr<QuadTree> &root)
{
    if (!root)
        return 0.0;
    return root->totalMass;
}
void updateCenterMass(std::unique_ptr<QuadTree> &root)
{
    if (!root)
        return;
    if (root->b)
    {
        root->totalMass = root->b->mass;
        root->centerMass = root->b->position;
        return;
    }

    updateCenterMass(root->topLeftTree);
    updateCenterMass(root->topRightTree);
    updateCenterMass(root->botLeftTree);
    updateCenterMass(root->botRightTree);

    double totalChildMass = getTotalMass(root->topLeftTree) + getTotalMass(root->topRightTree) + getTotalMass(root->botLeftTree) + getTotalMass(root->botRightTree);

    double totalCenterMassX = 0.0, totalCenterMassY = 0.0;
    if (root->topLeftTree)
    {
        totalCenterMassX += root->topLeftTree->centerMass.x * root->topLeftTree->totalMass;
        totalCenterMassY += root->topLeftTree->centerMass.y * root->topLeftTree->totalMass;
    }
    if (root->topRightTree)
    {
        totalCenterMassX += root->topRightTree->centerMass.x * root->topRightTree->totalMass;
        totalCenterMassY += root->topRightTree->centerMass.y * root->topRightTree->totalMass;
    }
    if (root->botLeftTree)
    {
        totalCenterMassX += root->botLeftTree->centerMass.x * root->botLeftTree->totalMass;
        totalCenterMassY += root->botLeftTree->centerMass.y * root->botLeftTree->totalMass;
    }
    if (root->botRightTree)
    {
        totalCenterMassX += root->botRightTree->centerMass.x * root->botRightTree->totalMass;
        totalCenterMassY += root->botRightTree->centerMass.y * root->botRightTree->totalMass;
    }
    root->totalMass = totalChildMass;
    root->centerMass = Vector(totalCenterMassX / totalChildMass, totalCenterMassY / totalChildMass);
}

bool QuadTree::inBoundary(Vector p)
{
    return (p.x >= topLeft.x && p.x <= botRight.x && p.y <= topLeft.y && p.y >= botRight.y);
}
double QuadTree::getWidth()
{
    return botRight.x - topLeft.x;
}
Vector QuadTree::getCenter()
{
    return centerMass;
}
void traverse(std::unique_ptr<QuadTree> &root)
{
    if (!root)
        return;

    if (root->b)
        std::cout << root->topLeft << " " << root->botRight << " " << root->totalMass << " " << root->centerMass << " " << *root->b << std::endl;
    else
        std::cout << root->topLeft << " " << root->botRight << " " << root->totalMass << " " << root->centerMass << " " << std::endl;

    traverse(root->topLeftTree);
    traverse(root->topRightTree);
    traverse(root->botLeftTree);
    traverse(root->botRightTree);
}
