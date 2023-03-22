#include "quadtree.h"
#include <math.h>
#include <iostream>

QuadTree::QuadTree()
{
    topLeft = Vector(0, 0);
    botRight = Vector(0, 0);
    centerMass = Vector(0, 0);
    totalMass = 0.0;
    b = nullptr;
    topLeftTree = nullptr;
    topRightTree = nullptr;
    botLeftTree = nullptr;
    botRightTree = nullptr;
}
QuadTree::QuadTree(Vector topL, Vector botR)
{
    topLeft = topL;
    botRight = botR;
    centerMass = Vector(0, 0);
    totalMass = 0.0;
    b = nullptr;
    topLeftTree = nullptr;
    topRightTree = nullptr;
    botLeftTree = nullptr;
    botRightTree = nullptr;
}
QuadTree::~QuadTree()
{
}
void QuadTree::insert(std::shared_ptr<Body> body)
{
    if (body == nullptr)
        return;

    if (!inBoundary(body->position))
        return;
    // If node x does not contain a body, put the new body here.
    if (b == nullptr && topLeftTree == nullptr && topRightTree == nullptr && botLeftTree == nullptr && botRightTree == nullptr)
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
    if ((topLeft.x + botRight.x) / 2 >= body->position.x)
    {
        // Indicates topLeftTree
        if ((topLeft.y + botRight.y) / 2 <= body->position.y)
        {
            if (topLeftTree == nullptr)
                topLeftTree = std::make_unique<QuadTree>(Vector(topLeft.x, topLeft.y), Vector((topLeft.x + botRight.x) / 2,
                                                                                              (topLeft.y + botRight.y) / 2));

            topLeftTree->insert(body);
        }

        // Indicates botLeftTree
        else
        {
            if (botLeftTree == nullptr)
                botLeftTree = std::make_unique<QuadTree>(Vector(topLeft.x,
                                                                (topLeft.y + botRight.y) / 2),
                                                         Vector((topLeft.x + botRight.x) / 2,
                                                                botRight.y));

            botLeftTree->insert(body);
        }
    }
    else
    {
        // Indicates topRightTree
        if ((topLeft.y + botRight.y) / 2 <= body->position.y)
        {
            if (topRightTree == nullptr)
                topRightTree = std::make_unique<QuadTree>(Vector((topLeft.x + botRight.x) / 2,
                                                                 topLeft.y),
                                                          Vector(botRight.x,
                                                                 (topLeft.y + botRight.y) / 2));

            topRightTree->insert(body);
        }

        // Indicates botRightTree
        else
        {
            if (botRightTree == nullptr)
                botRightTree = std::make_unique<QuadTree>(Vector((topLeft.x + botRight.x) / 2,
                                                                 (topLeft.y + botRight.y) / 2),
                                                          Vector(botRight.x, botRight.y));

            botRightTree->insert(body);
        }
    }
}
std::shared_ptr<Body> QuadTree::search(Vector p)
{

    // Current quad cannot contain it
    if (!inBoundary(p))
        return nullptr;

    // We are at a quad of unit length
    // We cannot subdivide this quad further
    if (b != nullptr)
        return b;

    if ((topLeft.x + botRight.x) / 2 >= p.x)
    {
        // Indicates topLeftTree
        if ((topLeft.y + botRight.y) / 2 <= p.y)
        {
            if (topLeftTree == nullptr)
                return nullptr;
            return topLeftTree->search(p);
        }

        // Indicates botLeftTree
        else
        {
            if (botLeftTree == nullptr)
                return nullptr;
            return botLeftTree->search(p);
        }
    }
    else
    {
        // Indicates topRightTree
        if ((topLeft.y + botRight.y) / 2 <= p.y)
        {
            if (topRightTree == nullptr)
                return nullptr;
            return topRightTree->search(p);
        }

        // Indicates botRightTree
        else
        {
            if (botRightTree == nullptr)
                return nullptr;
            return botRightTree->search(p);
        }
    }
}
double getTotalMass(std::unique_ptr<QuadTree> &root)
{
    if (!root)
        return 0.0;
    return root->totalMass;
}
void updateCM(std::unique_ptr<QuadTree> &root)
{
    if (!root)
        return;
    if (root->b)
    {
        root->totalMass = root->b->mass;
        root->centerMass = root->b->position;
        return;
    }

    updateCM(root->topLeftTree);
    updateCM(root->topRightTree);
    updateCM(root->botLeftTree);
    updateCM(root->botRightTree);

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
        std::cout << root->topLeft << " " << root->botRight << " " << root->totalMass << " " << root->centerMass << " " <<  *root->b << std::endl;
    else
        std::cout << root->topLeft << " " << root->botRight << " " << root->totalMass << " " << root->centerMass << " " << std::endl;

    traverse(root->topLeftTree);
    traverse(root->topRightTree);
    traverse(root->botLeftTree);
    traverse(root->botRightTree);
}
