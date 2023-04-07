#include <iostream>
#include "barnesHut_kernel.cuh"
#include "constants.h"

dim3 gridSize = GRID_SIZE;
dim3 blockSize = BLOCK_SIZE;

BarnesHutCuda::BarnesHutCuda(int n) : nBodies(n)
{
    nNodes = MAX_NODES;
    leaf_limit = MAX_NODES - N_LEAF;
    h_b = new Body[n];
    h_node = new Node[nNodes];

    // cudaMalloc((void **)&d_b, sizeof(Body) * n);
    // cudaMalloc((void **)&d_node, sizeof(Node) * nNodes);
    // cudaMalloc((void **)&d_mutex, sizeof(int));
}

BarnesHutCuda::~BarnesHutCuda()
{

    delete[] h_b;
    delete[] h_node;
    // cudaFree(d_b);
    // cudaFree(d_node);
    // cudaFree(d_mutex);
}

void BarnesHutCuda::computeBoundingBox()
{
    Vector topLeft = {INFINITY, -INFINITY}, botRight = {-INFINITY, INFINITY};
    for (int i = 0; i < nBodies; ++i)
    {
        topLeft.x = fminf(topLeft.x, h_b[i].position.x);
        topLeft.y = fmaxf(topLeft.y, h_b[i].position.y);
        botRight.x = fmaxf(botRight.x, h_b[i].position.x);
        botRight.y = fminf(botRight.y, h_b[i].position.y);
    }
    h_node[0].topLeft = topLeft;
    h_node[0].botRight = botRight;
}

int BarnesHutCuda::getQuadrant(Vector topLeft, Vector botRight, double x, double y)
{

    if ((topLeft.x + botRight.x) / 2 >= x)
    {
        // Indicates topLeftTree
        if ((topLeft.y + botRight.y) / 2 <= y)
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
        if ((topLeft.y + botRight.y) / 2 <= y)
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

bool BarnesHutCuda::inBoundary(Vector topLeft, Vector botRight, Vector p)
{
    return (p.x >= topLeft.x && p.x <= botRight.x && p.y <= topLeft.y && p.y >= botRight.y);
}

void BarnesHutCuda::insertQuadTree(int nodeIndex, int b)
{

    if (nodeIndex >= nNodes)
    {
        std::cout << "index out of bound" << std::endl;
        return;
    }

    Node &curNode = h_node[nodeIndex];
    Body &body = h_b[b];

    if (!inBoundary(curNode.topLeft, curNode.botRight, body.position))
    {
        std::cout << "wrong quat" << std::endl;
        return;
    }

    if (nodeIndex >= leaf_limit)
    {

        if (curNode.centerMass.x != -1)
        {

            double M = curNode.totalMass + body.mass;
            double Rx = (curNode.totalMass * curNode.centerMass.x + body.mass * curNode.centerMass.x) / M;
            double Ry = (curNode.totalMass * curNode.centerMass.y + body.mass * curNode.centerMass.y) / M;
            curNode.totalMass = M;
            curNode.centerMass = {Rx, Ry};
        }
        else
        {
            curNode.centerMass = body.position;
            curNode.totalMass = body.mass;
        }

        return;
    }

    // If node x does not contain a body, put the new body here.
    if (curNode.isLeaf)
    {

        if (curNode.centerMass.x != -1)
        {

            int quadrant = getQuadrant(curNode.topLeft, curNode.botRight, curNode.centerMass.x, curNode.centerMass.y);
            Node &childNode = h_node[((nodeIndex * 4) + quadrant)];

            if (quadrant == 1)
            {
                childNode.topLeft = {(curNode.topLeft.x + curNode.botRight.x) / 2, curNode.topLeft.y};
                childNode.botRight = {curNode.botRight.x, (curNode.topLeft.y + curNode.botRight.y) / 2};
            }
            else if (quadrant == 2)
            {
                childNode.topLeft = {curNode.topLeft.x, curNode.topLeft.y};
                childNode.botRight = {(curNode.topLeft.x + curNode.botRight.x) / 2, (curNode.topLeft.y + curNode.botRight.y) / 2};
            }
            else if (quadrant == 3)
            {
                childNode.topLeft = {curNode.topLeft.x, (curNode.topLeft.y + curNode.botRight.y) / 2};
                childNode.botRight = {(curNode.topLeft.x + curNode.botRight.x) / 2, curNode.botRight.y};
            }
            else
            {
                childNode.topLeft = {(curNode.topLeft.x + curNode.botRight.x) / 2, (curNode.topLeft.y + curNode.botRight.y) / 2};
                childNode.botRight = {curNode.botRight.x, curNode.botRight.y};
            }

            childNode.centerMass = curNode.centerMass;
            childNode.totalMass = curNode.totalMass;

            curNode.centerMass = {-1, -1};
            curNode.totalMass = 0.0;
            curNode.isLeaf = false;
        }
        else
        {
            curNode.centerMass = body.position;
            curNode.totalMass = body.mass;
            return;
        }
    }

    int quadrant = getQuadrant(curNode.topLeft, curNode.botRight, body.position.x, body.position.y);
    Node &childNode = h_node[((nodeIndex * 4) + quadrant)];
    if (quadrant == 1)
    {
        childNode.topLeft = {(curNode.topLeft.x + curNode.botRight.x) / 2, curNode.topLeft.y};
        childNode.botRight = {curNode.botRight.x, (curNode.topLeft.y + curNode.botRight.y) / 2};
    }
    else if (quadrant == 2)
    {
        childNode.topLeft = {curNode.topLeft.x, curNode.topLeft.y};
        childNode.botRight = {(curNode.topLeft.x + curNode.botRight.x) / 2, (curNode.topLeft.y + curNode.botRight.y) / 2};
    }
    else if (quadrant == 3)
    {
        childNode.topLeft = {curNode.topLeft.x, (curNode.topLeft.y + curNode.botRight.y) / 2};
        childNode.botRight = {(curNode.topLeft.x + curNode.botRight.x) / 2, curNode.botRight.y};
    }
    else
    {
        childNode.topLeft = {(curNode.topLeft.x + curNode.botRight.x) / 2, (curNode.topLeft.y + curNode.botRight.y) / 2};
        childNode.botRight = {curNode.botRight.x, curNode.botRight.y};
    }

    // std::cout << "Cur Node TL(" << curNode.topLeft.x << "," << curNode.topLeft.y << ") ";
    // std::cout << "BR(" << curNode.botRight.x << "," << curNode.botRight.y << ") ";
    // std::cout << "Body Pos (" << body.position.x << "," << body.position.y << ") ";
    // std::cout << "Quat " << quadrant << std::endl;

    insertQuadTree((nodeIndex * 4) + quadrant, b);
}

void BarnesHutCuda::constructQuadTree()
{

    for (int i = 0; i < nBodies; ++i)
    {
        insertQuadTree(0, i);
    }
}
double BarnesHutCuda::getTotalMass(int nodeIndex)
{

    return h_node[nodeIndex].totalMass;
}

void BarnesHutCuda::computeCenterMass(int nodeIndex)
{

    if (nodeIndex >= nNodes)
        return;

    Node &curNode = h_node[nodeIndex];

    if (curNode.isLeaf)
    {
        return;
    }

    computeCenterMass((nodeIndex * 4) + 1);
    computeCenterMass((nodeIndex * 4) + 2);
    computeCenterMass((nodeIndex * 4) + 3);
    computeCenterMass((nodeIndex * 4) + 4);

    double totalChildMass = getTotalMass((nodeIndex * 4) + 1) + getTotalMass((nodeIndex * 4) + 2) + getTotalMass((nodeIndex * 4) + 3) + getTotalMass((nodeIndex * 4) + 4);
    double totalCenterMassX = 0.0, totalCenterMassY = 0.0;

    Node &topLNode = h_node[(nodeIndex * 4) + 2], &topRNode = h_node[(nodeIndex * 4) + 1], &botLNode = h_node[(nodeIndex * 4) + 3], &botRNode = h_node[(nodeIndex * 4) + 4];

    totalCenterMassX += topLNode.centerMass.x * topLNode.totalMass;
    totalCenterMassY += topLNode.centerMass.y * topLNode.totalMass;

    totalCenterMassX += topRNode.centerMass.x * topRNode.totalMass;
    totalCenterMassY += topRNode.centerMass.y * topRNode.totalMass;

    totalCenterMassX += botLNode.centerMass.x * botLNode.totalMass;
    totalCenterMassY += botLNode.centerMass.y * botLNode.totalMass;

    totalCenterMassX += botRNode.centerMass.x * botRNode.totalMass;
    totalCenterMassY += botRNode.centerMass.y * botRNode.totalMass;

    curNode.totalMass = totalChildMass;
    curNode.centerMass = {totalCenterMassX / totalChildMass, totalCenterMassY / totalChildMass};
}

void BarnesHutCuda::randomInitBodies()
{
    srand(time(NULL));
    for (int i = 0; i < nBodies - 1; ++i)
    {
        int randPx = rand() % (1000 - 1 + 1) + 500;
        int randPy = rand() % (1000 - 1 + 1) + 500;
        // int randVx = rand() % (500 - 1 + 1) + 1;
        // int randVy = rand() % (500 - 1 + 1) + 1;
        h_b[i].isDynamic = true;
        h_b[i].mass = 1.0 / (double)nBodies;
        h_b[i].radius = 10;
        h_b[i].position = {(double)randPx, (double)randPy};
        h_b[i].velocity = {0.0, 0.0};
        h_b[i].acceleration = {0.0, 0.0};
    }

    h_b[nBodies - 1].isDynamic = false;
    h_b[nBodies - 1].mass = 200.0 / (double)nBodies;
    h_b[nBodies - 1].radius = 10;
    h_b[nBodies - 1].position = {CENTERX, CENTERY};
    h_b[nBodies - 1].velocity = {0.0, 0.0};
    h_b[nBodies - 1].acceleration = {0.0, 0.0};
}

Body *BarnesHutCuda::getBodies()
{

    return h_b;
}

double BarnesHutCuda::getDistance(Vector pos1, Vector pos2)
{

    return sqrt(pow(pos1.x - pos2.x, 2) + pow(pos1.y - pos2.y, 2));
}

double BarnesHutCuda::getWidth(Node &root)
{

    return root.botRight.x - root.topLeft.x;
}

bool BarnesHutCuda::isCollide(Body &b1, Body &b2)
{
    return b1.radius + b2.radius > getDistance(b1.position, b2.position);
}

void BarnesHutCuda::computeForceHelper(int nodeIndex, int b)
{

    if (nodeIndex >= nNodes)
    {
        return;
    }

    Node &curNode = h_node[nodeIndex];
    Body &bi = h_b[b];

    if (curNode.isLeaf)
    {

        if (curNode.centerMass.x != -1)
        {

            if (bi.radius + bi.radius > getDistance(bi.position, curNode.centerMass))
                return;

            Vector rij = {curNode.centerMass.x - bi.position.x, curNode.centerMass.y - bi.position.y};
            double inv_r3 = pow(rij.x * rij.x + rij.y * rij.y + E * E, -1.5);
            double f = (GRAVITY * curNode.totalMass) / inv_r3;
            Vector force = {rij.x * f, rij.y * f};
            bi.acceleration.x += (force.x / bi.mass);
            bi.acceleration.y += (force.y / bi.mass);
        }

        return;
    }

    double sd = getWidth(curNode) / getDistance(bi.position, curNode.centerMass);
    if (sd < THETA)
    {
        Vector rij = {curNode.centerMass.x - bi.position.x, curNode.centerMass.y - bi.position.y};
        if (bi.radius * 2 > getDistance(bi.position, curNode.centerMass))
            return;
        double inv_r3 = pow(rij.x * rij.x + rij.y * rij.y + E * E, -1.5);
        double f = (GRAVITY * curNode.totalMass) / inv_r3;
        Vector force = {rij.x * f, rij.y * f};
        bi.acceleration.x += (force.x / bi.mass);
        bi.acceleration.y += (force.y / bi.mass);
        return;
    }

    computeForceHelper((nodeIndex * 4) + 1, b);
    computeForceHelper((nodeIndex * 4) + 2, b);
    computeForceHelper((nodeIndex * 4) + 3, b);
    computeForceHelper((nodeIndex * 4) + 4, b);
}

void BarnesHutCuda::computeForce()
{

    for (int i = 0; i < nBodies; ++i)
    {

        Body &bi = h_b[i];
        if (bi.isDynamic)
        {

            bi.velocity.x += bi.acceleration.x * DT / 2.0;
            bi.velocity.y += bi.acceleration.y * DT / 2.0;

            bi.position.x += bi.velocity.x * DT;
            bi.position.y += bi.velocity.y * DT;

            bi.acceleration = {0.0, 0.0};
            computeForceHelper(0, i);
            bi.velocity.x += bi.acceleration.x * DT / 2.0;
            bi.velocity.y += bi.acceleration.y * DT / 2.0;
        }
    }
}
void BarnesHutCuda::reset()
{

    for (int b = 0; b < nNodes; ++b)
    {

        h_node[b].topLeft = {-1, -1};
        h_node[b].botRight = {-1, -1};
        h_node[b].centerMass = {-1, -1};
        h_node[b].totalMass = 0.0;
        h_node[b].isLeaf = true;
    }

    h_node[0].topLeft = {INFINITY, -INFINITY};
    h_node[0].botRight = {-INFINITY, INFINITY};
}

void BarnesHutCuda::setup()
{
    randomInitBodies();
    for (int i = 0; i < nBodies; ++i)
    {
        std::cout << "Set up Body " << h_b[i].position.x << " " << h_b[i].position.y << "\n"
                  << std::endl;
    }
    // cudaMemcpy(d_b, h_b, sizeof(Body) * nBodies, cudaMemcpyHostToDevice);
    // cudaMemcpy(d_node, h_node, sizeof(Node) * nNodes, cudaMemcpyHostToDevice);
}
void BarnesHutCuda::update()
{
    reset();
    computeBoundingBox();
    constructQuadTree();
    computeCenterMass(0);
    computeForce();
}
