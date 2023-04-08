#include <iostream>
#include "barnesHut_kernel.cuh"
#include "constants.h"

dim3 gridSize = GRID_SIZE;
dim3 blockSize = BLOCK_SIZE;

BarnesHutCuda::BarnesHutCuda(int n) : nBodies(n)
{
    nNodes = MAX_NODES;
    leafLimit = MAX_NODES - N_LEAF;
    h_b = new Body[n];
    h_node = new Node[nNodes];
    h_topLeft = new Vector;
    h_botRight = new Vector;

    cudaMalloc((void **)&d_b, sizeof(Body) * n);
    cudaMalloc((void **)&d_node, sizeof(Node) * nNodes);
    cudaMalloc((void **)&d_mutex, sizeof(int) * nNodes);
    cudaMalloc((void **)&d_topLeft, sizeof(Vector));
    cudaMalloc((void **)&d_botRight, sizeof(Vector));
}

BarnesHutCuda::~BarnesHutCuda()
{

    delete[] h_b;
    delete[] h_node;
    delete h_topLeft;
    delete h_botRight;
    cudaFree(d_b);
    cudaFree(d_node);
    cudaFree(d_mutex);
    cudaFree(d_topLeft);
    cudaFree(d_botRight);
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

    *h_topLeft = {topLeft.x - 1, topLeft.y + 1};
    *h_botRight = {botRight.x + 1, botRight.y - 1};
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
    // std::cout << "inBoundary TR After: " << topLeft.x << " " << topLeft.y << " ";
    // std::cout << "inBoundary BT After: " << botRight.x << " " << botRight.y << " ";
    // std::cout << "inBoundary Body: " << p.x << " " << p.y << " " << std::endl;

    return (p.x >= topLeft.x && p.x <= botRight.x && p.y <= topLeft.y && p.y >= botRight.y);
}

void BarnesHutCuda::updateBound(Vector &tl, Vector &br, int quadrant)
{
    if (quadrant == 1)
    {
        tl = {(tl.x + br.x) / 2, tl.y};
        br = {br.x, (tl.y + br.y) / 2};
    }
    else if (quadrant == 2)
    {
        tl = {tl.x, tl.y};
        br = {(tl.x + br.x) / 2, (tl.y + br.y) / 2};
    }
    else if (quadrant == 3)
    {
        tl = {tl.x, (tl.y + br.y) / 2};
        br = {(tl.x + br.x) / 2, br.y};
    }
    else
    {
        tl = {(tl.x + br.x) / 2, (tl.y + br.y) / 2};
        br = {br.x, br.y};
    }
}

void BarnesHutCuda::insertQuadTree(int nodeIndex, int b, Vector tl, Vector br)
{

    Body &body = h_b[b];
    while (nodeIndex < nNodes)
    {

        Node &curNode = h_node[nodeIndex];

        if (!inBoundary(tl, br, body.position))
        {
            std::cout << "wrong quat" << std::endl;
            break;
        }

        if (nodeIndex >= leafLimit)
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

            break;
        }

        // If node x does not contain a body, put the new body here.
        if (curNode.isLeaf)
        {

            if (curNode.centerMass.x != -1)
            {

                int quadrant = getQuadrant(tl, br, curNode.centerMass.x, curNode.centerMass.y);
                Node &childNode = h_node[(nodeIndex * 4) + quadrant];

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
                break;
            }
        }

        int quadrant = getQuadrant(tl, br, body.position.x, body.position.y);
        // std::cout << "TR Before: " << tl.x << " " << tl.y << " ";
        // std::cout << "BT Before: " << br.x << " " << br.y << " ";
        updateBound(tl, br, quadrant);
        nodeIndex = (nodeIndex * 4) + quadrant;
        // std::cout << "TR After: " << tl.x << " " << tl.y << " ";
        // std::cout << "BT After: " << br.x << " " << br.y << " ";
        // std::cout << "Body: " << body.position.x << " " << body.position.y << " ";
        // std::cout << "QD: " << quadrant << std::endl;
    }
}

void BarnesHutCuda::constructQuadTree()
{

    for (int i = 0; i < nBodies; ++i)
    {
        insertQuadTree(0, i, *h_topLeft, *h_botRight);
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

void BarnesHutCuda::resetCUDA()
{
    int blockSize = BLOCK_SIZE;
    dim3 gridSize = ceil((float)nNodes / blockSize);
    ResetKernel<<<gridSize, blockSize>>>(d_node, d_topLeft, d_botRight, d_mutex, nNodes);
}
void BarnesHutCuda::computeBoundingBoxCUDA()
{
    int blockSize = BLOCK_SIZE;
    dim3 gridSize = ceil((float)nBodies / blockSize);
    ComputeBoundingBox<<<gridSize, blockSize>>>(d_node, d_b, d_topLeft, d_botRight, d_mutex, nBodies);
}
void BarnesHutCuda::constructQuadTreeCUDA()
{
    int blockSize = BLOCK_SIZE;
    dim3 gridSize = ceil((float)nBodies / blockSize);
    ConstructQuadTreeKernel<<<gridSize, blockSize>>>(d_node, d_b, d_topLeft, d_botRight, d_mutex, nNodes, nBodies, leafLimit);
}
void BarnesHutCuda::computeCenterMassCUDA()
{
    int start = leafLimit, end = nNodes;
    int totalNodes = end - start, temp = 0;
    while (end > start)
    {
        // std::cout << start << "->" << end << std::endl;
        int blockSize = BLOCK_SIZE;
        dim3 gridSize = ceil((float)totalNodes / blockSize);
        ComputeCenterMass<<<gridSize, blockSize>>>(d_node, nNodes, start, end);
        totalNodes /= 4;
        temp = start;
        end = start;
        start = temp - totalNodes;
    }
}
void BarnesHutCuda::computeForceCUDA()
{
    int blockSize = BLOCK_SIZE;
    dim3 gridSize = ceil((float)nBodies / blockSize);
    ComputeForceKernel<<<gridSize, blockSize>>>(d_node, d_b, nNodes, nBodies, leafLimit, h_botRight->x - h_topLeft->x);
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

bool BarnesHutCuda::isCollide(Body &b1, Body &b2)
{
    return b1.radius + b2.radius > getDistance(b1.position, b2.position);
}

void BarnesHutCuda::computeForceHelper(int nodeIndex, int b, double width)
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

    double sd = width / getDistance(bi.position, curNode.centerMass);
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

    computeForceHelper((nodeIndex * 4) + 1, b, width / 2.0);
    computeForceHelper((nodeIndex * 4) + 2, b, width / 2.0);
    computeForceHelper((nodeIndex * 4) + 3, b, width / 2.0);
    computeForceHelper((nodeIndex * 4) + 4, b, width / 2.0);
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
            computeForceHelper(0, i, h_botRight->x - h_topLeft->x);
            bi.velocity.x += bi.acceleration.x * DT / 2.0;
            bi.velocity.y += bi.acceleration.y * DT / 2.0;
        }
    }
}
void BarnesHutCuda::reset()
{

    for (int b = 0; b < nNodes; ++b)
    {

        h_node[b].centerMass = {-1, -1};
        h_node[b].totalMass = 0.0;
        h_node[b].isLeaf = true;
    }

    *h_topLeft = {INFINITY, -INFINITY};
    *h_botRight = {-INFINITY, INFINITY};
}

void BarnesHutCuda::setup()
{
    randomInitBodies();
    cudaMemcpy(d_b, h_b, sizeof(Body) * nBodies, cudaMemcpyHostToDevice);
    cudaMemcpy(d_node, h_node, sizeof(Node) * nNodes, cudaMemcpyHostToDevice);
}
void BarnesHutCuda::update()
{
    // reset();
    resetCUDA();
    computeBoundingBoxCUDA();
    cudaMemcpy(h_node, d_node, sizeof(Node) * nNodes, cudaMemcpyDeviceToHost);
    cudaMemcpy(h_topLeft, d_topLeft, sizeof(Vector), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_botRight, d_botRight, sizeof(Vector), cudaMemcpyDeviceToHost);
    constructQuadTree();
    cudaMemcpy(d_node, h_node, sizeof(Node) * nNodes, cudaMemcpyHostToDevice);
    // constructQuadTreeCUDA();
    computeCenterMassCUDA();
    // computeBoundingBox();
    // std::cout << h_topLeft->x << " " << h_topLeft->y << " " << h_botRight->x << " " << h_botRight->y << std::endl;
    // cudaMemcpy(h_node, d_node, sizeof(Node) * nNodes, cudaMemcpyDeviceToHost);
    // std::cout << h_topLeft->x << " " << h_topLeft->y << " " << h_botRight->x << " " << h_botRight->y << std::endl;
    // constructQuadTree();
    // computeCenterMass(0);
    // computeForce();
    computeForceCUDA();
    // cudaMemcpy(d_b, h_b, sizeof(Body) * nBodies, cudaMemcpyHostToDevice);
    // cudaMemcpy(h_b, d_b, sizeof(Body) * nBodies, cudaMemcpyDeviceToHost);
}