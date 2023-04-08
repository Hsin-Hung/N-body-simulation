#ifndef BARNES_HUT_KERNEL_
#define BARNES_HUT_KERNEL_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <vector>
#include "constants.h"
#include "barnesHutCuda.cuh"
#include "barnesHut_kernel.cuh"

/*
----------------------------------------------------------------------------------------
RESET KERNEL
----------------------------------------------------------------------------------------
*/
__global__ void ResetKernel(Node *node, Vector *topLeft, Vector *botRight, int *mutex, int nNodes)
{
    int b = blockIdx.x * blockDim.x + threadIdx.x;

    if (b < nNodes)
    {

        node[b].centerMass = {-1, -1};
        node[b].totalMass = 0.0;
        node[b].isLeaf = true;
        mutex[b] = 0;
    }

    if (b == 0)
    {
        *topLeft = {INFINITY, -INFINITY};
        *botRight = {-INFINITY, INFINITY};
    }
}

/*
----------------------------------------------------------------------------------------
COMPUTE BOUNDING BOX
----------------------------------------------------------------------------------------
*/
__global__ void ComputeBoundingBox(Node *node, Body *bodies, Vector *topLeft, Vector *botRight, int *mutex, int nBodies)
{

    __shared__ double topLeftX[BLOCK_SIZE];
    __shared__ double topLeftY[BLOCK_SIZE];
    __shared__ double botRightX[BLOCK_SIZE];
    __shared__ double botRightY[BLOCK_SIZE];

    int tx = threadIdx.x;
    int b = blockIdx.x * blockDim.x + tx;

    topLeftX[tx] = INFINITY;
    topLeftY[tx] = -INFINITY;
    botRightX[tx] = -INFINITY;
    botRightY[tx] = INFINITY;

    __syncthreads();

    if (b < nBodies)
    {
        Body body = bodies[b];
        topLeftX[tx] = body.position.x;
        topLeftY[tx] = body.position.y;
        botRightX[tx] = body.position.x;
        botRightY[tx] = body.position.y;
    }

    for (int s = blockDim.x / 2; s > 0; s >>= 1)
    {
        __syncthreads();
        if (tx < s)
        {
            topLeftX[tx] = fminf(topLeftX[tx], topLeftX[tx + s]);
            topLeftY[tx] = fmaxf(topLeftY[tx], topLeftY[tx + s]);
            botRightX[tx] = fmaxf(botRightX[tx], botRightX[tx + s]);
            botRightY[tx] = fminf(botRightY[tx], botRightY[tx + s]);
        }
    }

    if (tx == 0 && b < nBodies)
    {
        while (atomicCAS(mutex, 0, 1) != 0)
            ;
        topLeft->x = fminf(topLeft->x, topLeftX[0] - 1);
        topLeft->y = fmaxf(topLeft->y, topLeftY[0] + 1);
        botRight->x = fmaxf(botRight->x, botRightX[0] + 1);
        botRight->y = fminf(botRight->y, botRightY[0] - 1);
        atomicExch(mutex, 0);
    }
}

/*
----------------------------------------------------------------------------------------
CONSTRUCT QUAD TREE
----------------------------------------------------------------------------------------
*/

__device__ int getQuadrant(Vector topLeft, Vector botRight, double x, double y)
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

__device__ bool inBoundary(Vector tl, Vector br, Vector p)
{
    return (p.x >= tl.x && p.x <= br.x && p.y <= tl.y && p.y >= br.y);
}

__device__ void updateBound(Vector &tl, Vector &br, int quadrant)
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

__device__ void ConstructQuadTreeHelper(Node *node, int nodeIndex, Body body, Vector *topLeft, Vector *botRight, int *mutex, int nNodes, int nBodies, int leafLimit)
{

    Vector tl = *topLeft;
    Vector br = *botRight;

    while (nodeIndex < nNodes)
    {

        Node &curNode = node[nodeIndex];

        if (!inBoundary(tl, br, body.position))
        {
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
                curNode.totalMass = body.mass;
                curNode.centerMass = body.position;
            }
            break;
        }

        // If node x does not contain a body, put the new body here.
        if (curNode.isLeaf)
        {
            if (curNode.centerMass.x != -1)
            {

                int quadrant = getQuadrant(tl, br, curNode.centerMass.x, curNode.centerMass.y);
                Node &childNode = node[(nodeIndex * 4) + quadrant];
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
        updateBound(tl, br, quadrant);
        nodeIndex = (nodeIndex * 4) + quadrant;
    }
}

__global__ void ConstructQuadTreeKernel(Node *node, Body *bodies, Vector *topLeft, Vector *botRight, int *mutex, int nNodes, int nBodies, int leafLimit)
{

    int b = blockIdx.x * blockDim.x + threadIdx.x;

    if (b < nBodies)
    {
        Body body = bodies[b];
        ConstructQuadTreeHelper(node, 0, body, topLeft, botRight, mutex, nNodes, nBodies, leafLimit);
    }
}

/*
----------------------------------------------------------------------------------------
COMPUTE CENTER MASS
----------------------------------------------------------------------------------------
*/

__global__ void ComputeCenterMass(Node *node, int nNodes, int start, int end)
{
    int nodeIndex = blockIdx.x * blockDim.x + threadIdx.x;
    nodeIndex += start;
    if (nodeIndex >= start && nodeIndex < end)
    {

        Node &curNode = node[nodeIndex];

        if (!curNode.isLeaf)
        {

            Node &topLNode = node[(nodeIndex * 4) + 2], &topRNode = node[(nodeIndex * 4) + 1], &botLNode = node[(nodeIndex * 4) + 3], &botRNode = node[(nodeIndex * 4) + 4];
            double totalChildMass = topLNode.totalMass + topRNode.totalMass + botLNode.totalMass + botRNode.totalMass;
            double totalCenterMassX = 0.0, totalCenterMassY = 0.0;

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
    }
}

// /*
// ----------------------------------------------------------------------------------------
// COMPUTE FORCE
// ----------------------------------------------------------------------------------------
// */

__device__ double getDistance(Vector pos1, Vector pos2)
{

    return sqrt(pow(pos1.x - pos2.x, 2) + pow(pos1.y - pos2.y, 2));
}

__device__ bool isCollide(Body &b1, Body &b2)
{
    return b1.radius + b2.radius > getDistance(b1.position, b2.position);
}

__device__ void ComputeForce(Node *node, Body *bodies, int bodyIndex, int nNodes, int nBodies, int leafLimit, double width)
{

    Body &bi = bodies[bodyIndex];
    int q_size = nNodes - leafLimit;
    int queue[MAX_NODES - N_LEAF];
    int front = 0, insert = 0;
    int size;
    queue[insert++] = 0;

    while (front != insert)
    {
        size = insert - front;

        for (int i = 0; i < size; ++i)
        {

            int nodeIndex = queue[front++];
            front %= q_size;
            Node &curNode = node[nodeIndex];

            if (curNode.isLeaf)
            {
                if (curNode.centerMass.x != -1)
                {
                    if (bi.radius * 2 > getDistance(bi.position, curNode.centerMass))
                        continue;

                    Vector rij = {curNode.centerMass.x - bi.position.x, curNode.centerMass.y - bi.position.y};
                    double inv_r3 = pow(rij.x * rij.x + rij.y * rij.y + E * E, -1.5);
                    double f = (GRAVITY * curNode.totalMass) / inv_r3;
                    Vector force = {rij.x * f, rij.y * f};
                    bi.acceleration.x += (force.x / bi.mass);
                    bi.acceleration.y += (force.y / bi.mass);
                }
                continue;
            }

            double sd = width / getDistance(bi.position, curNode.centerMass);
            if (sd < THETA)
            {
                Vector rij = {curNode.centerMass.x - bi.position.x, curNode.centerMass.y - bi.position.y};
                if (bi.radius * 2 > getDistance(bi.position, curNode.centerMass))
                    continue;
                double inv_r3 = pow(rij.x * rij.x + rij.y * rij.y + E * E, -1.5);
                double f = (GRAVITY * curNode.totalMass) / inv_r3;
                Vector force = {rij.x * f, rij.y * f};
                bi.acceleration.x += (force.x / bi.mass);
                bi.acceleration.y += (force.y / bi.mass);
                continue;
            }
            queue[insert++] = (nodeIndex * 4) + 1;
            queue[insert++] = (nodeIndex * 4) + 2;
            queue[insert++] = (nodeIndex * 4) + 3;
            queue[insert++] = (nodeIndex * 4) + 4;
            insert %= q_size;
        }

        width /= 2.0;
    }
}

__global__ void ComputeForceKernel(Node *node, Body *bodies, int nNodes, int nBodies, int leafLimit, double width)
{

    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i < nBodies)
    {
        Body &bi = bodies[i];
        if (bi.isDynamic)
        {

            bi.velocity.x += bi.acceleration.x * DT / 2.0;
            bi.velocity.y += bi.acceleration.y * DT / 2.0;

            bi.position.x += bi.velocity.x * DT;
            bi.position.y += bi.velocity.y * DT;

            bi.acceleration = {0.0, 0.0};
            ComputeForce(node, bodies, i, nNodes, nBodies, leafLimit, width);
            bi.velocity.x += bi.acceleration.x * DT / 2.0;
            bi.velocity.y += bi.acceleration.y * DT / 2.0;
                }
    }
}

#endif