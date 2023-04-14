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
__global__ void ResetKernel(Node *node, Vector *topLeft, Vector *botRight, int *mutex, int nNodes, int nBodies)
{
    int b = blockIdx.x * blockDim.x + threadIdx.x;

    if (b < nNodes)
    {
        node[b].topLeft = {INFINITY, -INFINITY};
        node[b].botRight = {-INFINITY, INFINITY};
        node[b].centerMass = {-1, -1};
        node[b].totalMass = 0.0;
        node[b].isLeaf = true;
        node[b].start = -1;
        node[b].end = -1;
        mutex[b] = 0;
    }

    if (b == 0)
    {
        node[b].start = 0;
        node[b].end = nBodies - 1;
        *topLeft = {INFINITY, -INFINITY};
        *botRight = {-INFINITY, INFINITY};
    }
}

/*
----------------------------------------------------------------------------------------
COMPUTE BOUNDING BOX
----------------------------------------------------------------------------------------
*/
__global__ void ComputeBoundingBoxKernel(Node *node, Body *bodies, Vector *topLeft, Vector *botRight, int *mutex, int nBodies)
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

    if (tx == 0)
    {
        while (atomicCAS(mutex, 0, 1) != 0)
            ;
        // topLeft->x = fminf(topLeft->x, topLeftX[0] - 1);
        // topLeft->y = fmaxf(topLeft->y, topLeftY[0] + 1);
        // botRight->x = fmaxf(botRight->x, botRightX[0] + 1);
        // botRight->y = fminf(botRight->y, botRightY[0] - 1);
        node[0].topLeft.x = fminf(node[0].topLeft.x, topLeftX[0] - 1);
        node[0].topLeft.y = fmaxf(node[0].topLeft.y, topLeftY[0] + 1);
        node[0].botRight.x = fmaxf(node[0].botRight.x, botRightX[0] + 1);
        node[0].botRight.y = fminf(node[0].botRight.y, botRightY[0] - 1);
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

__device__ void CountBodies(Body *bodies, Vector topLeft, Vector botRight, int *count, int start, int end, int nBodies)
{
    if (threadIdx.x < 4)
        count[threadIdx.x] = 0;
    __syncthreads();
    if (threadIdx.x == 0)
    {

        for (int i = start; i <= end; ++i)
        {
            Body &body = bodies[i];
            int quadrant = getQuadrant(topLeft, botRight, body.position.x, body.position.y);
            ++count[quadrant - 1];
        }
    }

    __syncthreads();
}

__device__ void GroupBodies(Body *bodies, Body *buffer, Vector topLeft, Vector botRight, int *count, int start, int end, int nBodies)
{
    int q1 = start, q2 = start + count[0], q3 = start + count[0] + count[1], q4 = start + count[0] + count[1] + count[2];
    if (threadIdx.x == 0)
    {

        for (int i = start; i <= end; ++i)
        {
            Body &body = bodies[i];
            int quadrant = getQuadrant(topLeft, botRight, body.position.x, body.position.y);
            if (quadrant == 1)
            {

                buffer[q1++] = body;
            }
            else if (quadrant == 2)
            {

                buffer[q2++] = body;
            }
            else if (quadrant == 3)
            {

                buffer[q3++] = body;
            }
            else
            {

                buffer[q4++] = body;
            }
        }
    }

    __syncthreads();
}

__global__ void ConstructQuadTreeDPKernel(Node *node, Body *bodies, Body *buffer, int nodeIndex, int nNodes, int nBodies, int leafLimit)
{
    __shared__ int count[4];
    int tx = threadIdx.x;
    nodeIndex += blockIdx.x;

    if (nodeIndex >= nNodes)
        return;

    Node &curNode = node[nodeIndex];
    int start = curNode.start, end = curNode.end;
    // if (tx == 0)
    // {
    //     printf("node: %d, start-end: %d -> %d , topleft x: %f y: %f -> botright x: %f y: %f\n", nodeIndex, start, end,
    //            curNode.topLeft.x, curNode.topLeft.y, curNode.botRight.x, curNode.botRight.y);
    // }

    Vector topLeft = curNode.topLeft, botRight = curNode.botRight;

    if (start == -1 && end == -1)
        return;

    double M = curNode.totalMass;
    double Rx = curNode.totalMass * curNode.centerMass.x;
    double Ry = curNode.totalMass * curNode.centerMass.y;
    for (int i = start; i <= end; ++i)
    {
        Body &body = bodies[i];
        M += body.mass;
        Rx += body.mass * body.position.x;
        Ry += body.mass * body.position.y;
    }
    Rx /= M;
    Ry /= M;
    curNode.totalMass = M;
    curNode.centerMass = {Rx, Ry};

    // if (tx == 0)
    //     printf("node: %d, total mass: %f, center mass: %f, %f\n", nodeIndex, M, Rx, Ry);

    if (nodeIndex >= leafLimit || start == end)
    {
        buffer[start] = bodies[start];
        return;
    }

    // if (threadIdx.x == 0)
    // {
    //     for (int i = start; i <= end; ++i)
    //     {
    //         printf("before node: %d, body id: %d \n", nodeIndex, bodies[i].id);
    //     }
    // }
    CountBodies(bodies, topLeft, botRight, count, start, end, nBodies);
    GroupBodies(bodies, buffer, topLeft, botRight, count, start, end, nBodies);
    // if (threadIdx.x == 0)
    // {
    //     for (int i = start; i <= end; ++i)
    //     {
    //         printf("after node: %d, body id: %d \n", nodeIndex, buffer[i].id);
    //     }
    // }
    Node &topLNode = node[(nodeIndex * 4) + 1],
         &topRNode = node[(nodeIndex * 4) + 2], &botLNode = node[(nodeIndex * 4) + 3], &botRNode = node[(nodeIndex * 4) + 4];

    if (tx == 0)
    {
        // printf("node: %d, count: %d %d %d %d\n", nodeIndex, count[0], count[1], count[2], count[3]);
        topLNode.topLeft = topLeft;
        topLNode.botRight = botRight;
        topRNode.topLeft = topLeft;
        topRNode.botRight = botRight;
        botLNode.topLeft = topLeft;
        botLNode.botRight = botRight;
        botRNode.topLeft = topLeft;
        botRNode.botRight = botRight;

        updateBound(topLNode.topLeft, topLNode.botRight, 1);
        updateBound(topRNode.topLeft, topRNode.botRight, 2);
        updateBound(botLNode.topLeft, botLNode.botRight, 3);
        updateBound(botRNode.topLeft, botRNode.botRight, 4);

        curNode.isLeaf = false;

        if (count[0] > 0)
        {
            topLNode.start = start;
            topLNode.end = start + count[0] - 1;
        }

        if (count[1] > 0)
        {
            topRNode.start = start + count[0];
            topRNode.end = start + count[0] + count[1] - 1;
        }

        if (count[2] > 0)
        {
            botLNode.start = start + count[0] + count[1];
            botLNode.end = start + count[0] + count[1] + count[2] - 1;
        }

        if (count[3] > 0)
        {
            botRNode.start = start + count[0] + count[1] + count[2];
            botRNode.end = end;
        }
        ConstructQuadTreeDPKernel<<<4, BLOCK_SIZE>>>(node, buffer, bodies, nodeIndex * 4 + 1, nNodes, nBodies, leafLimit);
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
__device__ void ComputeForceRecursive(Node *node, Body *bodies, int nodeIndex, int bodyIndex, int nNodes, int nBodies, int leafLimit, double width)
{

    if (nodeIndex >= nNodes)
    {
        return;
    }
    Node &curNode = node[nodeIndex];
    Body &bi = bodies[bodyIndex];
    if (curNode.isLeaf)
    {

        if (curNode.centerMass.x != -1)
        {
            if (bi.radius * 2 > getDistance(bi.position, curNode.centerMass))
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

    ComputeForceRecursive(node, bodies, (nodeIndex * 4) + 1, bodyIndex, nNodes, nBodies, leafLimit, width / 2);
    ComputeForceRecursive(node, bodies, (nodeIndex * 4) + 2, bodyIndex, nNodes, nBodies, leafLimit, width / 2);
    ComputeForceRecursive(node, bodies, (nodeIndex * 4) + 3, bodyIndex, nNodes, nBodies, leafLimit, width / 2);
    ComputeForceRecursive(node, bodies, (nodeIndex * 4) + 4, bodyIndex, nNodes, nBodies, leafLimit, width / 2);
}
__device__ void ComputeForce(Node *node, Body *bodies, int bodyIndex, int nNodes, int nBodies, int leafLimit, double width)
{
    Body &bi = bodies[bodyIndex];
    int q_size = 1024;
    __shared__ int queue[1024];
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

            if (node[(nodeIndex * 4) + 1].totalMass > 0)
                queue[insert++] = (nodeIndex * 4) + 1;
            if (node[(nodeIndex * 4) + 2].totalMass > 0)
                queue[insert++] = (nodeIndex * 4) + 2;
            if (node[(nodeIndex * 4) + 3].totalMass > 0)
                queue[insert++] = (nodeIndex * 4) + 3;
            if (node[(nodeIndex * 4) + 4].totalMass > 0)
                queue[insert++] = (nodeIndex * 4) + 4;
            insert %= q_size;
        }

        width /= 2.0;
    }
}

__global__ void ComputeForceKernel(Node *node, Body *bodies, int nNodes, int nBodies, int leafLimit)
{

    int i = blockIdx.x;
    double width = node[0].botRight.x - node[0].topLeft.x;

    if (i < nBodies)
    {
        Body &bi = bodies[i];
        if (bi.isDynamic)
        {
            if (threadIdx.x == 0)
            {
                bi.velocity.x += bi.acceleration.x * DT / 2.0;
                bi.velocity.y += bi.acceleration.y * DT / 2.0;
                // printf("velocity: x: %f  y: %f\n", bi.velocity.x, bi.velocity.y);
                bi.position.x += bi.velocity.x * DT;
                bi.position.y += bi.velocity.y * DT;

                bi.acceleration = {0.0, 0.0};

                ComputeForceRecursive(node, bodies, 0, i, nNodes, nBodies, leafLimit, width);
                bi.velocity.x += bi.acceleration.x * DT / 2.0;
                bi.velocity.y += bi.acceleration.y * DT / 2.0;
            }
        }
    }
}

#endif