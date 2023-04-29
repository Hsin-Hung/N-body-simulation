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
__global__ void ResetKernel(Node *node, int *mutex, int nNodes, int nBodies)
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
    }
}

/*
----------------------------------------------------------------------------------------
COMPUTE BOUNDING BOX
----------------------------------------------------------------------------------------
*/
__global__ void ComputeBoundingBoxKernel(Node *node, Body *bodies, int *mutex, int nBodies)
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
        node[0].topLeft.x = fminf(node[0].topLeft.x, topLeftX[0] - 1.0e10);
        node[0].topLeft.y = fmaxf(node[0].topLeft.y, topLeftY[0] + 1.0e10);
        node[0].botRight.x = fmaxf(node[0].botRight.x, botRightX[0] + 1.0e10);
        node[0].botRight.y = fminf(node[0].botRight.y, botRightY[0] - 1.0e10);
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

__device__ void UpdateChildBound(Vector &tl, Vector &br, Node &childNode, int quadrant)
{

    if (quadrant == 1)
    {
        childNode.topLeft = {(tl.x + br.x) / 2, tl.y};
        childNode.botRight = {br.x, (tl.y + br.y) / 2};
    }
    else if (quadrant == 2)
    {
        childNode.topLeft = {tl.x, tl.y};
        childNode.botRight = {(tl.x + br.x) / 2, (tl.y + br.y) / 2};
    }
    else if (quadrant == 3)
    {
        childNode.topLeft = {tl.x, (tl.y + br.y) / 2};
        childNode.botRight = {(tl.x + br.x) / 2, br.y};
    }
    else
    {
        childNode.topLeft = {(tl.x + br.x) / 2, (tl.y + br.y) / 2};
        childNode.botRight = {br.x, br.y};
    }
}

__device__ void warpReduce(volatile double *totalMass, volatile double2 *centerMass, int tx)
{
    totalMass[tx] += totalMass[tx + 32];
    centerMass[tx].x += centerMass[tx + 32].x;
    centerMass[tx].y += centerMass[tx + 32].y;
    totalMass[tx] += totalMass[tx + 16];
    centerMass[tx].x += centerMass[tx + 16].x;
    centerMass[tx].y += centerMass[tx + 16].y;
    totalMass[tx] += totalMass[tx + 8];
    centerMass[tx].x += centerMass[tx + 8].x;
    centerMass[tx].y += centerMass[tx + 8].y;
    totalMass[tx] += totalMass[tx + 4];
    centerMass[tx].x += centerMass[tx + 4].x;
    centerMass[tx].y += centerMass[tx + 4].y;
    totalMass[tx] += totalMass[tx + 2];
    centerMass[tx].x += centerMass[tx + 2].x;
    centerMass[tx].y += centerMass[tx + 2].y;
    totalMass[tx] += totalMass[tx + 1];
    centerMass[tx].x += centerMass[tx + 1].x;
    centerMass[tx].y += centerMass[tx + 1].y;
}

__device__ void ComputeCenterMass(Node &curNode, Body *bodies, double *totalMass, double2 *centerMass, int start, int end)
{
    int tx = threadIdx.x;
    int total = end - start + 1;
    int sz = ceil((double)total / blockDim.x);
    int s = tx * sz + start;
    double M = 0.0;
    double2 R = make_double2(0.0, 0.0);

    for (int i = s; i < s + sz; ++i)
    {
        if (i <= end)
        {
            Body &body = bodies[i];
            M += body.mass;
            R.x += body.mass * body.position.x;
            R.y += body.mass * body.position.y;
        }
    }

    totalMass[tx] = M;
    centerMass[tx] = R;

    for (unsigned int stride = blockDim.x / 2; stride > 32; stride >>= 1)
    {
        __syncthreads();
        if (tx < stride)
        {
            totalMass[tx] += totalMass[tx + stride];
            centerMass[tx].x += centerMass[tx + stride].x;
            centerMass[tx].y += centerMass[tx + stride].y;
        }
    }

    if (tx < 32)
    {
        warpReduce(totalMass, centerMass, tx);
    }
    __syncthreads();

    if (tx == 0)
    {
        centerMass[0].x /= totalMass[0];
        centerMass[0].y /= totalMass[0];
        curNode.totalMass = totalMass[0];
        curNode.centerMass = {centerMass[0].x, centerMass[0].y};
    }
}

__device__ void CountBodies(Body *bodies, Vector topLeft, Vector botRight, int *count, int start, int end, int nBodies)
{
    int tx = threadIdx.x;
    if (tx < 4)
        count[tx] = 0;
    __syncthreads();

    for (int i = start + tx; i <= end; i += blockDim.x)
    {
        Body body = bodies[i];
        int q = getQuadrant(topLeft, botRight, body.position.x, body.position.y);
        atomicAdd(&count[q - 1], 1);
    }

    __syncthreads();
}

__device__ void ComputeOffset(int *count, int start)
{
    int tx = threadIdx.x;
    if (tx < 4)
    {
        int offset = start;
        for (int i = 0; i < tx; ++i)
        {
            offset += count[i];
        }
        count[tx + 4] = offset;
    }
    __syncthreads();
}

__device__ void GroupBodies(Body *bodies, Body *buffer, Vector topLeft, Vector botRight, int *count, int start, int end, int nBodies)
{
    int *count2 = &count[4];
    for (int i = start + threadIdx.x; i <= end; i += blockDim.x)
    {
        Body body = bodies[i];
        int q = getQuadrant(topLeft, botRight, body.position.x, body.position.y);
        int dest = atomicAdd(&count2[q - 1], 1);
        buffer[dest] = body;
    }
    __syncthreads();
}

__global__ void ConstructQuadTreeKernel(Node *node, Body *bodies, Body *buffer, int nodeIndex, int nNodes, int nBodies, int leafLimit)
{
    __shared__ int count[8];
    __shared__ double totalMass[BLOCK_SIZE];
    __shared__ double2 centerMass[BLOCK_SIZE];
    int tx = threadIdx.x;
    nodeIndex += blockIdx.x;

    if (nodeIndex >= nNodes)
        return;

    Node &curNode = node[nodeIndex];
    int start = curNode.start, end = curNode.end;
    Vector topLeft = curNode.topLeft, botRight = curNode.botRight;

    if (start == -1 && end == -1)
        return;

    ComputeCenterMass(curNode, bodies, totalMass, centerMass, start, end);
    if (nodeIndex >= leafLimit || start == end)
    {
        for (int i = start; i <= end; ++i)
        {
            buffer[i] = bodies[i];
        }

        return;
    }

    CountBodies(bodies, topLeft, botRight, count, start, end, nBodies);
    ComputeOffset(count, start);
    GroupBodies(bodies, buffer, topLeft, botRight, count, start, end, nBodies);

    if (tx == 0)
    {
        Node &topLNode = node[(nodeIndex * 4) + 2],
             &topRNode = node[(nodeIndex * 4) + 1], &botLNode = node[(nodeIndex * 4) + 3], &botRNode = node[(nodeIndex * 4) + 4];

        UpdateChildBound(topLeft, botRight, topLNode, 2);
        UpdateChildBound(topLeft, botRight, topRNode, 1);
        UpdateChildBound(topLeft, botRight, botLNode, 3);
        UpdateChildBound(topLeft, botRight, botRNode, 4);

        curNode.isLeaf = false;

        if (count[0] > 0)
        {
            topRNode.start = start;
            topRNode.end = start + count[0] - 1;
        }

        if (count[1] > 0)
        {
            topLNode.start = start + count[0];
            topLNode.end = start + count[0] + count[1] - 1;
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
        ConstructQuadTreeKernel<<<4, BLOCK_SIZE>>>(node, buffer, bodies, nodeIndex * 4 + 1, nNodes, nBodies, leafLimit);
    }
}

/*
----------------------------------------------------------------------------------------
COMPUTE FORCE
----------------------------------------------------------------------------------------
*/
__device__ double getDistance(Vector pos1, Vector pos2)
{

    return sqrt(pow(pos1.x - pos2.x, 2) + pow(pos1.y - pos2.y, 2));
}

__device__ bool isCollide(Body &b1, Vector cm)
{
    return b1.radius * 2 + COLLISION_TH > getDistance(b1.position, cm);
}

__device__ void ComputeForce(Node *node, Body *bodies, int nodeIndex, int bodyIndex, int nNodes, int nBodies, int leafLimit, double width)
{

    if (nodeIndex >= nNodes)
    {
        return;
    }
    Node curNode = node[nodeIndex];
    Body bi = bodies[bodyIndex];
    if (curNode.isLeaf)
    {
        if (curNode.centerMass.x != -1 && !isCollide(bi, curNode.centerMass))
        {
            Vector rij = {curNode.centerMass.x - bi.position.x, curNode.centerMass.y - bi.position.y};
            double r = sqrt((rij.x * rij.x) + (rij.y * rij.y) + (E * E));
            double f = (GRAVITY * bi.mass * curNode.totalMass) / (r * r * r + (E * E));
            Vector force = {rij.x * f, rij.y * f};

            bodies[bodyIndex].acceleration.x += (force.x / bi.mass);
            bodies[bodyIndex].acceleration.y += (force.y / bi.mass);
        }
        return;
    }

    double sd = width / getDistance(bi.position, curNode.centerMass);
    if (sd < THETA)
    {
        if (!isCollide(bi, curNode.centerMass))
        {
            Vector rij = {curNode.centerMass.x - bi.position.x, curNode.centerMass.y - bi.position.y};
            double r = sqrt((rij.x * rij.x) + (rij.y * rij.y) + (E * E));
            double f = (GRAVITY * bi.mass * curNode.totalMass) / (r * r * r + (E * E));
            Vector force = {rij.x * f, rij.y * f};

            bodies[bodyIndex].acceleration.x += (force.x / bi.mass);
            bodies[bodyIndex].acceleration.y += (force.y / bi.mass);
        }

        return;
    }

    ComputeForce(node, bodies, (nodeIndex * 4) + 1, bodyIndex, nNodes, nBodies, leafLimit, width / 2);
    ComputeForce(node, bodies, (nodeIndex * 4) + 2, bodyIndex, nNodes, nBodies, leafLimit, width / 2);
    ComputeForce(node, bodies, (nodeIndex * 4) + 3, bodyIndex, nNodes, nBodies, leafLimit, width / 2);
    ComputeForce(node, bodies, (nodeIndex * 4) + 4, bodyIndex, nNodes, nBodies, leafLimit, width / 2);
}

__global__ void ComputeForceKernel(Node *node, Body *bodies, int nNodes, int nBodies, int leafLimit)
{

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    double width = node[0].botRight.x - node[0].topLeft.x;

    if (i < nBodies)
    {
        Body &bi = bodies[i];
        if (bi.isDynamic)
        {
            bi.acceleration = {0.0, 0.0};
            ComputeForce(node, bodies, 0, i, nNodes, nBodies, leafLimit, width);
            bi.velocity.x += bi.acceleration.x * DT;
            bi.velocity.y += bi.acceleration.y * DT;
            bi.position.x += bi.velocity.x * DT;
            bi.position.y += bi.velocity.y * DT;
        }
    }
}

#endif