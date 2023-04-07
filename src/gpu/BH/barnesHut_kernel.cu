#ifndef BARNES_HUT_KERNEL_H_
#define BARNES_HUT_KERNEL_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <vector>
#include "constants.h"
#include "barnesHutCuda.cuh"
#include "barnesHut_kernel.cuh"

// __global__ void ResetKernel(Node *node, int *mutex, int nNodes)
// {
//     int b = blockIdx.x * blockDim.x + threadIdx.x;

//     if (b < nNodes)
//     {

//         node[b].topLeft = {-1, -1};
//         node[b].botRight = {-1, -1};
//         node[b].centerMass = {-1, -1};
//         node[b].totalMass = 0.0;
//         node[b].bi = -1;
//         node[b].isNull = true;
//     }

//     if (b == 0)
//     {
//         node[0].topLeft = {INFINITY, -INFINITY};
//         node[0].botRight = {-INFINITY, INFINITY};
//         *mutex = 0;
//     }
// }

// __global__ void ComputeBoundingBox(Node *node, Body *bodies, int *mutex, int nBodies)
// {

//     __shared__ double topLeftX[BLOCK_SIZE];
//     __shared__ double topLeftY[BLOCK_SIZE];
//     __shared__ double botRightX[BLOCK_SIZE];
//     __shared__ double botRightY[BLOCK_SIZE];

//     int tx = threadIdx.x;
//     int b = blockIdx.x * blockDim.x + tx;

//     topLeftX[tx] = INFINITY;
//     topLeftY[tx] = -INFINITY;
//     botRightX[tx] = -INFINITY;
//     botRightY[tx] = INFINITY;

//     __syncthreads();

//     if (b < nBodies)
//     {

//         topLeftX[tx] = bodies[b].position.x;
//         topLeftY[tx] = bodies[b].position.y;
//         botRightX[tx] = bodies[b].position.x;
//         botRightY[tx] = bodies[b].position.y;
//     }

//     for (int s = blockDim.x / 2; s > 0; s >>= 1)
//     {
//         __syncthreads();
//         if (tx < s)
//         {
//             topLeftX[tx] = fminf(topLeftX[tx], topLeftX[tx + s]);
//             topLeftY[tx] = fmaxf(topLeftY[tx], topLeftY[tx + s]);
//             botRightX[tx] = fmaxf(botRightX[tx], botRightX[tx + s]);
//             botRightY[tx] = fminf(botRightY[tx], botRightY[tx + s]);
//         }
//     }

//     if (tx == 0)
//     {
//         while (atomicCAS(mutex, 0, 1) != 0)
//             ;
//         node[0].topLeft.x = fminf(node[0].topLeft.x, topLeftX[0]);
//         node[0].topLeft.y = fmaxf(node[0].topLeft.y, topLeftY[0]);
//         node[0].botRight.x = fmaxf(node[0].botRight.x, botRightX[0]);
//         node[0].botRight.y = fminf(node[0].botRight.y, botRightY[0]);
//         atomicExch(mutex, 0);
//     }
// }

// __device__ double getDistance(Vector pos1, Vector pos2)
// {

//     return sqrt(pow(pos1.x - pos2.x, 2) + pow(pos1.y - pos2.y, 2));
// }

// __device__ double getWidth(Node &root)
// {

//     return root.botRight.x - root.topLeft.x;
// }

// __device__ bool isCollide(Body &b1, Body &b2)
// {
//     return b1.radius + b2.radius > getDistance(b1.position, b2.position);
// }

// __device__ void ComputeForce(Node *node, int nodeIndex, Body *bodies, int bodyIndex, int nNodes, int nBodies)
// {
//     if (nodeIndex >= nNodes)
//     {
//         return;
//     }
//     Node &curNode = node[nodeIndex];
//     Body &bi = bodies[bodyIndex];
//     if (curNode.bi != -1)
//     {

//         Body &bj = bodies[curNode.bi];
//         if (isCollide(bi, bj))
//             return;

//         Vector rij = {bj.position.x - bi.position.x, bj.position.y - bi.position.y};
//         double inv_r3 = pow(rij.x * rij.x + rij.y * rij.y + E * E, -1.5);
//         double f = (GRAVITY * bj.mass) / inv_r3;
//         Vector force = {rij.x * f, rij.y * f};
//         bi.acceleration.x += (force.x / bi.mass);
//         bi.acceleration.y += (force.y / bi.mass);
//         return;
//     }

//     double sd = getWidth(curNode) / getDistance(bi.position, curNode.centerMass);
//     if (sd < THETA)
//     {
//         Vector rij = {curNode.centerMass.x - bi.position.x, curNode.centerMass.y - bi.position.y};
//         if (bi.radius * 2 > getDistance(bi.position, curNode.centerMass))
//             return;
//         double inv_r3 = pow(rij.x * rij.x + rij.y * rij.y + E * E, -1.5);
//         double f = (GRAVITY * curNode.totalMass) / inv_r3;
//         Vector force = {rij.x * f, rij.y * f};
//         bi.acceleration.x += (force.x / bi.mass);
//         bi.acceleration.y += (force.y / bi.mass);
//         return;
//     }

//     ComputeForce(node, (nodeIndex * 4) + 1, bodies, bodyIndex, nNodes, nBodies);
//     ComputeForce(node, (nodeIndex * 4) + 2, bodies, bodyIndex, nNodes, nBodies);
//     ComputeForce(node, (nodeIndex * 4) + 3, bodies, bodyIndex, nNodes, nBodies);
//     ComputeForce(node, (nodeIndex * 4) + 4, bodies, bodyIndex, nNodes, nBodies);
// }

// __global__ void ComputeForceKernel(Node *node, Body *bodies, int nNodes, int nBodies)
// {

//     int i = blockIdx.x * blockDim.x + threadIdx.x;

//     if (i < nBodies)
//     {
//         Body &bi = bodies[i];
//         if (bi.isDynamic)
//         {

//             bi.velocity.x += bi.acceleration.x * DT / 2.0;
//             bi.velocity.y += bi.acceleration.y * DT / 2.0;

//             bi.position.x += bi.velocity.x * DT;
//             bi.position.y += bi.velocity.y * DT;

//             bi.acceleration = {0.0, 0.0};
//             ComputeForce(node, 0, bodies, i, nNodes, nBodies);
//             bi.velocity.x += bi.acceleration.x * DT / 2.0;
//             bi.velocity.y += bi.acceleration.y * DT / 2.0;
//         }
//     }
// }

#endif