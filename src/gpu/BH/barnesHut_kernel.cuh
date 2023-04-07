#ifndef BARNES_HUT_KERNEL_H_
#define BARNES_HUT_KERNEL_H_
#include <stdio.h>
#include <stdlib.h>
#include "barnesHutCuda.cuh"
__global__ void ResetKernel(Node *node, Vector *topLeft, Vector *botRight, int *mutex, int nNodes);
__global__ void ComputeBoundingBox(Node *node, Body *bodies, Vector *topLeft, Vector *botRight, int *mutex, int nBodies);
__global__ void ConstructQuadTreeKernel(Node *node, Body *bodies, Vector *topLeft, Vector *botRight, int *mutex, int nNodes, int nBodies, int leafLimit);
// __global__ void ComputeForceKernel(Node *node, Body *bodies, int nNodes, int nBodies);
// __device__ void ComputeForce(Node *node, int nodeIndex, Body *bodies, int bodyIndex, int nNodes, int nBodies);

#endif