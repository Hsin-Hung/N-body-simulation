#ifndef DIRECT_SUM_KERNEL_H_
#define DIRECT_SUM_KERNEL_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>

#define N 100

__global__ void DirectSumKernel(Body *bodies, int n)
{
}

void randomizeBodies(Body *bodies)
{
}

int main(int argc, char **argv)
{

    int nBodies = N;
    const unsigned int iters = 20;

    int bytes = nBodies * sizeof(Body);
    Body *bodies = (Body *)malloc(bytes);

    randomizeBodies(bodies);

    Body *d_bodies;
    cudaMalloc((void **)&d_bodies, bytes);
    cudaMemcpy(d_bodies, bodies, bytes, cudaMemcpyHostToDevice);

    int blockSize = nBodies;
    int gridSize = 1;

    DirectSumKernel<<<gridSize, blockSize>>>(d_bodies);

    cudaMemcpy(bodies, d_bodies, bytes, cudaMemcpyDeviceToHost);

    cudaFree(d_bodies);
    free(bodies);
    return 0;
}

#endif