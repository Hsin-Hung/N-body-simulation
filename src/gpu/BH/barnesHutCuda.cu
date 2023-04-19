#include <iostream>
#include <cmath>
#include "barnesHut_kernel.cuh"
#include "constants.h"
#include "err.h"

BarnesHutCuda::BarnesHutCuda(int n) : nBodies(n)
{
    nNodes = MAX_NODES;
    leafLimit = MAX_NODES - N_LEAF;
    h_b = new Body[n];
    h_node = new Node[nNodes];
    h_topLeft = new Vector;
    h_botRight = new Vector;

    CHECK_CUDA_ERROR(cudaMalloc((void **)&d_b, sizeof(Body) * n));
    CHECK_CUDA_ERROR(cudaMalloc((void **)&d_node, sizeof(Node) * nNodes));
    CHECK_CUDA_ERROR(cudaMalloc((void **)&d_mutex, sizeof(int) * nNodes));
    CHECK_CUDA_ERROR(cudaMalloc((void **)&d_topLeft, sizeof(Vector)));
    CHECK_CUDA_ERROR(cudaMalloc((void **)&d_botRight, sizeof(Vector)));
    CHECK_CUDA_ERROR(cudaMalloc((void **)&d_b_buffer, sizeof(Body) * n));
}

BarnesHutCuda::~BarnesHutCuda()
{

    delete[] h_b;
    delete[] h_node;
    delete h_topLeft;
    delete h_botRight;
    CHECK_CUDA_ERROR(cudaFree(d_b));
    CHECK_CUDA_ERROR(cudaFree(d_node));
    CHECK_CUDA_ERROR(cudaFree(d_mutex));
    CHECK_CUDA_ERROR(cudaFree(d_topLeft));
    CHECK_CUDA_ERROR(cudaFree(d_botRight));
    CHECK_CUDA_ERROR(cudaFree(d_b_buffer));
}

void BarnesHutCuda::resetCUDA()
{
    int blockSize = BLOCK_SIZE;
    dim3 gridSize = ceil((float)nNodes / blockSize);
    ResetKernel<<<gridSize, blockSize>>>(d_node, d_topLeft, d_botRight, d_mutex, nNodes, nBodies);
}
void BarnesHutCuda::computeBoundingBoxCUDA()
{
    int blockSize = BLOCK_SIZE;
    dim3 gridSize = ceil((float)nBodies / blockSize);
    ComputeBoundingBoxKernel<<<gridSize, blockSize>>>(d_node, d_b, d_topLeft, d_botRight, d_mutex, nBodies);
}
void BarnesHutCuda::constructQuadTreeCUDA()
{
    int blockSize = BLOCK_SIZE;
    dim3 gridSize = ceil((float)nBodies / blockSize);
    ConstructQuadTreeDPKernel<<<1, blockSize>>>(d_node, d_b, d_b_buffer, 0, nNodes, nBodies, leafLimit);
}
void BarnesHutCuda::computeCenterMassCUDA()
{
    int start = leafLimit, end = nNodes;
    int totalNodes = end - start, temp = 0;
    while (end > start)
    {
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
    int blockSize = 256;
    dim3 gridSize = ceil((float)nBodies / blockSize);
    ComputeForceKernel<<<gridSize, blockSize>>>(d_node, d_b, nNodes, nBodies, leafLimit);
}

void BarnesHutCuda::initRandomBodies()
{
    srand(time(NULL));
    double maxDistance = 2.2790e11;
    double minDistance = 5.0e10;
    Vector centerPos = {CENTERX, CENTERY};
    for (int i = 0; i < nBodies - 1; ++i)
    {

        double angle = 2 * M_PI * (rand() / (double)RAND_MAX);
        // Generate random distance from center within the given max distance
        double radius = (maxDistance - minDistance) * (rand() / (double)RAND_MAX) + minDistance;

        // Calculate coordinates of the point
        double x = CENTERX + radius * std::cos(angle);
        double y = CENTERY + radius * std::sin(angle);
        Vector position = {x, y};
        h_b[i].isDynamic = true;
        h_b[i].mass = 5.974e24;
        h_b[i].radius = 1.3927e6;
        h_b[i].position = position;
        h_b[i].velocity = {0.0, 0.0};
        h_b[i].acceleration = {0.0, 0.0};
    }
    h_b[nBodies - 1].isDynamic = false;
    h_b[nBodies - 1].mass = 1.9890e30;
    h_b[nBodies - 1].radius = 1.3927e6;
    h_b[nBodies - 1].position = centerPos;
    h_b[nBodies - 1].velocity = {0.0, 0.0};
    h_b[nBodies - 1].acceleration = {0.0, 0.0};
}

void BarnesHutCuda::initSpiralBodies()
{

    srand(time(NULL));
    double maxDistance = 2.2790e11;
    double minDistance = 5.0e10;
    Vector centerPos = {CENTERX, CENTERY};
    for (int i = 0; i < nBodies - 1; ++i)
    {

        double angle = 2 * M_PI * (rand() / (double)RAND_MAX);
        // Generate random distance from center within the given max distance
        double radius = (maxDistance - minDistance) * (rand() / (double)RAND_MAX) + minDistance;

        // Calculate coordinates of the point
        double x = CENTERX + radius * std::cos(angle);
        double y = CENTERY + radius * std::sin(angle);
        Vector position = {x, y};

        double distance = sqrt(pow(x - centerPos.x, 2) + pow(y - centerPos.y, 2));
        Vector r = {position.x - centerPos.x, position.y - centerPos.y};
        Vector a = {r.x / distance, r.y / distance};

        // Calculate velocity vector components
        double esc = sqrt((GRAVITY * 1.9891e30) / (distance));
        Vector velocity = {-a.y * esc, a.x * esc};

        h_b[i].isDynamic = true;
        h_b[i].mass = 5.974e24;
        h_b[i].radius = 1.3927e6;
        h_b[i].position = position;
        h_b[i].velocity = velocity;
        h_b[i].acceleration = {0.0, 0.0};
    }
    h_b[nBodies - 1].isDynamic = false;
    h_b[nBodies - 1].mass = 1.9890e30;
    h_b[nBodies - 1].radius = 1.3927e6;
    h_b[nBodies - 1].position = centerPos;
    h_b[nBodies - 1].velocity = {0.0, 0.0};
    h_b[nBodies - 1].acceleration = {0.0, 0.0};
}

void BarnesHutCuda::setBody(int i, bool isDynamic, double mass, double radius, Vector position, Vector velocity, Vector acceleration)
{
    h_b[i].isDynamic = isDynamic;
    h_b[i].mass = mass;
    h_b[i].radius = radius;
    h_b[i].position = position;
    h_b[i].velocity = velocity;
    h_b[i].acceleration = acceleration;
}

void BarnesHutCuda::initSolarSystem()
{
    setBody(0, true, 5.9740e24, 1.3927e6, {1.4960e11, 0}, {0, 2.9800e4}, {0, 0});
    setBody(1, true, 6.4190e23, 1.3927e6, {2.2790e11, 0}, {0, 2.4100e4}, {0, 0});
    setBody(2, true, 3.3020e23, 1.3927e6, {5.7900e10, 0}, {0, 4.7900e4}, {0, 0});
    setBody(3, true, 4.8690e24, 1.3927e6, {1.0820e11, 0}, {0, 3.5000e4}, {0, 0});
    setBody(4, false, 1.9890e30, 1.3927e6, {CENTERX, CENTERY}, {0, 0}, {0, 0});
}

Body *BarnesHutCuda::getBodies()
{

    return h_b;
}

void BarnesHutCuda::readDeviceBodies()
{
    CHECK_CUDA_ERROR(cudaMemcpy(h_b, d_b, sizeof(Body) * nBodies, cudaMemcpyDeviceToHost));
}

void BarnesHutCuda::setup()
{
    initSpiralBodies();
    CHECK_CUDA_ERROR(cudaMemcpy(d_b, h_b, sizeof(Body) * nBodies, cudaMemcpyHostToDevice));
    CHECK_CUDA_ERROR(cudaMemcpy(d_node, h_node, sizeof(Node) * nNodes, cudaMemcpyHostToDevice));
}
void BarnesHutCuda::update()
{
    resetCUDA();
    computeBoundingBoxCUDA();
    constructQuadTreeCUDA();
    computeForceCUDA();
    CHECK_LAST_CUDA_ERROR();
}
