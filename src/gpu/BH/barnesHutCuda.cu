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

    CHECK_CUDA_ERROR(cudaMalloc((void **)&d_b, sizeof(Body) * n));
    CHECK_CUDA_ERROR(cudaMalloc((void **)&d_node, sizeof(Node) * nNodes));
    CHECK_CUDA_ERROR(cudaMalloc((void **)&d_mutex, sizeof(int) * nNodes));
    CHECK_CUDA_ERROR(cudaMalloc((void **)&d_b_buffer, sizeof(Body) * n));
}

BarnesHutCuda::~BarnesHutCuda()
{

    delete[] h_b;
    delete[] h_node;
    CHECK_CUDA_ERROR(cudaFree(d_b));
    CHECK_CUDA_ERROR(cudaFree(d_node));
    CHECK_CUDA_ERROR(cudaFree(d_mutex));
    CHECK_CUDA_ERROR(cudaFree(d_b_buffer));
}

void BarnesHutCuda::resetCUDA()
{
    int blockSize = BLOCK_SIZE;
    dim3 gridSize = ceil((float)nNodes / blockSize);
    ResetKernel<<<gridSize, blockSize>>>(d_node, d_mutex, nNodes, nBodies);
}
void BarnesHutCuda::computeBoundingBoxCUDA()
{
    int blockSize = BLOCK_SIZE;
    dim3 gridSize = ceil((float)nBodies / blockSize);
    ComputeBoundingBoxKernel<<<gridSize, blockSize>>>(d_node, d_b, d_mutex, nBodies);
}
void BarnesHutCuda::constructQuadTreeCUDA()
{
    int blockSize = BLOCK_SIZE;
    dim3 gridSize = ceil((float)nBodies / blockSize);
    ConstructQuadTreeKernel<<<1, blockSize>>>(d_node, d_b, d_b_buffer, 0, nNodes, nBodies, leafLimit);
}

void BarnesHutCuda::computeForceCUDA()
{
    int blockSize = 32;
    dim3 gridSize = ceil((float)nBodies / blockSize);
    ComputeForceKernel<<<gridSize, blockSize>>>(d_node, d_b, nNodes, nBodies, leafLimit);
}

void BarnesHutCuda::initRandomBodies()
{
    srand(time(NULL));
    double maxDistance = MAX_DIST;
    double minDistance = MIN_DIST;
    Vector centerPos = {CENTERX, CENTERY};
    for (int i = 0; i < nBodies - 1; ++i)
    {

        double angle = 2 * M_PI * (rand() / (double)RAND_MAX);
        // Generate random distance from center within the given max distance
        double radius = (maxDistance - minDistance) * (rand() / (double)RAND_MAX) + minDistance;

        // Calculate coordinates of the point
        double x = centerPos.x + radius * std::cos(angle);
        double y = centerPos.y + radius * std::sin(angle);
        Vector position = {x, y};
        h_b[i].isDynamic = true;
        h_b[i].mass = EARTH_MASS;
        h_b[i].radius = EARTH_DIA;
        h_b[i].position = position;
        h_b[i].velocity = {0.0, 0.0};
        h_b[i].acceleration = {0.0, 0.0};
    }
    h_b[nBodies - 1].isDynamic = false;
    h_b[nBodies - 1].mass = SUN_MASS;
    h_b[nBodies - 1].radius = SUN_DIA;
    h_b[nBodies - 1].position = centerPos;
    h_b[nBodies - 1].velocity = {0.0, 0.0};
    h_b[nBodies - 1].acceleration = {0.0, 0.0};
}

void BarnesHutCuda::initSpiralBodies()
{

    srand(time(NULL));
    double maxDistance = MAX_DIST;
    double minDistance = MIN_DIST;
    Vector centerPos = {CENTERX, CENTERY};
    for (int i = 0; i < nBodies - 1; ++i)
    {

        double angle = 2 * M_PI * (rand() / (double)RAND_MAX);
        // Generate random distance from center within the given max distance
        double radius = (maxDistance - minDistance) * (rand() / (double)RAND_MAX) + minDistance;

        // Calculate coordinates of the point
        double x = centerPos.x + radius * std::cos(angle);
        double y = centerPos.y + radius * std::sin(angle);
        Vector position = {x, y};

        double distance = sqrt(pow(x - centerPos.x, 2) + pow(y - centerPos.y, 2));
        Vector r = {position.x - centerPos.x, position.y - centerPos.y};
        Vector a = {r.x / distance, r.y / distance};

        // Calculate velocity vector components
        double esc = sqrt((GRAVITY * SUN_MASS) / (distance));
        Vector velocity = {-a.y * esc, a.x * esc};

        h_b[i].isDynamic = true;
        h_b[i].mass = EARTH_MASS;
        h_b[i].radius = EARTH_DIA;
        h_b[i].position = position;
        h_b[i].velocity = velocity;
        h_b[i].acceleration = {0.0, 0.0};
    }
    h_b[nBodies - 1].isDynamic = false;
    h_b[nBodies - 1].mass = SUN_MASS;
    h_b[nBodies - 1].radius = SUN_DIA;
    h_b[nBodies - 1].position = centerPos;
    h_b[nBodies - 1].velocity = {0.0, 0.0};
    h_b[nBodies - 1].acceleration = {0.0, 0.0};
}

void BarnesHutCuda::initCollideGalaxy()
{

    srand(time(NULL));
    double maxDistance = MAX_DIST / 4.0;
    double minDistance = MIN_DIST;
    Vector centerPos = {-NBODY_WIDTH / 6.0, CENTERY};

    int galaxy1 = nBodies / 2;

    for (int i = 0; i < galaxy1 - 1; ++i)
    {

        double angle = 2 * M_PI * (rand() / (double)RAND_MAX);
        // Generate random distance from center within the given max distance
        double radius = (maxDistance - minDistance) * (rand() / (double)RAND_MAX) + minDistance;

        // Calculate coordinates of the point
        double x = centerPos.x + radius * std::cos(angle);
        double y = centerPos.y + radius * std::sin(angle);
        Vector position = {x, y};

        double distance = sqrt(pow(x - centerPos.x, 2) + pow(y - centerPos.y, 2));
        Vector r = {position.x - centerPos.x, position.y - centerPos.y};
        Vector a = {r.x / distance, r.y / distance};

        // Calculate velocity vector components
        double esc = sqrt((GRAVITY * SUN_MASS) / (distance));
        Vector velocity = {-a.y * esc, a.x * esc};

        h_b[i].isDynamic = true;
        h_b[i].mass = EARTH_MASS;
        h_b[i].radius = EARTH_DIA;
        h_b[i].position = position;
        h_b[i].velocity = velocity;
        h_b[i].acceleration = {0.0, 0.0};
    }
    h_b[galaxy1 - 1].isDynamic = true;
    h_b[galaxy1 - 1].mass = SUN_MASS;
    h_b[galaxy1 - 1].radius = SUN_DIA;
    h_b[galaxy1 - 1].position = centerPos;
    h_b[galaxy1 - 1].velocity = {0.0, 0.0};
    h_b[galaxy1 - 1].acceleration = {0.0, 0.0};

    centerPos = {NBODY_WIDTH / 6.0, CENTERY};

    for (int i = galaxy1; i < nBodies - 1; ++i)
    {

        double angle = 2 * M_PI * (rand() / (double)RAND_MAX);
        // Generate random distance from center within the given max distance
        double radius = (maxDistance - minDistance) * (rand() / (double)RAND_MAX) + minDistance;

        // Calculate coordinates of the point
        double x = centerPos.x + radius * std::cos(angle);
        double y = centerPos.y + radius * std::sin(angle);
        Vector position = {x, y};

        double distance = sqrt(pow(x - centerPos.x, 2) + pow(y - centerPos.y, 2));
        Vector r = {position.x - centerPos.x, position.y - centerPos.y};
        Vector a = {r.x / distance, r.y / distance};

        // Calculate velocity vector components
        double esc = sqrt((GRAVITY * SUN_MASS) / (distance));
        Vector velocity = {-a.y * esc, a.x * esc};

        h_b[i].isDynamic = true;
        h_b[i].mass = EARTH_MASS;
        h_b[i].radius = EARTH_DIA;
        h_b[i].position = position;
        h_b[i].velocity = velocity;
        h_b[i].acceleration = {0.0, 0.0};
    }
    h_b[nBodies - 1].isDynamic = true;
    h_b[nBodies - 1].mass = SUN_MASS;
    h_b[nBodies - 1].radius = SUN_DIA;
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

void BarnesHutCuda::setup(int sim)
{
    if (sim == 0)
    {
        initSpiralBodies();
    }
    else if (sim == 1)
    {
        initRandomBodies();
    }
    else if (sim == 2)
    {
        initCollideGalaxy();
    }
    else
    {
        initSolarSystem();
    }

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
