#ifndef DIRECT_SUM_KERNEL_H_
#define DIRECT_SUM_KERNEL_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <vector>
#include <opencv2/opencv.hpp>
#include "constants.h"
#include "err.h"

#define BLOCK_SIZE 256

cv::VideoWriter video("nbody.avi", cv::VideoWriter::fourcc('M', 'J', 'P', 'G'), 30, cv::Size(WINDOW_WIDTH, WINDOW_HEIGHT));

typedef struct
{
    double x;
    double y;
} Vector;

typedef struct
{
    bool isDynamic;
    double mass;
    double radius;
    Vector position;
    Vector velocity;
    Vector acceleration;

} Body;

__device__ double getDistance(Vector pos1, Vector pos2)
{

    return sqrt(pow(pos1.x - pos2.x, 2) + pow(pos1.y - pos2.y, 2));
}

__device__ bool isCollide(Body &b1, Body &b2)
{
    return b1.radius + b2.radius > getDistance(b1.position, b2.position);
}

__device__ void calculateAcceleration(Body &bi, int i, Body *bodies, int n)
{

    bi.acceleration = {0.0, 0.0};
    for (int b = 0; b < n; ++b)
    {
        if (b != i)
        {
            Body &bj = bodies[b];
            if (!isCollide(bi, bj))
            {
                Vector rij = {bj.position.x - bi.position.x, bj.position.y - bi.position.y};
                double inv_r3 = pow(rij.x * rij.x + rij.y * rij.y + E * E, -1.5);
                double f = (GRAVITY * bj.mass) / inv_r3;
                Vector force = {rij.x * f, rij.y * f};
                bi.acceleration.x += (force.x / bi.mass);
                bi.acceleration.y += (force.y / bi.mass);
            }
        }
    }
}

__device__ void calculateVelocity(Body &bi)
{

    bi.velocity.x += bi.acceleration.x * DT / 2.0;
    bi.velocity.y += bi.acceleration.y * DT / 2.0;
}

__device__ void calculatePosition(Body &bi)
{
    int boundaryWidth = WINDOW_WIDTH, boundaryHeight = WINDOW_HEIGHT;
    bi.position.x += bi.velocity.x * DT;
    bi.position.y += bi.velocity.y * DT;
    if (bi.position.x < 0)
    {
        bi.position.x = 0;
        bi.velocity.x *= -1.0f;
    }
    if (bi.position.x > boundaryWidth)
    {
        bi.position.x = boundaryWidth;
        bi.velocity.x *= -1.0f;
    }
    if (bi.position.y < 0)
    {
        bi.position.y = 0;
        bi.velocity.y *= -1.0f;
    }
    if (bi.position.y > boundaryHeight)
    {
        bi.position.y = boundaryHeight;
        bi.velocity.y *= -1.0f;
    }
}

__global__ void DirectSumKernel(Body *bodies, int n)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i < n)
    {
        Body &bi = bodies[i];
        if (bi.isDynamic)
        {

            bi.velocity.x += bi.acceleration.x * DT / 2.0;
            bi.velocity.y += bi.acceleration.y * DT / 2.0;

            bi.position.x += bi.velocity.x * DT;
            bi.position.y += bi.velocity.y * DT;

            bi.acceleration = {0.0, 0.0};
            for (int b = 0; b < n; ++b)
            {
                if (b != i)
                {
                    Body &bj = bodies[b];
                    if (!isCollide(bi, bj))
                    {
                        Vector rij = {bj.position.x - bi.position.x, bj.position.y - bi.position.y};
                        double inv_r3 = pow(rij.x * rij.x + rij.y * rij.y + E * E, -1.5);
                        double f = (GRAVITY * bj.mass) / inv_r3;
                        Vector force = {rij.x * f, rij.y * f};
                        bi.acceleration.x += (force.x / bi.mass);
                        bi.acceleration.y += (force.y / bi.mass);
                    }
                }
            }
            bi.velocity.x += bi.acceleration.x * DT / 2.0;
            bi.velocity.y += bi.acceleration.y * DT / 2.0;
        }
    }
}

__global__ void DirectSumTiledKernel(Body *bodies, int n)
{
    __shared__ Body Bds[BLOCK_SIZE];

    int bx = blockIdx.x;
    int tx = threadIdx.x;
    int i = bx * blockDim.x + tx;

    if (i < n)
    {
        Body &bi = bodies[i];
        double fx = 0.0, fy = 0.0;
        bi.velocity.x += bi.acceleration.x * DT / 2.0;
        bi.velocity.y += bi.acceleration.y * DT / 2.0;
        bi.position.x += bi.velocity.x * DT;
        bi.position.y += bi.velocity.y * DT;
        bi.acceleration = {0.0, 0.0};
        for (int tile = 0; tile < gridDim.x; ++tile)
        {

            Bds[tx] = bodies[tile * blockDim.x + tx];
            __syncthreads();

            for (int b = 0; b < BLOCK_SIZE; ++b)
            {
                int j = tile * blockDim.x + b;
                if (j < n)
                {
                    Body bj = Bds[b];
                    if (!isCollide(bi, bj))
                    {
                        Vector rij = {bj.position.x - bi.position.x, bj.position.y - bi.position.y};
                        double inv_r3 = pow(rij.x * rij.x + rij.y * rij.y + E * E, -1.5);
                        double f = (GRAVITY * bj.mass) / inv_r3;
                        Vector force = {rij.x * f, rij.y * f};
                        fx += (force.x / bi.mass);
                        fy += (force.y / bi.mass);
                    }
                }
            }
            __syncthreads();
        }
        bi.acceleration.x += fx;
        bi.acceleration.y += fy;
        bi.velocity.x += bi.acceleration.x * DT / 2.0;
        bi.velocity.y += bi.acceleration.y * DT / 2.0;
    }
}

void display(Body *bodies)
{

    std::cout << bodies[0].position.x << " " << bodies[0].position.y << std::endl;
}

void storeFrame(Body *bodies, int n, int id)
{
    cv::Mat image = cv::Mat::zeros(WINDOW_HEIGHT, WINDOW_WIDTH, CV_8UC3);
    // Generate random points
    std::vector<cv::Point2f> points;
    cv::RNG rng;
    for (int i = 0; i < n; i++)
    {
        points.push_back(cv::Point2f(bodies[i].position.x, bodies[i].position.y));
    }

    cv::Scalar color = cv::Scalar(255, 255, 255); // White color
    for (int i = 0; i < n; i++)
    {
        cv::Point center(points[i].x, points[i].y);
        cv::circle(image, center, bodies[i].radius, color, -1);
    }
    video.write(image);
    // cv::imwrite("frame" + std::to_string(id) + ".jpg", image);
}

Body *initRandomBodies(int n)
{

    Body *bodies = new Body[n];
    srand(time(NULL));
    int maxDistance = 200;
    for (int i = 0; i < n - 1; ++i)
    {
        double angle = 2 * M_PI * (rand() / (double)RAND_MAX);
        // Generate random distance from center within the given max distance
        double distance = maxDistance * (rand() / (double)RAND_MAX);

        // Calculate coordinates of the point
        double x = CENTERX + distance * std::cos(angle);
        double y = CENTERY + distance * std::sin(angle);
        Vector position = {x, y};
        bodies[i].isDynamic = true;
        bodies[i].mass = 1.0 / (double)n;
        bodies[i].radius = 1;
        bodies[i].position = position;
        bodies[i].velocity = {0.0, 0.0};
        bodies[i].acceleration = {0.0, 0.0};
    }

    bodies[n - 1].isDynamic = false;
    bodies[n - 1].mass = 200.0 / (double)n;
    bodies[n - 1].radius = 1;
    bodies[n - 1].position = {CENTERX, CENTERY};
    bodies[n - 1].velocity = {0.0, 0.0};
    bodies[n - 1].acceleration = {0.0, 0.0};

    return bodies;
}

Body *initSpiralBodies(int n)
{

    Body *bodies = new Body[n];
    srand(time(NULL));
    int maxDistance = 200;
    for (int i = 0; i < n; ++i)
    {

        double angle = 2 * M_PI * (rand() / (double)RAND_MAX);
        // Generate random distance from center within the given max distance
        double distance = maxDistance * (rand() / (double)RAND_MAX);

        // Calculate coordinates of the point
        double x = CENTERX + distance * std::cos(angle);
        double y = CENTERY + distance * std::sin(angle);

        Vector position = {x, y};
        Vector r = {x - CENTERX, y - CENTERY};
        Vector velocity = {r.y, -r.x};

        bodies[i].isDynamic = true;
        bodies[i].mass = 1.0 / (double)n;
        bodies[i].radius = 1;
        bodies[i].position = position;
        bodies[i].velocity = velocity;
        bodies[i].acceleration = {0.0, 0.0};
    }

    return bodies;
}

int main(int argc, char **argv)
{
    int nBodies = NUM_BODIES;
    int iters = 300;
    if (argc == 3)
    {
        nBodies = atoi(argv[1]);
        iters = atoi(argv[2]);
    }

    Body *h_bodies = initRandomBodies(nBodies);

    int bytes = nBodies * sizeof(Body);

    Body *d_bodies;
    CHECK_CUDA_ERROR(cudaMalloc((void **)&d_bodies, bytes));
    CHECK_CUDA_ERROR(cudaMemcpy(d_bodies, h_bodies, bytes, cudaMemcpyHostToDevice));

    int blockSize = BLOCK_SIZE;
    int gridSize = ceil((double)nBodies / blockSize);
    int it = 0;
    while (it < iters) // main loop
    {

        DirectSumTiledKernel<<<gridSize, blockSize>>>(d_bodies, nBodies);
        CHECK_LAST_CUDA_ERROR();
        CHECK_CUDA_ERROR(cudaMemcpy(h_bodies, d_bodies, bytes, cudaMemcpyDeviceToHost));
        storeFrame(h_bodies, nBodies, ++it);
        // display(bodies);
    }
    video.release();

    // free memories
    CHECK_CUDA_ERROR(cudaFree(d_bodies));
    free(h_bodies);

    CHECK_LAST_CUDA_ERROR();
    return 0;
}

#endif