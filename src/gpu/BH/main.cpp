#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <vector>
#include <opencv2/opencv.hpp>
#include "barnesHutCuda.cuh"
#include "constants.h"

cv::VideoWriter video("nbody.avi", cv::VideoWriter::fourcc('M', 'J', 'P', 'G'), 30, cv::Size(WINDOW_WIDTH, WINDOW_HEIGHT));

void display(Body *bodies, int nBodies)
{
    for (int i = 0; i < nBodies; ++i)
    {
        std::cout << "Body: id:" << bodies[i].id << "  " << bodies[i].position.x << " " << bodies[i].position.y << std::endl;
    }

    std::cout << std::endl;
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

    int radius = 2;
    cv::Scalar color = cv::Scalar(255, 255, 255); // White color
    for (int i = 0; i < n; i++)
    {
        cv::Point center(points[i].x, points[i].y);
        cv::circle(image, center, bodies[i].radius, color, -1);
    }
    video.write(image);
    // cv::imwrite("frame" + std::to_string(id) + ".jpg", image);
}

int main(int argc, char **argv)
{
    int nBodies = NUM_BODIES;
    int iters = 1;
    if (argc == 3)
    {
        nBodies = atoi(argv[1]);
        iters = atoi(argv[2]);
    }
    BarnesHutCuda *bh = new BarnesHutCuda(nBodies);
    bh->setup();
    // display(bh->getBodies(), nBodies);
    for (int i = 0; i < iters; ++i)
    {
        bh->update();
        // display(bh->getBodies(), nBodies);
        storeFrame(bh->getBodies(), nBodies, i);
    }
    video.release();
    delete bh;
    return 0;
}
