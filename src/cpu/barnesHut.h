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

#ifndef BARNES_HUT_H
#define BARNES_HUT_H
#include <vector>
#include <memory>
#include "body.h"
#include "quadtree.h"
#include "constants.h"
#include "algorithm.h"

class BarnesHut : public Algorithm
{
    const double epsilon = 0.5;
    const double dt = 25000.0;
    const double theta = 0.5;

    std::unique_ptr<QuadTree> quadTree;
    void constructQuadTree();
    void computeCenterMass();
    void calculateForceHelper(std::unique_ptr<QuadTree> &root, std::shared_ptr<Body> body);
    void computeBoundingBox();
    void calculateForce(std::shared_ptr<Body> b);
    void calculateAcceleration();
    void calculateVelocity();
    void calculatePosition();
    bool isCollide(Body b1, Body b2);
    bool isCollide(Body b, Vector cm);

public:
    BarnesHut(std::vector<std::shared_ptr<Body>> &bs);
    void update() override;
};

#endif