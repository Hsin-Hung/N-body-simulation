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

#ifndef DIRECT_SUM_H
#define DIRECT_SUM_H
#include <vector>
#include <memory>
#include "body.h"
#include "constants.h"
#include "algorithm.h"

class DirectSum : public Algorithm
{
    const double epsilon = 0.5;
    const double dt = 25000.0;

    void calculateAcceleration();
    void calculateVelocity();
    void calculatePosition();
    bool isCollide(Body b1, Body b2);

public:
    DirectSum(std::vector<std::shared_ptr<Body>> &bs);
    void update() override;
};

#endif