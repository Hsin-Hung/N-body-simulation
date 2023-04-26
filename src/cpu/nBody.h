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

#ifndef NBODY_H
#define NBODY_H
#include <vector>
#include <memory>
#include "vector.h"
#include "body.h"
#include "algorithm.h"

class NBody
{
    std::unique_ptr<Algorithm> alg;
    int nBodies;
    void initRandomBodies();
    void initSpiralBodies();
    void initSolarSystem();

public:
    std::vector<std::shared_ptr<Body>> bodies;
    NBody(int n, int a, int s);
    void update();
};

#endif