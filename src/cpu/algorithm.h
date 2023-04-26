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

#ifndef ALGORITHM_H
#define ALGORITHM_H
#include <vector>
#include "body.h"

class Algorithm
{
protected:
    std::vector<std::shared_ptr<Body>> &bodies;
    int nBodies;

public:
    Algorithm(std::vector<std::shared_ptr<Body>> &bs, int n);
    virtual ~Algorithm() = default;
    virtual void update() = 0;
};

#endif