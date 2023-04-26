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

#ifndef BODY_H
#define BODY_H
#include <iostream>
#include "vector.h"

class Body
{
public:
    bool isDynamic;
    double mass;
    double radius;
    Vector position;
    Vector velocity;
    Vector acceleration;

    Body(double m, double r, Vector p, Vector v, Vector a, bool d = true) : mass(m), radius(r), position(p), velocity(v), acceleration(a), isDynamic(d) {}
    friend std::ostream &operator<<(std::ostream &os, const Body &b);
};

#endif