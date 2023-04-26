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

#ifndef QUADTREE_H
#define QUADTREE_H
#include <memory>
#include "body.h"
#include "vector.h"

class BarnesHut;

class QuadTree
{

    Vector topLeft;
    Vector botRight;
    Vector centerMass;
    double totalMass;
    bool isLeaf;
    std::shared_ptr<Body> b;
    std::unique_ptr<QuadTree> topLeftTree;
    std::unique_ptr<QuadTree> topRightTree;
    std::unique_ptr<QuadTree> botLeftTree;
    std::unique_ptr<QuadTree> botRightTree;

public:
    QuadTree();
    QuadTree(Vector topL, Vector botR);
    ~QuadTree();
    void insert(std::shared_ptr<Body> n);
    std::shared_ptr<Body> search(Vector point);
    bool inBoundary(Vector point);
    double getWidth();
    Vector getCenter();
    int getQuadrant(Vector pos);
    friend void updateCenterMass(std::unique_ptr<QuadTree> &root);
    friend void traverse(std::unique_ptr<QuadTree> &root);
    friend double getTotalMass(std::unique_ptr<QuadTree> &root);

    friend class BarnesHut;
};

#endif