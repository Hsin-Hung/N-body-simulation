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