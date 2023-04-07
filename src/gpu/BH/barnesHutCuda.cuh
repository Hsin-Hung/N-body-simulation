#ifndef BARNES_HUT_CUDA_H_
#define BARNES_HUT_CUDA_H_

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

typedef struct
{

    Vector topLeft;
    Vector botRight;
    Vector centerMass;
    double totalMass;
    bool isLeaf;

} Node;

class BarnesHutCuda
{
    int nBodies;
    int nNodes;
    int leaf_limit;

    Body *h_b;
    Node *h_node;

    Body *d_b;
    Node *d_node;
    int *d_mutex;

    bool inBoundary(Vector topLeft, Vector botRight, Vector p);
    int getQuadrant(Vector topLeft, Vector botRight, double x, double y);
    void randomInitBodies();
    double getTotalMass(int nodeIndex);
    void insertQuadTree(int nodeIndex, int b);
    void computeForceHelper(int nodeIndex, int b);
    double getDistance(Vector pos1, Vector pos2);
    double getWidth(Node &root);
    bool isCollide(Body &b1, Body &b2);

public:
    BarnesHutCuda(int n);
    ~BarnesHutCuda();
    void reset();
    void constructQuadTree();
    void computeBoundingBox();
    void computeCenterMass(int nodeIndex);
    void computeForce();
    void update();
    void setup();
    Body *getBodies();
};

#endif