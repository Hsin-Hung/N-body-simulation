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
    Vector centerMass;
    double totalMass;
    bool isLeaf;

} Node;

class BarnesHutCuda
{
    int nBodies;
    int nNodes;
    int leafLimit;

    Body *h_b;
    Node *h_node;
    Vector *h_topLeft;
    Vector *h_botRight;

    Body *d_b;
    Node *d_node;
    int *d_mutex;
    Vector *d_topLeft;
    Vector *d_botRight;

    bool inBoundary(Vector topLeft, Vector botRight, Vector p);
    int getQuadrant(Vector topLeft, Vector botRight, double x, double y);
    void randomInitBodies();
    double getTotalMass(int nodeIndex);
    void insertQuadTree(int nodeIndex, int b, Vector tl, Vector br);
    void computeForceHelper(int nodeIndex, int b, double width);
    void updateBound(Vector &tl, Vector &br, int quadrant);
    double getDistance(Vector pos1, Vector pos2);
    bool isCollide(Body &b1, Body &b2);

public:
    BarnesHutCuda(int n);
    ~BarnesHutCuda();
    void reset();
    void computeBoundingBox();
    void constructQuadTree();
    void computeCenterMass(int nodeIndex);
    void computeForce();
    void resetCUDA();
    void computeBoundingBoxCUDA();
    void constructQuadTreeCUDA();
    void computeCenterMassCUDA();
    void computeForceCUDA();
    void update();
    void setup();
    Body *getBodies();
};

#endif