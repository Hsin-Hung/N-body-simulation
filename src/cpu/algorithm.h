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