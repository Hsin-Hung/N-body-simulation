#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "nBody.h"
#include "directSum.h"
#include "barnesHut.h"
#include "quadtree.h"
#define NUM_BODIES 5
using namespace std;

int main()
{
     NBody nb;
     BarnesHut alg(nb.bodies);
     for (int i = 0; i < 10; ++i)
     {
          alg.update();
          nb.display();
     }

     return 0;
}