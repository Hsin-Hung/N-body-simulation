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

#ifndef BARNES_HUT_KERNEL_H_
#define BARNES_HUT_KERNEL_H_
#include <stdio.h>
#include <stdlib.h>
#include "barnesHutCuda.cuh"

__global__ void ResetKernel(Node *node, int *mutex, int nNodes, int nBodies);
__global__ void ComputeBoundingBoxKernel(Node *node, Body *bodies, int *mutex, int nBodies);
__global__ void ConstructQuadTreeKernel(Node *node, Body *bodies, Body *buffer, int nodeIndex, int nNodes, int nBodies, int leafLimit);
__global__ void ComputeForceKernel(Node *node, Body *bodies, int nNodes, int nBodies, int leafLimit);

#endif