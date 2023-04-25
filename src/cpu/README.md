# N-body-simulation CPU

## Dependencies
* OpenGL
* GLFW
  
## Building

```shell
mkdir build; cd build
cmake ..
cmake --build .
./nBodyCPU <N> <Alg> <Sim> 
```
**N** represents the number of bodies  
**Alg** represents the algorithm to use (0: Direct Sum, 1: Barnes-Hut)  
**Sim** represents the simulation (0: Spiral Galaxy, 1: Random, 2: Our Solar System)
