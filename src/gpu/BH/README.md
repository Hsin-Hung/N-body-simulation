# Barnes-Hut CUDA

## Dependencies

* NVIDIA GPU with CUDA support
* OpenCV
  
## Building

```shell
mkdir build; cd build
cmake ..
cmake --build .
./BarnesHut <N> <Sim> <Iter>
```
**N**: the number of bodies  
**Sim**: the simulation (0: Spiral Galaxy, 1: Random, 2: Colliding Galaxies, 3: Our Solar System)  
**Iter**: number of iterations/frames

## Visualization
An AVI file of the simulation with `<Iter>` frames will be produced.