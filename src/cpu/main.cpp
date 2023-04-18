#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <GLFW/glfw3.h>
#include "nBody.h"
#include "directSum.h"
#include "barnesHut.h"
#include "quadtree.h"
#include "constants.h"

using namespace std;

Vector scaleToWindow(Vector pos)
{

     double scaleX = WINDOW_HEIGHT / NBODY_HEIGHT;
     double scaleY = WINDOW_WIDTH / NBODY_WIDTH;
     return Vector((pos.x - 0) * scaleX + 800, (pos.y - 0) * scaleY + 800);
}

void drawDots(NBody &nb)
{

     glColor3f(1.0, 1.0, 1.0); // set drawing color to white

     for (auto &body : nb.bodies)
     {
          glPointSize(5);     // set point size to 5 pixels
          glBegin(GL_POINTS); // start drawing points
          Vector pos = scaleToWindow(body->position);
          // std::cout << "org: " << body->position << std::endl;
          glVertex2f(pos.x, pos.y);
          glEnd(); // end drawing points
     }
}

int main(int argc, char **argv)
{
     int nBodies = NUM_BODIES;
     int i = 0;
     if (argc == 3)
     {
          nBodies = atoi(argv[1]);
          i = atoi(argv[2]);
     }

     NBody nb(nBodies, i);

     // initialize GLFW
     if (!glfwInit())
          return -1;
     GLFWwindow *window = glfwCreateWindow(WINDOW_WIDTH, WINDOW_HEIGHT, "nBody", NULL, NULL); // create window
     if (!window)
     {
          glfwTerminate();
          return -1;
     }
     glfwMakeContextCurrent(window); // set context to current window

     glClearColor(0.0, 0.0, 0.0, 1.0); // set background color to black
     glMatrixMode(GL_PROJECTION);      // set up projection matrix
     glLoadIdentity();
     glOrtho(0.0f, WINDOW_WIDTH, WINDOW_HEIGHT, 0.0f, -1.0f, 1.0f);
     while (!glfwWindowShouldClose(window)) // main loop
     {
          glClear(GL_COLOR_BUFFER_BIT); // clear the screen

          nb.update();
          drawDots(nb);

          glfwSwapBuffers(window); // swap front and back buffers
          glfwPollEvents();        // poll for events
     }

     glfwTerminate(); // terminate GLFW

     return 0;
}