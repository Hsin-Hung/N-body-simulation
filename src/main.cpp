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

NBody nb(NUM_BODIES);
DirectSum alg(nb.bodies);

void updateDots()
{
     alg.update();
}

void drawDots()
{

     glColor3f(1.0, 1.0, 1.0); // set drawing color to white

     for (auto &body : nb.bodies)
     {
          glPointSize(body->radius/5.0);  // set point size to 5 pixels
          glBegin(GL_POINTS); // start drawing points
          glVertex2f(body->position.x, body->position.y);
          glEnd(); // end drawing points
     }
}

void display()
{
     glClear(GL_COLOR_BUFFER_BIT); // clear the screen
     updateDots();
     drawDots();
     // nb.display();
}

int main()
{

     // initialize GLFW
     if (!glfwInit())
          return -1;
     GLFWwindow *window = glfwCreateWindow(WINDOW_WIDTH, WINDOW_HEIGHT, "OpenGL Point", NULL, NULL); // create window
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
     int n = 10; // Run code once every 60 frames
     double last_time = 0;
     int frame_count = 0;
     while (!glfwWindowShouldClose(window)) // main loop
     {

          display();
          glfwSwapBuffers(window); // swap front and back buffers
          glfwPollEvents();        // poll for events
     }

     glfwTerminate(); // terminate GLFW

     return 0;
}