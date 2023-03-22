#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <GLFW/glfw3.h>
#include "nBody.h"
#include "directSum.h"
#include "barnesHut.h"
#include "quadtree.h"
#define NUM_BODIES 5
#define WINDOW_WIDTH 1600
#define WINDOW_HEIGHT 1600
using namespace std;

const int FRAME_INTERVAL = 10;
NBody nb;
DirectSum alg(nb.bodies);
void updateDots()
{
     static int frame_count = 0;

     // if (++frame_count < FRAME_INTERVAL)
     //      return;

     frame_count = 0;
     alg.update();
}

void drawDots()
{

     glColor3f(1.0, 1.0, 1.0); // set drawing color to white
     glPointSize(10.0);         // set point size to 5 pixels
     glBegin(GL_POINTS);       // start drawing points
     for (auto &body : nb.bodies)
     {
          glVertex2f(body->position.x, body->position.y);
     }
     glEnd(); // end drawing points
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
     // for (int i = 0; i < 10; ++i)
     // {
     //      alg.update();
     //      nb.display();
     // }
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