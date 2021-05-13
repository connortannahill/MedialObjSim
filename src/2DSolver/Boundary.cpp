#include "Boundary.h"
#include "../Utils/SimUtilities.h"
#include <iostream>

Boundary::Boundary(double xmin, double xmax, double ymin, double ymax) {
    this->xmin = xmin;
    this->xmax = xmax;
    this->ymin = ymin;
    this->ymax = ymax;
}

// Generates mesh of n subintervals withing the bounding region.
double *Boundary::generateXMesh(int nx) {
    double *x = new double[nx+1];
    simutils::linspace(this->xmin, this->xmax, nx, x);

    return x;
}

double *Boundary::generateYMesh(int ny) {
    double *y = new double[ny+1];
    simutils::linspace(this->ymin, this->ymax, ny, y);
    return y;
}

double Boundary::getXmin() {
    return this->xmin;
}

double Boundary::getXmax() {
    return this->ymax;
}

double Boundary::getYmin() {
    return this->ymin;
}

double Boundary::getYmax() {
    return this->ymax;
}

