#include "Boundary3D.h"
#include "../Utils/SimUtilities.h"
#include <iostream>

Boundary3D::Boundary3D(double xmin, double xmax, double ymin, double ymax,
                        double zmin, double zmax) {
    this->xmin = xmin;
    this->xmax = xmax;
    this->ymin = ymin;
    this->ymax = ymax;
    this->zmin = zmin;
    this->zmax = zmax;
}

// Generates mesh of n subintervals withing the bounding region.
double *Boundary3D::generateXMesh(int nx) {
    double *x = new double[nx+1];
    simutils::linspace(this->xmin, this->xmax, nx, x);

    return x;
}

double *Boundary3D::generateYMesh(int ny) {
    double *y = new double[ny+1];
    simutils::linspace(this->ymin, this->ymax, ny, y);
    return y;
}

double *Boundary3D::generateZMesh(int nz) {
    double *z = new double[nz+1];
    simutils::linspace(this->zmin, this->zmax, nz, z);
    return z;

}

double Boundary3D::getXmin() {
    return this->xmin;
}

double Boundary3D::getXmax() {
    return this->ymax;
}

double Boundary3D::getYmin() {
    return this->ymin;
}

double Boundary3D::getYmax() {
    return this->ymax;
}

double Boundary3D::getZmin() {
    return this->zmin;
}

double Boundary3D::getZmax() {
    return this->zmax;
}

