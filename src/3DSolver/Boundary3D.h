#ifndef BOUNDARY_3D_H
#define BOUNDARY_3D_H

/**
 * Very much a work in progress. Will probably be taking a callback function
 * or array at some point.
*/
class Boundary3D {
private:
    double xmin, xmax;
    double ymin, ymax;
    double zmin, zmax;
public:
    Boundary3D(double xmin, double xmax, double ymin, double ymax,
            double zmin, double zmax);
    Boundary3D() = default;
    double *generateXMesh(int n);
    double *generateYMesh(int n);
    double *generateZMesh(int n);
    double getXmin();
    double getXmax();
    double getYmin();
    double getYmax();
    double getZmin();
    double getZmax();
};

#endif