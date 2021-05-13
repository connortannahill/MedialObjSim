#ifndef BOUNDARY_H
#define BOUNDARY_H

/**
 * Very much a work in progress. Will probably be taking a callback function
 * or array at some point.
*/
class Boundary {
private:
    double xmin;
    double xmax;
    double ymin;
    double ymax;
public:
    Boundary(double xmin, double xmax, double ymin, double ymax);
    Boundary() = default;
    double *generateXMesh(int n);
    double *generateYMesh(int n);
    double getXmin();
    double getXmax();
    double getYmin();
    double getYmax();
};

#endif