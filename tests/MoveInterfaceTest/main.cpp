#include <iostream>
#include "./src/SimUtilities.h"
#include "./src/Boundary.h"

/* Testing the Pool class for shape embedding and other things. */

#include <iostream>
#include <unordered_map>
#include <string>
#include <math.h>
#include "./src/SimUtilities.h"
#include "./src/Boundary.h"
#include "./src/SolidObject.h"
#include "./src/NSSolver.h"

using namespace std;

double circleShapeFun(double x, double y, unordered_map<string, double> &params) {
    double cx = params["cx"];
    double cy = params["cy"];
    double r = params["r"];

    return simutils::square(x - cx) + simutils::square(y - cy) - simutils::square(r);
}

int main() {
    int nx = 128;
    int ny = 128;

    double xmin = 0;
    double xmax = 1;
    double ymin = 0;
    double ymax = 1;

    // Create the boundary object
    Boundary boundary(xmin, xmax, ymin, ymax);

    // Parameters of this circle
    unordered_map<string, double> params;
    params["cx"] = 0.5;
    params["cy"] = 0.5;
    params["r"] = 0.2;

    // Create the solid structure
    SolidObject circle(circleShapeFun, params);

    // Set the embedding tolerance
    double eps = sqrt(simutils::square(1.0/(simutils::square((double)nx)))
                + simutils::square(1.0/simutils::square((double)ny)));
    // double eps = 0.01;

    // Create the pool
    Pool2D pool(nx, ny, eps, boundary, 1, &circle);
    pool.printPool();

    double t = 0;
    double tEnd = 0.2;
    double dt = 0.0001; 

    // Test velocities for the level set methods
    double **u = simutils::new_constant(ny, nx, 0);
    double **v = simutils::new_constant(ny, nx, 0);
    
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            u[j][i] = 1;
            v[j][i] = 0;
        }
    }

    while (t < tEnd) {
        pool.updatePool(dt, u, v, 0, true); // Advance the level set function
        t+=dt;
    }

    pool.enumeratePool(eps);   // Enumerate based on structure and interface position
    pool.printPool();          // Print out the current pool state
    pool.outputPool("out.txt");
}