#include <iostream>
#include "./src/SimUtilities.h"
#include "./src/Boundary.h"
#include "./src/SolidObject.h"
#include "./src/NSSolver.h"

using namespace std;

void initialConditions(int nx, int ny, int nGhost, double *x, double *y, double **u, double **v) {
    int i, j;

    // Set up the initial conditions on the staggered grid.
    for (j = 1; j <= ny; j++) {
        for (i = 1; i < nx; i++) {
            u[j][i] = 0.0;
        }
    }

    // Set values on known boundaries
    for (j = 1; j <= ny; j++) {
        u[j][0] = 0.0;
        u[j][nx] = 0.0;
    }

    for (j = 1; j < ny; j++) {
        for (i = 1; i <= nx; i++) {
            v[j][i] = 0.0;
        }
    }

    for (i = 1; i <= nx; i++) {
        v[0][i] = 0.0;
        v[ny][i] = 0.0;
    }
}


int main() {
    // Testing the fluid solver
    double tEnd = 1.4;

    // The boundaries
    double xa = 0, xb = 1;
    double ya = 0, yb = 1;

    // Number of x, y points
    int nx = 164;
    int ny = 164;
    // int nx = 16;
    // int ny = 16;

    // Boundary object
    Boundary boundary(xa, xb, ya, yb);

    // Solid object null pointer
    SolidObject *solidObjects = NULL;

    // Hashmap of any required params
    std::unordered_map<std::string, double> params;
    params["Re"] = 1000.0;

    // Create the Solver object
    NSSolver solver(nx, ny, boundary, 0, solidObjects, params, initialConditions);

    // Current time
    double t = 0;
    double safetyFactor = 0.5;

    t = solver.step(tEnd, safetyFactor);
    while (t < tEnd) {
        t = solver.step(tEnd, safetyFactor);
        std::cout << "t = " << t;
        std::cout << std::endl;
    }

    // Write to the output file.
    solver.writeToFile("out.txt");

}