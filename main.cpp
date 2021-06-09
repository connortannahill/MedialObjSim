#include <iostream>
#include "./src/Utils/SimUtilities.h"
#include "./src/2DSolver/Boundary.h"
#include "./src/2DSolver/SolidObject.h"
#include "./src/2DSolver/NSSolver.h"
#include <math.h>
#include "./src/Utils/Discretizations.h"
#include <cassert>
#include <fenv.h>
#include "./src/2DSolver/ObjectSeeder.h"
#include "./src/2DSolver/SimParams.h"
#include "./src/Utils/TestFormatter.h"

using namespace std;
const double EPS = 1e-16;

double circleShapeFun(double x, double y, SolidParams &ps) {
    double cx, cy, r;
    ps.getParam("cx", cx);
    ps.getParam("cy", cy);
    ps.getParam("r", r);
    return simutils::square(x-cx)+simutils::square(y-cy)-simutils::square(r);
}

double coneShapeFun(double x, double y, SolidParams &ps) {
    double cx, cy, r;
    ps.getParam("cx", cx);
    ps.getParam("cy", cy);
    ps.getParam("r", r);
    return sqrt(simutils::square(x-cx)+simutils::square(y-cy))-r;
}

void initialConditions(int nx, int ny, int nGhost, double *x, double *y, double **u, double **v) {
    int i, j;
    double cons_u = 0.0;
    double cons_v = 0.0;

    int wg = nx + 2*nGhost; // "width"
    int hg = ny + 2*nGhost; // "height"

    simutils::set_constant(hg, wg-1, cons_u, u);
    simutils::set_constant(hg-1, wg, cons_v, v);
}

// Problem-dependent boundary condition
void boundaryConditions(int nx, int ny, double **u, double **v) {
    // Set all the boundary conditions to 0.
    for (int j = 1; j <= ny; j++) {
        u[j][0] = 0.0;
        u[j][nx] = 0.0;

        v[j][0] = -v[j][1];
        v[j][nx+1] = -v[j][nx];
    }

    for (int i = 1; i <= nx; i++) {
        u[0][i] = -u[1][i];
        u[ny+1][i] = -u[ny][i];

        v[0][i] = 0.0;
        v[ny][i] = 0.0;
    }

    // Lid driven cavity example
    double ubar = 1;
    for (int i = 1; i <= nx; i++) {
        u[ny+1][i] = 2*ubar - u[ny][i];
    }

    // // Flow past obstacle
    // for (int j = 1; j <= ny; j++) {
    //     // Inflow condition
    //     u[j][0] = 0.1;

    //     // Outflow condition
    //     u[j][nx] = u[j][nx-1];
    // }
}

int main(int argc, char **argv) {
    // Testing the fluid solver
    double tEnd = 1; 

    // The boundaries
    double xa = 0, xb = 1;
    double ya = 0, yb = 1;

    // Number of x, y points
    // int nx = 128;
    // int ny = 128;
    int nx = 20;
    int ny = 20;

    /* Creation of the solid objects */
    ///////////////////////////////////

    // Parameters of this circle
    SolidParams circParams1;
    SolidParams circParams2;

    double mass = 1.0;
    double density = 1.0;

    double E = 10;

    // circParams1.addParam("cx", 0.25);
    circParams1.addParam("cx", 0.5);
    circParams1.addParam("cy", 0.5);
    circParams1.addParam("r", 0.15);
    circParams1.addParam("mass", mass);
    circParams1.addParam("density", density);
    circParams1.addParam("E", E);
    circParams1.addParam("eta", 0.0);

    circParams2.addParam("cx", 0.75);
    circParams2.addParam("cy", 0.5);
    circParams2.addParam("r", 0.15);
    circParams2.addParam("mass", mass);
    circParams2.addParam("density", density);
    circParams2.addParam("E", E);
    circParams2.addParam("eta", 0.0);

    double u0 = 0.0;
    double v0 = 0.0;

    bool deformableBody = true;

    // Boundary object
    Boundary boundary(xa, xb, ya, yb);

    // ObjectSeeder seeder;
    // int nStructs = 5;
    // double r = 0.10;
    // double epsLoc = 2*sqrt(simutils::square(1/((double)nx)) + simutils::square(1/((double)ny)));
    // std::vector<SolidObject> shapes(seeder.randomInitialize(nStructs, deformableBody,
    //         r, epsLoc, boundary, circleShapeFun, circParams1));
    SolidObject::ObjectType objType = SolidObject::ObjectType::DEFORMABLE;
    SolidObject circle1(u0, v0, objType, coneShapeFun, circParams1);
    SolidObject circle2(u0, v0, objType, coneShapeFun, circParams2);
    // SolidObject circle1(u0, v0, deformableBody, circleShapeFun, circParams1);
    // SolidObject circle2(u0, v0, deformableBody, circleShapeFun, circParams2);

    // Create circle object array for input
    // SolidObject shapes[2] = {circle1, circle2};
    std::vector<SolidObject> shapes;
    shapes.push_back(circle1);
    // shapes.push_back(circle2);
    ///////////////////////////////////

    /* Create the solver object with appropriate parameters for the fluid and domain */
    ///////////////////////////////////////////////////////////////////////////////////

    SimParams params;
    params.setRe(1000);
    params.setNx(nx);
    params.setNy(ny);
    params.setUseEno(true);
    params.setMu(1.0/params.Re);
    params.setRepulseMode(2); // This turns on the KD tree error checking
    // simParams.setRepulseDist(5*sqrt(simutils::square(1/((double)nx)) + simutils::square(1/((double)ny))) );
    // params.setRepulseDist(0.1); // Actually need 0.1
    // params.setCollisionStiffness(2.0);
    // params.setCollisionDist(0.25);
    double h = sqrt(simutils::square(1.0/((double) nx)
        + simutils::square(1.0/((double) ny))));
    params.setRepulseDist(3*h); // Actually need 0.1
    params.setCollisionStiffness(2.0);
    params.setCollisionDist(3*h);
    params.setUpdateMode(1);
    double dt = 0.5/((double)nx) + 0.5/((double)ny);
    // params.setDtFix(dt);

    // Create the Solver object
    NSSolver solver(boundary, shapes, params, initialConditions, boundaryConditions);

    ///////////////////////////////////////////////////////////////////////////////////
    // Current time
    double t = 0;
    double safetyFactor = 0.5;

    int nsteps = 0;
    int max_steps = (argc == 1) ? 1 : atoi(argv[1]);
    while (t + EPS < tEnd && nsteps < max_steps) {
        t = solver.step(tEnd, safetyFactor);

        nsteps++;

        std::cout << "t = " << t << std::endl;
        std::cout << "step = " << nsteps << std::endl;
        std::cout << std::endl;
    }
    
    std::cout << "Outputting data" << std::endl;

    /* Output all of the relevant data */

    TestFormatter testFormatter("SimpleTest2D");
    string outStr;
    string outStr2;

    testFormatter.genOutStr("out", outStr);
    solver.writeToFile(outStr.c_str());
    cout << outStr << endl;

    testFormatter.genOutStr("poolOut", outStr);
    testFormatter.genOutStr("poolVel", outStr2);
    solver.writePoolToFile(outStr.c_str(), outStr2.c_str());
    cout << outStr << endl;
    cout << outStr2 << endl;

    testFormatter.genOutStr("MSSEdges", outStr);
    solver.outputAllStructures(outStr.c_str());
    cout << outStr << endl;

    testFormatter.genOutStr("MSSNodes", outStr);
    solver.outputAllStructureNodes(outStr.c_str());
    cout << outStr << endl;

    testFormatter.genOutStr("MSSVels", outStr);
    solver.outputAllStructureVels(outStr.c_str());
    cout << outStr << endl;

    testFormatter.genOutStr("MSSTracers", outStr);
    solver.outputTracers(outStr.c_str());

    testFormatter.genOutStr("medialAxis", outStr);
    solver.outputMedialAxis(outStr.c_str());
}
