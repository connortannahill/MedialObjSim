#include <iostream>
#include "./src/Utils/SimUtilities.h"
#include "./src/2DSolver/Boundary.h"
#include "./src/2DSolver/SolidObject.h"
#include "./src/2DSolver/NSSolver.h"
#include <math.h>
#include "./src/Utils/Discretizations.h"
#include <cassert>
#include <fenv.h>
#include <fstream>
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
    //     u[j][0] = 0.1; //simutils::dmin(t, 1.0)*((-6*simutils::square(y[j-1]) + 6*y[j-1])) + simutils::dmax(1.0 - t, 0);

    //     // Outflow condition
    //     u[j][nx] = u[j][nx-1];
    // }
}

// argv: main input_file max_steps
int main(int argc, char **argv) {
    map<string, double (*)(double, double, SolidParams&)> shapeFunctions;
    shapeFunctions["circleShapeFun"] = circleShapeFun;
    shapeFunctions["coneShapeFun"] = coneShapeFun;

    if (argc < 2) {
      std::cout << "need more args to run this one" << endl;
      return 0;
    }

    string input_file_name = argv[1];
    int max_steps = (argc == 2) ? 1 : atoi(argv[2]);

    /* Parse from input file */
    ///////////////////////////////////

    double xa, xb, ya, yb, tEnd;
    int nx, ny, re;
    int num_objects;
    vector<SolidObject> shapes;
    SimParams simParams;

    ifstream input_file(input_file_name);

    // get simulation params
    input_file >> xa >> xb >> ya >> yb;
    input_file >> nx >> ny >> re >> tEnd;

    // std::cout << xa << " " << xb << " " << ya << " " << yb << endl;
    // std::cout << nx << " " << ny << " " << re << " " << tEnd << endl;

    // get objects in simulation
    input_file >> num_objects;

    double cx, cy, r, mass, density, E, eta, u0, v0;
    int objectType;
    string objectFunc;
    for (int i = 0; i < num_objects; i++) {
      input_file >> cx >> cy >> r >> mass >> density >> E >> eta;
      input_file >> u0 >> v0 >> objectType >> objectFunc;
      // std::cout << cx << " " << cy << " " << r << " " << mass << " " << density << " " << E << " " << eta <<  endl;
      // std::cout << u0 << " " << v0 << " " << objectType << " " << objectFunc << endl;

      SolidParams params;
      params.addParam("cx", cx);
      params.addParam("cy", cy);
      params.addParam("r", r);
      params.addParam("mass", mass);
      params.addParam("density", density);
      params.addParam("E", E);
      params.addParam("eta", eta);

      SolidObject object(u0, v0, (SolidObject::ObjectType)objectType, shapeFunctions[objectFunc], params);
      shapes.push_back(object);
    }

    // actually set up simParams
    double h = sqrt(simutils::square(1.0/((double) nx)
        + simutils::square(1.0/((double) ny))));
    double dt = 0.5/((double)nx) + 0.5/((double)ny); // TODO: compute this more generally, perhaps make the time step computation method static.

    simParams.setRe(re);
    simParams.setNx(nx);
    simParams.setNy(ny);
    simParams.setMu(1.0/simParams.Re);
    simParams.setRepulseMode(2); // This turns on the KD tree error checking
    // simParams.setRepulseDist(5*sqrt(simutils::square(1/((double)nx)) + simutils::square(1/((double)ny))) );
    // simParams.setRepulseDist(0.1); // Actually need 0.1
    // simParams.setCollisionStiffness(2.0);
    // simParams.setCollisionDist(0.25);
    simParams.setRepulseDist(3*h); // Actually need 0.1
    simParams.setCollisionStiffness(2.0);
    simParams.setCollisionDist(3*h);
    simParams.setUpdateMode(2);
    simParams.setDtFix(dt);

    // Boundary object
    Boundary boundary(xa, xb, ya, yb);

    // Create the Solver object
    NSSolver solver(boundary, shapes, simParams, initialConditions, boundaryConditions);

    ///////////////////////////////////////////////////////////////////////////////////
    // Current time
    double t = 0;
    double safetyFactor = 1;

    // assert(false); // Think there is an issue with the boundary conditions for the obstacle domain

    int nsteps = 0;
    while (t+EPS < tEnd && nsteps < max_steps) {
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
    std::cout << outStr << endl;

    testFormatter.genOutStr("poolOut", outStr);
    testFormatter.genOutStr("poolVel", outStr2);
    solver.writePoolToFile(outStr.c_str(), outStr2.c_str());
    std::cout << outStr << endl;
    std::cout << outStr2 << endl;

    testFormatter.genOutStr("MSSEdges", outStr);
    solver.outputAllStructures(outStr.c_str());
    std::cout << outStr << endl;

    testFormatter.genOutStr("MSSNodes", outStr);
    solver.outputAllStructureNodes(outStr.c_str());
    std::cout << outStr << endl;

    testFormatter.genOutStr("MSSVels", outStr);
    solver.outputAllStructureVels(outStr.c_str());
    std::cout << outStr << endl;

    testFormatter.genOutStr("MSSTracers", outStr);
    solver.outputTracers(outStr.c_str());

    testFormatter.genOutStr("medialAxis", outStr);
    solver.outputMedialAxis(outStr.c_str());
}
