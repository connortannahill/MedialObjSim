#include <iostream>
#include <fstream>
#include "./src/3DSolver/Boundary3D.h"
#include "./src/3DSolver/SolidObject3D.h"
#include "./src/3DSolver/NSSolver3D.h"

#include "./src/Utils/SimUtilities.h"
#include "./src/Utils/Discretizations.h"
#include "./src/Utils/TestFormatter.h"
#include "./src/3DSolver/SimParams3D.h"

#include <math.h>
#include <cassert>

using namespace std;
const double EPS = 1e-12;

double circleShapeFun(double x, double y, double z, SolidParams &ps) {
    double cx, cy, cz, r;
    ps.getParam("cx", cx);
    ps.getParam("cy", cy);
    ps.getParam("cz", cz);
    ps.getParam("r", r);
    return simutils::square(x-cx)+simutils::square(y-cy)+simutils::square(z-cz)-simutils::square(r);
}

double coneShapeFun(double x, double y, double z, SolidParams &ps) {
    double cx, cy, cz, r;
    ps.getParam("cx", cx);
    ps.getParam("cy", cy);
    ps.getParam("cz", cz);
    ps.getParam("r", r);
    return sqrt(simutils::square(x-cx)+simutils::square(y-cy)+simutils::square(z-cz))-r;
}

// based on Cassini oval
double bloodCellShapeFun(double x, double y, double z, SolidParams &ps) {
    double cx, cy, cz, a, c, r;
    ps.getParam("cx", cx);
    ps.getParam("cy", cy);
    ps.getParam("cz", cz);
    ps.getParam("a", a);
    ps.getParam("c", c);
    ps.getParam("r", r);

    double b = 1.75 * r;

    double x_sqr = simutils::square((x-cx)/b);
    double y_sqr = simutils::square((y-cy)/b);
    double z_sqr = simutils::square((z-cz)/b);
    double a_sqr = simutils::square(a);
    double c_sqr = simutils::square(c);

    return simutils::square(x_sqr + y_sqr + z_sqr + a_sqr) - 4*a_sqr*(x_sqr + y_sqr) - c_sqr;
}

void initialConditions(int nx, int ny, int nz, int nGhost, double *x,
                        double *y, double *z, double ***u, double ***v,
                        double ***w) {
    double cons_u = 0.0;
    double cons_v = 0.0;
    double cons_w = 0.0;

    int wg = nx + 2*nGhost; // "width"
    int hg = ny + 2*nGhost; // "height"
    int vg = nz + 2*nGhost; // "vert"

    simutils::set_constant(vg, hg, wg-1, cons_u, u);
    simutils::set_constant(vg, hg-1, wg, cons_v, v);
    simutils::set_constant(vg-1, hg, wg, cons_w, w); 
}

void lidDrivenCavityBC(int nx, int ny, int nz, double ***u) {
    double ubar = 1.0;
    for (int j = 1; j <= ny; j++) {
        for (int i = 1; i <= nx; i++) {
            u[nz+1][j][i] = 2.0*ubar - u[nz][j][i];
        }
    }
}

void directionalFlowBC(int nx, int ny, int nz, double ***u) {
    for (int k = 1; k <= nz; k++) {
        for (int j = 1; j <= ny; j++) {
            // Inflow condition
            u[k][j][0] = 0.1; //simutils::dmin(t, 1.0)*((-6*simutils::square(y[j-1]) + 6*y[j-1])) + simutils::dmax(1.0 - t, 0);

            // Outflow condition
            u[k][j][nx] = u[k][j][nx-1];
        }
    }
}

void outputData(string f_name, NSSolver3D &solver) {
    TestFormatter testFormatter(f_name.c_str());
    std::cout << "Outputting data" << std::endl;

    string outStr;
    string outStr2;

    testFormatter.genOutStr("out3D", outStr);
    solver.writeToFile(outStr.c_str());
    cout << outStr << endl;

    testFormatter.genOutStr("pool3DOut", outStr);
    testFormatter.genOutStr("pool3DVel", outStr2);
    solver.writePoolToFile(outStr.c_str(), outStr2.c_str());
    cout << outStr << endl;
    cout << outStr2 << endl;

    testFormatter.genOutStr("MSS3DEdges", outStr);
    solver.outputAllStructures(outStr.c_str());
    cout << outStr << endl;

    testFormatter.genOutStr("MSS3DNodes", outStr);
    solver.outputAllStructureNodes(outStr.c_str());
    cout << outStr << endl;

    testFormatter.genOutStr("MSS3DVels", outStr);
    solver.outputAllStructureVels(outStr.c_str());
    cout << outStr << endl;

    // testFormatter.genOutStr("MSSTracers", outStr);
    // solver.outputTracers(outStr.c_str());

    // testFormatter.genOutStr("medialAxis", outStr);
    // solver.outputMedialAxis(outStr.c_str());
}

// argv: main input_file max_steps save_snapshots
int main(int argc, char **argv) {
    // set up dictionary of functions for input file
    map<string, double (*)(double, double, double, SolidParams&)> shapeFunctions;
    shapeFunctions["circleShapeFun"] = circleShapeFun;
    shapeFunctions["coneShapeFun"] = coneShapeFun;
    shapeFunctions["bloodCellShapeFun"] = bloodCellShapeFun;

    map<string, void (*)(int, int, int, double***)> boundaryConditionFunctions;
    boundaryConditionFunctions["lidDrivenCavityBC"] = lidDrivenCavityBC;
    boundaryConditionFunctions["directionalFlowBC"] = directionalFlowBC;

    if (argc < 2) {
      std::cout << "need more args to run this one" << endl;
      return 0;
    }

    string input_file_name = argv[1];
    int max_steps = (argc == 2) ? 1 : atoi(argv[2]);
    // bool save_snapshots = (argc <= 3) ? false : strcmp(argv[3], "1") == 0;
    bool save_snapshots = 0; // not implemented in 3D yet

    /* Parse from input file */
    ///////////////////////////////////

    string desc, boundaryConditionType;
    string testName;
    double xa, xb, ya, yb, za, zb, tEnd;
    int nx, ny, nz, re, num_objects;
    bool useEno;
    vector<SolidObject3D> shapes;
    SimParams3D simParams;

    ifstream input_file(input_file_name);

    // ignore first line, which is a description
    getline(input_file, desc);

    // Get the name of the test
    getline(input_file, testName);

    // get simulation params
    input_file >> xa >> xb >> ya >> yb >> za >> zb >> nx >> ny >> nz;
    input_file >> re >> useEno >> boundaryConditionType >> tEnd;

    // std::cout << xa << " " << xb << " " << ya << " " << yb << endl;
    // std::cout << nx << " " << ny << " " << re << " " << tEnd << endl;

    // get objects in simulation
    input_file >> num_objects;

    string objectFunc, paramName;
    int objectType;
    double u0, v0, w0, paramValue;
    for (int i = 0; i < num_objects; i++) {
        input_file >> objectFunc >> objectType >> u0 >> v0 >> w0;
        
        SolidParams params;
        input_file >> paramName;
        while (paramName != ".") {
            input_file >> paramValue;
            params.addParam(paramName, paramValue);
            // cout << paramName << " " << paramValue << endl;

            input_file >> paramName;
        }
        // std::cout << u0 << " " << v0 << " " << objectType << " " << objectFunc << endl;

        SolidObject3D object(u0, v0, w0, (SolidObject3D::ObjectType)objectType, shapeFunctions[objectFunc], params);
        shapes.push_back(object);
    }

    // actually set up simParams
    double h = sqrt(simutils::square(1.0/((double) nx)
        + simutils::square(1.0/((double) ny)
        + simutils::square(1.0/((double) nz)))));
    // double dt = 0.5/((double)nx) + 0.5/((double)ny); // TODO: compute this more generally, perhaps make the time step computation method static.

    simParams.setRe(re);
    simParams.setMu(1.0/simParams.Re);
    simParams.setNx(nx);
    simParams.setNy(ny);
    simParams.setNz(nz);
    simParams.setMssNx(nx);
    simParams.setMssNy(ny);
    simParams.setMssNz(nz);
    simParams.setUseEno(useEno);
    simParams.setRepulseMode(0);
    simParams.setRepulseDist(3*h); // Actually need 0.1
    simParams.setCollisionStiffness(2.0);
    simParams.setCollisionDist(3*h);
    simParams.setUpdateMode(1);

    // initial/boundary conditions and boundary object
    // auto initialConditions = getInitialConditionsFun(cons_u, cons_v);
    auto boundaryCondition = boundaryConditionFunctions[boundaryConditionType];
    Boundary3D boundary(xa, xb, ya, yb, za, zb);

    // Create the Solver object
    NSSolver3D solver(boundary, shapes, simParams, initialConditions, boundaryCondition);

    ///////////////////////////////////////////////////////////////////////////////////
    // Current time
    double t = 0;
    double safetyFactor = 0.5;

    // assert(false); // Think there is an issue with the boundary conditions for the obstacle domain

    int nsteps = 0;
    if (save_snapshots) {

        while (t+EPS < tEnd && nsteps < max_steps) {
            t = solver.step(tEnd, safetyFactor);

            if (nsteps % 20 == 0) {
                string f_name = testName + std::to_string(nsteps);
                outputData(f_name, solver);
            }

            nsteps++;

            std::cout << "t = " << t << std::endl;
            std::cout << "step = " << nsteps << std::endl;
            std::cout << std::endl;
        }

        string f_name = testName + std::to_string(nsteps);
        outputData(f_name, solver);

    } else {
        
        std::cout << "TAKING " << max_steps << " STEPS" << std::endl;
        while (t+EPS < tEnd && nsteps < max_steps) {
            t = solver.step(tEnd, safetyFactor);

            nsteps++;

            std::cout << "t = " << t << std::endl;
            std::cout << "step = " << nsteps << std::endl;
            std::cout << std::endl;
        }

        /* Output all of the relevant data */
        std::cout << "Outputting data" << std::endl;
        outputData(testName, solver);

    }
}
