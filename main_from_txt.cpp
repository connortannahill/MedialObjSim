#include <cassert>
#include <fenv.h>
#include <fstream>
#include <iostream>
#include <math.h>

#include "./src/2DSolver/Boundary.h"
#include "./src/2DSolver/SolidObject.h"
#include "./src/2DSolver/NSSolver.h"
#include "./src/2DSolver/ObjectSeeder.h"
#include "./src/2DSolver/SimParams.h"
#include "./src/Utils/Discretizations.h"
#include "./src/Utils/SimUtilities.h"
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

// based on Cassini oval
double bloodCellShapeFun(double x, double y, SolidParams &ps) {
    double cx, cy, a, c, deg;
    ps.getParam("cx", cx);
    ps.getParam("cy", cy);
    ps.getParam("a", a);
    ps.getParam("c", c);
    ps.getParam("deg", deg); // degree of rotation

    double rad = deg * M_PI / 180;
    double rotcx = (x-cx) * cos(rad) - (y-cy) * sin(rad);
    double rotcy = (x-cx) * sin(rad) + (y-cy) * cos(rad);

    double x_sqr = simutils::square(rotcx);
    double y_sqr = simutils::square(rotcy);
    double a_sqr = simutils::square(a);
    double c_sqr = simutils::square(c);

    return simutils::square(x_sqr + y_sqr + a_sqr) - 4*a_sqr*x_sqr - simutils::square(c_sqr);
}

// kinda complicated, should make better
typedef std::function<void (int,int,int,double*,double*,double**,double**)> initialConditionsFunType;
initialConditionsFunType getInitialConditionsFun(double cons_u, double cons_v) {
    // return lambda function with same parameters for NSSolver but user-given velocities
    return [cons_u, cons_v](int nx, int ny, int nGhost, 
                            double *x, double *y, double **u, double **v) {

        int wg = nx + 2*nGhost; // "width"
        int hg = ny + 2*nGhost; // "height"

        simutils::set_constant(hg, wg-1, cons_u, u);
        simutils::set_constant(hg-1, wg, cons_v, v);
    };
}

// Problem-dependent boundary conditions
void lidDrivenCavityBC(int nx, int ny, double **u, double **v) {
    double ubar = 1;
    for (int i = 1; i <= nx; i++) {
        u[ny+1][i] = 2*ubar - u[ny][i];
    }
}

void directionalFlowBC(int nx, int ny, double **u, double **v) {
    for (int j = 1; j <= ny; j++) {
        // Inflow condition
        u[j][0] = 0.1;
        //simutils::dmin(t, 1.0)*((-6*simutils::square(y[j-1]) + 6*y[j-1])) + simutils::dmax(1.0 - t, 0);

        // Outflow condition
        u[j][nx] = u[j][nx-1];
    }
}

void outputData(string f_name, NSSolver &solver) {
    TestFormatter testFormatter(f_name.c_str());
    std::cout << "Outputting data" << std::endl;

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

// argv: main input_file max_steps save_snapshots
int main(int argc, char **argv) {
    // set up dictionary of functions for input file
    map<string, double (*)(double, double, SolidParams&)> shapeFunctions;
    shapeFunctions["circleShapeFun"] = circleShapeFun;
    shapeFunctions["coneShapeFun"] = coneShapeFun;
    shapeFunctions["bloodCellShapeFun"] = bloodCellShapeFun;

    map<string, void (*)(int, int, double**, double**)> boundaryConditionFunctions;
    boundaryConditionFunctions["lidDrivenCavityBC"] = lidDrivenCavityBC;
    boundaryConditionFunctions["directionalFlowBC"] = directionalFlowBC;

    if (argc < 2) {
      std::cout << "need more args to run this one" << endl;
      return 0;
    }

    string input_file_name = argv[1];
    int max_steps = (argc == 2) ? 1 : atoi(argv[2]);
    bool save_snapshots = (argc <= 3) ? false : strcmp(argv[3], "1") == 0;

    /* Parse from input file */
    ///////////////////////////////////

    string desc, boundaryConditionType;
    string outFileName;
    double xa, xb, ya, yb, cons_u, cons_v, tEnd;
    int nx, ny, re, num_objects;
    bool useEno;
    vector<SolidObject> shapes;
    SimParams simParams;

    ifstream input_file(input_file_name);

    // ignore first line, which is a description
    getline(input_file, desc);

    // Get the output file name
    getline(input_file, outFileName);

    // get simulation params
    input_file >> xa >> xb >> ya >> yb >> nx >> ny >> cons_u >> cons_v;
    input_file >> re >> useEno >> boundaryConditionType >> tEnd;

    // std::cout << xa << " " << xb << " " << ya << " " << yb << endl;
    // std::cout << nx << " " << ny << " " << re << " " << tEnd << endl;

    // get objects in simulation
    input_file >> num_objects;

    string objectFunc, paramName;
    int objectType;
    double u0, v0, paramValue;
    for (int i = 0; i < num_objects; i++) {
        input_file >> objectFunc >> objectType >> u0 >> v0;
        
        SolidParams params;
        input_file >> paramName;
        while (paramName != ".") {
            input_file >> paramValue;
            params.addParam(paramName, paramValue);
            // cout << paramName << " " << paramValue << endl;

            input_file >> paramName;
        }
        // std::cout << u0 << " " << v0 << " " << objectType << " " << objectFunc << endl;

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
    simParams.setMssNx(nx);
    simParams.setMssNy(ny);
    simParams.setUseEno(useEno);
    simParams.setMu(1.0/simParams.Re);
    simParams.setRepulseMode(2); // This turns on the KD tree error checking
    // simParams.setRepulseDist(5*sqrt(simutils::square(1/((double)nx)) + simutils::square(1/((double)ny))) );
    // simParams.setRepulseDist(0.1); // Actually need 0.1
    // simParams.setCollisionStiffness(2.0);
    // simParams.setCollisionDist(0.25);
    simParams.setRepulseDist(0.1); // Actually need 0.1
    simParams.setCollisionStiffness(5.0);
    simParams.setCollisionDist(4*h);
    simParams.setUpdateMode(1);
    simParams.setAdmmTol(1e-10);
    simParams.setGx(0.0);
    simParams.setGy(0.0);
    // simParams.setDtFix(dt);

    // initial/boundary conditions and boundary object
    auto initialConditions = getInitialConditionsFun(cons_u, cons_v);
    auto boundaryCondition = boundaryConditionFunctions[boundaryConditionType];
    Boundary boundary(xa, xb, ya, yb);

    // Create the Solver object
    NSSolver solver(boundary, shapes, simParams, initialConditions, boundaryCondition);

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
                string f_name = outFileName + std::to_string(nsteps);
                outputData(f_name, solver);
            }

            nsteps++;

            std::cout << "t = " << t << std::endl;
            std::cout << "step = " << nsteps << std::endl;
            std::cout << std::endl;
        }

        string f_name = outFileName + std::to_string(nsteps);
        outputData(f_name, solver);

    } else {

        while (t+EPS < tEnd && nsteps < max_steps) {
            t = solver.step(tEnd, safetyFactor);

            nsteps++;

            std::cout << "t = " << t << std::endl;
            std::cout << "step = " << nsteps << std::endl;
            std::cout << std::endl;
        }

        /* Output all of the relevant data */
        std::cout << "Outputting data" << std::endl;
        outputData(outFileName, solver);

    }
}
