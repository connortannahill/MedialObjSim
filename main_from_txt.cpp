#include <cassert>
#include <fenv.h>
#include <fstream>
#include <iostream>
#include <math.h>
#include <chrono>

#include "./src/2DSolver/Boundary.h"
#include "./src/2DSolver/SolidObject.h"
#include "./src/2DSolver/NSSolver.h"
#include "./src/2DSolver/ObjectSeeder.h"
#include "./src/2DSolver/SimParams.h"
#include "./src/Utils/Discretizations.h"
#include "./src/Utils/SimUtilities.h"
#include "./src/Utils/TestFormatter.h"

using namespace std::chrono;
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
    double cx, cy, a, c, r, deg;
    ps.getParam("cx", cx);
    ps.getParam("cy", cy);
    ps.getParam("a", a);
    ps.getParam("c", c);
    ps.getParam("r", r);
    ps.getParam("deg", deg); // degree of rotation

    double b = 2.25 * r;
    double rad = deg * M_PI / 180;
    double rotcx = (x-cx) / (b) * cos(rad) - (y-cy) / (b) * sin(rad);
    double rotcy = (x-cx) / (b) * sin(rad) + (y-cy) / (b) * cos(rad);

    double x_sqr = simutils::square(rotcx);
    double y_sqr = simutils::square(rotcy);
    double a_sqr = simutils::square(a);
    double c_sqr = simutils::square(c);

    return simutils::square(x_sqr + y_sqr + a_sqr) - 4*a_sqr*x_sqr - c_sqr;
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
    double ubar = 0.5;
    for (int i = 1; i <= nx; i++) {
        u[ny+1][i] = 2*ubar - u[ny][i];
    }
}

void directionalFlowBC(int nx, int ny, double **u, double **v) {
    for (int j = 1; j <= ny; j++) {
        // Inflow condition
        u[j][0] = 0.5;
        //simutils::dmin(t, 1.0)*((-6*simutils::square(y[j-1]) + 6*y[j-1])) + simutils::dmax(1.0 - t, 0);

        // Outflow condition
        u[j][nx] = u[j][nx-1];
    }
}

void downDirFlowBC(int nx, int ny, double **u, double **v) {
    // cout << "in bc" << endl;
    for (int i = 1; i <= nx; i++) {
        // Inflow condition
        v[0][i] = v[1][i];
        //simutils::dmin(t, 1.0)*((-6*simutils::square(y[j-1]) + 6*y[j-1])) + simutils::dmax(1.0 - t, 0);

        // Outflow condition
        v[ny][i] = -0.1;
    }
    // cout << "out bc" << endl;
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
    cout << outStr << endl;

    testFormatter.genOutStr("medialAxis", outStr);
    solver.outputMedialAxis(outStr.c_str());
    cout << outStr << endl;

    testFormatter.genOutStr("domain", outStr);
    solver.outputMedialAxis(outStr.c_str());
    cout << outStr << endl;
}

void outputData(string f_name, NSSolver &solver, int testNum) {
    TestFormatter testFormatter(f_name.c_str(), testNum);
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
    cout << outStr << endl;

    testFormatter.genOutStr("medialAxis", outStr);
    solver.outputMedialAxis(outStr.c_str());
    cout << outStr << endl;

    testFormatter.genOutStr("domain", outStr);
    solver.outputDomain(outStr.c_str());
    cout << outStr << endl;
}

double getH(double xa, double xb, int nx, double ya, double yb, int ny) {
    return sqrt(simutils::square(abs(xb - xa)/((double) nx)
        + simutils::square(abs(yb - ya)/((double) ny))));
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
    boundaryConditionFunctions["downDirFlowBC"] = downDirFlowBC;

    if (argc < 2) {
      std::cout << "need more args to run this one" << endl;
      return 0;
    }

    string testInputDir = "./TestDrivers/2DDrivers/";
    string testTitle = argv[1];
    cout << "TestTitle = " << testTitle << endl;
    string input_file_name =  testInputDir + testTitle;
    int max_steps = (argc == 2) ? 1 : atoi(argv[2]);
    bool save_snapshots = (argc <= 3) ? false : strcmp(argv[3], "1") == 0;

    /* Parse from input file */
    ///////////////////////////////////

    string desc, boundaryConditionType;
    string outFileName;
    double xa, xb, ya, yb, cons_u, cons_v, tEnd, re, g_x, g_y;
    int nx, ny, num_objects;
    bool useEno;
    vector<SolidObject> shapes;
    SimParams simParams;

    ifstream input_file(input_file_name);

    // ignore first line, which is a description
    getline(input_file, desc);

    // Get the output file name
    getline(input_file, outFileName);

    // get simulation params
    input_file >> xa >> xb >> ya >> yb >> nx >> ny >> cons_u >> cons_v >> g_x >> g_y;
    input_file >> re >> useEno >> boundaryConditionType >> tEnd;

    cout << "xa = " << xa << ", xb = " << xb << endl;
    cout << "ya = " << ya << ", yb = " << yb << endl;
    cout << "nx = " << nx << ", ny = " << ny << endl;
    cout << "u = " << cons_u << ", v = " << cons_v << endl;
    cout << "re = " << re << endl;
    cout << "useEno = " << useEno << endl;
    cout << "bctype = " << boundaryConditionType << endl;
    cout << "tEnd = " << tEnd << endl;

    // get objects in simulation
    input_file >> num_objects;
    cout << "num_objects " << num_objects << endl;

    string objectFunc, paramName;
    int objectType;
    double u0, v0, paramValue;
    for (int i = 0; i < num_objects; i++) {
        input_file >> objectFunc >> objectType >> u0 >> v0;
        cout << "objFunc " << objectFunc << endl;
        cout << "objType " << objectType << endl;
        cout << "u0 " << u0 << endl;
        double velAdd = 0.0; //(2*((double) rand() / (RAND_MAX))-1)/10.0;
        cout << "v0 " << (v0 + velAdd) << endl;
        
        SolidParams params;
        input_file >> paramName;
        while (paramName != ".") {
            input_file >> paramValue;
            params.addParam(paramName, paramValue);

            input_file >> paramName;
        }

        SolidObject object(u0, v0 + velAdd, (SolidObject::ObjectType)objectType, shapeFunctions[objectFunc], params);
        shapes.push_back(object);
    }

    cout << "Number of objects = " << shapes.size() << endl;

    // actually set up simParams
    double h = getH(xa, xb, nx, ya, yb, ny);

    simParams.setRe(re);
    simParams.setNx(nx);
    simParams.setNy(ny);
    simParams.setUseEno(useEno);
    simParams.setMu(1.0/simParams.Re);
    simParams.setRepulseMode(2); // This turns on the KD tree error checking
    simParams.setCollisionStiffness(20);
    simParams.setCollisionDist(4.0*h);
    cout << "collisionDist = " << 4.0*h << endl;
    int updateMode = 1;
    simParams.setUpdateMode(updateMode);
    cout << "UpdateMode = " << updateMode << endl;
    simParams.setAdmmTol(1e-10);
    simParams.setGx(g_x);
    simParams.setGy(g_y);
    simParams.setU0(cons_u);
    simParams.setV0(cons_v);
    // double umax = 0.5;
    // double vmax = 0.5;
    // double dx = abs((xb- xa)/nx);
    // double dy = abs((yb- ya)/ny);
    // double dtFix = (re/2.0)*(1/(simutils::square(1.0/dx) + simutils::square(1.0/dy)));
    // simParams.setDtFix(dtFix);

    // initial/boundary conditions and boundary object
    auto initialConditions = getInitialConditionsFun(cons_u, cons_v);
    auto boundaryCondition = boundaryConditionFunctions[boundaryConditionType];
    Boundary boundary(xa, xb, ya, yb);

    // Create the Solver object
    cout << "Creating the solver" << endl;
    NSSolver solver(boundary, shapes, simParams, initialConditions, boundaryCondition);

    ///////////////////////////////////////////////////////////////////////////////////
    // Current time
    double t = 0;
    double safetyFactor = 0.75;

    // Start recording the time. We do not include the seeding as this is technically
    // not a part of the algorithm per se.
    auto start = high_resolution_clock::now();
    cout << "hi" << endl;

    int nsteps = 0;
    if (save_snapshots) {

        while (t+EPS < tEnd && nsteps < max_steps) {
            cout << "hi 274" << endl;
            if (nsteps % 25 == 0) {
                string f_name = outFileName;
                outputData(f_name, solver, nsteps);
            }
            
            t = solver.step(tEnd, safetyFactor);


            nsteps++;

            std::cout << "t = " << t << std::endl;
            std::cout << "step = " << nsteps << std::endl;
            auto stop = high_resolution_clock::now();

            auto duration = duration_cast<seconds>(stop - start);
            std::cout << std::endl;

            cout << "=========================================================" << endl;
            cout << "testName = " << testTitle << endl;
            cout << "The total run time of the algorithm: " << duration.count() << endl;
            cout << "The average run time per step of the algorithm: " << duration.count()/((double)nsteps) << endl;
            cout << "=========================================================" << endl;
        }

        string f_name = outFileName + "/" + std::to_string(nsteps);
        outputData(f_name, solver);

    } else {

        while (t+EPS < tEnd && nsteps < max_steps) {
            t = solver.step(tEnd, safetyFactor);

            nsteps++;

            std::cout << "t = " << t << std::endl;
            std::cout << "step = " << nsteps << std::endl;
            std::cout << std::endl;

            auto stop = high_resolution_clock::now();

            auto duration = duration_cast<seconds>(stop - start);

            cout << "=========================================================" << endl;
            cout << "The total run time of the algorithm: " << duration.count() << endl;
            cout << "The average run time per step of the algorithm: " << duration.count()/((double)nsteps) << endl;
            cout << "=========================================================" << endl;
        }

        auto stop = high_resolution_clock::now();

        auto duration = duration_cast<seconds>(stop - start);

        cout << "=========================================================" << endl;
        cout << "The total run time of the algorithm: " << duration.count() << endl;
        cout << "The average run time per step of the algorithm: " << duration.count()/((double)nsteps) << endl;
        cout << "=========================================================" << endl;

        /* Output all of the relevant data */
        std::cout << "Outputting data" << std::endl;
        outputData(outFileName, solver);

    }
}
