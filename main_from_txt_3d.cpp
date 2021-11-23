#include <cassert>
#include <fenv.h>
#include <fstream>
#include <iostream>
#include <math.h>
#include <chrono>

#include "./src/Utils/SimUtilities.h"
#include "./src/Utils/Discretizations.h"
#include "./src/Utils/TestFormatter.h"
#include "./src/3DSolver/SimParams3D.h"
#include "./src/3DSolver/Boundary3D.h"
#include "./src/3DSolver/SolidObject3D.h"
#include "./src/3DSolver/NSSolver3D.h"


using namespace std::chrono;
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
    double cx, cy, cz, a, c, r, deg;
    ps.getParam("cx", cx);
    ps.getParam("cy", cy);
    ps.getParam("cz", cz);
    ps.getParam("a", a);
    ps.getParam("c", c);
    ps.getParam("r", r);
    ps.getParam("deg", deg);

    double b = 2.25 * r;
    double rad = deg * M_PI / 180;
    double rotcy = (y-cy) / (b) * cos(rad) - (z-cz) / (b) * sin(rad);
    double rotcz = (y-cy) / (b) * sin(rad) + (z-cz) / (b) * cos(rad);

    double x_sqr = simutils::square((x-cx)/b);
    double y_sqr = simutils::square(rotcy);
    double z_sqr = simutils::square(rotcz);
    double a_sqr = simutils::square(a);
    double c_sqr = simutils::square(c);

    return simutils::square(x_sqr + y_sqr + z_sqr + a_sqr) - 4*a_sqr*(x_sqr + y_sqr) - c_sqr;
}

typedef std::function<void (int,int,int,int,double*,double*,double*,double***,double***,double***)> initialConditions3DFunType;
initialConditions3DFunType getInitialConditionsFun(double cons_u, double cons_v, double cons_w) {
    // return lambda function with same parameters for NSSolver but user-given velocities
    return [cons_u, cons_v, cons_w](int nx, int ny, int nz, int nGhost, 
                            double *x, double *y, double *z, double ***u, double ***v, double ***w) {

        int wg = nx + 2*nGhost; // "width"
        int hg = ny + 2*nGhost; // "height"
        int vg = nz + 2*nGhost; // "vert"

        cout << "consu = " << cons_u << endl;
        cout << "consv = " << cons_v << endl;
        cout << "consw = " << cons_w << endl;

        simutils::set_constant(vg, hg, wg-1, cons_u, u);
        simutils::set_constant(vg, hg-1, wg, cons_v, v);
        simutils::set_constant(vg-1, hg, wg, cons_w, w); 
    };
}

void lidDrivenCavityBC(int nx, int ny, int nz, double ***u, double ***v, double ***w) {
    double ubar = 1.0;
    for (int j = 1; j <= ny; j++) {
        for (int i = 1; i <= nx; i++) {
            u[nz+1][j][i] = 2.0*ubar - u[nz][j][i];
        }
    }
}

void directionalFlowBC(int nx, int ny, int nz, double ***u, double ***v, double ***w) {
    cout << "In directional flow BC" << endl;
    for (int k = 1; k <= nz; k++) {
        for (int j = 1; j <= ny; j++) {
            // Inflow condition
            u[k][j][0] = 0.1; //simutils::dmin(t, 1.0)*((-6*simutils::square(y[j-1]) + 6*y[j-1])) + simutils::dmax(1.0 - t, 0);

            // Outflow condition
            u[k][j][nx] = u[k][j][nx-1];
        }
    }
    cout << "OUT directional flow BC" << endl;
}
// for (int i = 1; i <= nx; i++) {
//         // Inflow condition
//         v[0][i] = v[1][i];
//         //simutils::dmin(t, 1.0)*((-6*simutils::square(y[j-1]) + 6*y[j-1])) + simutils::dmax(1.0 - t, 0);

//         // Outflow condition
//         v[ny][i] = -0.1;
//     }

void downDirFlowBC(int nx, int ny, int nz, double ***u, double ***v, double ***w) {
    for (int j = 1; j <= ny; j++) {
        for (int i = 1; i <= nx; i++) {
            w[0][j][i] = w[1][j][i];

            w[nz][j][i] = -0.1;
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

    testFormatter.genOutStr("MSS3DTracers", outStr);
    solver.outputTracers(outStr.c_str());
    cout << outStr << endl;

    testFormatter.genOutStr("medialAxis3D", outStr);
    solver.outputMedialAxis(outStr.c_str());
    cout << outStr << endl;
}


void outputData(string f_name, NSSolver3D &solver, int testNum) {
    TestFormatter testFormatter(f_name.c_str(), testNum);
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

    testFormatter.genOutStr("MSS3DTracers", outStr);
    solver.outputTracers(outStr.c_str());
    cout << outStr << endl;

    testFormatter.genOutStr("medialAxis3D", outStr);
    solver.outputMedialAxis(outStr.c_str());
    cout << outStr << endl;
}

// argv: main input_file max_steps save_snapshots
int main(int argc, char **argv) {
    // set up dictionary of functions for input file
    map<string, double (*)(double, double, double, SolidParams&)> shapeFunctions;
    shapeFunctions["circleShapeFun"] = circleShapeFun;
    shapeFunctions["coneShapeFun"] = coneShapeFun;
    shapeFunctions["bloodCellShapeFun"] = bloodCellShapeFun;

    map<string, void (*)(int, int, int, double***, double***, double***)> boundaryConditionFunctions;
    boundaryConditionFunctions["lidDrivenCavityBC"] = lidDrivenCavityBC;
    boundaryConditionFunctions["directionalFlowBC"] = directionalFlowBC;
    boundaryConditionFunctions["downDirFlowBC"] = downDirFlowBC;

    if (argc < 2) {
      std::cout << "need more args to run this one" << endl;
      return 0;
    }

    string testInputDir = "./TestDrivers/3DDrivers/";
    string testTitle = argv[1];
    string input_file_name =  testInputDir + testTitle;
    int max_steps = (argc == 2) ? 1 : atoi(argv[2]);
    bool save_snapshots = (argc <= 3) ? false : strcmp(argv[3], "1") == 0;
    // bool save_snapshots = 0; // not implemented in 3D yet

    /* Parse from input file */
    ///////////////////////////////////

    string desc, boundaryConditionType;
    string testName;
    double xa, xb, ya, yb, za, zb, cons_u, cons_v, cons_w, g_x, g_y, g_z, tEnd;
    double re;
    int nx, ny, nz, num_objects;
    bool useEno;
    vector<SolidObject3D> shapes;
    SimParams3D simParams;

    ifstream input_file(input_file_name);

    // ignore first line, which is a description
    getline(input_file, desc);

    // Get the name of the test
    getline(input_file, testName);

    cout << "tsetName = " << testName << endl;
    // assert(false);

    // get simulation params
    input_file >> xa >> xb >> ya >> yb >> za >> zb >> nx >> ny >> nz;
    input_file >> cons_u >> cons_v >> cons_w >> g_x >> g_y >> g_z;
    input_file >> re >> useEno >> boundaryConditionType >> tEnd;

    cout << "Reading in everything" << endl;
    cout << "xa = " << xa << endl;
    cout << "xb = " << xb << endl;
    cout << "ya = " << ya << endl;
    cout << "yb = " << yb << endl;
    cout << "za = " << za << endl;
    cout << "zb = " << zb << endl;

    cout << "nx = " << nx << endl;
    cout << "ny = " << ny << endl;
    cout << "nz = " << nz << endl;

    cout << "cons_u = " << cons_u << endl;
    cout << "cons_v = " << cons_v << endl;
    cout << "cons_w = " << cons_w << endl;

    cout << "g_x = " << g_x << endl;
    cout << "g_y = " << g_y << endl;
    cout << "g_z = " << g_z << endl;

    cout << "Re = " << re << endl;

    cout << "useEno = " << useEno << endl;

    cout << "bcType = " << boundaryConditionType << endl;

    // get objects in simulation
    input_file >> num_objects;

    cout << "numObjects = " << num_objects << endl;

    string objectFunc, paramName;
    int objectType;
    double u0, v0, w0, paramValue;
    for (int i = 0; i < num_objects; i++) {
        input_file >> objectFunc >> objectType >> u0 >> v0 >> w0;

        cout << "objectFunc = " << objectFunc << endl;
        cout << "objectType = " << objectType << endl;
        cout << "objvel = " << u0 << ", " << v0 << "< " << w0 << endl;
        
        SolidParams params;
        input_file >> paramName;
        while (paramName != ".") {
            input_file >> paramValue;
            cout << "Adding paramName = " << paramName << " paramValue = " << paramValue << endl;
            params.addParam(paramName, paramValue);

            input_file >> paramName;
        }

        SolidObject3D object(u0, v0, w0, (SolidObject3D::ObjectType)objectType, shapeFunctions[objectFunc], params);
        shapes.push_back(object);
    }

    // actually set up simParams
    double h = sqrt(simutils::square(abs(xa - xb)/((double) nx)
        + simutils::square(abs(ya - yb)/((double) ny)
        + simutils::square(abs(za - zb)/((double) nz)))));

    simParams.setRe(re);
    simParams.setMu(1.0/simParams.Re);
    simParams.setNx(nx);
    simParams.setNy(ny);
    simParams.setNz(nz);
    simParams.setMssNx(nx);
    simParams.setMssNy(ny);
    simParams.setMssNz(nz);
    simParams.setUseEno(useEno);
    simParams.setRepulseMode(2);
    // simParams.setRepulseDist(4*h); // Actually need 0.1
    simParams.setCollisionStiffness(0.05);
    double collisionDist = 2*h;
    simParams.setCollisionDist(collisionDist);
    cout << "collision dist in simulation is " << collisionDist << endl;
    simParams.setUpdateMode(1);
    simParams.setGx(g_x);
    simParams.setGy(g_y);
    simParams.setGz(g_z);
    simParams.setU0(cons_u);
    simParams.setV0(cons_v);
    simParams.setW0(cons_w);
    // simParams.setAdmmTol(1e-10);
    // double dtFix = 0.1*(1/2.0)*(simutils::square(h));
    // simParams.setDtFix(dtFix);

    // initial/boundary conditions and boundary object
    auto initialConditions = getInitialConditionsFun(cons_u, cons_v, cons_w);
    auto boundaryCondition = boundaryConditionFunctions[boundaryConditionType];
    Boundary3D boundary(xa, xb, ya, yb, za, zb);

    // Create the Solver object
    cout << "Making solver" << endl;
    NSSolver3D solver(boundary, shapes, simParams, initialConditions, boundaryCondition);
    cout << "Finished Making solver" << endl;

    ///////////////////////////////////////////////////////////////////////////////////
    // Current time
    double t = 0;
    double safetyFactor = 0.15;

    // assert(false); // Think there is an issue with the boundary conditions for the obstacle domain
    auto start = high_resolution_clock::now();

    int nsteps = 0;
    if (save_snapshots) {

        while (t+EPS < tEnd && nsteps < max_steps) {
            if (nsteps % 50 == 0) {
                string f_name = testName;
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
            // assert(false);
        }

        string f_name = testName;
        outputData(f_name, solver);

    } else {
        
        std::cout << "TAKING " << max_steps << " STEPS" << std::endl;
        while (t+EPS < tEnd && nsteps < max_steps) {
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

        /* Output all of the relevant data */
        std::cout << "Outputting data" << std::endl;
        outputData(testName, solver);

    }
}
