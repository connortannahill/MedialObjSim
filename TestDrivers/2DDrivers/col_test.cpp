#include <iostream>
#include "../../src/2DSolver/Boundary.h"
#include "../../src/2DSolver/SolidObject.h"
#include "../../src/2DSolver/Pool2D.h"
#include "../../src/2DSolver/MassSpring2D.h"
#include "../../src/2DSolver/ObjectSeeder.h"
#include "../../src/2DSolver/SimParams.h"

#include "ctime"

#include "../../src/Utils/Discretizations.h"
#include "../../src/Utils/SimUtilities.h"
#include "../../src/Utils/SolidParams.h"
#include "../../src/Utils/TestFormatter.h"

#include <math.h>
#include <cassert>
#include <Eigen/Dense>

using namespace std;
using namespace mass_spring;

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

int main(int argc, char **argv) { 
    int nx = 20;
    int ny = 20;

    double xmin = 0;
    double xmax = 1;
    double ymin = 0;
    double ymax = 1;

    // Create the boundary object
    Boundary boundary(xmin, xmax, ymin, ymax);

    // Creating the pool objects
    SolidObject::ObjectType objType = SolidObject::ObjectType::DEFORMABLE;
    SolidParams params1;
    params1.addParam("cx",  0.25);
    params1.addParam("cy",  0.5);
    params1.addParam("r",  0.15);
    params1.addParam("mass",  1.0);
    params1.addParam("density",  1.0);
    params1.addParam("E",  1.0);
    params1.addParam("eta",  0.0);

    SolidParams params2(params1);
    params2.addParam("cx", 0.75);
    params2.addParam("cy", 0.5);

    double u0 = 0.0;
    double v0 = 0.0;
    SolidObject circle1(u0, v0, objType, coneShapeFun, params1);
    SolidObject circle2(-u0, v0, objType, coneShapeFun, params2);

    vector<SolidObject> shapes;
    shapes.push_back(circle1);
    shapes.push_back(circle2);

    // Create the simparams
    SimParams simParams;
    simParams.setNx(nx);
    simParams.setNy(ny);
    simParams.setMssNx(nx);
    simParams.setMssNy(ny);
    simParams.setRe(0.0);
    simParams.setMu(0.0);
    simParams.setRepulseMode(2); // This turns on the KD tree error checking
    // Compute the repulse distance required
    double h = sqrt(simutils::square(1.0/((double) simParams.nx)
        + simutils::square(1.0/((double) simParams.ny))));
    simParams.setRepulseDist(5.0*h); // Actually need 0.1
    // simParams.setRepulseDist(0.1); // Actually need 0.1
    simParams.setCollisionStiffness(2.0);
    // simParams.setCollisionDist(3*h);
    simParams.setCollisionDist(2.0*h);
    simParams.setAdmmTol(1e-10);
    simParams.setUpdateMode(1);

    double dt = 0.01;//0.5/(1.0/(1/((double)nx))+1.0/((1/(double)ny)));
    simParams.setDtFix(dt);

    // Create the pool
    Pool2D pool(boundary, shapes, simParams);

    double t = 0;
    double tEnd = 10;

    cout << "dt = " << dt << endl;

    // Test velocities for the level set methods
    double **u = simutils::new_constant(ny, nx, 0.0);
    double **v = simutils::new_constant(ny, nx, 0.0);
    double **p = simutils::new_constant(ny, nx, 0.0);

    int nstep = 0;
    int max_steps = (argc > 1) ? atoi(argv[1]) : 1;
    bool reinitialize = true;
    while (t < tEnd && nstep < max_steps) {
	cout << "Updating the pool" << endl;
        pool.updatePool(dt, u, v, p, 0, reinitialize); // Advance the level set function
	cout << "Finished the pool" << endl;
        t += dt;
        cout << "t = " << t << endl;
        nstep++;
    }

    /* Output all of the relevant data */

    TestFormatter testFormatter("ColTest2D");
    string outStr;

    testFormatter.genOutStr("poolOut", outStr);
    pool.outputPool(outStr.c_str());
    cout << outStr << endl;

    testFormatter.genOutStr("poolVel", outStr);
    pool.outputPoolVelocity(outStr.c_str());
    cout << outStr << endl;

    testFormatter.genOutStr("MSSEdges", outStr);
    pool.outputAllStructures(outStr.c_str());
    cout << outStr << endl;

    testFormatter.genOutStr("MSSNodes", outStr);
    pool.outputAllStructureNodes(outStr.c_str());
    cout << outStr << endl;

    testFormatter.genOutStr("MSSVels", outStr);
    pool.outputAllStructureVels(outStr.c_str());
    cout << outStr << endl;

    testFormatter.genOutStr("MSSTracers", outStr);
    pool.outputTracers(outStr.c_str());
    cout << outStr << endl;

    testFormatter.genOutStr("medialAxis", outStr);
    pool.outputMedialAxisApprox(outStr.c_str());
    cout << outStr << endl;
}