#include <iostream>
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

double circleShapeFun(double x, double y, double z, SolidParams &ps) {
    double cx, cy, cz, r;
    ps.getParam("cx", cx);
    ps.getParam("cy", cy);
    ps.getParam("cz", cz);
    ps.getParam("r", r);
    // return simutils::square(x-ps["cx"])+simutils::square(y-ps["cy"])+simutils::square(z-ps["cz"])-simutils::square(ps["r"]);
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

int main(int argc, char **argv) {
    // Testing the fluid solver
    double tEnd = 6.4;

    // The boundaries
    double xa = 0, xb = 1;
    double ya = 0, yb = 1;
    double za = 0, zb = 1;

    // Number of x, y points
    int nx = 20;
    int ny = 20;
    int nz = 20;

    /* Creation of the solid objects */
    ///////////////////////////////////
    int numStructs = 1;
    // int numStructs = 2;

    // unordered_map<string, double> circParams1;
    // unordered_map<string, double> circParams2;

    SolidParams circParams1;
    SolidParams circParams2;


    double mass = 100.0;
    double density = 1.0;
    double E = 5.0;
    double eta = 0.0;


    // circParams1.addParam("cx", 0.25);
    circParams1.addParam("cx", 0.5);
    circParams1.addParam("cy", 0.5);
    circParams1.addParam("cz", 0.5);
    circParams1.addParam("r", 0.15);
    circParams1.addParam("mass", mass);
    circParams1.addParam("density", density);
    circParams1.addParam("E", E);
    // params["E"] = 5.0;
    circParams1.addParam("eta", eta);

    circParams2.addParam("cx", 0.75);
    circParams2.addParam("cy", 0.5);
    circParams2.addParam("cz", 0.5);
    circParams2.addParam("r", 0.15);
    circParams2.addParam("mass", mass);
    circParams2.addParam("density", density);
    circParams2.addParam("E", E);
    // params["E"] = 5.0;
    circParams2.addParam("eta", eta);

    double u0 = 0.0;
    double v0 = 0.0;
    double w0 = 0.0;

    SolidObject3D::ObjectType deformableBody = SolidObject3D::ObjectType::DEFORMABLE;

    // The individual circles
    SolidObject3D circle1(u0, v0, w0, deformableBody, coneShapeFun, circParams1);
    SolidObject3D circle2(u0, v0, w0, deformableBody, coneShapeFun, circParams2);

    // Create circle object array for embedding in Pool
    std::vector<SolidObject3D> circles;
    // circles.push_back(circle2);
    circles.push_back(circle1);

    // int numStructs = 1;
    // unordered_map<string, double> circParams;

    // circParams["cx"] = 0.5;
    // circParams["cy"] = 0.5;
    // circParams["cz"] = 0.5;
    // circParams["r"] = 0.20;
    // // Solid object null pointer
    // double mass = 10;
    // double u0 = 0.0;
    // double v0 = 0.0;
    // double w0 = 0.0;
    // double density = 1.0;
    // SolidObject3D circle(u0, v0, w0, mass, density, coneShapeFun, circParams);
    ///////////////////////////////////

    // Parameters of this circle
    // Boundary object
    Boundary3D boundary(xa, xb, ya, yb, za, zb);


    // Hashmap of any required params
    SimParams3D params;
    params.setRe(1000.0);
    // params["mu"] = 1/params["Re"];
    params.setMu(1.0/params.Re);
    params.setNx(nx);
    params.setNy(ny);
    params.setNz(nz);
    params.setMssNx(nx);
    params.setMssNy(ny);
    params.setMssNz(nz);
    params.setUseEno(false);
    params.setRepulseMode(0);
    double h = sqrt(simutils::square(1.0/((double) nx)
        + simutils::square(1.0/((double) ny)
        + simutils::square(1.0/((double) nz)))));
    params.setRepulseDist(3.0*h);
    params.setCollisionDist(3.0*h);
    params.setCollisionStiffness(2.0);
    params.setRepulseDist(3.0*h);
    params.setUpdateMode(1);

    //   params.setRe(1000);
    // params.setNx(nx);
    // params.setNy(ny);
    // params.setUseEno(true);
    // params.setMu(1.0/params.Re);
    // params.setRepulseMode(2); // This turns on the KD tree error checking
    // // simParams.setRepulseDist(5*sqrt(simutils::square(1/((double)nx)) + simutils::square(1/((double)ny))) );
    // // params.setRepulseDist(0.1); // Actually need 0.1
    // // params.setCollisionStiffness(2.0);
    // // params.setCollisionDist(0.25);
    // double h = sqrt(simutils::square(1.0/((double) nx)
    //     + simutils::square(1.0/((double) ny))));
    // params.setRepulseDist(3*h); // Actually need 0.1
    // params.setCollisionStiffness(2.0);
    // params.setCollisionDist(3*h);
    // params.setUpdateMode(1);

    // Create the Solver object
    std::cout << "creating the solver" << std::endl;
    NSSolver3D solver(boundary, circles, params, initialConditions);
    std::cout << "created the solver" << std::endl;
    // NSSolver3D solver(nx, ny, nz, boundary, numStructs, true, &circle1, params, initialConditions);

    // Current time
    double t = 0;
    double safetyFactor = 0.5;

    std::cout << "Going to update the pool" << std::endl;

    int nsteps = 0;
    // int max_steps = 1;
    int max_steps = (argc == 1) ? 1 : atoi(argv[1]);
    std::cout << "TAKING " << max_steps << " STEPS" << std::endl;
    double eps = 1e-12;
    while (t+eps < tEnd && nsteps < max_steps) {
        std::cout << "Taking a step" << std::endl;
        t = solver.step(tEnd, safetyFactor);
        std::cout << "Finished a step" << std::endl;

        nsteps++;

        std::cout << "t = " << t << std::endl;
        std::cout << "step = " << nsteps << std::endl;
        std::cout << std::endl;
    }

    // Write to the output file.
    // std::cout << "Outputting data" << std::endl;
    // solver.writeToFile("out3D.txt");
    // solver.writePoolToFile("pool3DOut.txt", "pool3DVel.txt");
    // solver.outputStructure(0, "MSS3DEdges0.txt");
    // // solver.outputStructure(1, "MSS3DEdges1.txt");
    // solver.outputStructureNodes(0, "MSS3DNodes0.txt");
    // // solver.outputStructureNodes(1, "MSS3DNodes1.txt");
    // solver.outputStructureVels(0, "MSS3DVels0.txt");
    // // solver.outputStructureVels(1, "MSS3DVels1.txt");

    // std::cout << "end of main" << std::endl;
    TestFormatter testFormatter("SimpleTest3D");
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