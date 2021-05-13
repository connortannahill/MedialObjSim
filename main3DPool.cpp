#include <iostream>
#include "./src/3DSolver/Boundary3D.h"
#include "./src/3DSolver/SolidObject3D.h"
#include "./src/3DSolver/Pool3D.h"
#include "./src/3DSolver/MassSpring3D.h"
#include "./src/Utils/SolidParams.h"
#include <math.h>
#include "./src/Utils/SimUtilities.h"
#include "./src/Utils/Discretizations.h"
#include <cassert>
#include "./src/3DSolver/SimParams3D.h"

using namespace std;
using namespace mass_spring;

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

int main(int argc, char **argv) { 
    int nx = 20;
    int ny = 20;
    int nz = 20;

    double xmin = 0;
    double xmax = 1;
    double ymin = 0;
    double ymax = 1;
    double zmin = 0;
    double zmax = 1;

    // Create the boundary object
    Boundary3D boundary(xmin, xmax, ymin, ymax, zmin, zmax);

    // Parameters of this circle
    SolidParams params1;
    params1.addParam("cx", 0.25);
    // params["cx"] = 0.5;
    params1.addParam("cy", 0.5);
    params1.addParam("cz", 0.5);
    params1.addParam("r", 0.15);
    params1.addParam("mass", 1.0);
    params1.addParam("density", 1.0);
    params1.addParam("E", 15.0);
    // params["E"] = 5.0;
    params1.addParam("eta", 0.0);

    // Create the solid structure
    double u0 = 0.0;
    double v0 = 0.0;
    double w0 = 0.0;
    bool deformableBody = true;
    SolidObject3D circle(u0, v0, w0, deformableBody, coneShapeFun, params1);
    // SolidObject3D circle(u0, v0, w0, deformableBody, circleShapeFun, params);

    SolidParams params2;
    params2.addParam("cx", 0.75);
    params2.addParam("cy", 0.5);
    params2.addParam("cz", 0.5);
    params2.addParam("r", 0.15);
    params2.addParam("mass", 1.0);
    params2.addParam("density", 1.0);
    params2.addParam("E", 15.0);
    params2.addParam("eta", 0.0);
    SolidObject3D circle2(u0, v0, w0, deformableBody, coneShapeFun, params2); 

    // SolidObject3D objs[2] = {circle, circle2};
    std::vector<SolidObject3D> circles;
    circles.push_back(circle);
    circles.push_back(circle2);
    // SolidObject3D objs[1] = {circle};

    // Create the pool
    SimParams3D simParams;
    simParams.setNx(nx);
    simParams.setNy(ny);
    simParams.setNz(nz);

    simParams.setMssNx(nx);
    simParams.setMssNy(ny);
    simParams.setMssNz(nz);

    simParams.setRe(0.0);
    simParams.setMu(0.0);

    simParams.setRepulseMode(2); // This turns on the KD tree error checking
    // simParams.setRepulseDist(5*sqrt(simutils::square(1/((double)nx)) + simutils::square(1/((double)ny))) );
    simParams.setRepulseDist(0.1); // Actually need 0.1
    simParams.setCollisionStiffness(2.0);
    simParams.setCollisionDist(0.25);

    double dt = 0.5/(1.0/(1/((double)nx))+1.0/((1/(double)ny))+1.0/((1/(double)nz)));

    simParams.setDtFix(dt);

    cout << "creating the pool" << endl;
    Pool3D pool(boundary, circles, simParams);
    cout << "FINSIEHD creating the pool" << endl;
    // assert(false);

    double t = 0;
    double tEnd = 10;
    // double dt = 0.5/(1.0/(1.0/((double)nx))+1.0/((1.0/(double)ny))+1.0/(1.0/((double)nz)));

    std::cout << "dt = " << dt << std::endl;

    // Test velocities for the level set methods
    double ***u = simutils::new_constant(nz, ny, nx, 0.0);
    double ***v = simutils::new_constant(nz, ny, nx, 0.0);
    double ***w = simutils::new_constant(nz, ny, nx, 0.0);
    double ***p = simutils::new_constant(nz, ny, nx, 0.0);

    int nstep = 0;
    std::cout << "Hi" << std::endl;
    int max_steps = (argc == 1) ? 1 : atoi(argv[1]);
    std::cout << "Bye" << std::endl;
    bool reinitialize = true;
    while (t < tEnd && nstep < max_steps) {
        pool.updatePool(dt, u, v, w, p, 0, reinitialize); // Advance the level set function

        t+=dt;
        std::cout << "t = " << t << std::endl;
        std::cout << "nStep = " << nstep << std::endl;
        nstep++;
    }

    string testTitle("PoolTests2D");

    std::cout << "Outputting all of the data" << std::endl;

    pool.outputPool("pool3DOut.txt");
    pool.outputPoolVelocity("pool3DVel.txt");
    pool.outputStructure(0, "MSS3DEdges0.txt");
    // pool.outputStructure(1, "MSS3DEdges1.txt");
    pool.outputStructureNodes(0, "MSS3DNodes0.txt");
    // pool.outputStructureNodes(1, "MSS3DNodes1.txt");
    pool.outputStructureVels(0, "MSS3DVels0.txt");
    // pool.outputStructureVels(1, "MSS3DVels1.txt");
    pool.outputSurfaceCentroids(0, "MSS3DCentroids0.txt");

    std::cout << "Finished updating all of the data" << std::endl;
}