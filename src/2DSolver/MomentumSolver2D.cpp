#include <iostream>
#include <unordered_map>
#include <functional>
#include "Pool2D.h"
#include "SolidObject.h"
#include "MomentumSolver2D.h"
#include "../Utils/SimUtilities.h"
#include <iostream>
#include "PressureSolver.h"
#include <fstream>
#include <string>
#include <cassert>
#include <vector>
#include "SimParams.h"

using namespace std;

MomentumSolver2D::MomentumSolver2D(Boundary &boundary,
    vector<SolidObject> &solidObjects, SimParams &params,
    void (*initialConditions)(int,int,int,double*,double*,double**,double**),
    void (*boundaryConditions)(int,int,double**,double**))
{

    // Check that the parameters object is set
    if (!params.checkParams()) {
        cout << "ERROR: 2D Sim params not set correctly!" << endl;
        assert(false);
    }

    // Copy basic info
    this->nx = params.nx;
    this->ny = params.nx;
    this->methodOrd = params.methodOrd;

    // Copy the constant body forces (gravity)
    this->gx = params.gx;
    this->gy = params.gy;

    // Create the uniform 1D meshes using the boundary object
    this->x = boundary.generateXMesh(this->nx);
    this->y = boundary.generateYMesh(this->ny);

    // The mesh spacing
    this->dx = this->x[1] - this->x[0];
    this->dy = this->y[1] - this->y[0];

    // Copy the input parameters using the default copy ctor
    this->params = &params;

    // Create the output arrays with required ghost cells for the staggered grid
    // configuration.
    // TODO: make sure that this is correct
    int hg = this->ny + 2*this->methodOrd;
    int wg = this->nx + 2*this->methodOrd;

    this->u = simutils::new_constant(hg, wg-1, 0.0);
    this->iu = simutils::new_constant(this->ny, this->nx, 0.0);
    this->FU = simutils::new_constant(hg, wg-1, 0.0);
    this->v = simutils::new_constant(hg-1, wg, 0.0);
    this->iv = simutils::new_constant(this->ny, this->nx, 0.0);
    this->FV = simutils::new_constant(hg-1, wg, 0.0);
    this->p = simutils::new_constant(hg, wg, 0.0);

    // Create pressure solver and its parameters. Parameters are default for now
    // TODO: parameters should be passed as an argument.
    this->pressureSolver = new PressureSolver(this->nx, this->x, this->ny, this->y);
    this->pressureParams = new ParamIter();

    // Set the default params
    const int ILU_LEVEL = 5;

    (this->pressureParams)->order = 0; // = 0 natural;
                    //  =1, RCM

    (this->pressureParams)->level = ILU_LEVEL; // level of ilu

    (this->pressureParams)->drop_ilu = 0; // = 0 level ilu
                        // = 1 drop tol

    (this->pressureParams)->iscal = 0; // = 0 no scaling by inverse of diag
                    // = 1 scale by inverse of diag

    (this->pressureParams)->nitmax = 10000; // max number of iterations

    (this->pressureParams)->ipiv = 0; // 0 no pivoting
                                    // 1 pivoting
                                    // only if drop_ilu = 1

    (this->pressureParams)->resid_reduc = 1.e-6; // residual reduction toleranc
                                // iterations stop if any of the following true
                                //   l^2 residual reduced by ctol
                                //   number of itns = nitmax
                                //   update in each variable x[i]<toler[i]

    (this->pressureParams)->info = 0; // = 0 no iteration info written to "output"
                    // = 1 write iteration info

    (this->pressureParams)->drop_tol = 1.e-3; // drop tolerance used if drop_ilu = 1

    (this->pressureParams)->new_rhat = 0; // = 0 use r^0 as rhat
                        // = 1 use (LU)^{-1^0 for rhat

    // (this->pressureParams)->iaccel = -1; // = 0 cgstab
    (this->pressureParams)->iaccel = 0; // = 0 cgstab
                        // = 1 orthomin
                        // = -1 conj gradiend (SPD only!)

    (this->pressureParams)->north = 10; // number of orthogs for orthomin

    // Use initial condition callback function to assign u, v
    // TODO: make initial conditions more general and perhaps more
    //       able to handle irregular boundaries.
    initialConditions(this->nx, this->ny, this->methodOrd, this->x, this->y, this->u, this->v);

    applyFluidBCs = boundaryConditions;

    // Create the pool object for this solver
    this->nStructs = solidObjects.size();
    this->pool = new Pool2D(boundary, solidObjects, params);

    // Pre-compute the Barycentric weights
    // TODO: ensure that this is correct for the FV methods used later
    this->baryWeights = new double[this->methodOrd+1];
    simutils::barycentricInterp(this->methodOrd, this->baryWeights);

    this->stepTaken = false;
}

/**
 * Set the boundary conditions within the structures to ensure things are OK on the next time step
 * 
*/
void MomentumSolver2D::test_setInternal() {
    int i, j;
    int mo = methodOrd;
    for (j = 0; j < this->ny; j++) {
        for (i = 0; i < this->nx-1; i++) {
            int xi = i + mo;
            int yi = j + mo;

            if (pool->objAtIndex(i, j)!=objects::FLUID_C) {
                this->u[yi][xi] = (this->pool)->getObjU(i, j);
                this->u[yi][xi-1] = (this->pool)->getObjU(i, j);
            }
        }
    }

    for (j = 0; j < this->ny-1; j++) {
        for (i = 0; i < this->nx; i++) {
            int xi = i + mo;
            int yi = j + mo;

            if (pool->objAtIndex(i, j)!=objects::FLUID_C) {
                this->v[yi][xi] = (this->pool)->getObjV(i, j);
                this->v[yi-1][xi] = (this->pool)->getObjV(i, j);
            }
        }
    }
}

/**
 * Method to update U using F and P. Note only the internal points are updated
 * as we only ever apply a relatively simple Dirichlet condition at the boundaries
 * for the fluids.
 * 
 * // TODO: generalize how the pressure derivatives are being evaluated to at least second order.
 */
void MomentumSolver2D::updateU() {
    int i, j;
    int mo = this->methodOrd;

    // Update u
    for (j = 0; j < this->ny; j++) {
        for (i = 0; i < this->nx-1; i++) {
            int xi = i + mo;
            int yi = j + mo;
            if (pool->objAtIndex(i, j) == objects::FLUID_C && pool->isUpdateableU(i, j)) {
                this->u[yi][xi] = this->FU[yi][xi] - (this->dt/this->dx)*(this->p[yi][xi+1] - this->p[yi][xi]);
            }
        }
    }

    // Update v
    for (j = 0; j < this->ny-1; j++) {
        for (i = 0; i < this->nx; i++) {
            int xi = i + mo;
            int yi = j + mo;
            if (pool->objAtIndex(i, j) == objects::FLUID_C && pool->isUpdateableV(i, j)) {
                this->v[yi][xi] = this->FV[yi][xi] - (this->dt/this->dy)*(this->p[yi+1][xi] - this->p[yi][xi]);
            }
        }
    }
}

/**
 * Method which interpolate fluid velocities to the cell centers.
 * 
*/
void MomentumSolver2D::interpolateVelocities() {
    int i, j;
    double x, y;

    double uvals[2];
    double vvals[2];

    // Interpolate u
    for (j = 0; j < this->ny; j++) {
        y = simutils::midpoint(this->y[j], this->y[j+1]);
        for (i = 0; i < this->nx; i++) {
            x = simutils::midpoint(this->x[i], this->x[i+1]);

            // The values being interpolated
            uvals[0] = this->u[this->methodOrd+j][i];
            uvals[1] = this->u[this->methodOrd+j][i+1];

            vvals[0] = this->v[j][this->methodOrd+i];
            vvals[1] = this->v[j+1][this->methodOrd+i];

            // Interpolate the value at the cell centers
            this->iu[j][i] = simutils::barycentricInterv(x, this->x[i], this->x[i+1],
                                this->methodOrd, this->baryWeights, uvals);
            this->iv[j][i] = simutils::barycentricInterv(y, this->y[j], this->y[j+1],
                                this->methodOrd, this->baryWeights, vvals);
        }
    }
}

/**
 * Method to update the Pressure equation using a second order method.
 * The linear system is solved using ILU preconditioned PCG.
*/
void MomentumSolver2D::updateP() {
    (this->pressureSolver)->solvePCG(this->dt, this->pool, this->FU, this->FV,
                                     this->p, true, *(this->pressureParams));
}

/**
 * General implementation of the explicit-implicit scheme.
*/
double MomentumSolver2D::step(double tEnd, double safetyFactor) {

    // Compute the time step
    if (params->dtFixSet) {
        cout << "SETTING POOL DTFIX" << endl;
        this->dt = params->dtFix;
        cout << "SETTING POOL DTFIX" << endl;
    } else {
        this->dt = safetyFactor*this->getDt();

    }
    if (!stepTaken) {
        this->dtPrev = this->dt;
    }

    // If the time step would exceed tEnd, truncate it
    if (this->t + this->dt > tEnd) {
        this->dt = tEnd - this->t;
    }

    this->applyInterfaceBCs();

    // Explicit step: advance the velocity-dependent terms using explicit time
    // discretization.
    this->updateF(pool);

    // Implicit step: update the pressure by solving Poisson equation.
    this->updateP();

    // Combine these together to update the fluid velocities.
    cout << "updating U" << endl;
    this->updateU();
    cout << "FINISHED updating U" << endl;

    // Interpolate the velocities to the cell centers
    cout << "Interpolating U" << endl;
    this->interpolateVelocities();
    cout << "FINISHED Interpolating U" << endl;

    if (this->nStructs > 0) {
        // Update the location of the interfaces
        cout << "Updating Pool" << endl;
        (this->pool)->updatePool(dt, iu, iv, p, methodOrd, true);
        cout << "FINISHED Updating Pool" << endl;

        // Apply object velocities to boundary points
        this->test_setInternal();
    }

    // Update the point in time
    this->t += this->dt;

    // Save the current time step
    this->dtPrev = this->dt;

    // The fist step has now been taken
    this->stepTaken = true;

    return this->t;
}

void MomentumSolver2D::writeToFile(const char *fname) {
    double x, y;
    int mo = methodOrd;

    // Open the output file
    ofstream outFile;
    outFile.open(fname);

    for (int j = 0; j < this->ny; j++) {
        for (int i = 0; i < this->nx; i++) {
            int xi = i + mo;
            int yi = j + mo;
            if (pool->objAtIndex(i, j) != objects::STRUCTURE) {
                x = simutils::midpoint(this->x[i], this->x[i+1]);
                y = simutils::midpoint(this->y[j], this->y[j+1]);
                outFile << x << ", ";
                outFile << y << ", ";
                outFile << this->iu[j][i] << ", ";
                outFile << this->iv[j][i] << ", ";
                outFile << this->p[yi][xi] << "\n";
            } else {
                outFile << 0.0 << ", ";
                outFile << 0.0 << ", ";
                outFile << 0.0 << ", ";
                outFile << 0.0 << ", ";
                outFile << 0.0 << "\n";
            }
        }
    }

    outFile.close();
}

void MomentumSolver2D::writePoolToFile(const char *poolFName, const char *poolVelFName) {
    pool->outputPool(poolFName);
    pool->outputPoolVelocity(poolVelFName);
}

MomentumSolver2D::~MomentumSolver2D() {
    delete this->pool;
    delete[] this->x;
    delete[] this->y;
    delete[] this->baryWeights;

    delete pressureSolver;
    delete pressureParams;

    // Free multidimensional arrays using helper
    int hg = this->ny + 2*this->methodOrd;
    int wg = this->nx + 2*this->methodOrd;

    // TODO: how necissary are the derivative arrays? My guess is not very? Will keep
    //       them around for ease of implementation
    simutils::free_double(hg, wg-1, this->u);
    simutils::free_double(this->ny, this->nx, this->iu);
    simutils::free_double(hg, wg-1, this->FU);
    simutils::free_double(hg-1, wg, this->v);
    simutils::free_double(this->ny, this->nx, this->iv);
    simutils::free_double(hg-1, wg, this->FV);
    simutils::free_double(hg, wg, this->p);
}