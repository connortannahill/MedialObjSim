#include <iostream>
#include <unordered_map>
#include <functional>
#include "Pool3D.h"
#include "SolidObject3D.h"
#include "MomentumSolver3D.h"
#include "../Utils/SimUtilities.h"
#include <iostream>
#include "PressureSolver3D.h"
#include <fstream>
#include <string>
#include <cassert>
#include "SimParams3D.h"

using namespace std;

MomentumSolver3D::MomentumSolver3D(
    Boundary3D &boundary, vector<SolidObject3D> &solidObjects,
    SimParams3D &params,
    std::function<void (int,int,int,int,double*,double*,double*,double***,double***,double***)> initialConditions,
    void (*boundaryConditions)(int,int,int,double***,double***,double***))
{

    if (!params.checkParams()) {
        cout << "3D Model Parameters not set correctly!" << endl;
        assert(false);
    }

    // Copy basic info
    this->nx = params.nx;
    this->ny = params.ny;
    this->nz = params.nz;
    this->methodOrd = params.methodOrd;

    // Body forces
    this->gx = params.gx;
    this->gy = params.gy;
    this->gz = params.gz;

    // Create the uniform 1D meshes using the boundary object
    this->x = boundary.generateXMesh(this->nx);
    this->y = boundary.generateYMesh(this->ny);
    this->z = boundary.generateZMesh(this->nz);

    // The mesh spacing
    this->dx = this->x[1] - this->x[0];
    this->dy = this->y[1] - this->y[0];
    this->dz = this->z[1] - this->z[0];

    // Create the pool object for this solver
    this->nStructs = solidObjects.size();
    // cout << "Creating the pool" << endl;
    this->pool = new Pool3D(boundary, solidObjects, params);
    // cout << "Finished Creating the pool" << endl;

    // Copy the input parameters using the default copy ctor
    this->params = &params;

    // Create the output arrays with required ghost cells for the staggered grid
    // configuration.
    int hg = this->ny + 2*this->methodOrd; // "height"
    int wg = this->nx + 2*this->methodOrd; // "width"
    int vg = this->nz + 2*this->methodOrd; // "vert"

    this->u = simutils::new_constant(vg, hg, wg-1, 0.0);
    this->iu = simutils::new_constant(this->nz, this->ny, this->nx, 0.0);
    this->FU = simutils::new_constant(vg, hg, wg-1, 0.0);
    this->v = simutils::new_constant(vg, hg-1, wg, 0.0);
    this->iv = simutils::new_constant(this->nz, this->ny, this->nx, 0.0);
    this->FV = simutils::new_constant(vg, hg-1, wg, 0.0);
    this->w = simutils::new_constant(vg-1, hg, wg, 0.0);
    this->iw = simutils::new_constant(this->nz, this->ny, this->nz, 0.0);
    this->FW = simutils::new_constant(vg-1, hg, wg, 0.0);
    this->p = simutils::new_constant(vg, hg, wg, 0.0);

    // Create pressure solver and its parameters. Parameters are default for now
    // TODO: parameters should be passed as an argument.
    this->pressureSolver = new PressureSolver3D(this->nx, this->x, this->ny, this->y, this->nz, this->z);
    this->pressureParams = new ParamIter();

    // Set the default params
    const int ILU_LEVEL = 1;

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
    initialConditions(this->nx, this->ny, this->nz, this->methodOrd, this->x,
                        this->y, this->z, this->u, this->v, this->w);

    applyFluidBCs = boundaryConditions;

    // Pre-compute the Barycentric weights
    this->baryWeights = new double[this->methodOrd+1];
    simutils::barycentricInterp(this->methodOrd, this->baryWeights);

    this->stepTaken = false;
}

/**
 * Set the boundary conditions within the structures to ensure things are OK on the next time step
*/
void MomentumSolver3D::test_setInternal() {
    int i, j, k;
    int mo = methodOrd;
    int xi, yi, zi;
    for (k = 0; k < this->nz; k++) {
        for (j = 0; j < this->ny; j++) {
            for (i = 0; i < this->nx-1; i++) {
                xi = i + mo;
                yi = j + mo;
                zi = k + mo;

                if (pool->objAtIndex(i, j, k) != objects::FLUID_C) {
                    this->u[zi][yi][xi] = (this->pool)->getObjU(i, j, k);
                    this->u[zi][yi][xi-1] = (this->pool)->getObjU(i, j, k);
                }
            }
        }
    }

    for (k = 0; k < this->nz; k++) {
        for (j = 0; j < this->ny-1; j++) {
            for (i = 0; i < this->nx; i++) {
                xi = i + mo;
                yi = j + mo;
                zi = k + mo;

                if (pool->objAtIndex(i, j, k) != objects::FLUID_C) {
                    this->v[zi][yi][xi] = (this->pool)->getObjV(i, j, k);
                    this->v[zi][yi-1][xi] = (this->pool)->getObjV(i, j, k);
                }
            }
        }
    }

    for (k = 0; k < this->nz-1; k++) {
        for (j = 0; j < this->ny; j++) {
            for (i = 0; i < this->nz; i++) {
                xi = i + mo;
                yi = j + mo;
                zi = k + mo;

                if (pool->objAtIndex(i, j, k) != objects::FLUID_C) {
                    this->w[zi][yi][xi] = (this->pool)->getObjW(i, j, k);
                    this->w[zi-1][yi][xi] = (this->pool)->getObjW(i, j, k);
                }
            }
        }
    }
}

/**
 * Method to update U using F and P. Note only the internal points are updated
 * as we only ever apply a relatively simple Dirichlet condition at the boundaries
 * for the fluids.
 */
void MomentumSolver3D::updateU() {
    int i, j, k;
    int xi, yi, zi;
    int mo = this->methodOrd;

    // Update u
    for (k = 0; k < this->nz; k++) {
        for (j = 0; j < this->ny; j++) {
            for (i = 0; i < this->nx-1; i++) {
                xi = i + mo;
                yi = j + mo;
                zi = k + mo;

                if (pool->objAtIndex(i, j, k) == objects::FLUID_C && pool->isUpdateableU(i, j, k)) {
                    this->u[zi][yi][xi] = this->FU[zi][yi][xi]
                                - (this->dt/this->dx)*(this->p[zi][yi][xi+1] - this->p[zi][yi][xi]);
                }
            }
        }
    }

    // Update v
    for (k = 0; k < this->nz; k++) {
        for (j = 0; j < this->ny-1; j++) {
            for (i = 0; i < this->nx; i++) {
                xi = i + mo;
                yi = j + mo;
                zi = k + mo;

                if (pool->objAtIndex(i, j, k) == objects::FLUID_C && pool->isUpdateableV(i, j, k)) {
                    this->v[zi][yi][xi] = this->FV[zi][yi][xi]
                                - (this->dt/this->dy)*(this->p[zi][yi+1][xi] - this->p[zi][yi][xi]);
                }
            }
        }
    }

    // Finally, update w
    for (k = 0; k < this->nz-1; k++) {
        for (j = 0; j < this->ny; j++) {
            for (i = 0; i < this->nx; i++) {
                xi = i + mo;
                yi = j + mo;
                zi = k + mo;

                if (pool->objAtIndex(i, j, k) == objects::FLUID_C && pool->isUpdateableW(i, j, k)) {
                    this->w[zi][yi][xi] = this->FW[zi][yi][xi]
                                - (this->dt/this->dz)*(this->p[zi+1][yi][xi] - this->p[zi][yi][xi]);
                }
            }
        }
    }
}

/**
 * Method which interpolate fluid velocities to the cell centers.
 * 
 * TODO: make sure that this implementation makes sense
*/
void MomentumSolver3D::interpolateVelocities() {
    int i, j, k;
    double x, y, z;
    int mo = this->methodOrd;

    double uvals[2];
    double vvals[2];
    double wvals[2];

    // Interpolate u
    for (k = 0; k < this->nz; k++) {
        z = simutils::midpoint(this->z[k], this->z[k+1]);
        for (j = 0; j < this->ny; j++) {
            y = simutils::midpoint(this->y[j], this->y[j+1]);
            for (i = 0; i < this->nx; i++) {
                x = simutils::midpoint(this->x[i], this->x[i+1]);

                // The values being interpolated
                uvals[0] = this->u[mo+k][mo+j][i];
                uvals[1] = this->u[mo+k][mo+j][i+1];

                vvals[0] = this->v[mo+k][j][mo+i];
                vvals[1] = this->v[mo+k][j+1][mo+i];

                wvals[0] = this->w[k][mo+j][mo+i];
                wvals[0] = this->w[k+1][mo+j][mo+i];

                // Interpolate the value at the cell centers
                this->iu[k][j][i] = simutils::barycentricInterv(x, this->x[i], this->x[i+1],
                                    this->methodOrd, this->baryWeights, uvals);
                this->iv[k][j][i] = simutils::barycentricInterv(y, this->y[j], this->y[j+1],
                                    this->methodOrd, this->baryWeights, vvals);
                this->iw[k][j][i] = simutils::barycentricInterv(z, this->z[k], this->z[k+1],
                                    this->methodOrd, baryWeights, wvals);
            }
        }
    }
}

/**
 * Method to update the Pressure equation using a second order method.
 * The linear system is solved using ILU preconditioned PCG.
*/
void MomentumSolver3D::updateP() {
    (this->pressureSolver)->solvePCG(this->dt, this->pool, this->FU, this->FV, this->FW,
                                    this->p, true, *(this->pressureParams));
}

/**
 * General implementation of the explicit-implicit scheme.
*/
double MomentumSolver3D::step(double tEnd, double safetyFactor) {
    // Compute the time step
    this->dt = safetyFactor*this->getDt();
    if (!stepTaken) {
        this->dtPrev = this->dt;
    }

    // If the time step would exceed tEnd, truncate it
    if (this->t + this->dt > tEnd) {
        this->dt = tEnd - this->t;
    }

    // Apply the fluid boundary conditions on the initial step using pure virtual function.
    this->applyInterfaceBCs();

    // Explicit step: advance the velocity-dependent terms using explicit time
    // discretization.
    this->updateF(pool);

    // Implicit step: update the pressure by solving Poisson equation.
    this->updateP();

    // Combine these together to update the fluid velocities.
    // cout << "Updating u" << endl;
    this->updateU();
    // cout << "FINISHED updating u" << endl;

    // Interpolate the velocities to the cell centers
    // cout << "interping velocities" << endl;
    this->interpolateVelocities();
    // cout << "FINISHED interping velocities" << endl;

    // If there are any structures in the pool, update the pool location.
    if (this->nStructs > 0) {
        // cout << "Updating pool" << endl;
        // Update the location of the interfaces
        (this->pool)->updatePool(dt, iu, iv, w, p, methodOrd, true);
        // cout << "FINISHED Updating pool" << endl;

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

void MomentumSolver3D::writePoolToFile(const char *poolFName, const char *poolVelFName) {
    pool->outputPool(poolFName);
    pool->outputPoolVelocity(poolVelFName);
}

void MomentumSolver3D::writeToFile(const char *fname) {
    double x, y, z;
    int mo = methodOrd;

    // Open the output file
    ofstream outFile;
    outFile.open(fname);

    int xi, yi, zi;

    for (int k = 0; k < this->nz; k++) {
        for (int j = 0; j < this->ny; j++) {
            for (int i = 0; i < this->nx; i++) {
                xi = i + mo;
                yi = j + mo;
                zi = k + mo;
                x = simutils::midpoint(this->x[i], this->x[i+1]);
                y = simutils::midpoint(this->y[j], this->y[j+1]);
                z = simutils::midpoint(this->z[k], this->z[k+1]);

                outFile << x << ", ";
                outFile << y << ", ";
                outFile << z << ", ";

                if (pool->objAtIndex(i, j, k) == objects::FLUID_C) {
                    outFile << this->iu[k][j][i] << ", ";
                    outFile << this->iv[k][j][i] << ", ";
                    outFile << this->iw[k][j][i] << ", ";
                    outFile << this->p[zi][yi][xi] << "\n";

                } else {
                    outFile << 0.0 << ", ";
                    outFile << 0.0 << ", ";
                    outFile << 0.0 << ", ";
                    outFile << 0.0 << "\n";
                }
            }
        }
    }

    outFile.close();
}

MomentumSolver3D::~MomentumSolver3D() {
    delete this->pool;
    delete[] this->x;
    delete[] this->y;
    delete[] this->z;
    delete[] this->baryWeights;

    delete pressureSolver;
    delete pressureParams;

    // Free multidimensional arrays using helper
    int hg = this->ny + 2*this->methodOrd;
    int wg = this->nx + 2*this->methodOrd;
    int vg = this->nz + 2*this->methodOrd;

    simutils::free_double(vg, hg, wg-1, this->u);
    simutils::free_double(this->nz, this->ny, this->nx, this->iu);
    simutils::free_double(vg, hg, wg-1, this->FU);

    simutils::free_double(vg, hg-1, wg, this->v);
    simutils::free_double(this->nz, this->ny, this->nx, this->iv);
    simutils::free_double(vg, hg-1, wg, this->FV);

    simutils::free_double(vg-1, hg, wg, this->w);
    simutils::free_double(this->nz, this->ny, this->nz, this->iw);
    simutils::free_double(vg-1, hg, wg, this->FW);

    simutils::free_double(vg, hg, wg, this->p);
}