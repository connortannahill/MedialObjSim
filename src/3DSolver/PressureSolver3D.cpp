#include "PressureSolver3D.h"
#include "../Utils/SimUtilities.h"
#include <math.h>
#include "../LASolver/SparseItObj.h"
#include "Boundary3D.h"
#include <cassert>

using namespace std;

/** 
 * Method which creates a new matrix object.
*/
void PressureSolver3D::setUpMatrixStruct(int nx, int ny, int nz) {
    int n = (nx+2)*(ny+2)*(nz+2);
    int i, j, k;
    int row = 0;

    // Set up the system matrix. Matrix is set up to have a 5-stencil at every point
    // except the sides & corners to permit future generality.
    MatrixStruc *matrixBuilder = NULL;
    try {
        matrixBuilder = new MatrixStruc(n /*Number of unknowns*/, 0 /*Add non-zeros to diagonals*/);

        /* Create the data structure */
        {
            // BOTTOM (in height)
            // Apply Dirichlet boundary condition to all cells on the bottom and top
            for (j = 0; j < ny+2; j++) {
                for (i = 0; i < nx+2; i++) {
                    // TODO: this is actually already handled
                    row++;
                }
            }

            // INTERNAL (in height)
            for (k = 1; k < nz+1; k++) {
                for (j = 0; j < ny+2; j++) {
                    for (i = 0; i < nx+2; i++) {
                        // Dirichlet cond at boundaries
                        if (j == 0 || j == ny+1 || i == 0 || i == nx+1) {
                        } else {
                            // Centered finite difference stencil
                            
                            matrixBuilder->set_entry(row, row+1); // (i+1, j, k)
                            matrixBuilder->set_entry(row, row-1); // (i-1, j, k)

                            matrixBuilder->set_entry(row, row+(nx+2)); // (i, j+1, k)
                            matrixBuilder->set_entry(row, row-(nx+2)); // (i, j-1, k)

                            matrixBuilder->set_entry(row, row+((nx+2)*(ny+2))); // (i, j, k+1)
                            matrixBuilder->set_entry(row, row-((nx+2)*(ny+2))); // (i, j, k-1)
                        }
                        row++;
                    }
                }
            }

            // TOP (in height)
            for (j = 0; j < ny+2; j++) {
                for (i = 0; i < nx+2; i++) {
                    row++;
                }
            }
        }
        assert(row == this->n);

        // Pack the data constructor: convert the LinkedList array to csr format
        matrixBuilder->pack();

    } catch(bad_alloc) {
        delete matrixBuilder; matrixBuilder = NULL;
        exit(1);
    } catch(General_Exception excep) {
        cout << "General exception" << endl;
        cout << excep.p << endl;
        delete matrixBuilder; matrixBuilder = NULL;
        exit(1);
    }

    // Create the actual matrix
    this->matrix = NULL;

    try {
        matrix = new MatrixIter(*matrixBuilder);
        assert(matrix != NULL);
        delete matrixBuilder; 
        matrixBuilder = NULL;
    } catch(bad_alloc) {
        delete matrixBuilder; matrixBuilder = NULL;
        exit(1);
    } catch(General_Exception excep) {
        delete matrixBuilder; matrixBuilder= NULL;
        exit(1);
    } 
}

/**
 * Base constructor. Initialize the PCG pressure solver.
*/
PressureSolver3D::PressureSolver3D(int nx, double *x, int ny, double *y, int nz, double *z) {
    // Counters and increments
    int i;

    this->nx = nx;
    this->ny = ny;
    this->nz = nz;

    this->n = (nx+2)*(ny+2)*(nz+2);

    this->x = new double[nx+1];
    this->y = new double[ny+1];
    this->z = new double[nz+1];

    for (i = 0; i < nx+1; i++)
        this->x[i] = x[i];

    for (i = 0; i < ny+1; i++)
        this->y[i] = y[i];
    
    for (i = 0; i < nz+1; i++)
        this->z[i] = z[i];
    
    this->pFlat = new double[this->n];
    this->tol = new double[this->n];

    // Set the update tolerances to 0, forces use of risidual reduction criteria
    for (i = 0; i < n; i++) {
        this->tol[i] = 0.0;
    }

    this->hx = this->x[1] - this->x[0];
    this->hy = this->y[1] - this->y[0];
    this->hz = this->z[1] - this->z[0];

    // Set up the memory for the PCG solver
    this->setUpMatrixStruct(nx, ny, nz);
}

void PressureSolver3D::discC(int row) {
    int colIndex;
    double hxsq = simutils::square(this->hx);
    double hysq = simutils::square(this->hy);
    double hzsq = simutils::square(this->hz);

    for (int i = matrix->rowBegin(row); i < matrix->rowEndPlusOne(row); i++) {
        colIndex = matrix->getColIndex(i);

        if (colIndex == row) {
            matrix->aValue(i) = -2*(1/hxsq + 1/hysq + 1/hzsq);
        } else if (colIndex == row-1 || colIndex == row+1) {
            matrix->aValue(i) = 1/hxsq;
        } else if (colIndex == row-(this->nx+2) || colIndex == row+(this->nx+2)) {
            matrix->aValue(i) = 1/hysq;
        } else if (colIndex == row-((this->nx+2)*(this->ny+2)) || colIndex == row+((this->nx+2)*(this->ny+2))) {
            matrix->aValue(i) = 1/hzsq;
        } else {
            matrix->aValue(i) = 0;
        }
    }
}

void PressureSolver3D::discNeumannBC(int row) {
    int colIndex;
    for (int i = matrix->rowBegin(row); i < matrix->rowEndPlusOne(row); i++) {
        colIndex = matrix->getColIndex(i);

        // Just avoiding a singular matrix, everything else is set to 0.
        if (colIndex == row) {
            matrix->aValue(i) = 1;
        } else {
            matrix->aValue(i) = 0;
        }
    }
}

/**
 * Set up the pressure matrix.
 * 
 * TODO: include boundary information. Using the helper functions with appropriately set boundary conditions
 * should make it very easy.
*/
void PressureSolver3D::setUpPressureMatrix(Pool3D *pool) {
    int i, j, k;
    int row = 0;
    objects::FSIObject obj;
    int mo = 1;

    // BOTTOM
    for (j = 0; j < ny+2; j++) {
        for (i = 0; i < nx+2; i++) {
            this->discNeumannBC(row);
            row++;
        }
    }

    // INTERNAL
    for (k = 1; k < nz+1; k++) {
        for (j = 0; j < ny+2; j++) {
            for (i = 0; i < nx+2; i++) {
                // Neumann BC in x, y
                if (j == 0 || j == ny+1 || i == 0 || i == nx+1) {
                    this->discNeumannBC(row);
                } else {
                    obj = pool->objAtIndex(i-mo, j-mo, k-mo);
                    if (pool->isInterface(obj)) {
                        this->discNeumannBC(row);
                    } else if (obj == objects::FLUID_C) {
                        this->discC(row);
                    } else {
                        this->discNeumannBC(row);
                    }
                }
                row++;
            }
        }
    }

    // TOP
    for (j = 0; j < ny+2; j++) {
        for (i = 0; i < nx+2; i++) {
            this->discNeumannBC(row);
            row++;
        }
    }

    assert(row == this->n);
}

/**
 * Set up the pressure matrix.
 * 
 * We only need to compute this for the fluid cells.
 * 
*/
void PressureSolver3D::setUpPressureRHS(double dt, double ***FU, double ***FV,
                                        double ***FW, double ***p, Pool3D *pool) {
    int i, j, k;

    int mo = 1;
    int row;
    int nDir[3];

    objects::FSIObject obj;

    try {

        row = 0;

        // Set the boundary points
        for (j = 0; j < this->ny+2; j++) {
            for (i = 0; i < this->nx+2; i++) {
                (this->matrix)->bValue(row) = p[1][j][i];
                row++;
            }
        }

        double nDirNorm[3];
        double velVec[3];
        double rhs, h;

        // Set the internal points
        for (k = 1; k < this->nz+1; k++) {
            for (j = 0; j < this->ny+2; j++) {
                for (i = 0; i < this->nx+2; i++) {

                    if (i == 0) {
                        (this->matrix)->bValue(row) = p[k][j][1];
                    } else if (i == nx+1) {
                        (this->matrix)->bValue(row) = p[k][j][nx];
                    } else if (j == 0) {
                        (this->matrix)->bValue(row) = p[k][1][i];
                    } else if (j == ny+1) {
                        (this->matrix)->bValue(row) = p[k][ny][i];
                    } else {
                        // Note: this code assumes that the boundary strip cannot be an interface (one possible source of weirdness)
                        obj = pool->objAtIndex(i-mo, j-mo, k-mo);

                        if (pool->isInterface(obj)) {
                            // Use negative pressure value in normal direction. (old)
                            pool->getNormalDir(obj, nDir);
                            // Get the normal unit vector.
                            nDirNorm[0] = nDir[0];
                            nDirNorm[1] = nDir[1];
                            nDirNorm[2] = nDir[2];
                            simutils::normalize3D(nDirNorm);

                            // Take the velocity of the solid on this interface point.
                            velVec[0] = FU[k+nDir[2]][j+nDir[1]][i+nDir[0]] - pool->getObjU(i-mo, j-mo, k-mo);
                            velVec[1] = FV[k+nDir[2]][j+nDir[1]][i+nDir[0]] - pool->getObjV(i-mo, j-mo, k-mo);
                            velVec[2] = FW[k+nDir[2]][j+nDir[1]][i+nDir[0]] - pool->getObjW(i-mo, j-mo, k-mo);

                            // Take dot product of structural velocity and the unit normal to apply the condition
                            rhs = simutils::ddot3d(velVec, nDirNorm);

                            h = sqrt(simutils::square(hx*nDir[0]) + simutils::square(hy*nDir[1]) + simutils::square(hz*nDir[2]));
                            (this->matrix)->bValue(row) = p[k+nDir[2]][j+nDir[1]][i+nDir[0]] - (h/dt)*rhs;
                        } else if (obj == objects::FLUID_C) {
                            (this->matrix)->bValue(row) = (
                                (FU[k][j][i] - FU[k][j][i-1])/this->hx
                                + (FV[k][j][i] - FV[k][j-1][i])/this->hy
                                + (FW[k][j][i] - FW[k-1][j][i])/this->hz)/dt;
                        } else {
                            (this->matrix)->bValue(row) = 0.0;
                        }
                    }
                    row++;
                }
            }
        }

        for (j = 0; j < this->ny+2; j++) {
            for (i = 0; i < this->nx+2; i++) {
                (this->matrix)->bValue(row) = p[nz][j][i];
                row++;
            }
        }
        assert(row == this->n);

    } catch(bad_alloc) {
        exit(1);
    } catch(General_Exception excep) {
        cout << excep.p << endl;
        exit(1);
    }
}

/**
 * Solve the pressure equation using PCG.
*/
void PressureSolver3D::solvePCG(double dt, Pool3D *pool, double ***FU, double ***FV,
                                double ***FW, double ***p, bool matrixRecompute,
                                ParamIter &params) {
    // Compute the matrix. Right now the iluRecompute flag
    // is being used to indicate setting up the pressure matrix at all
    if (matrixRecompute) {
        this->setUpPressureMatrix(pool);

        // Perform the ILU factorization
        if (params.drop_ilu == 0) {
            matrix->sfac(params);
        }
    }

    // Set up the rhs of the linear system
    this->setUpPressureRHS(dt, FU, FV, FW, p, pool);
    // assert(false);

    // Copy the pressure to the flat matrix
    int row = 0;
    for (int k = 0; k < nz+2; k++) {
        for (int j = 0; j < ny+2; j++) {
            for (int i = 0; i < nx+2; i++) {
                this->pFlat[row] = p[k][j][i];
                row++;
            }
        }
    }

    try {
        for (int i = 0; i < n; i++) {
            this->tol[i] = 0; // set update tolerance to 0. Forces the use of risidual reduction criteria.
                              // TODO: find out what this means.
                              // TODO: do this in the constructor.
        }
        matrix->set_toler(this->tol); // Set the update tolerance

        int niter; // The number of iterations done in the solve

        // Solve the linear system
        matrix->solve(params, this->pFlat, niter);

        cout << "Solved equation in " << niter << endl;

        assert(niter >= 0);
    }
    catch(bad_alloc) {
        cout << "line 376" << endl;
        exit(1);
    }
    catch(General_Exception excep) {
        cout << "line 380" << endl;
        exit(1);
    }

    // Copy the solution to the pressure matrix
    row = 0;
    double accum = 0.0;
    double tot = 0.0;
    for (int k = 0; k < nz+2; k++) {
        for (int j = 0; j < ny+2; j++) {
            for (int i = 0; i < nx+2; i++) {
                p[k][j][i] = this->pFlat[row];
                accum += p[k][j][i];
                tot += 1.0;
                row++;
            }
        }
    }
    cout << "Mean pressure = " << accum/tot << endl;
}

PressureSolver3D::~PressureSolver3D() {
    delete[] this->x;
    delete[] this->y;
    delete[] this->z;
    delete[] this->pFlat;
    delete[] this->tol;
    delete matrix;
}
