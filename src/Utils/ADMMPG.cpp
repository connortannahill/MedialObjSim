#include <Eigen/Sparse>
#include "ADMMPG.h"
#include <vector>
#include <unordered_map>
#include <string>
#include "Assembly.h"
#include <iostream>

using namespace std;

ADMMPG::ADMMPG(double dt, Assembly &a) {
    this->nSteps = 0;
    this->a = &a;
    this->dt = dt;

    // Allocate and assign initial values
    x = new Eigen::VectorXd(a.getD()*a.getNPnts());
    xPrev = new Eigen::VectorXd(a.getD()*a.getNPnts());
    xBar = new Eigen::VectorXd(a.getD()*a.getNPnts());

    // Assign initial values to x and z (xBar must be assigned at each step)
    (this->a)->copyX(*x);
    (this->a)->copyX(*xPrev);
    (this->a)->copyX(*xBar);

    // Compute the initial z vlaue
    z = new Eigen::VectorXd(*((this->a)->D) * (*this->x));

    uBar = new Eigen::VectorXd(Eigen::VectorXd::Constant(z->size(), 0.0));

    // Prefactor the matrix (assuming constant matrix for now)
    double dtsq = dt*dt;
    Eigen::SparseMatrix<double> t(*(this->a)->M);
    WD_T = new Eigen::SparseMatrix<double>(( *((this->a)->W) * *((this->a)->D)).transpose());

    t = t + dtsq*((*WD_T * (*(this->a)->W) * (*(this->a)->D)));

    // Compute the sparse Cholseky factorization.
    cgSol = new Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>>(t);

    DXpU = new Eigen::VectorXd(*z);
    vec = new Eigen::VectorXd(*z);
}

/**
 * Assumes a constant mass matrix for now
*/
void ADMMPG::step(int nIters, double tol) {
    // Get xBar, the predicted (explicit) location of the nodes independent of constraints
    double dtsq = dt*dt;

    // cout << "difference between x and xPrev " << (*x - *xPrev).norm() << endl;
    // assert(false);

    if (nSteps == 0) {
        a->copyX(*this->x);
        *this->xPrev = *this->x;
    }

    // Make prediction for next value of x (sorta)
    (this->a)->predictX(dt, *this->xPrev, *this->x, *this->xBar);

    double w = (*a->W).coeffRef(0,0);

    *xPrev = *x;
    *x = *xBar;
    *z = (*a->D) * *x;

    uBar->setZero();

    int i;
    double primalRes = 0.0;
    for (i = 0; i < nIters; i++) {
        // Update z_{n+1} using the assembly prox algorithm
        *DXpU = (*(a->D))*(*x) + (*uBar);
        a->prox(dt, *xPrev, *DXpU, *z);

        // Update the Lagrange multiplier uBar^{n+1}
        *uBar = *DXpU - *z;

        // Update the solution x^{n+1}
        *vec =  ((*a->M) * (*xBar)) + dtsq*(( *WD_T * ((*a->W) * (*z - *uBar))));
        *x = cgSol->solve(*vec);

        // Compute the primal residual. If it is beneath the tolerance, exit
        primalRes = w*((*a->D)*(*x) - *z).norm();
        if (primalRes < tol) {
            cout << "Converged after " << i+1 << "steps" << endl;
            break;
        }
    }

    if (!(primalRes < tol)) {
        cout << "failed to converge" << endl;
    }

    // Update the assembly using the new locations
    a->updateAfterStep(dt, *xPrev, *x);

    nSteps++;
}

ADMMPG::~ADMMPG() {
    delete x;
    delete xBar;
    delete z;
    delete xPrev;
    delete cgSol;
    delete WD_T;
    delete DXpU;
    delete vec;
}