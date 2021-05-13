#include "Assembly.h"
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <unordered_map>
#include <iostream>

using namespace std;

Assembly::Assembly() {
    // No arrays have been allocated
    dAlloc = false;
    wAlloc = false;
    mAlloc = false;
    nPntsSet = false;
    dSet = false;
}

int Assembly::getD() {
    if (!dSet) {
        cout << "d not set!" << endl;
        assert(false);
    }
    return this->d;
}

int Assembly::getNPnts() {
    if (!nPntsSet) {
        cout << "nPnts not set!" << endl;
        assert(false);
    }
    return this->nPnts;
}
void Assembly::setNPnts(int nPnts) {
    nPntsSet = true;
    this->nPnts = nPnts;
}
void Assembly::setD(int d) {
    dSet = true;
    this->d = d;
}

void Assembly::buildMassMatrix(Eigen::VectorXd &m) {
    if (!nPntsSet) {
        cout << "nPnts not set!" << endl;
        assert(false);
    }
    if (!mAlloc) {
        M = new Eigen::SparseMatrix<double>(m.size(), m.size());
        mAlloc = true;
    } else {
        assert(false);
    }
    M->reserve(Eigen::VectorXd::Constant(m.size(), 1));

    for (int i = 0; i < m.size(); i++) {
        M->insert(i, i) = m[i];
    }
}

/**
 * Default implementation, D is the identity matrix
*/
void Assembly::buildDMatrix() {
    if (!dSet) {
        cout << "d not set!" << endl;
        assert(false);
    }

    if (!nPntsSet) {
        cout << "nPnts not set!" << endl;
        assert(false);
    }
    if (!dAlloc) {
        D = new Eigen::SparseMatrix<double>(d*nPnts, d*nPnts);
        dAlloc = true;
    } else {
        assert(false);
    }

    D->reserve(Eigen::VectorXd::Constant(d*nPnts, 1));

    for (int i = 0; i < D->rows(); i++) {
        D->insert(i, i) = 1.0;
    }
}

/**
 * Default implementation, W is the identity matrix
*/
void Assembly::buildWMatrix(double w) {
    this->w = w;
    if (!nPntsSet) {
        cout << "nPnts not set!" << endl;
        assert(false);
    }

    if (!dSet) {
        cout << "d not set!" << endl;
        assert(false);
    }

    if (!wAlloc) {
        W = new Eigen::SparseMatrix<double>(d*nPnts, d*nPnts);
        wAlloc = true;
    } else {
        assert(false);
    }

    W->reserve(Eigen::VectorXd::Constant(d*nPnts, 1));

    for (int i = 0; i < W->rows(); i++) {
        W->insert(i, i) = w;
    }
}

Assembly::~Assembly() {
    if (mAlloc) {
        delete M;
    }
    if (dAlloc) {
        delete D;
    }
    if (wAlloc) {
        delete W;
    }
}