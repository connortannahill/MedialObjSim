#include "SimParams3D.h"

void SimParams3D::setRe(double Re) {
    this->Re = Re;
    reSet = true;
}

void SimParams3D::setMssNx(int mssNx) {
    this->mssNx = mssNx;
    mssNxSet = true;
}

void SimParams3D::setMssNy(int mssNy) {
    this->mssNy = mssNy;
    mssNySet = true;
}

void SimParams3D::setMssNz(int mssNy) {
    this->mssNy = mssNy;
    mssNySet = true;
}

void SimParams3D::setNx(int nx) {
    this->nx = nx;

    if (!mssNxSet) {
        mssNx = nx;
    }

    nxSet = true;
}

void SimParams3D::setNy(int ny) {
    this->ny = ny;

    if (!mssNySet) {
        mssNy = ny;
    }

    nySet = true;
}

void SimParams3D::setNz(int nz) {
    this->nz = nz;

    if (!mssNzSet) {
        mssNz = nz;
    }

    nzSet = true;
}

void SimParams3D::setDtFix(double dtFix) {
    this->dtFix = dtFix;
    dtFixSet = true;
}

void SimParams3D::setMu(double mu) {
    this->mu = mu;
    muSet = true;
}

bool SimParams3D::checkParams() {
    return nxSet && nySet && nzSet && reSet && muSet && repulseModeSet && repulseDistSet
            && collisionStiffnessSet && collisionDistSet;
}

void SimParams3D::setRepulseMode(int mode) {
    this->repulseMode = mode;
    this->repulseModeSet = true;
}

void SimParams3D::setRepulseDist(double dist) {
    this->repulseDist = dist;
    this->repulseDistSet = true;
}

void SimParams3D::setCollisionStiffness(double stiffness) {
    this->collisionDist = stiffness;
    this->collisionStiffnessSet = true;
}

void SimParams3D::setCollisionDist(double colDist) {
    this->collisionDist = colDist;
    this->collisionDistSet = true;
}

void SimParams3D::setAdmmTol(double admmTol) {
    this->admmTol = admmTol;
}

void SimParams3D::setUpdateMode(int mode) {
    this->updateMode = mode;
}

void SimParams3D::setElementMode(int mode) {
    this->elementMode = mode;
}