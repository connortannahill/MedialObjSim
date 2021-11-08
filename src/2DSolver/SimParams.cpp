#include "SimParams.h"
#include <math.h>

void SimParams::setRe(double Re) {
    this->Re = Re;
    reSet = true;
}

void SimParams::setMssNx(int mssNx) {
    this->mssNx = mssNx;
    mssNxSet = true;
}

void SimParams::setMssNy(int mssNy) {
    this->mssNy = mssNy;
    mssNySet = true;
}

void SimParams::setNx(int nx) {
    this->nx = nx;

    if (!mssNxSet) {
        mssNx = nx;
    }

    nxSet = true;
}

void SimParams::setUseEno(bool use) {
    this->useEno = use;
}

void SimParams::setNy(int ny) {
    this->ny = ny;

    if (!mssNySet) {
        mssNy = ny;
    }

    nySet = true;
}

void SimParams::setDtFix(double dtFix) {
    this->dtFix = dtFix;
    dtFixSet = true;
}

void SimParams::setMu(double mu) {
    this->mu = mu;
    muSet = true;
}

bool SimParams::checkParams() {
    return nxSet && nySet && reSet && muSet && repulseModeSet
            && collisionStiffnessSet
            && collisionDistSet;
}

void SimParams::setRepulseMode(int mode) {
    this->repulseMode = mode;
    this->repulseModeSet = true;
}


void SimParams::setGx(double gx) {
    this->gx = gx;
}
void SimParams::setGy(double gy) {
    this->gy = gy;
}

double SimParams::getGx() {
    return this->gx;
}
double SimParams::getGy() {
    return this->gy;
}

void SimParams::setU0(double u0) {
    this->u0 = u0;
}
void SimParams::setV0(double v0) {
    this->v0 = v0;
}

double SimParams::getU0() {
    return this->u0;
}
double SimParams::getV0() {
    return this->v0;
}

void SimParams::setCollisionStiffness(double stiffness) {
    this->collisionStiffness = stiffness;
    this->collisionStiffnessSet = true;
}

void SimParams::setCollisionDist(double colDist) {
    this->collisionDist = colDist;
    this->collisionDistSet = true;
}

void SimParams::setAdmmTol(double admmTol) {
    this->admmTol = admmTol;

}

void SimParams::setUpdateMode(int mode) {
    this->updateMode = mode;

}

void SimParams::setElementMode(int mode) {
    this->elementMode = mode;
}
