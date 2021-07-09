#include "SolidParams.h"
#include <unordered_map>
#include <string>

SolidParams::SolidParams() {
    intParamMap.clear();
    doubleParamMap.clear();
}

// SolidParams::SolidParams(const SolidParams &s) {
//     this->intParamMap = s.intParamMap;
//     this->doubleParamMap = s.doubleParamMap;
// }

void SolidParams::addParam(string &paramName, double val) {
    doubleParamMap[paramName] = val;
}

void SolidParams::addParam(string &paramName, int val) {
    // doubleParamMap[paramName] = val;
    intParamMap[paramName] = val;
}

void SolidParams::getParam(string &paramName, int &paramVal) {
    paramVal = intParamMap[paramName];
}

void SolidParams::getParam(string &paramName, double &paramVal) {
    paramVal = doubleParamMap[paramName];
}

void SolidParams::addParam(const char *paramName, double val) {
    string temp(paramName);
    this->addParam(temp, val);

}

void SolidParams::addParam(const char *paramName, int val) {
    string temp(paramName);
    this->addParam(temp, val);

}

void SolidParams::getParam(const char *paramName, double &paramVal) {
    string temp(paramName);
    this->getParam(temp, paramVal);

}

void SolidParams::getParam(const char *paramName, int &paramVal) {
    string temp(paramName);
    this->getParam(temp, paramVal);
}

