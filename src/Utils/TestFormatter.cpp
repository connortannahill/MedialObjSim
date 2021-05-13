#include "TestFormatter.h"
#include <iostream>
#include <filesystem>

using namespace std;

namespace fs = std::__fs::filesystem;

TestFormatter::TestFormatter(const char *testName) {

    this->testDir = new string("./output/");
    testDir->append(testName);
    testDir->append("/");

    // Make the test directory if it does not already exist
    if (!fs::exists(*testDir)) {
        fs::create_directory(*testDir);
    }
}

// void TestFormatter::outputData(void (*outFun)(const char *), const char *fname) {
//     (outFun)(fname);
// }

void TestFormatter::genOutStr(const char *fName, string &out) {
    out.clear();

    out += *testDir;
    out += fName;
}

TestFormatter::~TestFormatter() {
    delete testDir;
}