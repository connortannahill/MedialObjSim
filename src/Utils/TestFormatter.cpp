#include "TestFormatter.h"
#include <iostream>
#include <filesystem>
#include <assert.h>

using namespace std;

#if __APPLE__
namespace fs = std::__fs::filesystem;
#else /* */
namespace fs = std::filesystem;
#endif

TestFormatter::TestFormatter(const char *testName) {

    this->testDir = new string("./output/");
    testDir->append(testName);
    testDir->append("/");

    // Make the test directory if it does not already exist
    cout << "testdir = " << *testDir << endl;
    if (!fs::exists(*testDir)) {
        fs::create_directory(*testDir);
    }
    // assert(false);
}

TestFormatter::TestFormatter(const char *testName, int testNum) {

    this->testDir = new string("./output/");
    testDir->append(testName);
    testDir->append("/");
    cout << "testDir first = " << *this->testDir << endl;

    // Make the test directory if it does not already exist
    cout << "testdir = " << *testDir << endl;
    if (!fs::exists(*testDir)) {
        fs::create_directory(*testDir);
    }

    // string num(testNum);
    testDir->append(to_string(testNum));
    testDir->append("/");
    cout << "testDir second = " << *this->testDir << endl;
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