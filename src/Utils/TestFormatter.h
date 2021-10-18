#ifndef TEST_FORMATTER_H
#define TEST_FORMATTER_H

#include <string>
#include <iostream>

using namespace std;

class TestFormatter {
private:
    string *testDir;
public:
    TestFormatter(const char *testName);
    TestFormatter(const char *testName, int testNum);
    // void outputData(void (*outFun)(const char *), const char *fname);
    void genOutStr(const char *fname, string &out);
    ~TestFormatter();
};

#endif