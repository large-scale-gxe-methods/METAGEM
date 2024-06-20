#ifndef READPARAMETERS_H
#define READPARAMETERS_H

#include <string>
#include <stdio.h>
#include <stdlib.h>

class CommandLine {
public:

    std::vector<std::string> fileNames;
    std::vector<std::string> intNames;
    std::vector<std::string> lcIntNames;
    std::vector<std::string> intNames2;
    std::vector<std::string> lcIntNames2;
    std::vector<std::string> intNames3;
    std::vector<std::string> lcIntNames3;
    std::vector<std::string> additionalJointInfo;
    std::vector<std::string> additionalInteractionInfo;

    std::string outFile;
    std::string outFile2;
    std::string outFile3;
    std::string metaFileList;

    size_t nInt;
    size_t nInt2 = 1;
    size_t nInt3 = 1;
    size_t nCov;
    int metaOpt;

    bool mb = false;
    bool rb = false;
    bool additionalJoint = false;
    bool additionalInteraction = false;
    void processCommandLine(int argc, char* argv[]);
};


#endif
