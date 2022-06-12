#ifndef PRINT_H
#define PRINT_H

void printWelcome();

void printHeaderMissingError(std::string fileName, std::string column);

void printNColumnError(int z, int n, std::string fileName, int nheader);

void printZeroVariantsError(std::string fileName);

void printOpenFileError(std::string fileName);

void printProcessingFiles();

void printProcessedFiles(std::string fileName, int n);

void printMetaBegin(int nFiles, int nvars);

void printDone(int nbs);

void printOutputLocation(std::string outFile);

void printTimeCompleted(double wall0, double wall1, double cpu0, double cpu1);

#endif
